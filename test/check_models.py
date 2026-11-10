"""
Check gridgen output by running MODFLOW 6 and MODFLOW-USG on it.

Builds the grids of check_models_grids.dfn, writes their DISU data, and runs
each case with heads fixed on the cells along the edges of the active grid, so
that a head that varies linearly in x and y is the exact solution (in the
layered cases, on the cells next to them too; see FIX_CELLS_NEXT_TO_EDGE). A
case passes when the program terminates normally and the largest head error in
the other cells, which start at HEAD_AT_ORIGIN, is less than HEAD_TOLERANCE.
The largest flow between layers in those cells is reported but not checked. The
programs are installed with FloPy's get-modflow, from the latest release of
each, into bin/modflow.

Cases:
  - layered_gncv    : layers refined differently and an active domain in layer
                      2, with horizontal and vertical ghost nodes (DISU).
  - layered_vp_gncv : the same, with vertical pass-through cells (DISU).
  - same            : the same refinement in every layer (DISU, and DISV in
                      MODFLOW 6).
"""

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

import flopy
import numpy as np
from flopy.utils.get_modflow import run_main as get_modflow
from flopy.utils.gnc import get_gnc_dtype, get_gridprops_gnc5

TEST_DIRECTORY = Path(__file__).resolve().parent
WORKSPACE = TEST_DIRECTORY / "check_models"

# repository that get-modflow installs each MODFLOW program from
PROGRAMS = {"mf6": "modflow6", "mfusg": "mfusg", "mfusgt": "mfusgt"}

# gridgen blocks, in the order they are run, and the definition file of each
GRIDGEN_BLOCKS = (
    ("check_models_build_layered", "check_models_grids.dfn"),
    ("check_models_build_same", "check_models_grids.dfn"),
    ("check_models_layered_gncv", "check_models_layered_data.dfn"),
    ("check_models_layered_vp_gncv", "check_models_layered_data.dfn"),
    ("check_models_same_usgdata", "check_models_same_data.dfn"),
)

# case name and the discretizations it is run with
CASES = {
    "layered_gncv": ("disu",),
    "layered_vp_gncv": ("disu",),
    "same": ("disu", "disv"),
}

# ghost nodes at the edge of an active domain are missing or inexact (#170), so
# in these cases the cells connected to an edge cell have fixed heads too
FIX_CELLS_NEXT_TO_EDGE = ("layered_gncv", "layered_vp_gncv")

# origin of the grid and the head gradient in x and y; the heads vary by 86
# length units over the grid
X_ORIGIN = 543750.0
Y_ORIGIN = 2792250.0
X_GRADIENT = 0.001
Y_GRADIENT = -0.0005
HEAD_AT_ORIGIN = 100.0

# largest head error allowed, in length units
HEAD_TOLERANCE = 1.0e-4

HORIZONTAL_CONDUCTIVITY = 10.0
VERTICAL_CONDUCTIVITY = 1.0


def linear_head(x, y):
    return HEAD_AT_ORIGIN + X_GRADIENT * (x - X_ORIGIN) + Y_GRADIENT * (y - Y_ORIGIN)


def read_numbers(path, dtype=float):
    return np.array(Path(path).read_text().split(), dtype=dtype)


class UsgData:
    """The DISU data gridgen writes with one USG_DATA_PREFIX"""

    def __init__(self, prefix):
        self.prefix = prefix
        nod = np.loadtxt(f"{prefix}.nod", skiprows=1)
        self.x, self.y, z, self.dx, self.dy, self.thickness = (
            nod[:, i] for i in range(2, 8)
        )
        self.top = z + self.thickness / 2.0
        self.bottom = z - self.thickness / 2.0
        self.nodes = len(nod)
        self.ia = read_numbers(f"{prefix}.ia.dat", int) - 1
        self.ja = read_numbers(f"{prefix}.ja.dat", int) - 1
        self.iac = np.diff(self.ia)
        self.cl12 = read_numbers(f"{prefix}.c1.dat")
        self.fahl = read_numbers(f"{prefix}.fahl.dat")
        self.fldr = read_numbers(f"{prefix}.fldr.dat", int)
        self.area = read_numbers(f"{prefix}.area.dat")
        self.nodes_per_layer = read_numbers(f"{prefix}.nodesperlay.dat", int)
        self.ghost_nodes = np.atleast_1d(
            np.genfromtxt(f"{prefix}.gnc.dat", dtype=get_gnc_dtype(2))
        )
        for name in ("n", "m", "j0", "j1"):
            self.ghost_nodes[name] -= 1
        self.vertices = [
            [int(row[0]) - 1, row[1], row[2]]
            for row in np.loadtxt(f"{prefix}.vertices.dat")
        ]
        self.cell2d = []
        with open(f"{prefix}.cell2d.dat") as file:
            for line in file:
                values = line.split()
                self.cell2d.append(
                    [int(values[0]) - 1, float(values[1]), float(values[2])]
                    + [int(values[3])]
                    + [int(vertex) - 1 for vertex in values[4:]]
                )

        # node of each entry of ja, and the connections other than a cell to itself
        self.row = np.repeat(np.arange(self.nodes), self.iac)
        self.connection = self.ja != self.row
        self.vertical = self.connection & (np.abs(self.fldr) == 3)
        self.horizontal = self.connection & ~self.vertical

        # face width of a horizontal connection; the layers are flat, so the
        # face area is the width times the average thickness of the two cells
        self.hwva = self.fahl.copy()
        average_thickness = 0.5 * (self.thickness[self.row] + self.thickness[self.ja])
        self.hwva[self.horizontal] = (
            self.fahl[self.horizontal] / average_thickness[self.horizontal]
        )

        # a cell whose horizontal faces are not all shared with active cells is
        # on the edge of the grid or of an active domain
        shared_width = np.bincount(
            self.row[self.horizontal],
            weights=self.hwva[self.horizontal],
            minlength=self.nodes,
        )
        perimeter = 2.0 * (self.dx + self.dy)
        self.fixed = shared_width < perimeter * (1.0 - 1.0e-9)
        self.exact_head = linear_head(self.x, self.y)


def install_programs(bindir, programs, update):
    bindir.mkdir(parents=True, exist_ok=True)
    for program in programs:
        if update or shutil.which(program, path=str(bindir)) is None:
            get_modflow(str(bindir), repo=PROGRAMS[program], subset=program, quiet=True)
        if shutil.which(program, path=str(bindir)) is None:
            raise RuntimeError(f"get-modflow did not install {program} in {bindir}")


def run_gridgen(gridgen):
    WORKSPACE.mkdir(exist_ok=True)
    for block, definition_file in GRIDGEN_BLOCKS:
        result = subprocess.run(
            [str(gridgen), block, definition_file],
            cwd=TEST_DIRECTORY,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"gridgen {block} {definition_file} failed:\n{result.stderr}"
            )


def run_mf6(data, workspace, executable, discretization):
    name = "check"
    simulation = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=workspace, exe_name=executable
    )
    flopy.mf6.ModflowTdis(simulation)
    flopy.mf6.ModflowIms(
        simulation,
        complexity="simple",
        linear_acceleration="bicgstab",
        inner_dvclose=1.0e-9,
        outer_dvclose=1.0e-9,
        inner_maximum=5000,
        outer_maximum=50,
        rcloserecord=[1.0e-6, "strict"],
    )
    model = flopy.mf6.ModflowGwf(simulation, modelname=name, save_flows=True)
    gnc = data.ghost_nodes
    if discretization == "disu":
        flopy.mf6.ModflowGwfdisu(
            model,
            nodes=data.nodes,
            nja=len(data.ja),
            top=data.top,
            bot=data.bottom,
            area=data.area,
            iac=data.iac,
            ja=data.ja,
            ihc=np.where(data.vertical, 0, 1),
            cl12=data.cl12,
            hwva=data.hwva,
            nvert=len(data.vertices),
            vertices=data.vertices,
            cell2d=data.cell2d,
        )

        def cellid(node):
            return (int(node),)

        fixed = [cellid(node) for node in np.flatnonzero(data.fixed)]
    else:
        # every layer has the same cells, so the vertices and the CELL2D
        # records of the first layer describe each layer
        ncpl = int(data.nodes_per_layer[0])
        nlay = len(data.nodes_per_layer)
        flopy.mf6.ModflowGwfdisv(
            model,
            nlay=nlay,
            ncpl=ncpl,
            nvert=len(data.vertices),
            top=data.top[:ncpl],
            botm=data.bottom.reshape(nlay, ncpl),
            vertices=data.vertices,
            cell2d=data.cell2d[:ncpl],
        )

        def cellid(node):
            return (int(node) // ncpl, int(node) % ncpl)

        fixed = [cellid(node) for node in np.flatnonzero(data.fixed)]
    flopy.mf6.ModflowGwfnpf(
        model, k=HORIZONTAL_CONDUCTIVITY, k33=VERTICAL_CONDUCTIVITY, icelltype=0
    )
    flopy.mf6.ModflowGwfic(model, strt=data.starting_head)
    flopy.mf6.ModflowGwfchd(
        model,
        stress_period_data=[
            [cell, data.exact_head[node]]
            for cell, node in zip(fixed, np.flatnonzero(data.fixed))
        ],
    )
    flopy.mf6.ModflowGwfgnc(
        model,
        numgnc=len(gnc),
        numalphaj=2,
        gncdata=[
            [
                cellid(record["n"]),
                cellid(record["m"]),
                cellid(record["j0"]),
                cellid(record["j1"]),
                record["alpha0"],
                record["alpha1"],
            ]
            for record in gnc
        ],
    )
    flopy.mf6.ModflowGwfoc(
        model,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    simulation.write_simulation(silent=True)
    success, _ = simulation.run_simulation(silent=True)
    if not success:
        return False, None, None
    heads = model.output.head().get_data().ravel()
    flow = model.output.budget().get_data(text="FLOW-JA-FACE")[0].ravel()
    if discretization == "disu":
        row, column, vertical = data.row, data.ja, data.vertical
    else:
        # the connections are numbered by MODFLOW 6; a vertical one joins
        # the same cell in two layers
        grid = flopy.mf6.utils.MfGrdFile(Path(workspace) / f"{name}.disv.grb")
        row = np.repeat(np.arange(grid.nodes), np.diff(grid.ia))
        column = grid.ja
        vertical = np.abs(column - row) == ncpl
    return True, heads, largest_vertical_flow(data, row, column, vertical, flow)


def run_mfusg(data, workspace, executable):
    name = "check"
    model = flopy.mfusg.MfUsg(
        modelname=name, model_ws=workspace, exe_name=executable, structured=False
    )
    flopy.mfusg.MfUsgDisU(
        model,
        nodes=data.nodes,
        nlay=len(data.nodes_per_layer),
        njag=len(data.ja),
        nodelay=data.nodes_per_layer,
        top=data.top,
        bot=data.bottom,
        area=data.area,
        iac=data.iac.copy(),
        ja=data.ja.copy(),
        ivc=data.vertical.astype(int),
        cl12=data.cl12.copy(),
        fahl=data.fahl.copy(),
        idsymrd=0,
        ivsd=1,
    )
    flopy.mfusg.MfUsgBas(
        model, ibound=np.where(data.fixed, -1, 1), strt=data.starting_head
    )
    flopy.mfusg.MfUsgLpf(
        model,
        hk=HORIZONTAL_CONDUCTIVITY,
        vka=VERTICAL_CONDUCTIVITY,
        laytyp=0,
        ipakcb=50,
    )
    flopy.mfusg.MfUsgSms(
        model,
        hclose=1.0e-9,
        hiclose=1.0e-10,
        mxiter=50,
        iter1=5000,
        linmeth=2,
        iacl=2,
        norder=0,
        level=7,
        north=7,
        rrctol=0.0,
    )
    flopy.mfusg.MfUsgGnc(
        model, **get_gridprops_gnc5(data.ghost_nodes, ia=data.ia, ja=data.ja)
    )
    flopy.modflow.ModflowOc(
        model, stress_period_data={(0, 0): ["save head", "save budget"]}
    )
    model.write_input()
    success, _ = model.run_model(silent=True)
    if not success:
        return False, None, None
    heads = flopy.utils.HeadUFile(Path(workspace) / f"{name}.hds").get_data()
    heads = np.concatenate([np.asarray(layer).ravel() for layer in heads])
    flow = (
        flopy.utils.CellBudgetFile(Path(workspace) / f"{name}.cbc")
        .get_data(text="FLOW JA FACE")[0]
        .ravel()
    )
    return (
        True,
        heads,
        largest_vertical_flow(data, data.row, data.ja, data.vertical, flow),
    )


def largest_vertical_flow(data, row, column, vertical, flow):
    """Largest flow between two cells in different layers that are not fixed"""
    interior = vertical & ~data.fixed[row] & ~data.fixed[column]
    return np.abs(flow[interior]).max() if interior.any() else 0.0


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument(
        "--gridgen",
        default=str(TEST_DIRECTORY.parent / "bin" / "gridgen"),
        help="gridgen executable (default: bin/gridgen)",
    )
    parser.add_argument(
        "--bindir",
        default=str(TEST_DIRECTORY.parent / "bin" / "modflow"),
        help="folder for the MODFLOW programs (default: bin/modflow)",
    )
    parser.add_argument(
        "--update",
        action="store_true",
        help="install the latest release of each program even if it is present",
    )
    parser.add_argument(
        "--programs",
        nargs="+",
        choices=list(PROGRAMS),
        default=list(PROGRAMS),
        help="programs to run (default: all)",
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        choices=list(CASES),
        default=list(CASES),
        help="cases to run (default: all)",
    )
    args = parser.parse_args()

    bindir = Path(args.bindir).resolve()
    install_programs(bindir, args.programs, args.update)
    run_gridgen(Path(args.gridgen).resolve())

    failures = 0
    print(
        f"{'case':16s} {'program':8s} {'grid':5s} {'nodes':>7s} {'fixed':>6s} "
        f"{'seconds':>8s} {'head error':>11s} {'vertical flow':>14s}  result"
    )
    for case in args.cases:
        data = UsgData(WORKSPACE / case)
        if case in FIX_CELLS_NEXT_TO_EDGE:
            next_to_edge = data.ja[data.connection & data.fixed[data.row]]
            data.fixed[next_to_edge] = True
        # the cells that are not fixed start away from the solution, so that the
        # solver has to reach it
        data.starting_head = np.where(data.fixed, data.exact_head, HEAD_AT_ORIGIN)
        for program in args.programs:
            executable = shutil.which(program, path=str(bindir))
            for discretization in CASES[case]:
                if program != "mf6" and discretization == "disv":
                    continue
                workspace = WORKSPACE / f"{case}_{program}_{discretization}"
                start = time.perf_counter()
                if program == "mf6":
                    success, heads, flow = run_mf6(
                        data, workspace, executable, discretization
                    )
                else:
                    success, heads, flow = run_mfusg(data, workspace, executable)
                seconds = time.perf_counter() - start
                label = (
                    f"{case:16s} {program:8s} {discretization:5s} {data.nodes:7d} "
                    f"{int(data.fixed.sum()):6d} {seconds:8.1f}"
                )
                if not success:
                    failures += 1
                    print(
                        f"{label} {'':>11s} {'':>14s}  FAIL: did not terminate normally"
                    )
                    continue
                error = np.abs(heads - data.exact_head)[~data.fixed].max()
                passed = error < HEAD_TOLERANCE
                failures += not passed
                result = (
                    "pass" if passed else f"FAIL: error of {HEAD_TOLERANCE:g} or more"
                )
                print(f"{label} {error:11.3e} {flow:14.3e}  {result}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
