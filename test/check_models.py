"""
Check gridgen output by running MODFLOW 6 and MODFLOW-USG on it.

Builds the grids of check_models_grids.dfn, writes their DISU data, and runs
each case with heads fixed on the cells along the edges of the active grid, so
that a head that varies linearly in x and y is the exact solution (in the
layered cases, on the cells next to them too; see FIX_CELLS_NEXT_TO_EDGE). A
case passes when the program terminates normally and the largest head error in
the other cells, which start at HEAD_AT_ORIGIN, is less than HEAD_TOLERANCE,
and, in MODFLOW 6, their specific discharge differs from the exact one by less
than DISCHARGE_TOLERANCE. The DISU input is gridgen's arrays, including IHC,
HWVA, ANGLDEGX, and IVC. The largest flow between layers in those cells is
reported but not checked. The programs are installed with FloPy's get-modflow,
from the latest release of each, into bin/modflow.

Cases:
  - layered_gncv    : layers refined differently and an active domain in layer
                      2, with horizontal and vertical ghost nodes (DISU).
  - layered_vp_gncv : the same, with vertical pass-through cells (DISU).
  - same            : the same refinement in every layer (DISU, and DISV in
                      MODFLOW 6).
"""

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import flopy
import numpy as np
from flopy.utils.get_modflow import flopy_appdata_path
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

# extent of the base grid, the thickness of each layer from the elevations of
# check_models_mfg, and the layers of each case with an active domain
X_EXTENT = 58 * 800.0
Y_EXTENT = 100 * 800.0
LAYER_THICKNESS = (20.0, 25.0, 20.0)
ACTIVE_DOMAIN_LAYERS = {"layered_gncv": (2,), "layered_vp_gncv": (2,), "same": ()}

HORIZONTAL_CONDUCTIVITY = 10.0
VERTICAL_CONDUCTIVITY = 1.0

# specific discharge in x, y, and z of the linear head field, and the largest
# error allowed in the specific discharge MODFLOW 6 saves
EXACT_DISCHARGE = (
    -HORIZONTAL_CONDUCTIVITY * X_GRADIENT,
    -HORIZONTAL_CONDUCTIVITY * Y_GRADIENT,
    0.0,
)
DISCHARGE_TOLERANCE = 1.0e-6


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
        self.ihc = read_numbers(f"{prefix}.ihc.dat", int)
        self.ivc = read_numbers(f"{prefix}.ivc.dat", int)
        self.hwva = read_numbers(f"{prefix}.hwva.dat")
        self.angldegx = read_numbers(f"{prefix}.angldegx.dat")
        self.area = read_numbers(f"{prefix}.area.dat")
        self.nodes_per_layer = read_numbers(f"{prefix}.nodesperlay.dat", int)
        self.layer = np.repeat(
            np.arange(1, len(self.nodes_per_layer) + 1), self.nodes_per_layer
        )
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


def check_geometry(data, case):
    """Messages for the data that differs from the definition of the grid"""
    messages = []
    if data.nodes_per_layer.sum() != data.nodes or len(data.layer) != data.nodes:
        return [f"{data.nodes_per_layer.sum()} cells per layer but {data.nodes} cells"]
    if not np.allclose(data.area, data.dx * data.dy, rtol=1.0e-12):
        messages.append("area is not dx times dy")
    for layer, thickness in enumerate(LAYER_THICKNESS, start=1):
        in_layer = data.layer == layer
        if not np.allclose(data.thickness[in_layer], thickness):
            messages.append(f"layer {layer} is not {thickness:g} thick")
        if layer in ACTIVE_DOMAIN_LAYERS[case]:
            continue
        # without an active domain the cells fill the base grid, and the cells
        # fixed are those on its edge
        if not np.isclose(data.area[in_layer].sum(), X_EXTENT * Y_EXTENT, rtol=1e-12):
            messages.append(f"the cells of layer {layer} do not fill the grid")
        x, y, dx, dy = (
            data.x - X_ORIGIN,
            data.y - Y_ORIGIN,
            data.dx / 2.0,
            data.dy / 2.0,
        )
        on_edge = (
            np.isclose(x - dx, 0.0)
            | np.isclose(x + dx, X_EXTENT)
            | np.isclose(y - dy, 0.0)
            | np.isclose(y + dy, Y_EXTENT)
        )
        if np.any(data.fixed[in_layer] != on_edge[in_layer]):
            messages.append(f"the fixed cells of layer {layer} are not its edge")

    # each connection is listed by both cells, with the same area and the
    # opposite direction
    key = data.row * data.nodes + data.ja
    order = np.argsort(key)
    reverse_key = data.ja * data.nodes + data.row
    position = np.minimum(np.searchsorted(key, reverse_key, sorter=order), len(key) - 1)
    reverse = order[position]
    connection = data.connection
    if not np.array_equal(key[reverse][connection], reverse_key[connection]):
        return [*messages, "a connection is listed by only one of its cells"]
    if not np.allclose(data.fahl[connection], data.fahl[reverse][connection]):
        messages.append("the two cells of a connection give different areas")
    if np.any(data.fldr[connection] != -data.fldr[reverse][connection]):
        messages.append("the two cells of a connection give the same direction")

    # IHC and IVC give the type of each connection; HWVA is the same from both
    # cells, and, since the layers are flat, a horizontal face area is the
    # width times the average thickness of the two cells
    if not np.array_equal(data.ihc, data.horizontal.astype(int)):
        messages.append("ihc is not 1 for exactly the horizontal connections")
    if not np.array_equal(data.ivc, data.vertical.astype(int)):
        messages.append("ivc is not 1 for exactly the vertical connections")
    if not np.array_equal(data.hwva[connection], data.hwva[reverse][connection]):
        messages.append("the two cells of a connection give different hwva")
    horizontal = data.horizontal
    average_thickness = 0.5 * (data.thickness[data.row] + data.thickness[data.ja])
    if not np.allclose(
        (data.hwva * average_thickness)[horizontal], data.fahl[horizontal]
    ) or not np.array_equal(data.hwva[data.vertical], data.fahl[data.vertical]):
        messages.append("hwva is not the face width or area")

    # the face normal points toward the center of the connected cell
    toward = np.arctan2(
        data.y[data.ja] - data.y[data.row], data.x[data.ja] - data.x[data.row]
    )
    if np.any(np.cos(np.radians(data.angldegx) - toward)[horizontal] <= 0.0):
        messages.append("angldegx points away from a connected cell")

    # a horizontal connection is half a cell from each center, and one between
    # adjacent layers is the area the cells share at the distance between their
    # centers; one across a vertical pass-through cell is not checked
    half_size = (
        np.where(np.abs(data.fldr) == 1, data.dx[data.row], data.dy[data.row]) / 2
    )
    if not np.allclose(data.cl12[horizontal], half_size[horizontal]):
        messages.append("a horizontal distance is not half the cell size")
    adjacent = data.vertical & (np.abs(data.layer[data.row] - data.layer[data.ja]) == 1)
    shared_area = np.minimum(data.area[data.row], data.area[data.ja])
    if not np.allclose(data.fahl[adjacent], shared_area[adjacent]):
        messages.append("a vertical area is not the area the cells share")
    center = (data.top + data.bottom) / 2.0
    distance = np.abs(center[data.row] - center[data.ja])
    if not np.allclose((data.cl12 + data.cl12[reverse])[adjacent], distance[adjacent]):
        messages.append("a vertical distance is not the distance between centers")
    return messages


def print_versions(bindir, programs):
    """Print the release get-modflow installed each program from"""
    # the file lists every get-modflow install of this user; one made without
    # a subset installed every program of its release
    try:
        installs = json.loads((flopy_appdata_path / "get_modflow.json").read_text())
    except (OSError, json.JSONDecodeError):
        installs = []
    for program in programs:
        releases = [
            install.get("release_id", "unknown release")
            for install in installs
            if Path(install.get("bindir", "")) == bindir
            and install.get("repo") == PROGRAMS[program]
            and program in install.get("subset", [program])
        ]
        release = releases[-1] if releases else "unknown release"
        print(f"{program} from {PROGRAMS[program]} {release}")


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
            ihc=data.ihc,
            cl12=data.cl12,
            hwva=data.hwva,
            angldegx=data.angldegx,
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
        model,
        k=HORIZONTAL_CONDUCTIVITY,
        k33=VERTICAL_CONDUCTIVITY,
        icelltype=0,
        save_specific_discharge=True,
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
        return False, None, None, None
    heads = model.output.head().get_data().ravel()
    budget = model.output.budget()
    flow = budget.get_data(text="FLOW-JA-FACE")[0].ravel()
    discharge = budget.get_data(text="DATA-SPDIS")[0]
    discharge_error = np.zeros(data.nodes)
    for component, exact in zip(("qx", "qy", "qz"), EXACT_DISCHARGE):
        discharge_error[discharge["node"] - 1] = np.maximum(
            discharge_error[discharge["node"] - 1], np.abs(discharge[component] - exact)
        )
    if discretization == "disu":
        row, column, vertical = data.row, data.ja, data.vertical
    else:
        # the connections are numbered by MODFLOW 6; a vertical one joins
        # the same cell in two layers
        grid = flopy.mf6.utils.MfGrdFile(Path(workspace) / f"{name}.disv.grb")
        row = np.repeat(np.arange(grid.nodes), np.diff(grid.ia))
        column = grid.ja
        vertical = np.abs(column - row) == ncpl
    return (
        True,
        heads,
        largest_vertical_flow(data, row, column, vertical, flow),
        discharge_error[~data.fixed].max(),
    )


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
        ivc=data.ivc,
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
        return False, None, None, None
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
        None,
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
    print_versions(bindir, args.programs)
    run_gridgen(Path(args.gridgen).resolve())

    failures = 0
    print(
        f"{'case':16s} {'program':8s} {'grid':5s} {'nodes':>7s} {'fixed':>6s} "
        f"{'seconds':>8s} {'head error':>11s} {'discharge error':>15s} "
        f"{'vertical flow':>14s}  result"
    )
    for case in args.cases:
        data = UsgData(WORKSPACE / case)
        messages = check_geometry(data, case)
        if messages:
            failures += 1
            print(f"{case:16s} FAIL: " + "; ".join(messages))
            continue
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
                    success, heads, flow, discharge_error = run_mf6(
                        data, workspace, executable, discretization
                    )
                else:
                    success, heads, flow, discharge_error = run_mfusg(
                        data, workspace, executable
                    )
                seconds = time.perf_counter() - start
                label = (
                    f"{case:16s} {program:8s} {discretization:5s} {data.nodes:7d} "
                    f"{int(data.fixed.sum()):6d} {seconds:8.1f}"
                )
                if not success:
                    failures += 1
                    print(
                        f"{label} {'':>11s} {'':>15s} {'':>14s}  "
                        "FAIL: did not terminate normally"
                    )
                    continue
                error = np.abs(heads - data.exact_head)[~data.fixed].max()
                failed = []
                if error >= HEAD_TOLERANCE:
                    failed.append(f"head error of {HEAD_TOLERANCE:g} or more")
                if (
                    discharge_error is not None
                    and discharge_error >= DISCHARGE_TOLERANCE
                ):
                    failed.append(f"discharge error of {DISCHARGE_TOLERANCE:g} or more")
                failures += bool(failed)
                result = "FAIL: " + "; ".join(failed) if failed else "pass"
                discharge = (
                    f"{discharge_error:15.3e}"
                    if discharge_error is not None
                    else f"{'':>15s}"
                )
                print(f"{label} {error:11.3e} {discharge} {flow:14.3e}  {result}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
