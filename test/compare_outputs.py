"""
Run the gridgen test cases and compare every output file byte for byte.

With --ref-exe the outputs of two executables (e.g., a branch build and a main
build) are compared; without it the same executable is run twice to detect
run-to-run differences. Both runs use the inputs in this directory.
"""

import argparse
import filecmp
import fnmatch
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent

# (action, definition file), run in order in one working directory
CASES = (
    ("buildqtg", "action01_buildqtg.dfn"),
    ("grid02qtg-to-usgdata", "action02_writeusgdata.dfn"),
    ("grid01mfg-to-polyshapefile", "action03_shapefile.dfn"),
    ("grid02qtg-to-polyshapefile", "action03_shapefile.dfn"),
    ("grid01mfg-to-pointshapefile", "action03_shapefile.dfn"),
    ("grid02qtg-to-pointshapefile", "action03_shapefile.dfn"),
    ("canal_grid02qtg_lay1_intersect", "action04_intersect.dfn"),
    ("chd_grid02qtg_lay1_intersect", "action04_intersect.dfn"),
    ("grid01mfg-to-vtkfile", "action05_vtkfile.dfn"),
    ("grid02qtg-to-vtkfile", "action05_vtkfile.dfn"),
    ("grid02qtg-to-vtkfilesv", "action05_vtkfile.dfn"),
    # three layers, point refinement by attribute, interpolated elevations
    ("buildqtg_multilayer", "action06_buildqtg_multilayer.dfn"),
    ("grid03qtg-to-usgdata", "action07_multilayer.dfn"),
    ("grid03qtg-to-usgdata-vp", "action07_multilayer.dfn"),
    ("grid03qtg-to-polyshapefile", "action07_multilayer.dfn"),
    ("wells_grid03qtg_lay3_intersect", "action07_multilayer.dfn"),
    ("buffer_grid03qtg_lay2_intersect", "action07_multilayer.dfn"),
    # rotated grid with cell sizes that are not exact in binary
    ("buildqtg_rotated", "action08_buildqtg_rotated.dfn"),
    ("grid04qtg-to-usgdata", "action09_rotated.dfn"),
    ("grid04qtg-to-polyshapefile", "action09_rotated.dfn"),
    ("grid04qtg-to-vtkfile", "action09_rotated.dfn"),
    ("canal_grid04qtg_lay1_intersect", "action09_rotated.dfn"),
)

# generated files and test infrastructure that are not copied as inputs
SKIP_PATTERNS = (
    "*.tsf",
    "grid02qtg*.dat",
    "grid03qtg*",
    "grid04qtg*",
    "*.py",
    "__pycache__",
    "compare",
    "unit",
    "builddir",
)
# output path prefixes that change from run to run; reported but not counted
# as failures
# - polygon intersection cell order depends on set<Box*> order (#51)
KNOWN_UNSTABLE = (
    "output_intersection/buffer_grid03qtg_lay2_intersect.",
)

OUTPUT_DIRS = (
    "output_intersection",
    "output_shapefiles",
    "output_usgdata",
    "output_vtkfiles",
)


def _ignore(src, names):
    skip = []
    in_output_dir = Path(src).name in OUTPUT_DIRS
    for name in names:
        if any(fnmatch.fnmatch(name, p) for p in SKIP_PATTERNS):
            skip.append(name)
        elif in_output_dir and name != "readme.txt":
            skip.append(name)
    return skip


def run_cases(exe, workdir):
    """Copy the inputs to workdir and run all cases with exe."""
    shutil.copytree(TEST_DIR, workdir, ignore=_ignore)
    for action, dfn in CASES:
        before = _snapshot(workdir)
        result = subprocess.run(
            [str(exe), action, dfn],
            cwd=workdir,
            capture_output=True,
            text=True,
            errors="replace",
        )
        if result.returncode != 0:
            print(result.stdout)
            print(result.stderr, file=sys.stderr)
            raise RuntimeError(f"{exe} {action} {dfn} failed ({result.returncode})")
        # an action that exits 0 without writing anything did not run
        if _snapshot(workdir) == before:
            print(result.stdout)
            raise RuntimeError(f"{exe} {action} {dfn} wrote no output")


def _snapshot(workdir):
    """Size and modification time of every file under workdir."""
    snap = {}
    for p in Path(workdir).rglob("*"):
        if p.is_file():
            st = p.stat()
            snap[p] = (st.st_size, st.st_mtime_ns)
    return snap


def compare_dirs(dir0, dir1):
    """Return a list of files that differ or exist in only one directory."""
    files0 = {p.relative_to(dir0).as_posix() for p in dir0.rglob("*") if p.is_file()}
    files1 = {p.relative_to(dir1).as_posix() for p in dir1.rglob("*") if p.is_file()}
    diffs = [f"only in {dir0.name}: {f}" for f in sorted(files0 - files1)]
    diffs += [f"only in {dir1.name}: {f}" for f in sorted(files1 - files0)]
    for f in sorted(files0 & files1):
        if not filecmp.cmp(dir0 / f, dir1 / f, shallow=False):
            diffs.append(f"differs: {f}")
    return diffs, len(files0 | files1)


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--exe", required=True, help="gridgen executable to test")
    parser.add_argument("--ref-exe", help="reference gridgen executable")
    parser.add_argument(
        "--keep",
        type=Path,
        help="directory to keep the outputs in (default: temporary, removed)",
    )
    args = parser.parse_args()

    exe = Path(args.exe).resolve()
    ref_exe = Path(args.ref_exe).resolve() if args.ref_exe else exe
    for e in {exe, ref_exe}:
        if not e.exists() and not e.with_suffix(".exe").exists():
            parser.error(f"executable not found: {e}")

    root = args.keep if args.keep else Path(tempfile.mkdtemp(prefix="gridgen_"))
    if args.keep:
        for sub in ("ref", "test"):
            shutil.rmtree(root / sub, ignore_errors=True)
    try:
        run_cases(ref_exe, root / "ref")
        run_cases(exe, root / "test")
        diffs, nfiles = compare_dirs(root / "ref", root / "test")
    finally:
        if not args.keep:
            shutil.rmtree(root, ignore_errors=True)

    label = f"{exe} vs {ref_exe}" if args.ref_exe else f"{exe} run twice"
    # a known unstable file may differ but must still exist in both runs
    known = [
        d
        for d in diffs
        if d.startswith("differs: ") and d.split(": ", 1)[1].startswith(KNOWN_UNSTABLE)
    ]
    diffs = [d for d in diffs if d not in known]
    for d in known:
        print(f"known unstable output {d}")
    if diffs:
        print(f"{len(diffs)} of {nfiles} files differ ({label}):")
        for d in diffs:
            print(f"  {d}")
        return 1
    nsame = nfiles - len(known)
    print(f"{nsame} of {nfiles} files identical, {len(known)} known unstable ({label})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
