"""
Check the names in src/ and test/unit/ against the convention in .clang-tidy.

clang-tidy reads the compile commands of a build directory, builddir-lint/,
which is configured on the first run. The vendored shapelib directory is skipped.
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "builddir-lint"
FILES = r"^(?!.*shapelib).*/(src|test/unit)/.*\.cpp$"


def run(cmd):
    subprocess.run(cmd, cwd=ROOT, check=True)


def main():
    if not (BUILD / "compile_commands.json").exists():
        run(["meson", "setup", "-Dunit_tests=enabled", str(BUILD)])
    # gridgenusg.cpp includes the generated version header
    run(["ninja", "-C", str(BUILD), "src/gridgen_version.h"])
    cmd = ["run-clang-tidy", "-p", str(BUILD), "-quiet", FILES]
    sys.exit(subprocess.run(cmd, cwd=ROOT).returncode)


if __name__ == "__main__":
    main()
