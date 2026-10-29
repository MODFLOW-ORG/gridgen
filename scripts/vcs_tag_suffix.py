"""
Print the suffix meson adds to the version of a development build.

The suffix is '+shortsha', '+shortsha.dirty' with uncommitted changes, or
empty for a release version (no 'dev'), when HEAD is a tag, or when git is
not available.
"""

import subprocess
import sys


def get_suffix():
    try:
        subprocess.run(["git", "status"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, OSError):
        return ""
    try:
        subprocess.check_output(
            ["git", "describe", "--exact-match", "--tags", "HEAD"],
            stderr=subprocess.DEVNULL,
        )
        return ""
    except subprocess.CalledProcessError:
        sha = (
            subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
            .decode()
            .strip()
        )
        dirty = (
            subprocess.run(["git", "diff", "--quiet", "HEAD"], capture_output=True)
            .returncode
            != 0
        )
        return f"+{sha}.dirty" if dirty else f"+{sha}"


if __name__ == "__main__":
    version = sys.argv[1] if len(sys.argv) > 1 else ""
    print(get_suffix() if "dev" in version else "", end="")
