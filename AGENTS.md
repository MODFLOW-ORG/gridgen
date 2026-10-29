# AGENTS.md

Guidance for AI agents working in the GRIDGEN repository. Keep it short and
current; when a command or convention changes, update this file in the same
pull request.

## Quick commands

Run from the project root. Pixi is the canonical workflow and matches CI.

- **Configure the build directory:** `pixi run setup`
- **Build** (installs to `bin/`): `pixi run build`
- **Run the tests:** `pixi run test`
- **Show the OS tag used for release archives:** `pixi run get-ostag`

## Repository orientation

- `src/` — C++ source: `basic/` (geometry), `grid/` (MODFLOW and quadtree
  grids, connectivity), `intersection/`, `tools/` (definition file parsing,
  exporters, ghost nodes, shapefiles, vertical pass-through), and the
  `gridgenusg.cpp` driver.
- `test/` — definition (`.dfn`) files and input shapefiles for the test cases;
  outputs are written in place and ignored by git.
- `doc/` — the GRIDGEN report (OFR 2014-1109), release notes, and release note
  items.
- `builddir/` — disposable Meson build output; never edit by hand.

## Documentation

- The GRIDGEN report, `doc/ofr2014-1109.pdf`, documents every block and
  record; check changes to input or output against it.
- `doc/ReleaseNotes.pdf` lists corrections to the report (e.g., the
  `ROTATION_ANGLE` definition) and known issues.

## Conventions and guardrails

- **Issues and milestones:** every change has an issue with a label (`bug`,
  `enhancement`, `refactor`, `ci`, or `documentation`) and a version milestone.
  A bug issue lists the tests to add with the fix.
- **Pull requests:** one issue per pull request; fill in the pull request
  template and close the issue from its checklist.
- **Release notes:** a pull request that changes gridgen adds a TOML item in
  `doc/ReleaseNotes/items/` (see the README there).
- **Commit messages:** conventional commits (`fix:`, `feat:`, `refactor:`,
  `test:`, `docs:`, `ci:`, `chore:`).
- **Generated files:** do not commit test outputs or `builddir/`.
- **Tone:** reserve compliments for genuinely exceptional finds or ideas.
