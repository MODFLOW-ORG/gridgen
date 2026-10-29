# AGENTS.md

Guidance for AI agents working in the GRIDGEN repository. Keep it short and
current; when a command or convention changes, update this file in the same
pull request.

## Quick commands

Run from the project root. Pixi is the canonical workflow and matches CI.

- **Configure the build directory:** `pixi run setup`
- **Build** (installs to `bin/`): `pixi run build`
- **Run the tests** (smoke, unit, and run-twice determinism): `pixi run test`
- **Run one unit test:** `builddir/test/unit/unit_tests --gtest_filter='QuadTree.*'`
- **Compare outputs with main** (needs a build of main at `../gridgen-main`,
  e.g. `git worktree add ../gridgen-main main`, then `pixi run setup-main`
  and `pixi run build-main`): `pixi run compare-main`
- **Show the OS tag used for release archives:** `pixi run get-ostag`

## Repository orientation

- `src/` — C++ source: `basic/` (geometry), `grid/` (MODFLOW and quadtree
  grids, connectivity), `intersection/`, `tools/` (definition file parsing,
  exporters, ghost nodes, shapefiles, vertical pass-through), and the
  `gridgenusg.cpp` driver.
- `test/` — definition (`.dfn`) files and input shapefiles for the test cases;
  outputs are written in place and ignored by git. `test/compare_outputs.py`
  runs the cases (its `CASES` list) and compares every output file;
  `test/unit/` holds the GoogleTest unit tests.
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
- **Pull request review:** after the first green CI run, have an agent on a
  smaller model review the diff, given the pull request's intent and the files
  but not the author's conclusions, and post its findings as a `## Review`
  comment. Verify each finding before acting: fix the ones that hold and
  answer the rest with the reason in a `## Response` comment. Repeat until a
  review has no actionable findings, at most three rounds; findings still in
  dispute are left for the maintainer to decide.
- **Commit messages:** conventional commits (`fix:`, `feat:`, `refactor:`,
  `test:`, `docs:`, `ci:`, `chore:`).
- **Output comparison:** CI requires every output file to match a build of
  main byte for byte. A pull request that changes output on purpose gets the
  `output-change` label and says what changed and why.
- **Known defects in tests:** a known defect is recorded as an entry in
  `KNOWN_UNSTABLE` (`test/compare_outputs.py`) or a `DISABLED_` unit test, with
  a comment naming its issue. The pull request that fixes the issue removes the
  entry or enables the test, and adds the tests listed in the issue.
- **Generated files:** do not commit test outputs, `builddir/`, or
  `subprojects/` downloads.
- **CI:** changes only to Markdown, `doc/`, issue templates, `dependabot.yml`,
  or `LICENSE` skip the build and tests (`paths-ignore` in
  `.github/workflows/ci.yml`); update that list when adding non-code files.
- **Tone:** reserve compliments for genuinely exceptional finds or ideas.
