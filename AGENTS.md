# AGENTS.md

Guidance for AI agents working in the GRIDGEN repository. Keep it short and
current; when a command or convention changes, update this file in the same
pull request. It points to the documents and files that hold a procedure
rather than repeating them; follow those instead of re-deriving a command.

## Quick commands

Run from the project root. Pixi is the canonical workflow and matches CI.

- **Configure the build directory:** `pixi run setup`
- **Build** (installs to `bin/`): `pixi run build`
- **Run the tests** (smoke, unit, and run-twice determinism): `pixi run test`;
  it does not rebuild, so run `pixi run build` first
- **Run one group of tests:** `pixi run test-unit` (also `test-initial`,
  `test-write`, `test-shapefile`, `test-intersect`, `test-vtk`, and
  `test-determinism`; see `pixi.toml`)
- **Run one unit test:** `builddir/test/unit/unit_tests --gtest_filter='QuadTree.*'`
- **Compare outputs with main** (needs a build of main at `../gridgen-main`,
  e.g. `git worktree add ../gridgen-main main`, then `pixi run setup-main`
  and `pixi run build-main`): `pixi run compare-main`
- **Check Python format, lint, and spelling:** `pixi run check-format`;
  `pixi run fix-format` applies the fixes
- **Show the OS tag used for release archives:** `pixi run get-ostag`
- **Build the release notes** (needs TeX Live and the USGS class from
  `MODFLOW-ORG/usgslatex`): `pixi run build-release-notes`, which first runs
  `doc/ReleaseNotes/mk_releasenotes.py` to check the items and write the
  History entry for the release in progress (`develop.tex`)
- **Set the version** (`version.txt` is the source; a version containing
  `dev` is a development version): `pixi run update-version --version <v>`;
  add `--date YYYY-MM-DD --releasemode` for a release

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
- The release notes, `doc/ReleaseNotes/ReleaseNotes.tex`, list corrections to
  the report (e.g., the `ROTATION_ANGLE` definition) and known issues. Each
  past release is in `doc/ReleaseNotes/previous/`; the version, date, and
  develop mode are in `doc/ReleaseNotes/version.tex`.
- **References:** cite with natbib (`\citet`, `\citep`) from
  `doc/ReleaseNotes/ReleaseNotes.bib`. Write entries in the style of the
  MODFLOW 6 bibliography (`doc/MODFLOW6References.bib` in
  MODFLOW-ORG/modflow6), and copy an entry from there when it has one: a USGS
  report is a `@book` with `author`, `series` (e.g., `{{U.S. Geological Survey
  Open-File Report 2014--1109, 26 p.}}`), `title`, `url`, `urldate`, and
  `year`, keyed like `lien2015gridgen`.

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
  review has no actionable findings, at most three rounds for code changes; a
  pull request that changes only documentation (Markdown, `doc/`, comments)
  gets one round, and its findings are addressed without another review.
  Findings still in dispute are left for the maintainer to decide.
- **Branching:** `origin` is a fork and `upstream` is `MODFLOW-ORG/gridgen`.
  Cut each branch from an up-to-date `upstream/main`, and read `upstream`
  refs, not `origin` refs, for the state of the repository.
- **Commit messages:** conventional commits (`fix:`, `feat:`, `refactor:`,
  `test:`, `docs:`, `ci:`, `chore:`), kept concise: an imperative, lowercase
  subject and a body that contrasts the previous and new behavior.
- **Testing scope:** rebuild, then run the unit tests or test group that
  covers the change before the full `pixi run test`. Run
  `pixi run compare-main` when a change can alter output.
- **Error output in unit tests:** do not assert on captured stderr.
  `testing::internal::CaptureStderr()` returns nothing on the
  `windows-2022 intel 2025.0` CI build, though it works on the others, so
  such a test passes locally and fails only in CI. Have the code under test
  record what it reports, for example with an accessor like
  `DefinitionParser::getWarnings()` (see `UnrecognizedRecordIsReported`),
  and check a death test by its exit code alone, with an empty message
  matcher: `EXPECT_EXIT(..., testing::ExitedWithCode(1), "")` (see
  `test_asciigrid.cpp`).
- **Formatting:** run `pixi run check-format` before committing Python;
  ruff (`ruff.toml`) formats and lints it, and codespell (`.codespellrc`)
  checks spelling in the text files other than the C++ source, data files,
  generated files, and the archived release notes. No formatter runs on the
  C++ source yet (#64); match the layout and naming of the surrounding code.
- **Comments and documentation:** keep code comments and docstrings short,
  stating what a function does; put the reason, when one is needed, in an
  inline comment. Write the release notes for groundwater modelers rather
  than programmers, in U.S. English, with one line per paragraph in the
  LaTeX source.
- **Output comparison:** on pull requests, CI requires every output file to
  match a build of the base branch byte for byte. A pull request that changes
  output on purpose gets the `output-change` label and says what changed and
  why.
- **Known issues in the release notes:** the Known Issues section of
  `doc/ReleaseNotes/ReleaseNotes.tex` lists the open defects a user can run
  into, each in user terms with any workaround and `\ghissue{N}`. Add an
  entry when such a bug issue is filed, and remove it in the pull request that
  fixes the bug, which also adds a release note item.
- **Known defects in tests:** a known defect is recorded as an entry in
  `KNOWN_UNSTABLE` (`test/compare_outputs.py`) or a `DISABLED_` unit test, with
  a comment naming its issue. The pull request that fixes the issue removes the
  entry or enables the test, and adds the tests listed in the issue.
- **Generated files:** do not commit test outputs, `builddir/`, or
  `subprojects/` downloads. Do not hand-edit a generated file; edit its
  source and regenerate: `develop.tex` comes from the release note items
  (`mk_releasenotes.py`), `gridgen_version.h` from `version.txt` (the
  build), and the versions in `pixi.toml` and `version.tex` from
  `pixi run update-version`.
- **CI:** changes only to Markdown, `doc/`, issue templates, `dependabot.yml`,
  `LICENSE`, or the lint configuration skip the build and tests
  (`paths-ignore` in `.github/workflows/ci.yml`); update that list when adding
  non-code files. Adding or removing a label on a pull request starts a new
  run that cancels the one in progress (`cancel-in-progress`), including the
  labels that `gh pr create --label` applies just after opening; judge a
  pull request by its latest run (`gh run list -b <branch>`), not by the
  cancelled jobs `gh pr checks` reports. `lint.yml` runs
  `pixi run check-format` on every push and on pull requests to `main`.
  `release-notes.yml` builds the release notes on Ubuntu when
  `doc/ReleaseNotes/` changes and attaches the PDF to each GitHub release.
- **Tone:** reserve compliments for genuinely exceptional finds or ideas.
