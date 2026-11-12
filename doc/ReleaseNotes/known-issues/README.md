# Known issues

Each open defect a user can run into is one TOML file here, named after the
defect rather than the issue number (for example
`active-domain-edge-ghost-nodes.toml`). Add the file when the bug issue is
filed, and delete it in the pull request that fixes the bug, which adds a
release note item in `../items/` instead.

```toml
issue = 170
description = "Along the edge of an active domain, some ghost nodes in the .gnc.dat file do not reproduce a head that varies linearly, ..."
```

- `issue` is the number of the GitHub issue that tracks the defect.
- `description` says, for someone using GRIDGEN, what goes wrong, when, and
  any workaround; the issue number is added after it.

`../mk_releasenotes.py` checks the files and renders them as the Known Issues
section (`knownissues.tex`); `pixi run build-release-notes` runs it. They are
kept at a release, since the defects remain open.
