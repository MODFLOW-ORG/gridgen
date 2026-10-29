# Release note items

Each pull request that changes gridgen adds one TOML file here describing the
change for the release notes. Name the file after the change, not the pull
request number (for example `active-domain-layer.toml`).

```toml
section = "fixes"
subsection = "input"
description = "ACTIVE_DOMAIN LAYER n was ignored unless every layer above it also had one. Each layer's active domain is now read independently."
```

- `section` and `subsection` must be keys in `../schema.toml`; `subsection` may
  be empty.
- `description` is one to four sentences for someone using gridgen: what was
  wrong or what is new, and what it means for an existing model. Detail
  belongs in the pull request.

`../mk_releasenotes.py` checks the items and renders them as the History entry
for the release in progress (`develop.tex`); `pixi run build-release-notes` runs
it. At a release the entry is kept in `../previous/` and the items are removed.
