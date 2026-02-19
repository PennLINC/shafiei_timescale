# R Environment Setup

This project uses `renv` for reproducible R package management.

## One-time setup

From repo root:

```bash
cd /Volumes/cbica/projects/developmental_gradients/gitrepo/shafiei_timescale/code
Rscript r_env/setup_renv.R
```

This script will:
- use a writable local library at `.r-user-lib/`,
- install `renv` if missing,
- initialize `renv` for this repo,
- hydrate required packages from your existing system R library (no source compile if already installed),
- install only core missing packages (required for analysis),
- prefer binary package installs to avoid local compiler/toolchain issues,
- for `ggseg`, `ggsegSchaefer`, and `ggseg3d`, try CRAN first and then GitHub fallback,
- treat `ggseg*` packages as optional (setup continues even if those fail),
- include `svglite` so `ggsave(..., ".svg")` works out of the box,
- install missing packages with `dependencies = FALSE` to avoid unnecessary upgrades / source rebuilds,
- write `renv.lock`.

## Recreate environment later

```bash
cd /Volumes/cbica/projects/developmental_gradients/gitrepo/shafiei_timescale/code
R
```

Then in R:

```r
renv::restore()
```

## Package versions source

Pinned versions come from:

`/Volumes/cbica/projects/developmental_gradients/gitrepo/shafiei_timescale/code/r_env/timescale_r_package_versions.csv`
