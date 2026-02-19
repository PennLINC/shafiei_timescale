#!/usr/bin/env Rscript

# Bootstrap a project-local renv environment for this repository.
# Run from repo root:
#   Rscript r_env/setup_renv.R

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
local_lib <- file.path(repo_root, ".r-user-lib")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"))
# Prefer prebuilt binaries on macOS to avoid local C/C++ compilation issues.
options(pkgType = "binary")

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", lib = local_lib)
}

required_core <- c(
  "stringr",
  "tidyr",
  "mgcv",
  "gratia",
  "tidyverse",
  "dplyr",
  "paletteer",
  "pals",
  "ggplot2",
  "svglite",
  "scales",
  "tibble",
  "purrr"
)

optional_plot <- c("ggseg", "ggsegSchaefer", "ggseg3d")

# Discover likely package libraries (system + user) and search all of them.
split_paths <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  strsplit(x, .Platform$path.sep, fixed = TRUE)[[1]]
}

ver <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
candidate_libs <- unique(c(
  .libPaths(),
  split_paths(Sys.getenv("R_LIBS_USER")),
  split_paths(Sys.getenv("R_LIBS_SITE")),
  .Library,
  .Library.site,
  path.expand(file.path("~", "Library", "R", ver, "library")),
  path.expand(file.path("~", "Library", "R", "arm64", ver, "library")),
  path.expand(file.path("~", "Library", "R", "x86_64", ver, "library"))
))
candidate_libs <- candidate_libs[nzchar(candidate_libs) & dir.exists(candidate_libs)]

# Ensure hydrate can see all discovered libraries.
.libPaths(unique(c(local_lib, candidate_libs, .libPaths())))

# Capture installed packages BEFORE renv project isolation.
installed_preinit <- unique(rownames(installed.packages(lib.loc = candidate_libs)))

if (!dir.exists(file.path(repo_root, "renv"))) {
  renv::init(bare = TRUE)
}

present_core <- intersect(required_core, installed_preinit)
missing_core <- setdiff(required_core, present_core)
present_optional <- intersect(optional_plot, installed_preinit)
missing_optional <- setdiff(optional_plot, present_optional)

install_with_fallback <- function(pkg, github_remote = NULL) {
  ok <- tryCatch({
    renv::install(pkg, dependencies = FALSE, type = "binary")
    TRUE
  }, error = function(e) FALSE)

  # If binary is unavailable for a package, try source once.
  if (!ok) {
    ok <- tryCatch({
      renv::install(pkg, dependencies = FALSE, type = "source")
      TRUE
    }, error = function(e) FALSE)
  }

  if (!ok && !is.null(github_remote)) {
    message(sprintf("CRAN install failed for %s; trying GitHub %s", pkg, github_remote))
    renv::install(github_remote, dependencies = FALSE)
    ok <- TRUE
  }

  ok
}

message("Hydrating any required packages that already exist in system library.")
if (length(present_core) > 0) {
  renv::hydrate(packages = present_core)
}
if (length(present_optional) > 0) {
  renv::hydrate(packages = present_optional)
}

if (length(missing_core) > 0) {
  message("Some core packages are missing; installing missing core packages.")

  for (pkg in missing_core) {
    if (!install_with_fallback(pkg, NULL)) {
      stop(sprintf("Failed to install required core package '%s'", pkg))
    }
  }
}

optional_failed <- character()
if (length(missing_optional) > 0) {
  message("Some plotting packages are missing; attempting optional installs.")

  # Some ggseg ecosystem packages are distributed via GitHub rather than CRAN.
  # Try CRAN first; if unavailable, fall back to canonical GitHub remotes.
  github_map <- list(
    ggseg = "ggseg/ggseg",
    ggsegSchaefer = "ggseg/ggsegSchaefer",
    ggseg3d = "ggseg/ggseg3d"
  )

  for (pkg in missing_optional) {
    remote <- github_map[[pkg]]
    if (!install_with_fallback(pkg, remote)) {
      optional_failed <- c(optional_failed, pkg)
    }
  }
}

renv::snapshot(prompt = FALSE)

if (length(optional_failed) > 0) {
  message(
    sprintf(
      "renv setup complete, but optional plotting packages failed: %s",
      paste(optional_failed, collapse = ", ")
    )
  )
  message("Core analysis packages are installed; ggseg plotting in scp_timescale.R may be unavailable.")
} else {
  message("renv setup complete. Created/updated renv.lock")
}
