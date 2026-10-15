# Use the renv library of the BExTE checkout next to this folder, so that
# BExTEsim builds and tests against the package versions pinned in
# ../BExTE/renv.lock. Install BExTE into that library first.
local({
  bexte <- normalizePath(file.path("..", "BExTE"), mustWork = FALSE)
  if (file.exists(file.path(bexte, "renv", "activate.R"))) {
    Sys.setenv(RENV_PROJECT = bexte)
    source(file.path(bexte, "renv", "activate.R"))
  }
})
