## A run started from BExTE-app happens in a background process, so the
## filesystem is the only channel back to the app's Run tab. Result rows are too
## coarse a signal on their own: a method that runs in parallel collects its
## scenarios in the master and writes them all at once when the method is
## over, so the readout would jump a whole method at a time. The tracker
## below is fed by the same per-scenario callback that drives the console
## progress bar - foreach's .options.snow$progress, which reports the number
## of tasks finished so far within one case study/method block, and the
## sequential loop's own index.

## run_progress_path() is relative to the working directory the simulation
## runs in, the same convention as ./results/<env>/ and ./logs/<env>/.
local_run_root <- function(env = parent.frame()) {
  root <- tempfile("bexte-run")
  dir.create(root, recursive = TRUE)
  old <- setwd(root)
  withr::defer(setwd(old), envir = env)
  root
}

test_that("a new tracker records no progress against the run's total", {
  local_run_root()

  run_progress_tracker("demo", total = 108)

  snapshot <- jsonlite::read_json(run_progress_path("demo"), simplifyVector = TRUE)
  expect_identical(snapshot$done, 0L)
  expect_identical(snapshot$total, 108L)
})

test_that("tick records the scenarios finished so far", {
  local_run_root()

  tracker <- run_progress_tracker("demo", total = 108)
  tracker$starting("botox", "separate")
  tracker$tick(1)
  tracker$tick(2)

  snapshot <- jsonlite::read_json(run_progress_path("demo"), simplifyVector = TRUE)
  expect_identical(snapshot$done, 2L)
})

test_that("tick counts on from the previous case study and method", {
  local_run_root()

  tracker <- run_progress_tracker("demo", total = 108)
  tracker$starting("botox", "separate")
  tracker$tick(54)
  tracker$starting("botox", "pooling")
  tracker$tick(1)

  snapshot <- jsonlite::read_json(run_progress_path("demo"), simplifyVector = TRUE)
  expect_identical(snapshot$done, 55L)
})

test_that("the snapshot names the case study and method being simulated", {
  local_run_root()

  tracker <- run_progress_tracker("demo", total = 108)
  tracker$starting("botox", "pooling")

  snapshot <- jsonlite::read_json(run_progress_path("demo"), simplifyVector = TRUE)
  expect_identical(snapshot$case_study, "botox")
  expect_identical(snapshot$method, "pooling")
})

test_that("writing a snapshot leaves no partial file behind for the app to read", {
  local_run_root()

  tracker <- run_progress_tracker("demo", total = 108)
  tracker$starting("botox", "separate")
  tracker$tick(1)

  written <- list.files(dirname(run_progress_path("demo")))
  expect_identical(written, "progress.json")
})

test_that("a tracker starts a run's progress over rather than resuming the last one", {
  local_run_root()

  tracker <- run_progress_tracker("demo", total = 108)
  tracker$starting("botox", "separate")
  tracker$tick(54)

  run_progress_tracker("demo", total = 108)

  snapshot <- jsonlite::read_json(run_progress_path("demo"), simplifyVector = TRUE)
  expect_identical(snapshot$done, 0L)
})
