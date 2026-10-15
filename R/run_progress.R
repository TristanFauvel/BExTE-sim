## Progress reporting for a simulation run.
##
## A run launched from BExTE-app happens in a background process, so the
## filesystem is the only channel back to the app's Run tab. Result rows are too
## coarse a signal on their own: a method that runs in parallel collects its
## scenarios in the master and writes them in one go once the method is over,
## which makes a row-counting readout jump a whole method at a time. The
## tracker here is fed by the same per-scenario callback that drives the
## console progress bar, so the app can follow a run scenario by scenario.

#' Path of the file a run reports its progress to.
#'
#' Relative to the working directory the simulation runs in, the same
#' convention as ./results/<env>/ and ./logs/<env>/.
#'
#' @param env Name of the environment being run.
#'
#' @return Path of the environment's progress file.
#'
#' @noRd
run_progress_path <- function(env) {
  return(file.path("logs", env, "progress.json"))
}

#' Write a snapshot of how far a run has got.
#'
#' The snapshot goes to a temporary file in the same directory and is renamed
#' into place. The app polls this file about once a second while the run
#' writes it, and a rename is atomic within a filesystem, so the app never
#' reads a half-written snapshot.
#'
#' @param path Path of the progress file, see [run_progress_path()].
#' @param done Number of scenarios simulated so far.
#' @param total Number of scenarios the run covers.
#' @param case_study,method What is being simulated right now.
#'
#' @return The path, invisibly.
#'
#' @noRd
write_run_progress <- function(path,
                               done,
                               total,
                               case_study = NA_character_,
                               method = NA_character_) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  snapshot <- list(
    done = as.integer(done),
    total = as.integer(total),
    case_study = as.character(case_study),
    method = as.character(method),
    updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  )

  tmp <- tempfile(pattern = "progress", tmpdir = dirname(path), fileext = ".json")
  jsonlite::write_json(snapshot, tmp, auto_unbox = TRUE, null = "null", na = "null")
  if (!file.rename(tmp, path)) {
    unlink(tmp)
  }
  return(invisible(path))
}

#' Track how many scenarios of a run have been simulated.
#'
#' `total` counts the scenarios of the whole run - every case study and
#' method it covers - whereas the callbacks that drive the console progress
#' bar count within the case study/method block being worked on: `foreach`
#' calls `.options.snow$progress` with the number of tasks finished so far,
#' and the sequential loop passes its own index. `starting()` rebases on a
#' block boundary so `tick()` can take those block-local counts and still
#' record a run-wide total.
#'
#' Creating a tracker resets the environment's progress file, so a relaunch
#' starts from zero rather than resuming whatever the previous run left.
#'
#' @param env Name of the environment being run.
#' @param total Number of scenarios the run covers.
#'
#' @return A list with `path`, `starting(case_study, method)` and `tick(n)`.
#' @export
run_progress_tracker <- function(env, total) {
  path <- run_progress_path(env)

  state <- new.env(parent = emptyenv())
  state$done <- 0L
  state$base <- 0L
  state$case_study <- NA_character_
  state$method <- NA_character_

  flush <- function() {
    return(write_run_progress(path, state$done, total, state$case_study, state$method))
  }
  flush()

  return(list(
    path = path,
    starting = function(case_study, method) {
      state$base <- state$done
      state$case_study <- case_study
      state$method <- method
      return(flush())
    },
    tick = function(n = 1L) {
      state$done <- state$base + as.integer(n)
      return(flush())
    }
  ))
}
