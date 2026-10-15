#' Determine a safe number of parallel workers
#'
#' Every worker holds its own copy of a scenario's replicates, so memory grows
#' with the worker count. The `BEXTE_MAX_WORKERS` environment variable caps it
#' for a whole run, for instance to leave room for other work on the machine.
#'
#' @param max_workers Maximum number of workers to use.
#' @param detected_cores Number of detected logical CPU cores.
#' @param env_max_workers Value of `BEXTE_MAX_WORKERS`; empty means no cap.
#'
#' @return A positive integer worker count.
#' @noRd
get_parallel_worker_count <- function(
  max_workers = Inf,
  detected_cores = parallel::detectCores(),
  env_max_workers = Sys.getenv("BEXTE_MAX_WORKERS")
) {
  if (nzchar(env_max_workers)) {
    env_cap <- suppressWarnings(as.integer(env_max_workers))
    if (is.na(env_cap) || env_cap < 1L) {
      stop("BEXTE_MAX_WORKERS must be a positive integer, but it is \"",
           env_max_workers, "\".", call. = FALSE)
    }
    max_workers <- min(max_workers, env_cap)
  }

  if (length(detected_cores) != 1L ||
      is.na(detected_cores) ||
      detected_cores < 2L) {
    return(1L)
  }

  return(as.integer(min(detected_cores - 1L, max_workers)))
}

#' Cap chain-level parallelism when scenarios already run in parallel
#'
#' When the scenario loop is parallelised, every worker would otherwise start
#' `parallel_chains` CmdStan processes of its own, oversubscribing the machine by
#' the product of the two. Running the chains sequentially inside each worker
#' keeps the total number of processes equal to the number of workers.
#'
#' @param mcmc_config MCMC configuration read from `mcmc_config.yml`.
#' @param parallelization Whether scenarios are simulated in parallel.
#'
#' @return The MCMC configuration, with `parallel_chains` capped at 1 when
#'   scenarios run in parallel.
#' @noRd
limit_mcmc_chain_parallelism <- function(mcmc_config, parallelization) {
  if (isTRUE(parallelization)) {
    mcmc_config$parallel_chains <- 1L
  }

  return(mcmc_config)
}

#' Check and normalise the `parallelization` setting
#'
#' @param parallelization The `parallelization` entry of the scenarios config.
#'
#' @return `FALSE` for `NULL`, a single logical as it is, or a character vector
#'   of method names.
#' @noRd
normalise_parallelization <- function(parallelization) {
  if (is.null(parallelization)) {
    return(FALSE)
  }

  if (is.logical(parallelization)) {
    if (length(parallelization) != 1L || is.na(parallelization)) {
      stop(
        "parallelization must be a single TRUE or FALSE, or a list of method names.",
        call. = FALSE
      )
    }
    return(parallelization)
  }

  methods <- unlist(parallelization, use.names = FALSE)
  if (!is.character(methods) || length(methods) == 0L) {
    stop(
      "parallelization must be a single TRUE or FALSE, or a list of method names.",
      call. = FALSE
    )
  }

  return(methods)
}

#' Whether one method's scenarios run in parallel
#'
#' `parallelization` in the scenarios configuration is either a single logical,
#' which applies to every method in the run, or a list of method names, which
#' parallelises only those. The per-method form exists because the cost profiles
#' differ: a method left on the replicate loop gains the full worker count from
#' scenario-level parallelism, while a vectorised one already holds a large
#' posterior mixture in memory and multiplies that footprint by every worker.
#'
#' @param parallelization The `parallelization` entry of the scenarios config.
#' @param method Method name.
#'
#' @return A single logical.
#' @noRd
method_runs_in_parallel <- function(parallelization, method) {
  parallelization <- normalise_parallelization(parallelization)
  if (is.logical(parallelization)) {
    return(parallelization)
  }

  return(method %in% parallelization)
}

#' Whether the post-processing analysis runs in parallel
#'
#' `parallelization` is either a single logical or a list of method names (see
#' [method_runs_in_parallel()]). The analysis loops over every result row at
#' once rather than one method at a time, so the per-method form collapses to
#' "run in parallel if any method does".
#'
#' @param parallelization The `parallelization` entry of the scenarios config.
#'
#' @return A single logical.
#' @noRd
analysis_runs_in_parallel <- function(parallelization) {
  return(!isFALSE(normalise_parallelization(parallelization)))
}

#' Result rows below which the analysis is not worth a cluster
#'
#' Standing the cluster up costs around 12 seconds on a source checkout,
#' because every worker loads the package with pkgload::load_all(), while a
#' row costs about 0.07 seconds to analyse. Parallelism only pays for itself
#' past roughly `startup / (row_cost * (1 - 1 / workers))` rows - about 185
#' on a 12-core machine - and below that the cluster is pure overhead: a
#' 108-row environment takes 7.6s sequentially against 27.5s on 11 workers.
#' @noRd
ANALYSIS_PARALLEL_MIN_ROWS <- 200L

#' Fewest simulated trials worth generating on a cluster
#'
#' The separate analysis's trials are simulated once per design, and their cost
#' scales with the replicate count: about 0.25 ms a trial on a time-to-event
#' endpoint, which fits a model to each, so 2.5 s a design at the paper's
#' 10,000 replicates. The cluster costs 15-30 s to stand up, so the threshold is
#' on designs times replicates rather than on designs: 24 designs of 2,000
#' trials took 6.7 s sequentially and 28 s on a cluster, while the paper's 264
#' designs of 10,000 are 2.6 million trials.
#' @noRd
ANALYSIS_PARALLEL_MIN_TRIALS <- 200000L

#' Whether to analyse this many rows on a cluster
#'
#' @param parallelization Whether the caller asked for parallelism.
#' @param n_rows Number of result rows to analyse.
#' @param min_rows Row count below which the cluster costs more than it saves.
#'
#' @return A single logical.
#' @noRd
analysis_uses_cluster <- function(parallelization, n_rows,
                                  min_rows = ANALYSIS_PARALLEL_MIN_ROWS) {
  return(isTRUE(parallelization) && n_rows >= min_rows)
}

#' Make BExTEsim (and any other packages) available in every worker of a cluster
#'
#' BExTEsim is not necessarily an installed package: BExTE-app's background
#' process `pkgload::load_all()`s it from source, and a worker that only
#' calls `library(BExTEsim)` fails outright there.
#' Each worker therefore falls back to loading the same source tree.
#'
#' @param cl A cluster from `parallel::makeCluster()`.
#' @param packages Additional package names to attach in each worker.
#'
#' @return `NULL`, invisibly.
#' @noRd
load_bexte_in_workers <- function(cl, packages = character()) {
  paths <- .libPaths()
  pkg_root <- find.package("BExTEsim")

  parallel::clusterExport(
    cl,
    varlist = c("paths", "pkg_root", "packages"),
    envir = environment()
  )

  parallel::clusterEvalQ(cl, {
    .libPaths(paths)
    if (requireNamespace("BExTEsim", quietly = TRUE)) {
      library(BExTEsim)
    } else {
      pkgload::load_all(pkg_root, quiet = TRUE)
    }
    for (package in packages) {
      library(package, character.only = TRUE)
    }
    NULL
  })

  return(invisible(NULL))
}

#' Register the parallel backend simulations run their scenarios on
#'
#' doSNOW rather than doParallel: the simulation loops pass a per-scenario
#' progress callback as `.options.snow$progress`, and only doSNOW calls it as
#' each result comes back. doParallel reads `preschedule` and
#' `attachExportEnv` out of those options and warns that it is ignoring the
#' rest, which leaves the console progress bar and the progress file
#' BExTE-app reads frozen until a whole case study/method block is finished.
#'
#' @param cl A cluster from `parallel::makeCluster()`.
#'
#' @return The backend registration's return value, invisibly.
#' @noRd
register_parallel_backend <- function(cl) {
  return(invisible(doSNOW::registerDoSNOW(cl)))
}

#' Check that the packages for parallel runs are installed
#'
#' foreach and doSNOW are in Suggests: only parallel runs need them. Call this
#' before starting a cluster, so that a missing package stops the run before
#' any worker starts.
#'
#' @return `NULL`, invisibly. Raises an error if a package is missing and the
#'   user does not install it.
#' @noRd
check_parallel_packages <- function() {
  return(rlang::check_installed(c("foreach", "doSNOW"), reason = "to run in parallel."))
}

#' Packages every analysis worker attaches
#'
#' The union of what the analysis steps' workers need, so that one cluster
#' can serve all of them.
#' @noRd
ANALYSIS_WORKER_PACKAGES <- c("dplyr", "yaml", "pwr", "BSDA")

#' A worker cluster the analysis steps share, started on first use
#'
#' Each analysis step used to stand up its own cluster and stop it on the
#' way out, at 15-30 seconds a time on a source checkout - most of the
#' parallel run time of a small environment. [simulation_analysis()] creates
#' one of these instead and hands it to every step; the cluster starts the
#' first time a step actually wants one, so a run whose steps all stay
#' sequential never pays for it.
#'
#' @return A list of two functions: `get()`, which returns the cluster,
#'   starting it if need be, with the foreach backend registered on it; and
#'   `stop()`, which stops it if it was started.
#' @noRd
new_analysis_cluster <- function() {
  state <- new.env(parent = emptyenv())

  return(list(
    get = function() {
      if (is.null(state$cl)) {
        check_parallel_packages()
        state$cl <- parallel::makeCluster(get_parallel_worker_count())
        load_bexte_in_workers(state$cl, packages = ANALYSIS_WORKER_PACKAGES)
      }
      # Registered on every use: something in between may have registered
      # another backend.
      register_parallel_backend(state$cl)
      return(state$cl)
    },
    stop = function() {
      if (!is.null(state$cl)) {
        parallel::stopCluster(state$cl)
        state$cl <- NULL
      }
      return(invisible(NULL))
    }
  ))
}
