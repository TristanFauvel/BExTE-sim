#' The analysis steps that follow a simulation
#'
#' An environment's `scenarios_config.yml` selects a subset of these by name
#' through its `analysis_steps` key. They are a constant so that a name which
#' does not exist is rejected where the config is read, rather than quietly
#' running nothing.
#' @noRd
ANALYSIS_STEPS <- c(
  "frequentist_power_at_equivalent_tie",
  "frequentist_power_at_nominal_tie",
  "sweet_spot",
  "bayesian_ocs"
)



#' Resolve which analysis steps to run for an environment
#'
#' @description A `scenarios_config.yml` without an `analysis_steps` key runs
#'   all of [ANALYSIS_STEPS], so nothing written before that key existed
#'   changes behaviour. Naming a subset skips the rest, which is worth doing
#'   because no step is cheap and not every run needs every output - the
#'   paper's figures, for one, read neither the sweet spot nor the Bayesian
#'   operating characteristics.
#'
#' @param scenarios_config The environment's scenarios configuration.
#' @param all_case_studies_shipped Whether every case study of the run ships
#'   with the package. Both power steps call `load_data()` without threading
#'   `case_studies_config_dir` through, so they only resolve shipped case
#'   studies and are dropped rather than left to crash.
#'
#' @return A character vector of step names, which may be empty.
#' @noRd
resolve_analysis_steps <- function(scenarios_config, all_case_studies_shipped) {
  steps <- if (is.null(scenarios_config$analysis_steps)) {
    ANALYSIS_STEPS
  } else {
    as.character(unlist(scenarios_config$analysis_steps, use.names = FALSE))
  }

  # A misspelled step would otherwise be skipped in silence, and the omission
  # only surfaces much later, as a figure with nothing to plot.
  unknown_steps <- setdiff(steps, ANALYSIS_STEPS)
  if (length(unknown_steps) > 0) {
    stop(
      "analysis_steps names steps that do not exist: ",
      paste(unknown_steps, collapse = ", "),
      ". The steps are: ", paste(ANALYSIS_STEPS, collapse = ", "),
      call. = FALSE
    )
  }

  if (!all_case_studies_shipped) {
    steps <- setdiff(steps, c("frequentist_power_at_equivalent_tie",
                              "frequentist_power_at_nominal_tie"))
  }

  return(steps)
}



#' Run a full simulation for one environment
#'
#' @description Runs the frequentist and/or Bayesian Monte Carlo operating
#'   characteristics simulation for a single environment. This mirrors the
#'   per-environment body of the driver loop in BExTEpaper's `scripts/main.R`, extracted
#'   into a callable function so it can be invoked directly (e.g. from a
#'   background process launched by BExTE-app).
#'
#' @param env The environment name (used to name the `./logs/<env>/` and
#'   `./results/<env>/` directories).
#' @param config_dir Directory containing `scenarios_config.yml`,
#'   `mcmc_config.yml` and `methods_config.R` for this environment (must end
#'   with a trailing slash).
#' @param case_studies_config_dir Directory containing the case study YAML
#'   files referenced by this environment (must end with a trailing slash).
#' @param simulation_config Simulation configuration list (as read from
#'   `simulation_config.yml`).
#' @param analysis_config Analysis configuration list (as read from
#'   `analysis_config.yml`).
#' @param frequentist_metrics List of frequentist metrics, as defined by
#'   sourcing `metrics_config.R`.
#' @param inference_metrics List of inference metrics, as defined by sourcing
#'   `metrics_config.R`.
#' @param results_dir Directory the results are written to. Defaults to
#'   `./results/<env>/`.
#' @param check_results_completeness Reserved for a future completeness
#'   check; currently unused, mirroring BExTEpaper's `scripts/main.R` where the
#'   corresponding call is commented out.
#'
#' @return `TRUE`, invisibly, once frequentist and/or Bayesian OCs have been
#'   computed and concatenated according to `simulation_config`.
#' @export
run_simulation_env <- function(env,
                               config_dir,
                               case_studies_config_dir,
                               simulation_config,
                               analysis_config,
                               frequentist_metrics,
                               inference_metrics,
                               results_dir = paste0("./results/", env, "/"),
                               check_results_completeness = TRUE) {
  check_decision_threshold(simulation_config, "simulation_config")

  # Set up the log file location
  dir.create(paste0("./logs/", env, "/checkpoints/"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0("./logs/", env, "/error_logs/"), showWarnings = FALSE, recursive = TRUE)
  # Empty the directory
  unlink(paste0("./logs/", env, "/checkpoints/*"), recursive = TRUE)
  unlink(paste0("./logs/", env, "/error_logs/*"), recursive = TRUE)

  LOGGING_FILE_PATH <- generate_log_filename(base_name = paste0("./logs/", env, "/error_logs/error_log.log"), suffix_type = "timestamp")

  scenarios_config <- read_config(paste0(config_dir, "scenarios_config.yml"), scenarios_config_schema)

  # Set up the new log file. Sequential and parallel runs use the same
  # logger: the workers' futile.logger calls go to their own console.
  futile.logger::flog.appender(futile.logger::appender.file(LOGGING_FILE_PATH))

  # Set the global error handler for this run only
  withr::local_options(error = global_error_handler)

  outputs_config <- yaml::read_yaml(system.file("conf/outputs_config.yml", package = "BExTEsim"))
  ocs_filename <- outputs_config$frequentist_ocs_results_filename

  if (simulation_config$compute_frequentist_ocs == TRUE) {
    if (simulation_config$delete_old_results) {
      unlink(results_dir, recursive = TRUE, force = TRUE)
      dir.create(results_dir, recursive = TRUE)  # This recreates the empty directory
    }

    simulation_frequentist_ocs(env = env,
    simulation_config = simulation_config,
    scenarios_config = scenarios_config,
    analysis_config = analysis_config,
    config_dir = config_dir,
    case_studies_config_dir = case_studies_config_dir,
    logging_file_path = LOGGING_FILE_PATH,
    frequentist_metrics = frequentist_metrics,
    inference_metrics = inference_metrics)

    # Concatenate case_study/method results into a single file for the environment.
    concatenate_simulation_results(results_dir = results_dir, ocs_filename = ocs_filename)

    # frequentist_power_at_equivalent_tie()/frequentist_power_at_nominal_tie()
    # (in R/analysis_operating_characteristics.R) call load_data() without
    # threading case_studies_config_dir through, so they only resolve
    # package-shipped case studies. Skip them for a case study outside the
    # package rather than crash; sweet_spot and bayesian_ocs (which do
    # thread case_studies_config_dir through) still run either way.
    all_case_studies_shipped <- all(vapply(scenarios_config$case_studies, function(cs) {
      return(nzchar(system.file(file.path("conf", "case_studies", paste0(cs, ".yml")), package = "BExTE")))
    }, logical(1)))
    analysis_to_compute <- resolve_analysis_steps(scenarios_config,
                                                  all_case_studies_shipped)

    if (length(analysis_to_compute) > 0) {
      # The analysis is performed on the results concatenated at the level of the environment.
      simulation_analysis(env, analysis_config, config_dir, frequentist_metrics,
                          to_compute = analysis_to_compute,
                          case_studies_config_dir = case_studies_config_dir,
                          parallelization = scenarios_config$parallelization)
    }
  }

  if (simulation_config$compute_bayesian_ocs_mc == TRUE) {
    if (simulation_config$delete_old_results) {
      folder_dir <- file.path(results_dir, "bayesian")
      if (file.exists(folder_dir)) {
        unlink(folder_dir, recursive = TRUE)
      }
    }

    simulation_bayesian_ocs(
      env = env,
      simulation_config = simulation_config,
      scenarios_config = scenarios_config,
      analysis_config = analysis_config,
      config_dir = config_dir,
      case_studies_config_dir = case_studies_config_dir
    )
    ocs_filename <- outputs_config$bayesian_ocs_mc_results_filename
    concatenate_simulation_results(results_dir = results_dir, ocs_filename = ocs_filename)
  }

  return(invisible(TRUE))
}
