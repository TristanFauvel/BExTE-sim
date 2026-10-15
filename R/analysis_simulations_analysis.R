#' Perform Simulation Analysis
#'
#' @description This function performs a simulation analysis for a given environment using the specified
#' analysis configuration and configuration directory. It reads the output configuration,
#' loads results, applies frequentist power analysis, and computes the sweet spot.
#'
#' @param env A character string specifying the environment for the analysis (e.g., "development", "production").
#' @param analysis_config A list containing the analysis configuration parameters.
#' @param config_dir A character string specifying the directory containing the configuration files.
#' @param frequentist_metrics List of frequentist metrics.
#' @param case_studies Case studies to analyse, or "all".
#' @param to_compute Which analyses to run: any of
#'   "frequentist_power_at_equivalent_tie", "frequentist_power_at_nominal_tie",
#'   "sweet_spot" and "bayesian_ocs".
#' @param methods Methods to analyse, or "all".
#' @param case_studies_config_dir Directory holding the case study YAMLs.
#' @param parallelization Whether the analysis runs in parallel: a single
#'   logical, or a list of method names, in the same form as the
#'   `parallelization` entry of `scenarios_config.yml`. Defaults to that
#'   entry, read from `config_dir`.
#' @param n_replicates Number of Monte Carlo replicates the simulated power
#'   estimates are built from. Defaults to the `n_replicates` entry of
#'   `scenarios_config.yml`, read from `config_dir`.
#'
#' @return This function does not return a value. It writes the results and sweet spot analysis
#' to CSV files in the specified results directory.
#'
#' @importFrom yaml read_yaml
#' @importFrom readr read_csv write_csv
#'
#' @export
simulation_analysis <- function(env,
                                analysis_config,
                                config_dir,
                                frequentist_metrics, case_studies = "all", to_compute = c("frequentist_power_at_equivalent_tie", "frequentist_power_at_nominal_tie", "sweet_spot", "bayesian_ocs"), methods = "all", case_studies_config_dir = NULL,
                                parallelization = NULL,
                                n_replicates = NULL) {
  futile.logger::flog.info("Starting the analysis of results in environment %s", env)

  # The caller usually already holds the scenarios config; fall back to
  # reading it so a direct call still honours the environment's settings. A
  # config_dir without one keeps the behaviour these arguments replaced:
  # analyse sequentially, and leave the replicate count at the power
  # functions' own default.
  if (is.null(parallelization) || is.null(n_replicates)) {
    scenarios_config_path <- paste0(config_dir, "scenarios_config.yml")
    scenarios_config <- if (file.exists(scenarios_config_path)) {
      read_config(scenarios_config_path, scenarios_config_schema)
    } else {
      list(parallelization = FALSE, n_replicates = 1000)
    }
    if (is.null(parallelization)) {
      parallelization <- scenarios_config$parallelization
    }
    if (is.null(n_replicates)) {
      n_replicates <- scenarios_config$n_replicates
    }
  }
  run_in_parallel <- analysis_runs_in_parallel(parallelization)

  results_dir <- paste0("./results/", env)
  outputs_config <- yaml::read_yaml(system.file("conf/outputs_config.yml", package = "BExTEsim"))
  simulation_config <- read_simulation_config(system.file("conf/simulation_config.yml", package = "BExTEsim"))

  freq_filename <- paste0(results_dir,
                          "/",
                          outputs_config$frequentist_ocs_results_filename)

  # Load the results and apply frequentist_power_at_equivalent_tie
  results_freq_df <- readr::read_csv(freq_filename, col_types = frequentist_col_types, col_names = TRUE)
  results_freq_df <- upgrade_target_size_columns(results_freq_df)

  if (length(case_studies) == 1 && case_studies == "all"){
    case_studies <- unique(results_freq_df$case_study)
  }

  results_freq_subset <- subset(results_freq_df, case_study %in% case_studies)

  if (length(methods) == 1 && methods == "all"){
    methods <- unique(results_freq_subset$method)
  }

  results_freq_subset <- subset(results_freq_subset, method %in% methods)

  if (nrow(results_freq_subset) == 0) {
    warning("Results dataframe is empty.")
    return()
  }

  matching_columns <- c(
    "method",
    "parameters",
    "control_drift",
    "source_denominator",
    "source_denominator_change_factor",
    "case_study",
    "target_to_source_std_ratio",
    "target_sample_size_control", "target_sample_size_treatment",
    "theta_0",
    "null_space",
    "sampling_approximation",
    "summary_measure_likelihood",
    "source_sample_size_treatment",
    "source_sample_size_control",
    "endpoint",
    "source_standard_error",
    "source_treatment_effect_estimate",
    "equivalent_source_sample_size_per_arm"
  )

  # Both power steps read the same simulated trials, which depend on the
  # design alone: the first simulates them, the second reads them back for its
  # separate and pooled powers instead of simulating every design again.
  trial_cache <- new.env(parent = emptyenv())

  # One worker cluster for every step, started the first time a step wants
  # one: standing it up costs 15-30 s, which each step used to pay again.
  cluster <- new_analysis_cluster()
  on.exit(cluster$stop(), add = TRUE)

  if ("frequentist_power_at_equivalent_tie" %in% to_compute){
    new_results_power_df <- frequentist_power_at_equivalent_tie(results = results_freq_subset, analysis_config = analysis_config, simulation_config = simulation_config, parallelization = run_in_parallel, n_replicates = n_replicates, trial_cache = trial_cache, cluster = cluster)

    updated_results_df <- results_freq_df %>%
      dplyr::anti_join(new_results_power_df, by = matching_columns)

    # Keep the updated results in memory for all subsequent analysis steps.
    results_freq_df <- dplyr::bind_rows(updated_results_df, new_results_power_df)
    results_freq_subset <- subset(
      results_freq_df,
      case_study %in% case_studies & method %in% methods
    )

    readr::write_csv(results_freq_df, freq_filename)
  }
  if ("frequentist_power_at_nominal_tie" %in% to_compute){
    new_results_power_df <- frequentist_power_at_nominal_tie(results = results_freq_subset, analysis_config = analysis_config, simulation_config = simulation_config, parallelization = run_in_parallel, n_replicates = n_replicates, trial_cache = trial_cache, cluster = cluster)

    updated_results_df <- results_freq_df %>%
      dplyr::anti_join(new_results_power_df, by = matching_columns)

    # Keep the updated results in memory for all subsequent analysis steps.
    results_freq_df <- dplyr::bind_rows(updated_results_df, new_results_power_df)
    results_freq_subset <- subset(
      results_freq_df,
      case_study %in% case_studies & method %in% methods
    )

    readr::write_csv(results_freq_df, freq_filename)
  }
  if ("sweet_spot" %in% to_compute){
    sweet_spot_df <- sweet_spot(
      results_freq_subset,
      frequentist_metrics,
      nominal_tie = analysis_config$nominal_tie
    )

    jsonlite::write_json(sweet_spot_df, paste0(results_dir, "/sweet_spot.json"))
  }
  if ("bayesian_ocs" %in% to_compute){
    # Compute Bayesian OCs (in a deterministic manner) based on the results
    results_bayesian_ocs <- compute_bayesian_ocs(results_freq_subset, env, config_dir = config_dir, case_studies_config_dir = case_studies_config_dir, parallelization = run_in_parallel, cluster = cluster)
    bayes_filename <- paste0(results_dir,
                             "/",
                             outputs_config$bayesian_ocs_deterministic_results_filename)

    # Append new Bayesian results to existing data
    if (file.exists(bayes_filename)) {
      existing_bayesian_df <- upgrade_target_size_columns(
        readr::read_csv(bayes_filename, col_names = TRUE)
      )

      # Remove rows in `existing_bayesian_df` that match the keys in `results_bayesian_ocs`
      updated_bayesian_df <- existing_bayesian_df %>%
        dplyr::anti_join(results_bayesian_ocs, by = c("case_study", "method"))

      # Append new Bayesian results
      combined_bayesian_df <- dplyr::bind_rows(updated_bayesian_df, results_bayesian_ocs)

      # Write updated results back to the CSV
      return(readr::write_csv(combined_bayesian_df, bayes_filename))
    } else {
      # Write new results if file does not exist
      return(readr::write_csv(results_bayesian_ocs, bayes_filename))
    }
  }
  return(invisible(NULL))
}
