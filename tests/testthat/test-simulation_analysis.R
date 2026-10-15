test_that("simulation analysis carries results through sequential steps", {
  original_directory <- getwd()
  temporary_directory <- tempfile("simulation-analysis-")
  dir.create(file.path(temporary_directory, "results", "test"), recursive = TRUE)
  on.exit({
    setwd(original_directory)
    unlink(temporary_directory, recursive = TRUE)
  }, add = TRUE)
  setwd(temporary_directory)

  matching_result <- data.frame(
    method = "separate",
    parameters = "{}",
    control_drift = 0,
    source_denominator = NA_real_,
    source_denominator_change_factor = 1,
    case_study = "example",
    target_to_source_std_ratio = 1,
    target_sample_size_control = 30,
    target_sample_size_treatment = 30,
    theta_0 = 0,
    null_space = "left",
    sampling_approximation = TRUE,
    summary_measure_likelihood = "normal",
    source_sample_size_treatment = 50,
    source_sample_size_control = 50,
    endpoint = "continuous",
    source_standard_error = 0.2,
    source_treatment_effect_estimate = 0.5,
    equivalent_source_sample_size_per_arm = 50
  )
  missing_columns <- setdiff(
    names(BExTEsim:::frequentist_col_types$cols),
    names(matching_result)
  )
  matching_result[missing_columns] <- NA
  results_path <- file.path("results", "test", "results_frequentist.csv")
  readr::write_csv(matching_result, results_path)

  analysis_inputs <- new.env(parent = emptyenv())
  mock_equivalent_power <- function(results, ...) {
    results$equivalent_power_computed <- TRUE
    results
  }
  mock_nominal_power <- function(results, ...) {
    expect_true(all(results$equivalent_power_computed))
    results$nominal_power_computed <- TRUE
    results
  }
  mock_sweet_spot <- function(results, ...) {
    analysis_inputs$sweet_spot <- results
    list(status = "computed")
  }
  mock_bayesian_ocs <- function(results_freq_df, ...) {
    analysis_inputs$bayesian_ocs <- results_freq_df
    data.frame(
      case_study = unique(results_freq_df$case_study),
      method = unique(results_freq_df$method)
    )
  }

  with_mocked_bindings(
    simulation_analysis(
      env = "test",
      analysis_config = list(
        frequentist_test = "t-test",
        nominal_tie = 0.025
      ),
      config_dir = tempfile(),
      frequentist_metrics = list(),
      to_compute = c(
        "frequentist_power_at_equivalent_tie",
        "frequentist_power_at_nominal_tie",
        "sweet_spot",
        "bayesian_ocs"
      )
    ),
    frequentist_power_at_equivalent_tie = mock_equivalent_power,
    frequentist_power_at_nominal_tie = mock_nominal_power,
    sweet_spot = mock_sweet_spot,
    compute_bayesian_ocs = mock_bayesian_ocs,
    .package = "BExTEsim"
  )

  final_results <- readr::read_csv(results_path, show_col_types = FALSE)
  expect_true(all(final_results$equivalent_power_computed))
  expect_true(all(final_results$nominal_power_computed))
  expect_true(all(analysis_inputs$sweet_spot$equivalent_power_computed))
  expect_true(all(analysis_inputs$sweet_spot$nominal_power_computed))
  expect_true(all(analysis_inputs$bayesian_ocs$equivalent_power_computed))
  expect_true(all(analysis_inputs$bayesian_ocs$nominal_power_computed))
})
