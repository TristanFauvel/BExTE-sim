## The Monte Carlo branch of compute_freq_power(), compute_freq_power_pooling()
## and compute_power_with_tie_ci() draws `n_replicates` trials per scenario. The
## analysis drivers used to call all three without that argument, so every
## simulated power estimate was built from the function default of 1000 draws
## however many replicates the run was configured with. On the case studies that
## take the Monte Carlo branch - the ones uses_analytical_power() rejects, such
## as belimumab - that left the nominal power series roughly three times noisier
## than the 10000-replicate success probability plotted beside it.

minimal_results <- data.frame(
  case_study = "example",
  theta_0 = 0,
  null_space = "left"
)

mock_load_data <- function(results_row, type, reload_data_objects = FALSE) {
  list(type = type)
}

skip_target_data_assertions <- function(target_data) {
  invisible(NULL)
}

test_that("frequentist_power_at_nominal_tie draws the replicates it is given", {
  drawn <- new.env(parent = emptyenv())

  record_separate <- function(..., n_replicates = 1000) {
    drawn$separate <- n_replicates
    list(power = 0.5, conf_int_power = c(0.4, 0.6))
  }
  record_pooling <- function(..., n_replicates = 1000) {
    drawn$pooling <- n_replicates
    list(power = 0.5, conf_int_power = c(0.4, 0.6))
  }

  with_mocked_bindings(
    frequentist_power_at_nominal_tie(
      results = minimal_results,
      analysis_config = list(frequentist_test = "t-test", nominal_tie = 0.025),
      simulation_config = list(),
      n_replicates = 10000
    ),
    load_data = mock_load_data,
    assert_target_data_numbers = skip_target_data_assertions,
    compute_freq_power = record_separate,
    compute_freq_power_pooling = record_pooling,
    .package = "BExTEsim"
  )

  expect_equal(drawn$separate, 10000)
  expect_equal(drawn$pooling, 10000)
})

test_that("frequentist_power_at_equivalent_tie draws the replicates it is given", {
  drawn <- new.env(parent = emptyenv())

  record_power_with_tie_ci <- function(..., n_replicates = 1000) {
    drawn$equivalent <- n_replicates
    list(power = 0.8, conf_int_power = c(0.75, 0.85))
  }

  tie_results <- data.frame(
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
    equivalent_source_sample_size_per_arm = 50,
    target_treatment_effect = c(0, 0.5),
    success_proba = c(0.025, 0.8),
    mcse_success_proba = c(0.005, 0.02),
    conf_int_success_proba_lower = c(0.015, 0.76),
    conf_int_success_proba_upper = c(0.035, 0.84)
  )

  with_mocked_bindings(
    frequentist_power_at_equivalent_tie(
      results = tie_results,
      analysis_config = list(frequentist_test = "t-test"),
      simulation_config = list(),
      parallelization = FALSE,
      n_replicates = 10000
    ),
    load_data = mock_load_data,
    compute_power_with_tie_ci = record_power_with_tie_ci,
    .package = "BExTEsim"
  )

  expect_equal(drawn$equivalent, 10000)
})

test_that("simulation_analysis forwards the configured n_replicates", {
  original_directory <- getwd()
  temporary_directory <- tempfile("n-replicates-")
  config_directory <- file.path(temporary_directory, "config")
  dir.create(file.path(temporary_directory, "results", "test"), recursive = TRUE)
  dir.create(config_directory, recursive = TRUE)
  on.exit({
    setwd(original_directory)
    unlink(temporary_directory, recursive = TRUE)
  }, add = TRUE)

  writeLines(
    c(
      "n_replicates: 4321",
      "ndrift: 30",
      "parallelization: FALSE",
      "denominator_change_factor:",
      "    - 1",
      "sample_size_factors:",
      "    - 1",
      "case_studies:",
      "    - example",
      "methods:",
      "    - separate"
    ),
    file.path(config_directory, "scenarios_config.yml")
  )

  setwd(temporary_directory)

  results_row <- data.frame(
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
    names(results_row)
  )
  results_row[missing_columns] <- NA
  readr::write_csv(results_row, file.path("results", "test", "results_frequentist.csv"))

  forwarded <- new.env(parent = emptyenv())
  record_equivalent <- function(results, ..., n_replicates = 1000) {
    forwarded$equivalent <- n_replicates
    results
  }
  record_nominal <- function(results, ..., n_replicates = 1000) {
    forwarded$nominal <- n_replicates
    results
  }

  with_mocked_bindings(
    simulation_analysis(
      env = "test",
      analysis_config = list(frequentist_test = "t-test", nominal_tie = 0.025),
      config_dir = paste0(config_directory, .Platform$file.sep),
      frequentist_metrics = list(),
      to_compute = c(
        "frequentist_power_at_equivalent_tie",
        "frequentist_power_at_nominal_tie"
      )
    ),
    frequentist_power_at_equivalent_tie = record_equivalent,
    frequentist_power_at_nominal_tie = record_nominal,
    .package = "BExTEsim"
  )

  expect_equal(forwarded$equivalent, 4321)
  expect_equal(forwarded$nominal, 4321)
})
