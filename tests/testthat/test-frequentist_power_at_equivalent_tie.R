sample_results <- data.frame(
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

mock_load_data <- function(results_row, type, reload_data_objects = FALSE) {
  list(type = type)
}

test_that("frequentist_power_at_equivalent_tie adds power estimates", {
  mock_power_with_tie_ci <- function(...) {
    list(power = 0.8, conf_int_power = c(0.75, 0.85))
  }

  with_mocked_bindings(
    {
      final_results <- frequentist_power_at_equivalent_tie(
        results = sample_results,
        analysis_config = list(frequentist_test = "t-test"),
        simulation_config = list(),
        parallelization = FALSE
      )

      expect_equal(final_results$tie, rep(0.025, 2))
      expect_equal(final_results$frequentist_power_at_equivalent_tie, rep(0.8, 2))
      expect_equal(final_results$frequentist_power_at_equivalent_tie_lower, rep(0.75, 2))
      expect_equal(final_results$frequentist_power_at_equivalent_tie_upper, rep(0.85, 2))
      expect_equal(final_results$frequentist_test, rep("t-test", 2))
    },
    load_data = mock_load_data,
    compute_power_with_tie_ci = mock_power_with_tie_ci,
    .package = "BExTEsim"
  )
})

test_that("frequentist_power_at_equivalent_tie rejects empty results", {
  expect_error(
    frequentist_power_at_equivalent_tie(
      results = data.frame(),
      analysis_config = list(frequentist_test = "t-test"),
      simulation_config = list()
    ),
    "The results dataframe is empty"
  )
})

test_that("frequentist_power_at_equivalent_tie requires a null scenario", {
  results_without_null <- data.frame(
    case_study = "example",
    target_treatment_effect = 0.5,
    theta_0 = 0
  )

  expect_error(
    frequentist_power_at_equivalent_tie(
      results = results_without_null,
      analysis_config = list(frequentist_test = "t-test"),
      simulation_config = list()
    ),
    "theta_0 not included among the target study treatment effects"
  )
})
