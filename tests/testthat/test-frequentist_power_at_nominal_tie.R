sample_results <- data.frame(
  case_study = "example",
  method = "separate",
  parameters = "{}",
  theta_0 = 0,
  null_space = "left",
  target_treatment_effect = c(0, 0.5)
)

mock_load_data <- function(results_row, type, reload_data_objects = FALSE) {
  list(
    type = type,
    treatment_effect = 0.5,
    standard_deviation = 1,
    sample_size_control = 30,
    sample_size_treatment = 30
  )
}

mock_freq_power <- function(...) {
  list(power = 0.8, conf_int_power = c(0.75, 0.85))
}

mock_freq_power_pooling <- function(...) {
  list(power = 0.9, conf_int_power = c(0.87, 0.93))
}

# The function drives a txtProgressBar, which writes straight to the console.
run_at_nominal_tie <- function(results) {
  final_results <- NULL
  capture.output(
    final_results <- frequentist_power_at_nominal_tie(
      results = results,
      analysis_config = list(nominal_tie = 0.025, frequentist_test = "t-test"),
      simulation_config = list(),
      n_replicates = 100
    )
  )
  final_results
}

test_that("frequentist_power_at_nominal_tie keeps the power confidence bounds", {
  with_mocked_bindings(
    {
      final_results <- run_at_nominal_tie(sample_results)

      expect_equal(final_results$nominal_frequentist_power_separate, rep(0.8, 2))
      expect_equal(final_results$nominal_frequentist_power_separate_lower, rep(0.75, 2))
      expect_equal(final_results$nominal_frequentist_power_separate_upper, rep(0.85, 2))

      expect_equal(final_results$nominal_frequentist_power_pooling, rep(0.9, 2))
      expect_equal(final_results$nominal_frequentist_power_pooling_lower, rep(0.87, 2))
      expect_equal(final_results$nominal_frequentist_power_pooling_upper, rep(0.93, 2))
    },
    load_data = mock_load_data,
    compute_freq_power = mock_freq_power,
    compute_freq_power_pooling = mock_freq_power_pooling,
    .package = "BExTEsim"
  )
})

test_that("frequentist_power_at_nominal_tie overwrites stale bounds rather than duplicating them", {
  stale_results <- sample_results
  stale_results$nominal_frequentist_power_separate <- 0.1
  stale_results$nominal_frequentist_power_separate_lower <- 0.05
  stale_results$nominal_frequentist_power_separate_upper <- 0.15
  stale_results$nominal_frequentist_power_pooling <- 0.2
  stale_results$nominal_frequentist_power_pooling_lower <- 0.15
  stale_results$nominal_frequentist_power_pooling_upper <- 0.25

  with_mocked_bindings(
    {
      final_results <- run_at_nominal_tie(stale_results)

      bound_columns <- c(
        "nominal_frequentist_power_separate_lower",
        "nominal_frequentist_power_separate_upper",
        "nominal_frequentist_power_pooling_lower",
        "nominal_frequentist_power_pooling_upper"
      )
      expect_equal(
        vapply(bound_columns, function(x) sum(names(final_results) == x), integer(1)),
        setNames(rep(1L, length(bound_columns)), bound_columns)
      )

      expect_equal(final_results$nominal_frequentist_power_separate_lower, rep(0.75, 2))
      expect_equal(final_results$nominal_frequentist_power_pooling_upper, rep(0.93, 2))
    },
    load_data = mock_load_data,
    compute_freq_power = mock_freq_power,
    compute_freq_power_pooling = mock_freq_power_pooling,
    .package = "BExTEsim"
  )
})

test_that("frequentist_power_at_nominal_tie leaves degenerate closed-form bounds degenerate", {
  # compute_freq_power() returns c(power, power) whenever it has a closed form,
  # which is what the plot layer reads to decide against drawing error bars.
  mock_analytical_power <- function(...) {
    list(power = 0.8, conf_int_power = c(0.8, 0.8))
  }

  with_mocked_bindings(
    {
      final_results <- run_at_nominal_tie(sample_results)

      expect_equal(
        final_results$nominal_frequentist_power_separate_lower,
        final_results$nominal_frequentist_power_separate_upper
      )
    },
    load_data = mock_load_data,
    compute_freq_power = mock_analytical_power,
    compute_freq_power_pooling = mock_analytical_power,
    .package = "BExTEsim"
  )
})

test_that("frequentist_power_at_nominal_tie rejects empty results", {
  expect_error(
    frequentist_power_at_nominal_tie(
      results = data.frame(),
      analysis_config = list(nominal_tie = 0.025, frequentist_test = "t-test"),
      simulation_config = list()
    ),
    "The results dataframe is empty"
  )
})
