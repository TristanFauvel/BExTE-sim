analytical_target_data <- list(
  sample_size_control = 30,
  sample_size_treatment = 30,
  treatment_effect = 0.5,
  standard_deviation = 1,
  summary_measure_likelihood = "normal",
  endpoint = "continuous"
)

simulated_target_data <- list(
  sample_size_control = 30,
  sample_size_treatment = 30,
  treatment_effect = 0.5,
  standard_deviation = 1,
  summary_measure_likelihood = "normal",
  endpoint = "recurrent_event",
  generate = function(n_replicates) {
    data.frame(
      treatment_effect_estimate = rnorm(n_replicates, 0.5, 1 / sqrt(30)),
      standard_deviation = 1,
      sample_size_control = 30,
      sample_size_treatment = 30
    )
  }
)

tie_estimate <- function(mean, lower, upper, mcse = NA_real_) {
  list(
    mean = mean,
    conf_int_lower = lower,
    conf_int_upper = upper,
    mcse = mcse,
    n_replicates = NA_real_
  )
}

test_that("compute_power_with_tie_ci propagates Monte Carlo error when the type I error is known exactly", {
  # With a degenerate type I error interval the only remaining source of
  # uncertainty is the Monte Carlo error of the simulated power, which the
  # reported interval must still reflect.
  result <- compute_power_with_tie_ci(
    alpha = tie_estimate(0.025, 0.025, 0.025, mcse = 0),
    target_data = simulated_target_data,
    frequentist_test = "t-test",
    theta_0 = 0,
    null_space = "left",
    simulation_config = list(seed = 42),
    case_study = "example",
    n_replicates = 200,
    n_samples = 50
  )

  expect_false(is.na(result$power))
  expect_gt(result$conf_int_power[2] - result$conf_int_power[1], 0)
  expect_lte(result$conf_int_power[1], result$power)
  expect_gte(result$conf_int_power[2], result$power)
})

test_that("compute_power_with_tie_ci widens the interval as the type I error interval widens", {
  narrow <- compute_power_with_tie_ci(
    alpha = tie_estimate(0.025, 0.020, 0.031, mcse = NA_real_),
    target_data = analytical_target_data,
    frequentist_test = "t-test",
    theta_0 = 0,
    null_space = "left",
    simulation_config = list(seed = 42),
    case_study = "example",
    n_samples = 4000
  )
  wide <- compute_power_with_tie_ci(
    alpha = tie_estimate(0.025, 0.004, 0.080, mcse = NA_real_),
    target_data = analytical_target_data,
    frequentist_test = "t-test",
    theta_0 = 0,
    null_space = "left",
    simulation_config = list(seed = 42),
    case_study = "example",
    n_samples = 4000
  )

  expect_gt(
    wide$conf_int_power[2] - wide$conf_int_power[1],
    narrow$conf_int_power[2] - narrow$conf_int_power[1]
  )
})

test_that("compute_power_with_tie_ci uses every draw when the type I error estimate sits against zero", {
  # A normal approximation centred on a near-zero type I error puts a large part
  # of its mass below zero; those draws used to be silently discarded, biasing
  # the retained sample upwards and shrinking the effective sample size.
  result <- compute_power_with_tie_ci(
    alpha = tie_estimate(0, 0, 0.0037, mcse = 0),
    target_data = analytical_target_data,
    frequentist_test = "t-test",
    theta_0 = 0,
    null_space = "left",
    simulation_config = list(seed = 42),
    case_study = "example",
    n_samples = 500
  )

  expect_false(is.na(result$power))
  expect_equal(result$n_effective_samples, 500)
})

test_that("compute_power_with_tie_ci reports the number of draws actually used", {
  result <- compute_power_with_tie_ci(
    alpha = tie_estimate(0.025, 0.016, 0.037, mcse = NA_real_),
    target_data = analytical_target_data,
    frequentist_test = "t-test",
    theta_0 = 0,
    null_space = "left",
    simulation_config = list(seed = 42),
    case_study = "example",
    n_samples = 137
  )

  expect_equal(result$n_effective_samples, 137)
})

test_that("compute_power_with_tie_ci returns NA when the type I error interval is missing", {
  result <- compute_power_with_tie_ci(
    alpha = tie_estimate(0.025, NA_real_, NA_real_),
    target_data = analytical_target_data,
    frequentist_test = "t-test",
    theta_0 = 0,
    null_space = "left",
    simulation_config = list(seed = 42),
    case_study = "example"
  )

  expect_true(is.na(result$power))
  expect_true(all(is.na(result$conf_int_power)))
})
