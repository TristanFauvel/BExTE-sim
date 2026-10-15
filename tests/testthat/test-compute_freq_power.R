target_data <- list(
  sample_size_control = 30,
  sample_size_treatment = 30,
  treatment_effect = 1,
  standard_deviation = 1,
  summary_measure_likelihood = "normal",
  endpoint = "continuous"
)

compute_test_power <- function(frequentist_test, null_space) {
  compute_freq_power(
    alpha = 0.05,
    target_data = target_data,
    frequentist_test = frequentist_test,
    theta_0 = 0,
    null_space = null_space,
    simulation_config = list()
  )
}

test_that("compute_freq_power works for t-test with left null space", {
  result <- compute_test_power("t-test", "left")
  expected_power <- pwr::pwr.t.test(
    d = 1,
    n = 30,
    sig.level = 0.05,
    type = "one.sample",
    alternative = "greater"
  )$power

  expect_equal(result$power, expected_power, tolerance = 1e-6)
  expect_equal(result$conf_int_power, rep(expected_power, 2))
})

test_that("compute_freq_power works for t-test with right null space", {
  result <- compute_test_power("t-test", "right")
  expected_power <- pwr::pwr.t.test(
    d = 1,
    n = 30,
    sig.level = 0.05,
    type = "one.sample",
    alternative = "less"
  )$power

  expect_equal(result$power, expected_power, tolerance = 1e-6)
  expect_equal(result$conf_int_power, rep(expected_power, 2))
})

test_that("compute_freq_power works for z-test with left null space", {
  result <- compute_test_power("z-test", "left")
  expected_power <- pwr::pwr.norm.test(
    d = 1,
    n = 30,
    sig.level = 0.05,
    alternative = "greater"
  )$power

  expect_equal(result$power, expected_power, tolerance = 1e-6)
  expect_equal(result$conf_int_power, rep(expected_power, 2))
})

test_that("compute_freq_power works for z-test with right null space", {
  result <- compute_test_power("z-test", "right")
  expected_power <- pwr::pwr.norm.test(
    d = 1,
    n = 30,
    sig.level = 0.05,
    alternative = "less"
  )$power

  expect_equal(result$power, expected_power, tolerance = 1e-6)
  expect_equal(result$conf_int_power, rep(expected_power, 2))
})

test_that("compute_freq_power rejects an invalid null space", {
  expect_error(
    compute_test_power("t-test", "middle"),
    "Null space must be either 'left' or 'right'"
  )
})

test_that("compute_freq_power preserves its return structure for missing alpha", {
  result <- compute_freq_power(
    alpha = NA_real_,
    target_data = target_data,
    frequentist_test = "t-test",
    theta_0 = 0,
    null_space = "left",
    simulation_config = list()
  )

  expect_named(result, c("power", "conf_int_power"))
  expect_true(is.na(result$power))
  expect_equal(result$conf_int_power, rep(NA_real_, 2))
})

test_that("compute_power_with_tie_ci preserves its return structure for a missing interval", {
  result <- compute_power_with_tie_ci(
    alpha = list(
      mean = 0.05,
      conf_int_lower = NA_real_,
      conf_int_upper = 0.06
    ),
    target_data = NULL,
    frequentist_test = "t-test",
    theta_0 = 0,
    null_space = "left",
    simulation_config = list()
  )

  expect_named(result, c("power", "conf_int_power", "n_effective_samples"))
  expect_true(is.na(result$power))
  expect_equal(result$conf_int_power, rep(NA_real_, 2))
  expect_equal(result$n_effective_samples, 0L)
})

test_that("compute_power_with_tie_ci handles a degenerate type I error of zero", {
  result <- compute_power_with_tie_ci(
    alpha = list(
      mean = 0,
      conf_int_lower = 0,
      conf_int_upper = 0
    ),
    target_data = NULL,
    frequentist_test = "t-test",
    theta_0 = 0,
    null_space = "left",
    simulation_config = list(),
    n_samples = 10
  )

  expect_named(result, c("power", "conf_int_power", "n_effective_samples"))
  expect_true(is.na(result$power))
  expect_equal(result$conf_int_power, rep(NA_real_, 2))
  expect_equal(result$n_effective_samples, 0L)
})
