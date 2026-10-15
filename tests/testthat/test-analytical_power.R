## compute_power_with_tie_ci() prices 1000 sampled type I errors per result
## row. Every pwr function used here takes a vector of significance levels, so
## those are computed in one call rather than 1000; the per-call overhead of
## pwr (match.arg, argument assembly, simplify2array) costs more than the power
## calculation itself. These tests pin the vectorised answers to the values the
## one-at-a-time loop produced.

normal_target <- list(
  summary_measure_likelihood = "normal",
  endpoint = "continuous",
  treatment_effect = 0.4,
  standard_deviation = 1.1,
  sample_size_control = 60,
  sample_size_treatment = 60
)

binomial_target <- list(
  summary_measure_likelihood = "binomial",
  endpoint = "binary",
  treatment_rate = 0.55,
  control_rate = 0.4,
  sample_size_control = 80,
  sample_size_treatment = 80
)

alphas <- seq(0.005, 0.15, length.out = 50)

test_that("analytical_power over a vector matches one call per alpha (t-test)", {
  vectorised <- analytical_power(alphas, normal_target, "t-test", 0, "greater")
  one_at_a_time <- vapply(alphas, function(a) {
    analytical_power(a, normal_target, "t-test", 0, "greater")
  }, numeric(1))

  expect_equal(vectorised, one_at_a_time)
  expect_length(vectorised, length(alphas))
})

test_that("analytical_power over a vector matches one call per alpha (z-test)", {
  vectorised <- analytical_power(alphas, normal_target, "z-test", 0, "less")
  one_at_a_time <- vapply(alphas, function(a) {
    analytical_power(a, normal_target, "z-test", 0, "less")
  }, numeric(1))

  expect_equal(vectorised, one_at_a_time)
})

test_that("analytical_power over a vector matches one call per alpha (binomial)", {
  vectorised <- analytical_power(alphas, binomial_target, "t-test", 0, "greater")
  one_at_a_time <- vapply(alphas, function(a) {
    analytical_power(a, binomial_target, "t-test", 0, "greater")
  }, numeric(1))

  expect_equal(vectorised, one_at_a_time)
})

test_that("the vectorised path returns what compute_freq_power returns per alpha", {
  vectorised <- analytical_power(alphas, normal_target, "t-test", 0, "greater")
  through_compute_freq_power <- vapply(alphas, function(a) {
    compute_freq_power(
      alpha = a, target_data = normal_target, frequentist_test = "t-test",
      theta_0 = 0, null_space = "left", simulation_config = list()
    )$power
  }, numeric(1))

  expect_equal(vectorised, through_compute_freq_power)
})

test_that("the null space still decides the direction, and is still validated", {
  expect_identical(alternative_from_null_space("left"), "greater")
  expect_identical(alternative_from_null_space("right"), "less")
  expect_error(alternative_from_null_space("middle"), "Null space must be either")

  # compute_freq_power() reports it through the same helper.
  expect_error(
    compute_freq_power(
      alpha = 0.05, target_data = normal_target, frequentist_test = "t-test",
      theta_0 = 0, null_space = "middle", simulation_config = list()
    ),
    "Null space must be either 'left' or 'right'"
  )
})

test_that("unsupported tests and likelihoods are still rejected", {
  expect_error(analytical_power(0.05, normal_target, "chi-squared", 0, "greater"),
               "Unsupported test type.")
  expect_error(
    analytical_power(0.05, utils::modifyList(normal_target,
      list(summary_measure_likelihood = "poisson")), "t-test", 0, "greater"),
    "Unsupported likelihood type."
  )
})
