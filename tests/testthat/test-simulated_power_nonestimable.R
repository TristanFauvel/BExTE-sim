## The frequentist baselines are simulated whenever the endpoint is not
## continuous. Mepolizumab used to be the exception, priced in closed form
## although its trials are generated patient by patient from a negative
## binomial, so its baselines assumed a test the Bayesian methods were never
## compared against.
##
## A simulated trial can also have no estimable summary measure - a
## recurrent-event arm with no event at all. It has no p-value, and before it
## was dropped a single such replicate turned the whole simulated power into
## NA. The Bayesian operating characteristics leave those replicates out of
## their denominator, and the baselines now do the same.

test_that("recurrent events take the simulated power branch", {
  expect_false(uses_analytical_power(list(endpoint = "recurrent_event")))
  expect_false(uses_analytical_power(list(endpoint = "time_to_event")))
  expect_false(uses_analytical_power(list(endpoint = "binary")))
  expect_true(uses_analytical_power(list(endpoint = "continuous")))
})

## Every third replicate is non-estimable, as negative_binomial_regression()
## reports an arm with no event: an NA estimate and an NA standard error.
estimates <- c(-1, -1, NA, -1, 0.5, NA)

recurrent_target_data <- list(
  summary_measure_likelihood = "normal",
  endpoint = "recurrent_event",
  treatment_effect = -0.5,
  standard_deviation = 1,
  sample_size_control = 30,
  sample_size_treatment = 30,
  generate = function(n_replicates) {
    estimate <- rep_len(estimates, n_replicates)
    data.frame(
      treatment_effect_estimate = estimate,
      treatment_effect_standard_error = ifelse(is.na(estimate), NA_real_, 0.2),
      sample_size_control = 30,
      sample_size_treatment = 30,
      standard_deviation = ifelse(is.na(estimate), NA_real_, 0.2 * sqrt(30))
    )
  }
)

## BSDA::tsum.test() notes "argument 'var.equal' ignored for one-sample test."
## on every one-sample call. Muffle just that, so a real warning still shows.
without_bsda_one_sample_note <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("var.equal", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

test_that("the separate analysis power leaves non-estimable replicates out", {
  result <- without_bsda_one_sample_note(compute_freq_power(
    alpha = 0.025,
    target_data = recurrent_target_data,
    frequentist_test = "t-test",
    theta_0 = 0,
    null_space = "right",
    simulation_config = list(seed = 1),
    n_replicates = 60
  ))

  ## Of the four estimable replicates in each block, the three at -1 reject
  ## and the one at 0.5 does not.
  expect_equal(result$power, 0.75)
  expect_true(all(is.finite(result$conf_int_power)))
})

test_that("the power at equivalent type I error leaves them out too", {
  result <- without_bsda_one_sample_note(compute_power_with_tie_ci(
    alpha = list(mean = 0.025, conf_int_lower = 0.02, conf_int_upper = 0.03,
                 mcse = NA_real_, n_replicates = NA_real_),
    target_data = recurrent_target_data,
    frequentist_test = "t-test",
    theta_0 = 0,
    null_space = "right",
    simulation_config = list(seed = 1),
    n_replicates = 60,
    n_samples = 50
  ))

  expect_true(is.finite(result$power))
  expect_equal(result$power, 0.75, tolerance = 0.1)
})

test_that("the pooled analysis power leaves them out too", {
  result <- without_bsda_one_sample_note(compute_freq_power_pooling(
    alpha = 0.025,
    target_data = recurrent_target_data,
    source_data = list(
      treatment_effect_estimate = -0.5,
      standard_error = 0.1,
      equivalent_source_sample_size_per_arm = 300
    ),
    frequentist_test = "t-test",
    theta_0 = 0,
    null_space = "right",
    simulation_config = list(seed = 1),
    case_study = "mepolizumab",
    n_replicates = 60
  ))

  expect_true(is.finite(result$power))
  expect_true(all(is.finite(result$conf_int_power)))
})

test_that("a scenario with no estimable replicate fails loudly", {
  all_missing <- recurrent_target_data
  all_missing$generate <- function(n_replicates) {
    data.frame(
      treatment_effect_estimate = rep(NA_real_, n_replicates),
      treatment_effect_standard_error = NA_real_,
      sample_size_control = 30,
      sample_size_treatment = 30,
      standard_deviation = NA_real_
    )
  }

  expect_error(
    compute_freq_power(
      alpha = 0.025,
      target_data = all_missing,
      frequentist_test = "t-test",
      theta_0 = 0,
      null_space = "right",
      simulation_config = list(seed = 1),
      n_replicates = 10
    ),
    "No simulated trial has an estimable summary measure"
  )
})
