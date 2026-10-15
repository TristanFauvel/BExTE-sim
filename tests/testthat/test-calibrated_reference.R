# The reference test every borrowing method's power is read against is
# calibrated on its actual type I error in the design's null scenario, not on
# its nominal level: by simulation where the operating characteristics are
# simulated, and exactly, with a randomised boundary, where they are
# enumerated. These tests pin both, and that calibration removes the spurious
# power gain a t-test on an approximately normal statistic gave every method.

test_that("calibrated levels reject the target share of the null trials", {
  set.seed(1)
  null_p_values <- runif(10000)

  for (alpha in c(0, 0.00005, 0.0238, 0.025, 0.1, 0.5, 1)) {
    threshold <- calibrated_levels(alpha, null_p_values)
    expect_equal(mean(null_p_values < threshold), floor(alpha * 10000 + 1e-9) / 10000,
                 info = alpha)
  }
})

test_that("calibration removes the gain a t-test on a normal statistic gives", {
  # The statistic is exactly standard normal under the null, as a Wald
  # statistic is approximately. A z-test at 1.96 - which is what a separate
  # Bayesian analysis with a flat prior amounts to - is the "method"; the
  # reference is the t-test on 122 degrees of freedom, which is conservative
  # on this statistic.
  set.seed(2)
  n <- 123
  standard_error <- 0.1
  trials <- function(mean) {
    data.frame(treatment_effect_estimate = rnorm(40000, mean, standard_error),
               standard_deviation = standard_error * sqrt(n), sample_size_control = n, sample_size_treatment = n)
  }
  null_trials <- trials(0)
  alternative_trials <- trials(2.2 * standard_error)

  method_tie <- mean(null_trials$treatment_effect_estimate / standard_error > 1.96)
  method_power <- mean(alternative_trials$treatment_effect_estimate / standard_error > 1.96)

  null_p <- trial_p_values(null_trials, theta_0 = 0, alternative = "greater")
  alternative_p <- trial_p_values(alternative_trials, theta_0 = 0, alternative = "greater")

  # At the method's true size, 0.025, the t-test's nominal level understates
  # its critical value - 1.980 against 1.960 - and so its power, by about
  # 0.008 here.
  uncalibrated_gain <- method_power - mean(alternative_p < 0.025)
  expect_gt(uncalibrated_gain, 0.005)

  # Calibrated on the null trials the method's type I error was estimated on,
  # the reference matches it in those trials, which also cancels the Monte
  # Carlo error of that estimate: the gain is the window between the two
  # critical values only, a fraction of a trial's worth.
  calibrated_gain <- method_power - mean(alternative_p < calibrated_levels(method_tie, null_p))
  expect_lt(abs(calibrated_gain), 0.002)
})

exact_aprepitant <- function(n, treatment_drift) {
  case_study_config <- yaml::read_yaml(
    system.file("conf", "case_studies", "aprepitant.yml", package = "BExTE")
  )
  source_data <- SourceData$new(case_study_config)
  TargetDataFactory$new()$create(
    source_data = source_data, case_study_config = case_study_config,
    target_sample_size_control = n,
    target_sample_size_treatment = n, control_drift = 0, treatment_drift = treatment_drift,
    summary_measure_likelihood = case_study_config$summary_measure_likelihood
  )
}

test_that("the exact reference rejects exactly alpha in the null scenario", {
  null_target <- exact_aprepitant(40, treatment_drift = -0.0796)
  null_support <- null_target$enumerate_support()

  for (alpha in c(0.001, 0.0238, 0.025, 0.0361, 0.2)) {
    # Priced on the null scenario itself, the power is the size.
    expect_equal(
      exact_reference_power(alpha, null_support, null_support, theta_0 = 0, alternative = "greater"),
      alpha, tolerance = 1e-12, info = alpha
    )
  }
})

test_that("the exact reference is non-decreasing in alpha and unrandomised at attainable levels", {
  null_support <- exact_aprepitant(40, treatment_drift = -0.0796)$enumerate_support()
  alternative_support <- exact_aprepitant(40, treatment_drift = 0.1)$enumerate_support()

  levels <- seq(0.001, 0.2, length.out = 60)
  power <- exact_reference_power(levels, alternative_support, null_support, 0, "greater")
  expect_false(is.unsorted(power))

  # At the null mass below some outcome's p-value, the test rejects exactly the
  # outcomes below it, with no randomisation.
  null <- reference_outcomes(null_support, 0, "greater")
  alternative <- reference_outcomes(alternative_support, 0, "greater")
  threshold <- sort(unique(null$p[is.finite(null$p)]))[25]
  attainable <- sum(null$w[null$p < threshold])
  expect_equal(
    exact_reference_power(attainable, alternative_support, null_support, 0, "greater"),
    sum(alternative$w[alternative$p < threshold]),
    tolerance = 1e-12
  )
})

test_that("an enumerated design's reference is priced exactly, without simulation", {
  null_target <- exact_aprepitant(40, treatment_drift = -0.0796)
  alternative_target <- exact_aprepitant(40, treatment_drift = 0)
  targets <- list(null_target, alternative_target)

  results <- data.frame(
    method = "separate", parameters = "{}", control_drift = 0,
    source_denominator = NA_real_, source_denominator_change_factor = 1,
    case_study = "aprepitant", target_to_source_std_ratio = NA_real_,
    target_sample_size_control = 40,
    target_sample_size_treatment = 40, theta_0 = 0, null_space = "left",
    sampling_approximation = FALSE, summary_measure_likelihood = "binomial",
    source_sample_size_treatment = 293, source_sample_size_control = 280,
    endpoint = "binary", source_standard_error = 0.04,
    source_treatment_effect_estimate = 0.0796,
    equivalent_source_sample_size_per_arm = 286,
    target_treatment_effect = c(0, 0.0796),
    success_proba = c(0.0361, 0.4), mcse_success_proba = 0,
    conf_int_success_proba_lower = c(0.0361, 0.4),
    conf_int_success_proba_upper = c(0.0361, 0.4),
    exact_ocs = TRUE
  )

  powers <- with_mocked_bindings(
    frequentist_power_at_equivalent_tie(
      results = results,
      analysis_config = list(frequentist_test = "t-test"),
      simulation_config = list(seed = 1),
      n_replicates = 100
    ),
    load_data = function(results_row, type, reload_data_objects = FALSE) {
      targets[[which(c(0, 0.0796) == results_row$target_treatment_effect)]]
    },
    simulated_trials = function(...) stop("an enumerated design was simulated"),
    .package = "BExTEsim"
  )

  expected <- exact_reference_power(0.0361, alternative_target$enumerate_support(),
                                    null_target$enumerate_support(), 0, "greater")
  alternative_row <- powers[powers$target_treatment_effect == 0.0796, ]
  expect_equal(alternative_row$frequentist_power_at_equivalent_tie, expected)
  expect_equal(alternative_row$frequentist_power_at_equivalent_tie_lower, expected)
  expect_equal(alternative_row$frequentist_power_at_equivalent_tie_upper, expected)
  expect_equal(alternative_row$frequentist_reference_calibration, "exact, randomised")
  # In the null scenario the reference's power is its size, the method's TIE.
  expect_equal(powers$frequentist_power_at_equivalent_tie[powers$target_treatment_effect == 0],
               0.0361, tolerance = 1e-12)
})
