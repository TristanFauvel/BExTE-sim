# The Monte Carlo branch of compute_freq_power_pooling() is exercised whenever
# the summary measure is normal but the endpoint is not continuous (belimumab,
# teriflunomide). It must run the *pooled* test: the pooled treatment effect
# compared against the pooled standard error, on the pooled sample size. Using
# the target-only standard error instead is not a pooled test at all, and at the
# teriflunomide design point it turns a power of 1 into a power of 0.

# Teriflunomide-like inputs: a very precise source study and a small target one.
source_data <- list(
  treatment_effect_estimate = -0.393,
  standard_error = 0.081,
  equivalent_source_sample_size_per_arm = 2 * 752 * 731 / (752 + 731)
)

target_sample_size_per_arm <- 83
target_standard_error <- 0.2668287130568081

# A target data stub whose replicates are all identical, so that the Monte Carlo
# power is a deterministic 0 or 1 and pins down which test was run.
make_target_data <- function(treatment_effect_estimate,
                             sample_size_per_arm = target_sample_size_per_arm,
                             standard_error = target_standard_error) {
  standard_deviation <- standard_error * sqrt(sample_size_per_arm)

  list(
    summary_measure_likelihood = "normal",
    endpoint = "time_to_event",
    treatment_effect = treatment_effect_estimate,
    standard_deviation = standard_deviation,
    sample_size_control = sample_size_per_arm,
    sample_size_treatment = sample_size_per_arm,
    generate = function(n_replicates) {
      data.frame(
        treatment_effect_estimate = rep(treatment_effect_estimate, n_replicates),
        treatment_effect_standard_error = rep(standard_error, n_replicates),
        sample_size_control = rep(sample_size_per_arm, n_replicates),
        sample_size_treatment = rep(sample_size_per_arm, n_replicates),
        standard_deviation = rep(standard_deviation, n_replicates)
      )
    }
  )
}

pooled_estimate <- function(target_treatment_effect, target_standard_error) {
  weight_source <- target_standard_error ^ 2 /
    (source_data$standard_error ^ 2 + target_standard_error ^ 2)

  weight_source * source_data$treatment_effect_estimate +
    (1 - weight_source) * target_treatment_effect
}

# BSDA::tsum.test() emits "argument 'var.equal' ignored for one-sample test."
# on every one-sample call, whether or not var.equal was supplied. Muffle just
# that message so a genuine warning would still surface.
without_bsda_one_sample_note <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("var.equal", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

run_pooling_power <- function(target_data, alpha) {
  compute_freq_power_pooling(
    alpha = alpha,
    target_data = target_data,
    source_data = source_data,
    frequentist_test = "t-test",
    theta_0 = 0,
    null_space = "right",
    case_study = "teriflunomide",
    simulation_config = list(seed = 1),
    n_replicates = 5
  )
}

test_that("the Monte Carlo branch tests the pooled effect against the pooled standard error", {
  target_treatment_effect <- -0.4155154439616658
  target_data <- make_target_data(target_treatment_effect)

  # The pooled test is decisive here (t = -8.1), while the same pooled estimate
  # compared against the target-only standard error is not (t = -1.5, p = 0.06).
  target_only_p_value <- stats::pt(
    pooled_estimate(target_treatment_effect, target_standard_error) / target_standard_error,
    df = target_sample_size_per_arm - 1
  )
  expect_gt(target_only_p_value, 0.025)

  result <- without_bsda_one_sample_note(run_pooling_power(target_data, alpha = 0.025))

  expect_equal(result$power, 1)
})

test_that("the Monte Carlo branch uses the pooled sample size for the test", {
  # A balanced source/target pair placing the pooled test near the boundary, so
  # that the reported decision pins both the pooled standard error and the
  # pooled degrees of freedom.
  balanced_source_standard_error <- 0.2
  balanced_target_standard_error <- 0.25
  balanced_source_data <- list(
    treatment_effect_estimate = -0.26,
    standard_error = balanced_source_standard_error,
    equivalent_source_sample_size_per_arm = 120
  )
  target_data <- make_target_data(-0.25, standard_error = balanced_target_standard_error)

  weight_source <- balanced_target_standard_error ^ 2 /
    (balanced_source_standard_error ^ 2 + balanced_target_standard_error ^ 2)
  estimate <- weight_source * balanced_source_data$treatment_effect_estimate +
    (1 - weight_source) * (-0.25)
  standard_error <- sqrt(
    1 / (1 / balanced_source_standard_error ^ 2 + 1 / balanced_target_standard_error ^ 2)
  )
  pooled_sample_size <- target_sample_size_per_arm +
    balanced_source_data$equivalent_source_sample_size_per_arm
  expected_p_value <- stats::pt(estimate / standard_error, df = pooled_sample_size - 1)

  run_balanced <- function(alpha) {
    compute_freq_power_pooling(
      alpha = alpha,
      target_data = target_data,
      source_data = balanced_source_data,
      frequentist_test = "t-test",
      theta_0 = 0,
      null_space = "right",
      case_study = "teriflunomide",
      simulation_config = list(seed = 1),
      n_replicates = 5
    )
  }

  expect_equal(without_bsda_one_sample_note(run_balanced(expected_p_value + 1e-6))$power, 1)
  expect_equal(without_bsda_one_sample_note(run_balanced(expected_p_value - 1e-6))$power, 0)
})
