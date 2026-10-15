# Three fixes from a review of the Bayesian and pooled operating
# characteristics:
# - under a binomial likelihood the design prior is taken given the control
#   rate the trials are simulated at, so it puts no mass on effects those trials
#   cannot have;
# - the analytic power of the pooled analysis holds the source data fixed, as
#   the simulation does, rather than treating them as resampled;
# - the Monte Carlo Bayesian operating characteristics average the decisions
#   of the replicates actually analysed.

review_case_study <- function(name) {
  yaml::read_yaml(system.file("conf", "case_studies", paste0(name, ".yml"), package = "BExTE"))
}

review_mcmc_config <- function() {
  yaml::read_yaml(system.file("conf", "full", "mcmc_config.yml", package = "BExTEsim"))
}

review_model <- function(method, parameters, config = review_case_study("aprepitant")) {
  Model$new()$create(
    case_study_config = config,
    method = method,
    method_parameters = lapply(parameters, list),
    source_data = SourceData$new(config, NA),
    mcmc_config = review_mcmc_config()
  )
}

review_design_prior <- function(type, model, config, control_rate) {
  DesignPrior$new()$create(
    design_prior_type = type,
    model = model,
    source_data = SourceData$new(config, NA),
    case_study_config = config,
    simulation_config = list(n_samples_mixture_approx = 500),
    mcmc_config = review_mcmc_config(),
    case_study = config$name,
    target_control_rate = control_rate
  )
}


test_that("the analytic pooled power holds the source estimate fixed", {
  # Equal standard errors, a null source estimate and a null target effect: a
  # one-sided z-test at 5% rejects when the target estimate exceeds
  # z sqrt(2) standard errors, which happens about 1% of the time.
  config <- review_case_study("botox")
  source_data <- SourceData$new(config, NA)
  source_data$treatment_effect_estimate <- config$theta_0
  target_data <- TargetDataFactory$new()$create(
    source_data = source_data, case_study_config = config,
    target_sample_size_control = as.integer(source_data$equivalent_source_sample_size_per_arm),
    target_sample_size_treatment = as.integer(source_data$equivalent_source_sample_size_per_arm),
    control_drift = 0, treatment_drift = config$theta_0 - source_data$treatment_effect_estimate,
    summary_measure_likelihood = "normal", target_to_source_std_ratio = 1
  )
  target_se <- target_data$standard_deviation / sqrt(target_effective_arm_size(target_data))
  source_data$standard_error <- target_se

  power <- compute_freq_power_pooling(
    alpha = 0.05, target_data = target_data, source_data = source_data,
    frequentist_test = "z-test", theta_0 = config$theta_0,
    null_space = config$null_space, simulation_config = list()
  )$power
  expect_equal(power, stats::pnorm(stats::qnorm(0.95) * sqrt(2), lower.tail = FALSE), tolerance = 1e-10)
})


test_that("the analytic pooled power agrees with simulating the target estimate", {
  config <- review_case_study("botox")
  source_data <- SourceData$new(config, NA)
  target_data <- TargetDataFactory$new()$create(
    source_data = source_data, case_study_config = config,
    target_sample_size_control = 58L,
    target_sample_size_treatment = 58L,
    control_drift = 0, treatment_drift = -source_data$treatment_effect_estimate / 2,
    summary_measure_likelihood = "normal", target_to_source_std_ratio = 1
  )
  alternative <- if (config$null_space == "left") "greater" else "less"

  power <- compute_freq_power_pooling(
    alpha = 0.025, target_data = target_data, source_data = source_data,
    frequentist_test = "z-test", theta_0 = config$theta_0,
    null_space = config$null_space, simulation_config = list()
  )$power

  target_se <- target_data$standard_deviation / sqrt(target_effective_arm_size(target_data))
  set.seed(3)
  estimates <- stats::rnorm(200000, target_data$treatment_effect, target_se)
  pooled_se <- sqrt(1 / (1 / source_data$standard_error^2 + 1 / target_se^2))
  pooled <- (source_data$treatment_effect_estimate / source_data$standard_error^2 +
               estimates / target_se^2) * pooled_se^2
  statistic <- (pooled - config$theta_0) / pooled_se
  simulated <- if (alternative == "greater") {
    mean(statistic > stats::qnorm(0.975))
  } else {
    mean(statistic < -stats::qnorm(0.975))
  }
  expect_equal(power, simulated, tolerance = 0.005)
})


test_that("the binomial pooled power sums over the target counts only", {
  config <- review_case_study("aprepitant")
  source_data <- SourceData$new(config, NA)
  target_data <- TargetDataFactory$new()$create(
    source_data = source_data, case_study_config = config,
    target_sample_size_control = 71L,
    target_sample_size_treatment = 71L, control_drift = 0,
    treatment_drift = -source_data$treatment_effect_estimate / 2,
    summary_measure_likelihood = "binomial", target_to_source_std_ratio = 1
  )
  power <- pooled_binomial_power(0.025, target_data, source_data, "greater")

  set.seed(4)
  n <- 71L
  y_t <- stats::rbinom(200000, n, target_data$treatment_rate)
  y_c <- stats::rbinom(200000, n, target_data$control_rate)
  s_t <- counts_from_rate(source_data$treatment_rate, source_data$sample_size_treatment)
  s_c <- counts_from_rate(source_data$control_rate, source_data$sample_size_control)
  n_t <- n + source_data$sample_size_treatment
  n_c <- n + source_data$sample_size_control
  statistic <- (2 * asin(sqrt((s_t + y_t) / n_t)) - 2 * asin(sqrt((s_c + y_c) / n_c))) /
    sqrt(1 / n_t + 1 / n_c)
  expect_equal(power, mean(statistic > stats::qnorm(0.975)), tolerance = 0.005)
})
