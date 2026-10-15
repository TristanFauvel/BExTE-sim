oc_estimates <- list(
  average_tie = 0.05,
  average_power = 0.8,
  prior_proba_no_benefit = 0.4,
  prepost_proba_FP = 0.02,
  prepost_proba_TP = 0.5,
  upper_bound_proba_FP = 0.03,
  prior_proba_success = 0.6
)

test_source_data <- list(
  summary_measure_likelihood = "normal",
  standard_error = 0.2,
  equivalent_source_sample_size_per_arm = 25,
  treatment_effect_estimate = 0.5,
  to_dict = function() {
    list(
      source_treatment_effect_estimate = 0.5,
      source_standard_error = 0.2,
      endpoint = "continuous",
      summary_measure_likelihood = "normal",
      source_sample_size_control = 25,
      source_sample_size_treatment = 25,
      equivalent_source_sample_size_per_arm = 25,
      source_control_rate = NA,
      source_treatment_rate = NA
    )
  }
)

test_scenario <- data.frame(
  case_study = "unit_test",
  method = "EB_PP",
  target_sample_size_control = 50,
  target_sample_size_treatment = 50,
  source_denominator = 1,
  target_to_source_std_ratio = 1,
  parameters = "{}",
  stringsAsFactors = FALSE
)

test_simulation_config <- list(
  confidence_level = 0.95,
  n_samples_design_prior = 10,
  n_samples_quantiles_estimation = 100,
  n_samples_mixture_approx = 100
)

run_estimate_bayesian_ocs <- function(model, design_prior_type) {
  BExTEsim:::estimate_bayesian_ocs(
    scenario = test_scenario,
    case_study_config = list(name = "unit_test", null_space = "left"),
    source_data = test_source_data,
    target_sample_size_control = 50,
    target_sample_size_treatment = 50,
    model = model,
    method_parameters = list(),
    theta_0 = 0,
    n_replicates = 10,
    critical_value = 0.975,
    computation_state = list(rng_state = "seed", computation_time = 0),
    design_prior_type = design_prior_type,
    json_parameters = "{}",
    simulation_config = test_simulation_config,
    target_to_source_std_ratio = 1,
    mcmc_config = NULL
  )
}

empirical_bayes_model <- list(
  empirical_bayes = TRUE,
  estimate_bayesian_operating_characteristics = function(...) {
    stop("The Monte Carlo estimator must not run for an empirical Bayes analysis prior.")
  }
)

test_that("an empirical Bayes analysis prior yields missing OCs instead of aborting", {
  result <- run_estimate_bayesian_ocs(empirical_bayes_model, "analysis_prior")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(result$design_prior, "analysis_prior")
  expect_true(all(is.na(result[, BExTEsim:::bayesian_ocs_metric_names()])))
})

test_that("the skipped analysis prior row matches the layout of a computed row", {
  computed_model <- list(
    empirical_bayes = TRUE,
    estimate_bayesian_operating_characteristics = function(...) oc_estimates
  )

  skipped <- run_estimate_bayesian_ocs(empirical_bayes_model, "analysis_prior")
  computed <- run_estimate_bayesian_ocs(computed_model, "ui_design_prior")

  expect_identical(names(skipped), names(computed))
  expect_equal(nrow(rbind(computed, skipped)), 2)
})

test_that("a non-empirical-Bayes model still computes its analysis prior OCs", {
  classical_model <- list(
    empirical_bayes = FALSE,
    prior_to_RBesT = function(n_samples) invisible(NULL),
    estimate_bayesian_operating_characteristics = function(...) oc_estimates
  )

  result <- run_estimate_bayesian_ocs(classical_model, "analysis_prior")

  expect_equal(result$design_prior, "analysis_prior")
  expect_equal(result$average_power, 0.8)
})
