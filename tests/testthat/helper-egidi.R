# Shared fixtures for the Egidi empirical mixture prior. The method reuses the
# robust mixture prior's components, so the fixtures build both from one case
# study configuration: a test that compares the two has to be sure they differ
# only in where the weight came from, and that is only true if the source data,
# the target data and the replicates are literally the same objects.

egidi_case_study_config <- function(summary_measure_likelihood = "normal") {
  list(
    name = "unit_test",
    endpoint = if (summary_measure_likelihood == "normal") "continuous" else "binary",
    summary_measure_likelihood = summary_measure_likelihood,
    sampling_approximation = TRUE,
    theta_0 = 0,
    null_space = "left",
    source = list(
      control = 100,
      treatment = 100,
      treatment_effect = 0.5,
      standard_error = 0.12
    )
  )
}

egidi_method_parameters <- function(alpha_pc = 0.05,
                                    weight_grid_step = 0.001,
                                    pvalue_method = "exact") {
  list(
    initial_prior = list("noninformative"),
    alpha_pc = list(alpha_pc),
    pvalue_method = list(pvalue_method),
    weight_grid_step = list(weight_grid_step),
    empirical_bayes = list(TRUE)
  )
}

#' The model, its target data and one draw of replicates.
#'
#' `method` selects between the empirical mixture prior and the robust mixture
#' prior it is built on, so both can be run over the identical replicates.
egidi_fixture <- function(method = "egidi_empirical_mixture",
                          method_parameters = egidi_method_parameters(),
                          treatment_drift = 0,
                          n_replicates = 200,
                          target_sample_size_per_arm = 60,
                          seed = 20260917) {
  case_study_config <- egidi_case_study_config()
  source_data <- SourceData$new(case_study_config)
  target_data <- TargetDataFactory$new()$create(
    source_data = source_data,
    case_study_config = case_study_config,
    target_sample_size_control = target_sample_size_per_arm,
    target_sample_size_treatment = target_sample_size_per_arm,
    control_drift = 0,
    treatment_drift = treatment_drift,
    summary_measure_likelihood = "normal",
    target_to_source_std_ratio = 1
  )
  model <- Model$new()$create(
    case_study_config = case_study_config,
    method = method,
    method_parameters = method_parameters,
    source_data = source_data
  )
  set.seed(seed)
  list(
    model = model,
    target_data = target_data,
    samples = target_data$generate(n_replicates)
  )
}

#' The conflict p-value by dense numerical integration, with no root-finding.
#'
#' Independent of the implementation, which splits the line at the turning points
#' and integrates between the crossings in closed form; this one evaluates the
#' indicator everywhere and sums.
#'
#' The nodes are the union of the two components' own standardised grids rather
#' than one uniform grid. A uniform grid has to span the wider component while
#' resolving the narrower, so when the two scales differ by a factor of ten its
#' remaining error - a step-function boundary term of order the spacing times the
#' peak density - is around 1e-4, which is larger than anything worth testing
#' for. The union grid keeps that term near 1e-6 at the same cost.
egidi_reference_conflict_pvalue <- function(t_obs, psi, mu_p, sigma_p,
                                            mu_q, sigma_q, n_nodes = 1e6) {
  density <- function(t) {
    (1 - psi) * stats::dnorm(t, mu_p, sigma_p) + psi * stats::dnorm(t, mu_q, sigma_q)
  }
  level <- density(t_obs)
  if (level <= 0) {
    return(0)
  }
  deviations <- seq(-42, 42, length.out = n_nodes)
  grid <- sort(c(mu_p + sigma_p * deviations, mu_q + sigma_q * deviations))
  values <- density(grid)
  widths <- c(diff(grid), 0)
  sum(ifelse(values <= level, values, 0) * widths)
}

#' The conflict p-value by prior-predictive simulation.
egidi_reference_monte_carlo <- function(t_obs, psi, mu_p, sigma_p, mu_q, sigma_q,
                                        draws = 2e6, seed = 1) {
  withr::local_seed(seed)
  from_weak <- stats::runif(draws) < psi
  simulated <- ifelse(
    from_weak,
    stats::rnorm(draws, mu_q, sigma_q),
    stats::rnorm(draws, mu_p, sigma_p)
  )
  density <- function(t) {
    (1 - psi) * stats::dnorm(t, mu_p, sigma_p) + psi * stats::dnorm(t, mu_q, sigma_q)
  }
  estimate <- mean(density(simulated) <= density(t_obs))
  list(pvalue = estimate, standard_error = sqrt(estimate * (1 - estimate) / draws))
}
