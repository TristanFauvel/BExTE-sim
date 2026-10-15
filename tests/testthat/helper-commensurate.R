# Shared fixtures for the commensurate power prior. These models compile their
# Stan program only when they first sample, but the compiler is still stubbed
# here so that no test can reach cmdstan by accident; keeping that in one place
# stops the stub and the configuration from drifting between the Stan-program
# tests and the fast-path tests.

commensurate_mcmc_config <- function() {
  list(
    num_chains = 1L,
    parallel_chains = 1L,
    tune = 1L,
    target_accept = 0.8,
    chain_length = 1L,
    max_chain_length = 2L,
    target_ess = 1L,
    rhat_threshold = 1.1,
    max_divergence_rate = 0.01
  )
}

commensurate_model <- function(prior, mcmc_config = commensurate_mcmc_config(),
                               generator = GaussianCommensuratePowerPrior) {
  model <- testthat::with_mocked_bindings(
    generator$new(prior = prior, mcmc_config = mcmc_config),
    compile_stan_model = function(...) NULL,
    .package = "BExTE"
  )
  # Model$create() normally installs the prior after constructing the subclass.
  model$prior <- prior
  model
}

# The same, for the plain commensurate prior - the gamma == 1 case.
commensurate_prior_model <- function(prior,
                                     mcmc_config = commensurate_mcmc_config()) {
  commensurate_model(prior, mcmc_config, generator = GaussianCommensuratePrior)
}

# Enough of a prior to read the Stan program off the constructed model.
commensurate_stan_model <- function() {
  commensurate_model(list(
    source = list(),
    method_parameters = list(
      initial_prior = list("noninformative"),
      heterogeneity_prior = list(family = "half_normal", std_dev = 1)
    )
  ))
}

# Enough of a prior to read the gamma-free Stan program off the model.
commensurate_prior_stan_model <- function() {
  commensurate_prior_model(list(
    source = list(),
    method_parameters = list(
      initial_prior = list("noninformative"),
      heterogeneity_prior = list(family = "half_normal", std_dev = 1)
    )
  ))
}

# A fully specified prior, so the quadrature mixture can be built from it.
commensurate_fast_path_prior <- function(heterogeneity_prior) {
  list(
    source = list(
      standard_error = 0.12,
      equivalent_source_sample_size_per_arm = 50,
      treatment_effect_estimate = 0.5
    ),
    method_parameters = list(
      initial_prior = list("noninformative"),
      heterogeneity_prior = heterogeneity_prior
    )
  )
}

commensurate_fast_path_model <- function(heterogeneity_prior) {
  commensurate_model(commensurate_fast_path_prior(heterogeneity_prior))
}

# The same source and target configuration, under the gamma == 1 model, so the
# two can be compared component by component.
commensurate_prior_fast_path_model <- function(heterogeneity_prior) {
  commensurate_prior_model(commensurate_fast_path_prior(heterogeneity_prior))
}

# The heterogeneity priors the shipped configurations run, as listed in
# inst/conf/full/methods_config.R; the last is the Cauchy(0, 30) on log(tau) of
# Hobbs et al. (2011). The Cauchy(0, 10) is kept as a second scale, because the
# quadrature sizes its log-Cauchy rule by the scale.
commensurate_configured_priors <- function() {
  list(
    list(family = "inverse_gamma", alpha = 1 / 3, beta = 1),
    list(family = "inverse_gamma", alpha = 1 / 7, beta = 1),
    list(family = "half_normal", std_dev = 1),
    list(family = "half_normal", std_dev = 5),
    list(family = "cauchy", location = 0, scale = 10),
    list(family = "cauchy", location = 0, scale = 30)
  )
}

# Posterior mean, standard deviation and quantiles of the treatment effect
# under a mixture built at a given resolution.
commensurate_posterior_summary <- function(model, n_tau, n_gamma,
                                           estimate = 0.9,
                                           standard_error = 0.1) {
  prior <- commensurate_prior_mixture(model, n_tau, n_gamma)
  posterior <- normal_mixture_posterior(
    weights = prior$weights,
    means = prior$means,
    sds = prior$sds,
    estimate = estimate,
    standard_error = standard_error
  )
  summary <- normal_mixture_summary(
    posterior$weights,
    posterior$means,
    posterior$sds
  )

  c(summary$mean, summary$sd, summary$quantiles)
}

# Three replicates of target data, and the fast path run over them. Both
# commensurate models take the same inputs here, which is what lets the gamma
# == 1 model be compared against the power prior component by component.
commensurate_replicate_samples <- function() {
  data.frame(
    treatment_effect_estimate = c(0.2, 0.5, 0.9),
    treatment_effect_standard_error = rep(0.1, 3),
    standard_deviation = rep(0.1 * sqrt(60), 3)
  )
}

commensurate_fast_path_run <- function(model, to_return,
                                       critical_value = 0.975,
                                       confidence_level = 0.95,
                                       samples = commensurate_replicate_samples()) {
  model$vectorised_replicate_inference(
    target_data = list(sample_size_control = 60, sample_size_treatment = 60),
    samples = samples,
    to_return = to_return,
    critical_value = critical_value,
    theta_0 = 0,
    confidence_level = confidence_level,
    null_space = "left"
  )
}

# The case study both Stan-equivalence tests fit, so the commensurate power
# prior and its gamma == 1 case are compared against Stan on the same data.
commensurate_equivalence_config <- function() {
  list(
    name = "unit_test",
    endpoint = "continuous",
    summary_measure_likelihood = "normal",
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

# Long chains, and thresholds loose enough that the replicate loop never
# restarts an iteration: the point here is to compare the two paths on the same
# data, not to exercise the chain-length adaptation.
commensurate_equivalence_mcmc_config <- function(chain_length = 4000L) {
  list(
    num_chains = 4L,
    parallel_chains = 4L,
    tune = 1000L,
    target_accept = 0.9,
    chain_length = chain_length,
    max_chain_length = chain_length,
    target_ess = 1L,
    rhat_threshold = 1.5,
    max_divergence_rate = 1,
    # inference() would otherwise use quadrature, and compare it with itself.
    engine = "stan"
  )
}
