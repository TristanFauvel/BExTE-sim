#' Read a YAML configuration and validate it against a schema
#'
#' @param path Path to the YAML file.
#' @param schema The schema to validate against, as taken by
#'   [validate_config()].
#'
#' @return The parsed configuration.
#'
#' @keywords internal
read_config <- function(path, schema) {
  config <- yaml::read_yaml(path)
  validate_config(config, schema, path)
  return(config)
}


#' Schema for `simulation_config.yml`
#'
#' @keywords internal
simulation_config_schema <- list(
  critical_value = "probability",
  confidence_level = "probability",
  seed = "number",
  n_samples_design_prior = "count",
  n_samples_mixture_approx = "count",
  n_samples_quantiles_estimation = "count",
  compute_bayesian_ocs_mc = "flag",
  compute_frequentist_ocs = "flag",
  delete_old_results = "flag"
)


#' Check that the decision threshold matches the credible interval level
#'
#' @description Models sampled by MCMC reject the null when the central
#'   credible interval at `confidence_level` excludes it, and stop on any
#'   other threshold. The other models reject it when the posterior
#'   probability of benefit exceeds `critical_value`, and report the interval
#'   separately. The two rules agree only when `critical_value` is
#'   `(1 + confidence_level) / 2`. Any other pair would make some methods
#'   report intervals that contradict their decisions, and would stop the
#'   MCMC methods partway through a run, so it is refused up front.
#'
#' @param config The simulation configuration.
#' @param context The name of the configuration, used in the error.
#'
#' @return No return value, called for side effects.
#'
#' @keywords internal
check_decision_threshold <- function(config, context) {
  implied_critical_value <- (1 + config$confidence_level) / 2
  if (!isTRUE(all.equal(config$critical_value, implied_critical_value))) {
    stop(paste0(
      context, " is not valid: critical_value should be ",
      "(1 + confidence_level) / 2 = ", format(implied_critical_value),
      " so that the test decision agrees with the credible interval, but is ",
      format(config$critical_value)
    ), call. = FALSE)
  }

  return(invisible(NULL))
}


#' Read `simulation_config.yml` and validate it
#'
#' @description Checks the configuration against [simulation_config_schema],
#'   then with [check_decision_threshold()].
#'
#' @param path Path to the YAML file.
#'
#' @return The parsed configuration.
#'
#' @keywords internal
read_simulation_config <- function(path) {
  config <- read_config(path, simulation_config_schema)
  check_decision_threshold(config, path)
  return(config)
}


#' Schema for `mcmc_config.yml`
#'
#' @keywords internal
mcmc_config_schema <- list(
  num_chains = "count",
  parallel_chains = "count",
  tune = "count",
  target_accept = "probability",
  chain_length = "count",
  max_chain_length = "count?",
  target_ess = "count",
  rhat_threshold = "number",
  max_divergence_rate = "probability",
  # Optional. "quadrature", the default, or "stan": how the models that can do
  # either compute their posterior - see MCMCModel.
  engine = "string?"
)


#' Schema for `scenarios_config.yml`
#'
#' @keywords internal
scenarios_config_schema <- list(
  n_replicates = "count",
  ndrift = "nonnegative_count",
  parallelization = "flag_or_list?",
  denominator_change_factor = "numeric_vector",
  sample_size_factors = "numeric_vector",
  target_to_source_std_ratio_range = "numeric_vector?",
  # Optional, and only simulated for the time-to-event endpoint. They default
  # to the primary design, so a config that omits them keeps its current grid -
  # see compute_control_drift_range() and compute_time_to_event_ranges().
  control_drift_range = "numeric_vector?",
  dropout_probability = "numeric_vector?",
  event_time_distribution = "character_vector?",
  treatment_delay = "numeric_vector?",
  # Optional. Restricts the sensitivity designs above to one sample size factor
  # and one source denominator change factor, instead of crossing them with the
  # whole grid - see time_to_event_sensitivity_reference().
  sensitivity_reference = "named_numeric_list?",
  # Optional. When TRUE, a sensitivity design moves one of the axes above off
  # its primary value at a time, rather than every combination of them - see
  # time_to_event_axes_off_primary().
  sensitivity_one_at_a_time = "flag?",
  case_studies = "character_vector",
  methods = "character_vector",
  # Optional. The analysis steps to run after the simulation, out of
  # ANALYSIS_STEPS. Absent means all of them, so a config written before
  # this key keeps its behaviour - see run_simulation_env().
  analysis_steps = "character_vector?",
  # Optional. Restricts the sample size factors simulated for a given
  # case study, so a run need not take the cross product of every case
  # study with every factor - see simulation_frequentist_ocs(). Case
  # studies absent from it keep the full sample_size_factors.
  case_study_sample_size_factors = "named_numeric_list?",
  # Optional. Ratio of the treatment to the control arm size of the target
  # study, which sample_size_factors sizes in total. Absent means 1, equal
  # arms - see target_allocation_ratio().
  target_allocation_ratio = "number?",
  # Optional. Overrides target_allocation_ratio for a given case study.
  case_study_target_allocation_ratio = "named_numeric_list?",
  # Optional. Overrides n_replicates for a given case study - see
  # case_study_n_replicates(). Case studies absent from it keep n_replicates.
  case_study_n_replicates = "named_numeric_list?",
  # Optional. The case studies whose operating characteristics are computed
  # exactly, by enumerating the trial outcomes, instead of by simulating
  # n_replicates trials - see case_study_exact_enumeration().
  exact_enumeration = "character_vector?"
)
