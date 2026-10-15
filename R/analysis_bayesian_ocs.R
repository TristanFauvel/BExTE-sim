# Preposterior probability of success over one region of the treatment-effect
# grid. The design prior mass beyond the grid edge is added by
# design_prior_tail_contribution(). Only the null region warns when it holds
# fewer than two grid points.
preposterior_proba_region <- function(conditional_proba_success,
                                      treatment_effect_values,
                                      theta_0,
                                      null_space,
                                      design_prior_pdf,
                                      design_prior_cdf,
                                      region = c("null", "alternative")) {
  region <- match.arg(region)

  input_lengths <- c(
    length(conditional_proba_success),
    length(treatment_effect_values),
    length(design_prior_pdf)
  )
  if (length(unique(input_lengths)) != 1L) {
    stop("Probability, treatment-effect, and prior-density vectors must have the same length.")
  }

  values_in_region <- space_mask(treatment_effect_values, theta_0, null_space, region)

  # The null region of a left null space extends to the left, and so on.
  outward_left <- (null_space == "left") == (region == "null")
  outward_side <- if (outward_left) "lower" else "upper"
  limit <- if (region == "null") 0 else 1

  if (sum(values_in_region) < 2) {
    if (region == "null") {
      warning(
        "The number of values in the null space is less than 2, the preposterior probability of FP cannot be computed."
      )
    }
    return(NA)
  }

  region_treatment_effects <- treatment_effect_values[values_in_region]
  region_proba_success <- conditional_proba_success[values_in_region]

  p <- Bolstad2::sintegral(
    x = region_treatment_effects,
    fx = design_prior_pdf[values_in_region] * region_proba_success,
    n.pts = sum(values_in_region)
  )$int

  edge <- if (outward_side == "lower") {
    which.min(region_treatment_effects)
  } else {
    which.max(region_treatment_effects)
  }

  p <- p + design_prior_tail_contribution(
    design_prior_cdf = design_prior_cdf,
    edge_value = region_treatment_effects[edge],
    side = outward_side,
    limit = limit,
    edge_proba_success = region_proba_success[edge],
    space_label = region
  )

  p <- sapply(p, check_probability_value)
  assert_single_number(p)
  return(p)
}

# Preposterior probability of a False Positive
preposterior_proba_FP <- function(conditional_proba_success,
                                  treatment_effect_values,
                                  theta_0,
                                  null_space,
                                  design_prior_pdf,
                                  design_prior_cdf) {
  return(preposterior_proba_region(
    conditional_proba_success, treatment_effect_values, theta_0, null_space,
    design_prior_pdf, design_prior_cdf, region = "null"
  ))
}

# Preposterior probability of a True Positive
preposterior_proba_TP <- function(conditional_proba_success,
                                  treatment_effect_values,
                                  theta_0,
                                  null_space,
                                  design_prior_pdf,
                                  design_prior_cdf) {
  return(preposterior_proba_region(
    conditional_proba_success, treatment_effect_values, theta_0, null_space,
    design_prior_pdf, design_prior_cdf, region = "alternative"
  ))
}

#' Bayesian operating characteristics of one method setting on one design
#'
#' @description Rebuilds the model of one method and parameter combination for
#'   one trial design, and integrates its probability of success curve against
#'   each design prior. [compute_bayesian_ocs()] collects one job per
#'   combination and hands them to this function, sequentially or on a
#'   cluster; the jobs share nothing, so they can run in any order.
#'
#' @param job A list with the combination's results rows (`results_df`), its
#'   `parameters_combination` row, and the case study, method, design and
#'   source data they belong to.
#' @param simulation_config The simulation configuration.
#' @param mcmc_config The MCMC configuration.
#' @param design_prior_types The design priors to integrate against.
#'
#' @return A data frame with one row per design prior.
#' @noRd
bayesian_ocs_for_combination <- function(job,
                                         simulation_config,
                                         mcmc_config,
                                         design_prior_types) {
  case_study <- job$case_study
  case_study_config <- job$case_study_config
  null_space <- case_study_config$null_space
  theta_0 <- case_study_config$theta_0
  source_data <- job$source_data
  source_data_df <- job$source_data_df
  source_denominator_change_factor <- job$source_denominator_change_factor
  target_to_source_std_ratio <- job$target_to_source_std_ratio
  method <- job$method
  target_sample_size_control <- job$target_sample_size_control
  target_sample_size_treatment <- job$target_sample_size_treatment
  parameters_combination <- job$parameters_combination
  results_df <- job$results_df

  # The mixture approximations of the priors are fitted to random draws. Each
  # job is seeded on its own, so its result does not depend on which jobs ran
  # before it, or on whether they ran on a cluster.
  withr::local_seed(simulation_config$seed)

  rows <- list()

  # We unpack the parameter inside this loop (and not inside the previous one), because for some methods such as the commensurate power prior, there is a nested parameters structure which implies that they cannot all be stored in a single dataframe.
  method_parameters <- as.list(get_parameters(
    parameters_combination[, "parameters", drop = FALSE]
  ))
  # convert strings to numeric or boolean if possible
  method_parameters <- data.frame(lapply(method_parameters, convert_if_possible))

  # The following applies to a single row dataframe, to recover the nested list structure
  method_params <- extract_nested_parameter(method_parameters)

  conditional_proba_success <- results_df$success_proba
  treatment_effect_values <- results_df$target_treatment_effect

  model <- Model$new()
  model <- model$create(
    case_study_config = case_study_config,
    method = method,
    method_parameters = method_params,
    source_data = source_data,
    mcmc_config = mcmc_config
  )

  # A method whose hyperparameters follow from the design has to be
  # given one here as well: this model is rebuilt from the results
  # file rather than handed over by the simulation, so it arrives
  # uncalibrated.
  #
  # The design used is the nominal one - the trial this row
  # describes, at zero treatment drift - rather than any particular
  # drift, because the model stands for the whole curve here: it
  # supplies the design prior, which is a property of the method and
  # the design and not of the treatment effect the curve is indexed
  # by. On a continuous endpoint that is the same design the
  # simulation calibrated against at every drift, since the standard
  # error the design implies does not move with the drift there.
  first_row <- results_df[1, , drop = FALSE]
  column_or <- function(name, fallback) {
    if (name %in% names(first_row) && !is.na(first_row[[name]])) {
      return(first_row[[name]])
    } else {
      return(fallback)
    }
  }
  nominal_target_data <- TargetDataFactory$new()$create(
    source_data = source_data,
    case_study_config = case_study_config,
    target_sample_size_control = target_sample_size_control,
    target_sample_size_treatment = target_sample_size_treatment,
    control_drift = column_or("control_drift", 0),
    treatment_drift = 0,
    summary_measure_likelihood =
      case_study_config$summary_measure_likelihood,
    target_to_source_std_ratio = target_to_source_std_ratio,
    dropout_probability = column_or("dropout_probability", 0),
    event_time_distribution =
      column_or("event_time_distribution", "exponential"),
    treatment_delay = column_or("treatment_delay", 0)
  )
  model$calibrate_for_design(nominal_target_data)

  for (design_prior_type in design_prior_types) {
    if (design_prior_type == "analysis_prior" &&
        model$empirical_bayes) {
      # We cannot compute the Bayesian OCs with an analysis design prior for empirical Bayes methods.
      results_to_add <- data.frame(
        case_study = case_study,
        method = method,
        target_sample_size_control = target_sample_size_control,
        target_sample_size_treatment = target_sample_size_treatment,
        parameters = parameters_combination[1, ],
        source_denominator_change_factor = source_denominator_change_factor,
        target_to_source_std_ratio = target_to_source_std_ratio,
        design_prior_type = design_prior_type,
        prior_proba_success = NA,
        prior_proba_no_benefit = NA,
        prior_proba_benefit = NA,
        prepost_proba_FP = NA,
        prepost_proba_TP = NA,
        average_tie = NA,
        average_power = NA,
        upper_bound_proba_FP = NA
      )
    } else {
      design_prior <- DesignPrior$new()
      design_prior <- design_prior$create(
        design_prior_type = design_prior_type,
        model = model,
        source_data = source_data,
        case_study_config = case_study_config,
        simulation_config = simulation_config,
        mcmc_config = mcmc_config,
        case_study = case_study,
        # The curve is simulated at this control rate, whatever the treatment
        # effect, so under a binomial likelihood the design prior is taken
        # given it - see DesignPrior$create().
        target_control_rate = nominal_target_data$control_rate
      )

      design_prior_pdf <- design_prior$pdf(treatment_effect_values)

      if (null_space == "left") {
        prior_proba_no_benefit <- design_prior$cdf(theta_0)
      } else if (null_space == "right") {
        prior_proba_no_benefit <- 1 - design_prior$cdf(theta_0)
      } else {
        stop("Null space must be either 'left' or 'right'.")
      }

      prior_proba_benefit <- 1 - prior_proba_no_benefit

      prior_proba_success_SI <- prior_proba_success(
        conditional_proba_success = conditional_proba_success,
        treatment_effect_values = treatment_effect_values,
        design_prior_pdf = design_prior_pdf,
        design_prior_cdf = design_prior$cdf,
        null_space = null_space
      )

      prepost_proba_FP_SI <- preposterior_proba_FP(
        conditional_proba_success = conditional_proba_success,
        treatment_effect_values = treatment_effect_values,
        theta_0 = theta_0,
        null_space = null_space,
        design_prior_pdf = design_prior_pdf,
        design_prior_cdf = design_prior$cdf
      )

      prepost_proba_TP_SI <- preposterior_proba_TP(
        conditional_proba_success = conditional_proba_success,
        treatment_effect_values = treatment_effect_values,
        theta_0 = theta_0,
        null_space = null_space,
        design_prior_pdf = design_prior_pdf,
        design_prior_cdf = design_prior$cdf
      )

      average_tie_SI <- average_tie(
        prepost_proba_FP = prepost_proba_FP_SI,
        prior_proba_no_benefit = prior_proba_no_benefit
      )

      average_power_SI <- average_power(
        prepost_proba_TP = prepost_proba_TP_SI,
        prior_proba_no_benefit = prior_proba_no_benefit
      )

      if(length(conditional_proba_success[treatment_effect_values == theta_0]) > 1){
        stop("Too many values corresponding to TIE.")
      }
      upper_bound_proba_FP_SI <- upper_bound_proba_FP(
        prior_proba_no_benefit = prior_proba_no_benefit,
        treatment_effect_values = treatment_effect_values,
        conditional_proba_success = conditional_proba_success,
        theta_0 = theta_0
      )

      results_to_add <- data.frame(
        case_study = case_study,
        method = method,
        target_sample_size_control = target_sample_size_control,
        target_sample_size_treatment = target_sample_size_treatment,
        parameters = parameters_combination[1, ],
        source_denominator_change_factor = source_denominator_change_factor,
        target_to_source_std_ratio = target_to_source_std_ratio,
        design_prior_type = design_prior_type,
        prior_proba_success = prior_proba_success_SI,
        prior_proba_no_benefit = prior_proba_no_benefit,
        prior_proba_benefit = prior_proba_benefit,
        prepost_proba_FP = prepost_proba_FP_SI,
        prepost_proba_TP = prepost_proba_TP_SI,
        average_tie = average_tie_SI,
        average_power = average_power_SI,
        upper_bound_proba_FP = upper_bound_proba_FP_SI
      )
    }
    rows[[design_prior_type]] <- cbind(results_to_add, source_data_df)
  }

  return(dplyr::bind_rows(rows))
}


compute_bayesian_ocs <- function(results_freq_df, env, config_dir = NULL, case_studies_config_dir = NULL, parallelization = FALSE, cluster = NULL) {
  jobs <- list()

  if (is.null(config_dir)) {
    config_dir <- paste0(system.file(paste0("conf/", env), package = "BExTEsim"), "/")
  }
  if (is.null(case_studies_config_dir)) {
    case_studies_config_dir <- paste0(system.file("conf/case_studies", package = "BExTE"), "/")
  }
  scenarios_config <- read_config(paste0(config_dir, "scenarios_config.yml"), scenarios_config_schema)
  simulation_config <- read_simulation_config(system.file("conf/simulation_config.yml", package = "BExTEsim"))

  if (is.null(simulation_config)) {
    stop("Simulation config is NULL")
  }

  design_prior_types <- c("ui_design_prior", "analysis_prior", "source_posterior")

  # Get the list of case studies
  case_studies <- unique(results_freq_df$case_study)

  mcmc_config <- read_config(paste0(config_dir, "/mcmc_config.yml"), mcmc_config_schema)

  for (case_study in case_studies) {
    case_study_config <- yaml::yaml.load_file(file.path(case_studies_config_dir, paste0(case_study, ".yml")))
    results_df_0 <- results_freq_df %>%
      dplyr::filter(case_study == !!case_study)
    null_space <- case_study_config$null_space
    theta_0 <- case_study_config$theta_0


    source_denominator_change_factors <- unique(results_df_0$source_denominator_change_factor)

    for (source_denominator_change_factor in source_denominator_change_factors) {
      results_df_1 <- results_df_0 %>%
        dplyr::filter(
          matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor)
        )

      source_data <- load_data(results_df_1[1, ],
                               type = "source",
                               reload_data_objects = TRUE,
                               case_studies_config_dir = case_studies_config_dir)

      source_data_df <- data.frame(source_data$to_dict())
      to_remove <- c("source_control_rate", "source_treatment_rate")

      source_data_df <- remove_columns_from_df(source_data_df, to_remove)


      target_to_source_std_ratio_range <- unique(results_df_1$target_to_source_std_ratio)
      for (target_to_source_std_ratio in target_to_source_std_ratio_range) {
        results_df_2 <- results_df_1 %>%
          dplyr::filter(matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio))


        if (nrow(results_df_2) == 0) {
          warning("Dataframe is empty")
          next
        }

        theta_0 <- case_study_config$theta_0

        # Get the list of methods
        methods <- unique(results_df_2$method)

        for (method in methods) {
          results_df_3 <- results_df_2 %>%
            dplyr::filter(method == !!method)

          target_sample_sizes <- unique(
            results_df_3[, c("target_sample_size_control", "target_sample_size_treatment")]
          )

          for (size_index in seq_len(nrow(target_sample_sizes))) {
            target_sample_size_control <- target_sample_sizes$target_sample_size_control[size_index]
            target_sample_size_treatment <- target_sample_sizes$target_sample_size_treatment[size_index]
            results_df_4 <- results_df_3 %>%
              dplyr::filter(
                target_sample_size_control == !!target_sample_size_control,
                target_sample_size_treatment == !!target_sample_size_treatment
              )

            # Get the different parameters combinations studies for this method,
            # once per trial design. Control-arm heterogeneity, loss to follow-up
            # and the event time distribution each describe a different trial, and
            # the preposterior quantities below interpolate the probability of
            # success over the treatment effect, so they have to be computed
            # within a single design rather than across a mixture of them.
            design_columns <- intersect(
              c("control_drift", time_to_event_design_columns),
              names(results_df_4)
            )
            parameters_combinations <- unique(
              results_df_4[, c("parameters", design_columns), drop = FALSE]
            )

            for (i in seq_len(nrow(parameters_combinations))) {
              # Joining against the one-row combination keeps this working
              # whichever of the design columns the results happen to carry.
              jobs[[length(jobs) + 1]] <- list(
                case_study = case_study,
                case_study_config = case_study_config,
                source_data = source_data,
                source_data_df = source_data_df,
                source_denominator_change_factor = source_denominator_change_factor,
                target_to_source_std_ratio = target_to_source_std_ratio,
                method = method,
                target_sample_size_control = target_sample_size_control,
                target_sample_size_treatment = target_sample_size_treatment,
                parameters_combination = parameters_combinations[i, , drop = FALSE],
                results_df = dplyr::inner_join(
                  results_df_4,
                  parameters_combinations[i, , drop = FALSE],
                  by = c("parameters", design_columns)
                )
              )
            }
          }
        }
      }
    }
  }

  if (length(jobs) == 0) {
    return(data.frame())
  }

  # Each combination costs ~0.1 s, most of it fitting a mixture to the
  # analysis prior, and the paper's grid has hundreds per case study.
  if (analysis_uses_cluster(parallelization, length(jobs))) {
    # The cluster simulation_analysis() shares between its steps, or one of
    # its own for a direct call.
    if (is.null(cluster)) {
      cluster <- new_analysis_cluster()
      on.exit(cluster$stop(), add = TRUE)
    }
    cluster$get()
    `%dopar%` <- foreach::`%dopar%`

    # The foreach body is evaluated outside the package namespace, where an
    # internal function is not visible; bound here, it travels to the workers
    # as an exported variable, with the namespace as its environment.
    run_job <- bayesian_ocs_for_combination
    results_list <- foreach::foreach(job = jobs, .packages = c("dplyr", "yaml")) %dopar% {
      run_job(job, simulation_config, mcmc_config, design_prior_types)
    }
  } else {
    results_list <- lapply(jobs, bayesian_ocs_for_combination,
                           simulation_config = simulation_config,
                           mcmc_config = mcmc_config,
                           design_prior_types = design_prior_types)
  }

  return(dplyr::bind_rows(results_list))
}
