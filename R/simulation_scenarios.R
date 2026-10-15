#' The design a time-to-event case study is primarily simulated under
#'
#' @description No control-arm heterogeneity, no loss to follow-up,
#'   exponential event times and proportional hazards (no delay before the
#'   treatment effect). It is the design a configuration that does not
#'   mention the sensitivity axes produces, and the one the whole scenario grid
#'   is simulated at.
#'
#' @param control_drift Control-arm heterogeneity, kappa.
#' @param dropout_probability Probability of loss to follow-up over the maximum
#'   follow-up time.
#' @param event_time_distribution Either "exponential" or "weibull".
#' @param treatment_delay Time before the treatment effect starts, in years.
#'
#' @return TRUE when the four describe the primary design.
is_primary_time_to_event_design <- function(control_drift,
                                            dropout_probability,
                                            event_time_distribution,
                                            treatment_delay = 0) {
  return(control_drift == 0 &&
    dropout_probability == 0 &&
    identical(as.character(event_time_distribution), "exponential") &&
    treatment_delay == 0)
}

#' The time-to-event design axes
#'
#' @description Each axis's value under the primary design, and how a figure
#'   that varies it titles its legend. The names are the results columns.
#'
#' @noRd
TIME_TO_EVENT_AXES <- list(
  control_drift = list(
    primary = 0,
    title = "Control-arm heterogeneity (log HR)"
  ),
  dropout_probability = list(
    primary = 0,
    title = "Probability of loss to follow-up"
  ),
  event_time_distribution = list(
    primary = "exponential",
    title = "Event time distribution"
  ),
  # Simulated in years like every other time; figures show it in weeks.
  treatment_delay = list(
    primary = 0,
    title = "Delay before the treatment effect (weeks)",
    format = function(x) round(52 * x)
  )
)

#' Rows at which every time-to-event axis but one is at its primary value
#'
#' @param df A results frame carrying the three design columns.
#' @param axis The axis left free, one of `names(TIME_TO_EVENT_AXES)`.
#'
#' @return A logical vector, one entry per row.
#'
#' @noRd
time_to_event_varies_only <- function(df, axis) {
  others <- setdiff(names(TIME_TO_EVENT_AXES), axis)
  return(Reduce(`&`, lapply(others, function(other) {
    # A column the results do not carry was only ever at its primary value.
    if (!other %in% names(df)) {
      return(rep(TRUE, nrow(df)))
    }
    return(df[[other]] == TIME_TO_EVENT_AXES[[other]]$primary)
  })))
}

#' How many time-to-event design axes are off their primary value
#'
#' @description A sensitivity analysis that varies one axis at a time only
#'   needs the designs for which this is at most one; crossing the axes
#'   simulates every combination of them as well.
#'
#' @inheritParams is_primary_time_to_event_design
#'
#' @return An integer vector, 0 for the primary design.
time_to_event_axes_off_primary <- function(control_drift,
                                           dropout_probability,
                                           event_time_distribution,
                                           treatment_delay = 0) {
  return(as.integer(control_drift != 0) +
    as.integer(dropout_probability != 0) +
    as.integer(as.character(event_time_distribution) != "exponential") +
    as.integer(treatment_delay != 0))
}

#' The point at which the time-to-event sensitivity designs are simulated
#'
#' @description Crossing control-arm heterogeneity, loss to follow-up and the
#'   event time distribution with the whole grid multiplies a run by the product
#'   of their lengths, and most of that is spent re-simulating designs nobody
#'   reads. When `sensitivity_reference` names a sample size factor and a source
#'   denominator change factor, a design other than the primary one is simulated
#'   only there, so the sensitivity analysis still varies the drift across its
#'   full range at a single trial size. Without the key the axes cross the whole
#'   grid, which is what a configuration written before it existed does.
#'
#' @param scenarios_config Configuration for the simulation.
#' @param case_study_config Configuration for the case study.
#'
#' @return A list with `sample_size_factor` and `denominator_change_factor`, or
#'   NULL when every design is simulated across the whole grid.
time_to_event_sensitivity_reference <- function(scenarios_config,
                                                case_study_config) {
  if (!identical(case_study_config$endpoint, "time_to_event")) {
    return(NULL)
  }

  reference <- scenarios_config$sensitivity_reference
  if (is.null(reference)) {
    return(NULL)
  }

  missing_keys <- setdiff(
    c("sample_size_factor", "denominator_change_factor"), names(reference)
  )
  if (length(missing_keys) > 0) {
    stop(paste0(
      "sensitivity_reference is missing ", paste(missing_keys, collapse = " and "),
      ". It names the one sample size factor and source denominator change",
      " factor the sensitivity designs are simulated at."
    ), call. = FALSE)
  }

  return(list(
    sample_size_factor = as.numeric(reference$sample_size_factor),
    denominator_change_factor = as.numeric(reference$denominator_change_factor)
  ))
}

#' Get important drift values for a given source treatment effect and case study configuration.
#'
#' @param source_treatment_effect The source treatment effect.
#' @param case_study_config Configuration for the case study.
#'
#' @return A vector of important drift values to consider.
important_drift_values <- function(source_treatment_effect,
                                   case_study_config) {
  vec <- c(-source_treatment_effect, -source_treatment_effect / 2, 0)
  if (case_study_config$summary_measure_likelihood == "normal") {
    return(vec)
  } else if (case_study_config$summary_measure_likelihood == "binomial") {
    lower_bound <- -1 - source_treatment_effect
    upper_bound <- 1 - source_treatment_effect
    vec <- vec[vec <= upper_bound & vec >= lower_bound]
    return(vec)
  } else {
    return(invisible(NULL))
  }
}

#' Compute the control drift range for a given source treatment effect.
#'
#' @description Control drift is the control-arm heterogeneity between the source
#'   and the target population, kappa = log(control rate in the target) - log(control
#'   rate in the source). Only the time-to-event endpoint simulates it; every other
#'   case study keeps a single scenario with no control drift. The no-heterogeneity
#'   scenario is always included, as it is the reference the others are read against.
#'
#' @param source_treatment_effect The source treatment effect.
#' @param scenarios_config Configuration for the simulation.
#' @param case_study_config Configuration for the case study.
#'
#' @return A vector representing the control drift range.
compute_control_drift_range <- function(source_treatment_effect,
                                        scenarios_config = NULL,
                                        case_study_config = NULL) {
  if (!identical(case_study_config$endpoint, "time_to_event")) {
    return(c(0))
  }

  control_drift_range <- unlist(scenarios_config$control_drift_range)
  if (is.null(control_drift_range)) {
    return(c(0))
  }

  return(sort(unique(c(0, control_drift_range))))
}

#' Compute the ranges of the time-to-event design axes.
#'
#' @description The probability of loss to follow-up, the distribution of the
#'   event times and the delay before the treatment effect only change how a
#'   time-to-event trial is simulated, so every other
#'   case study keeps a single scenario at the primary design rather than paying for
#'   a cross product that would generate identical data.
#'
#' @param scenarios_config Configuration for the simulation.
#' @param case_study_config Configuration for the case study.
#'
#' @return A list with the `dropout_probability`, `event_time_distribution` and
#'   `treatment_delay` ranges.
compute_time_to_event_ranges <- function(scenarios_config, case_study_config) {
  if (!identical(case_study_config$endpoint, "time_to_event")) {
    return(list(dropout_probability = 0, event_time_distribution = "exponential",
                treatment_delay = 0))
  }

  dropout_probability <- unlist(scenarios_config$dropout_probability)
  event_time_distribution <- unlist(scenarios_config$event_time_distribution)
  treatment_delay <- unlist(scenarios_config$treatment_delay)

  return(list(
    dropout_probability = if (is.null(dropout_probability)) {
      0
    } else {
      sort(unique(dropout_probability))
    },
    event_time_distribution = if (is.null(event_time_distribution)) {
      "exponential"
    } else {
      unique(event_time_distribution)
    },
    # Zero is always simulated: it is the proportional-hazards reference.
    treatment_delay = sort(unique(c(0, treatment_delay)))
  ))
}

#' Compute the source denominator range for a given source data, simulation configuration, and case study configuration.
#'
#' @param source_data The source data.
#' @param scenarios_config Configuration for the simulation.
#' @param case_study_config Configuration for the case study.
#'
#' @return A vector representing the source denominator range.
compute_source_denominator_range <- function(source_data,
                                             scenarios_config,
                                             case_study_config) {
  # Compute, for each value of the source denominator, the corresponding numerator so that the treatment effec stays the same
  # Default behavior, used when the endpoint is continuous
  source_denominator_range <- c(NA)
  source_treatment_rate_range <- NA

  denominator_change_factor <- unlist(scenarios_config$denominator_change_factor)
  if (source_data$endpoint == "binary" &&
      case_study_config$summary_measure_likelihood == "normal") {
    odds_control <- source_data$control_rate / (1 - source_data$control_rate)
    source_denominator_range <- denominator_change_factor * odds_control

    source_treatment_rate_range <- 1 / (1 + exp(-(
      source_data$treatment_effect_estimate + log(source_denominator_range)
    )))
  } else if (source_data$endpoint == "time_to_event" ||
             source_data$endpoint == "recurrent_event") {
    source_denominator_range <- source_data$control_rate * denominator_change_factor

    source_treatment_rate_range <- source_denominator_range * exp(source_data$treatment_effect_estimate)
  } else {
    denominator_change_factor <- c(1)
  }

  return(
    list(
      source_denominator = source_denominator_range,
      source_denominator_change_factor = denominator_change_factor
    )
  )
}

#' Combine parameters for a given method.
#'
#' @param method The method name.
#'
#' @return A list of combinations of parameters for the method.
combine_parameters <- function(method) {
  if (is.null(methods_dict[[method]])) {
    return(NULL)
  } else {
    param_info <- methods_dict[[method]]
    ranges <- lapply(param_info, function(x) {
      return(x$range)
    })
    combinations <- expand.grid(ranges)
    return(combinations)
  }
}

#' Generate simulation scenarios based on the configuration directory.
#'
#' @param config_dir The directory containing the configuration files.
#' @param scenarios_config The simulation configuration file.
#' @param case_studies_config_dir Directory containing the case study YAML
#'   files.
#'
#' @return A list of simulation scenarios.
#' @export
simulation_scenarios <- function(config_dir, scenarios_config, case_studies_config_dir) {
  source(paste0(config_dir, "methods_config.R"))
  case_studies <- scenarios_config$case_studies

  methods_configurations_df <- data.frame()

  # We store the parameters in a dataframe with two columns : method and parameters
  for (method in names(methods_dict)) {
    if (method %in% scenarios_config$methods) {
      parameters_combinations <- combine_parameters(method)

      converted_parameters_df <- data.frame(method = method, parameters = I(apply(parameters_combinations, 1, function(row) {
        params_list <- as.list(row)
        names(params_list) <- colnames(parameters_combinations)
        return(params_list)
      })))

      methods_configurations_df <- rbind(methods_configurations_df, converted_parameters_df)
    }
  }

  if (any(!(scenarios_config$methods %in% names(methods_dict)))) {
    stop("Some methods listed in the simulation configuration need configuration.")
  }

  scenarios_configurations <- list()

  for (case_study in case_studies) {
    case_study_config <- yaml::read_yaml(file.path(case_studies_config_dir, paste0(case_study, ".yml")))
    source_data <- ObservedSourceData$new(case_study_config)

    sampling_approximation <- case_study_config$sampling_approximation
    theta_0 <- case_study_config$theta_0
    null_space <- case_study_config$null_space

    source_denominator <- compute_source_denominator_range(source_data, scenarios_config, case_study_config)
    source_denominator <- data.frame(source_denominator)

    drift_range <- compute_drift_range(scenarios_config, case_study_config)

    source_treatment_effect_estimate <- source_data$treatment_effect_estimate

    case_study_config <- yaml::read_yaml(file.path(case_studies_config_dir, paste0(case_study, ".yml")))

    mandatory_drift_values <- important_drift_values(source_treatment_effect_estimate, case_study_config)
    drift_range <- sort(c(drift_range, mandatory_drift_values))

    control_drift_range <- compute_control_drift_range(
      source_treatment_effect_estimate,
      scenarios_config = scenarios_config,
      case_study_config = case_study_config
    )

    time_to_event_ranges <- compute_time_to_event_ranges(
      scenarios_config = scenarios_config,
      case_study_config = case_study_config
    )

    if (case_study_config$endpoint == "continuous") {
      target_to_source_std_ratio_range <- scenarios_config$target_to_source_std_ratio_range
    } else {
      target_to_source_std_ratio_range <- NA
    }

    n <- nrow(source_denominator)

    # Create drift/control drift pairs
    # stringsAsFactors = FALSE keeps the event time distribution a character
    # column: a factor would reach the results as its integer code.
    drift_combinations <- expand.grid(
      drift = drift_range,
      control_drift = control_drift_range,
      target_to_source_std_ratio = target_to_source_std_ratio_range,
      dropout_probability = time_to_event_ranges$dropout_probability,
      event_time_distribution = time_to_event_ranges$event_time_distribution,
      treatment_delay = time_to_event_ranges$treatment_delay,
      stringsAsFactors = FALSE
    )

    source_denominator_repeated <- source_denominator[rep(seq_len(n), each = nrow(drift_combinations)), ]
    drift_combinations_repeated <- drift_combinations[rep(seq_len(nrow(drift_combinations)), times = n), ]

    drift_combinations <- cbind(source_denominator_repeated,
                                drift_combinations_repeated)

    drift_combinations <- unique(drift_combinations)

    # Each sensitivity design moves one axis at a time, when asked to, so the
    # figures that vary one axis do not pay for every combination of them.
    if (isTRUE(scenarios_config$sensitivity_one_at_a_time)) {
      off_primary <- time_to_event_axes_off_primary(
        drift_combinations$control_drift,
        drift_combinations$dropout_probability,
        drift_combinations$event_time_distribution,
        drift_combinations$treatment_delay
      )
      drift_combinations <- drift_combinations[off_primary <= 1, , drop = FALSE]
    }

    # Exclude combinations that would result in rates that are outside [0,1]
    if (case_study_config$summary_measure_likelihood == "binomial") {
      lower_bound <- pmax(
        -source_data$treatment_rate - drift_combinations$control_drift,
        -1 - source_data$treatment_effect_estimate
      )
      upper_bound <- pmin(
        1 - source_data$treatment_rate - drift_combinations$control_drift,
        1 - source_data$treatment_effect_estimate
      )

      condition <- drift_combinations$drift >= lower_bound &
        drift_combinations$drift <= upper_bound
      drift_combinations <- subset(drift_combinations, condition)
    }

    sample_size_factors <- unlist(scenarios_config$sample_size_factors)
    total_target_sample_sizes <- (source_data$sample_size_control + source_data$sample_size_treatment) / sample_size_factors
    allocation_ratio <- target_allocation_ratio(scenarios_config, case_study)
    current_scenario_configurations <- list()

    sensitivity_reference <- time_to_event_sensitivity_reference(
      scenarios_config = scenarios_config,
      case_study_config = case_study_config
    )

    if (!is.null(sensitivity_reference) &&
        !sensitivity_reference$sample_size_factor %in% sample_size_factors) {
      stop(paste0(
        "sensitivity_reference$sample_size_factor (",
        sensitivity_reference$sample_size_factor,
        ") is not one of the sample_size_factors being simulated, so the",
        " sensitivity designs would never be simulated at all."
      ), call. = FALSE)
    }

    for (i in seq_len(nrow(drift_combinations))) {
      drift <- drift_combinations[i, "drift"]
      control_drift <- drift_combinations[i, "control_drift"]
      source_denominator <- drift_combinations[i, "source_denominator"]
      source_denominator_change_factor <- drift_combinations[i, "source_denominator_change_factor"]
      dropout_probability <- drift_combinations[i, "dropout_probability"]
      event_time_distribution <- drift_combinations[i, "event_time_distribution"]
      treatment_delay <- drift_combinations[i, "treatment_delay"]
      treatment_drift <- drift + control_drift
      target_treatment_effect <- source_treatment_effect_estimate + drift

      primary_design <- is_primary_time_to_event_design(
        control_drift, dropout_probability, event_time_distribution, treatment_delay
      )

      for (factor_index in seq_along(total_target_sample_sizes)) {
        total_target_sample_size <- total_target_sample_sizes[factor_index]
        target_arms <- split_target_sample_size(total_target_sample_size, allocation_ratio)

        # The whole grid is simulated under the primary design; the sensitivity
        # designs are simulated only at the reference point, so they still vary
        # the drift without multiplying the run by the sample sizes and source
        # denominators as well.
        if (!is.null(sensitivity_reference) && !primary_design &&
            !(sample_size_factors[factor_index] == sensitivity_reference$sample_size_factor &&
              isTRUE(all.equal(
                source_denominator_change_factor,
                sensitivity_reference$denominator_change_factor
              )))) {
          next
        }

        current_scenario_configurations[[length(current_scenario_configurations) + 1]] <- list(
          target_sample_size_control = target_arms$control,
          target_sample_size_treatment = target_arms$treatment,
          drift = drift,
          treatment_drift = treatment_drift,
          control_drift = control_drift,
          source_denominator = source_denominator,
          source_denominator_change_factor = source_denominator_change_factor,
          case_study = case_study,
          sampling_approximation = sampling_approximation,
          source_treatment_effect_estimate = source_treatment_effect_estimate,
          target_treatment_effect = target_treatment_effect,
          target_to_source_std_ratio = drift_combinations[i, "target_to_source_std_ratio"],
          dropout_probability = dropout_probability,
          event_time_distribution = event_time_distribution,
          treatment_delay = treatment_delay,
          theta_0 = theta_0,
          null_space = null_space
        )
      }
    }
    scenarios_configurations <- c(scenarios_configurations,
                                  current_scenario_configurations)
  }
  # Cartesian product between methods configurations and scenarios
  # rbind() over a list of lists gives a matrix of mode list, which turns
  # every column into a list of length-1 values - see
  # unwrap_scalar_list_columns().
  scenarios_configurations_df <- unwrap_scalar_list_columns(
    data.frame(do.call(rbind, scenarios_configurations))
  )

  cases_df <- merge(methods_configurations_df,
                    scenarios_configurations_df,
                    by = NULL)

  return(cases_df)
}

scenarios_table_ranges <- function(cases_df, results_path) {
  dir.create(results_path, showWarnings = FALSE, recursive = TRUE)
  file_path <- paste0(results_path, "/drift_ranges.tex")

  # Initialize an empty dataframe
  df <- data.frame()


  for (case_study in unique(cases_df$case_study)) {
    subdf <- cases_df[cases_df$case_study == case_study, ]
    subdf$theta_0 <- as.numeric(subdf$theta_0)
    subdf$source_treatment_effect_estimate <- as.numeric(subdf$source_treatment_effect_estimate)

    # Compute the ranges
    drift_range <- range(subdf$drift)
    treatment_range <- range(subdf$target_treatment_effect)

    # Bind rows to df
    df <- rbind(
      df,
      data.frame(
        Drug = paste0(toupper(substring(case_study, 1, 1)), substring(case_study, 2)),
        Drift_range_lower = signif(drift_range[1], 3),
        Drift_range_upper = signif(drift_range[2], 3),
        Drift_with_theta_T_equals_theta_0 = signif(
          unique(subdf$theta_0 - subdf$source_treatment_effect_estimate),
          3
        ),
        Treatment_effect_range_lower = signif(treatment_range[1], 3),
        Treatment_effect_range_upper = signif(treatment_range[2], 3)
      )
    )
  }

  # Function to generate the LaTeX table
  generate_latex_table <- function(df, file_path) {
    withr::local_output_sink(file_path)

    cat(
      "\\begin{table}[h]
\\centering
\\begin{tabular}{|l|c|c|c|}
\\textbf{Case study} & \\textbf{Drift range} & \\textbf{Drift with $\\theta_T^{(true)} = \\theta_0$} & \\textbf{Treatment effect range} \\\\
\\midrule
"
    )

    for (i in seq_len(nrow(df))) {
      # Concatenate ranges
      drift_range <- paste0("[",
                            df$Drift_range_lower[i],
                            ",",
                            df$Drift_range_upper[i],
                            "]")
      treatment_range <- paste0(
        "[",
        df$Treatment_effect_range_lower[i],
        ",",
        df$Treatment_effect_range_upper[i],
        "]"
      )

      cat(
        paste0(
          "\\hline \\text{",
          df$Drug[i],
          "} & ",
          drift_range,
          " & ",
          df$Drift_with_theta_T_equals_theta_0[i],
          " & ",
          treatment_range,
          " \\\\
"
        )
      )
    }

    cat(
      "\\hline
    \\end{tabular}
\\caption{Drift ranges considered for each case study.}
\\label{tab:drift_ranges}
\\end{table}"
    )

    # withr closes the sink when the function exits
    return(invisible(NULL))
  }

  # Generate the LaTeX table
  return(generate_latex_table(df, file_path))
}
