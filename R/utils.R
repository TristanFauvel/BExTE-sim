#' Match a scenario value, or keep rows where it does not apply
#'
#' Scenario columns such as `source_denominator_change_factor` hold `NA` for
#' methods and case studies the factor does not concern. Those rows belong to
#' every scenario. Use inside [dplyr::filter()].
#'
#' @param column A column of the data frame.
#' @param value The scenario value to keep.
#'
#' @return A logical vector.
#' @keywords internal
matches_or_na <- function(column, value) {
  column == value | is.na(column)
}

## Restrict results to the parameter combinations a figure or table plots.
##
## The merge on the shared parameter columns is the filter, and it has to be
## handed the data frame itself. `parameters_combinations[, ]` drops a
## single-column frame to a bare vector, leaving merge() with no column in
## common: it cross-joins - keeping every parameter value - or returns nothing
## when a name happens to collide. Figure S3 selects the conditional power
## prior at gamma = 0.25 and is the only manifest entry whose selection has
## exactly one parameter, so it drew all five gamma values as a single
## zig-zagging line while two-parameter figures were unaffected.
filter_to_parameter_combinations <- function(results_df, parameters_combinations) {
  return(merge(results_df, parameters_combinations))
}

frequentist_col_types <- readr::cols(
  method = readr::col_character(),
  parameters = readr::col_character(),
  target_sample_size_control = readr::col_double(),
  target_sample_size_treatment = readr::col_double(),
  drift = readr::col_double(),
  treatment_drift = readr::col_double(),
  control_drift = readr::col_double(),
  source_denominator = readr::col_double(),
  source_denominator_change_factor = readr::col_double(),
  case_study = readr::col_character(),
  parallelization = readr::col_logical(),
  sampling_approximation = readr::col_logical(),
  source_treatment_effect_estimate = readr::col_double(),
  target_treatment_effect = readr::col_double(),
  target_to_source_std_ratio = readr::col_double(),
  dropout_probability = readr::col_double(),
  event_time_distribution = readr::col_character(),
  treatment_delay = readr::col_double(),
  theta_0 = readr::col_double(),
  null_space = readr::col_character(),
  summary_measure_likelihood = readr::col_character(),
  target_standard_deviation = readr::col_double(),
  target_treatment_rate = readr::col_double(),
  target_control_rate = readr::col_double(),
  source_standard_error = readr::col_double(),
  source_sample_size_control = readr::col_double(),
  source_sample_size_treatment = readr::col_double(),
  equivalent_source_sample_size_per_arm = readr::col_double(),
  endpoint = readr::col_character(),
  source_control_rate = readr::col_double(),
  source_treatment_rate = readr::col_double(),
  rng_state = readr::col_character(),
  computation_time = readr::col_double(),
  success_proba = readr::col_double(),
  mcse_success_proba = readr::col_double(),
  conf_int_success_proba_lower = readr::col_double(),
  conf_int_success_proba_upper = readr::col_double(),
  coverage = readr::col_double(),
  conf_int_coverage_lower = readr::col_double(),
  conf_int_coverage_upper = readr::col_double(),
  mse = readr::col_double(),
  conf_int_mse_lower = readr::col_double(),
  conf_int_mse_upper = readr::col_double(),
  bias = readr::col_double(),
  conf_int_bias_lower = readr::col_double(),
  conf_int_bias_upper = readr::col_double(),
  posterior_mean = readr::col_double(),
  conf_int_posterior_mean_lower = readr::col_double(),
  conf_int_posterior_mean_upper = readr::col_double(),
  posterior_median = readr::col_double(),
  conf_int_posterior_median_lower = readr::col_double(),
  conf_int_posterior_median_upper = readr::col_double(),
  precision = readr::col_double(),
  conf_int_precision_lower = readr::col_double(),
  conf_int_precision_upper = readr::col_double(),
  interval_score = readr::col_double(),
  conf_int_interval_score_lower = readr::col_double(),
  conf_int_interval_score_upper = readr::col_double(),
  credible_interval_lower = readr::col_double(),
  credible_interval_upper = readr::col_double(),
  posterior_parameters = readr::col_character(),
  ess_moment = readr::col_double(),
  conf_int_ess_moment_lower = readr::col_double(),
  conf_int_ess_moment_upper = readr::col_double(),
  ess_precision = readr::col_double(),
  conf_int_ess_precision_lower = readr::col_double(),
  conf_int_ess_precision_upper = readr::col_double(),
  ess_elir = readr::col_double(),
  conf_int_ess_elir_lower = readr::col_double(),
  conf_int_ess_elir_upper = readr::col_double(),
  rhat = readr::col_double(),
  conf_int_rhat_lower = readr::col_double(),
  conf_int_rhat_upper = readr::col_double(),
  mcmc_ess = readr::col_double(),
  conf_int_mcmc_ess_lower = readr::col_double(),
  conf_int_mcmc_ess_upper = readr::col_double(),
  n_divergences = readr::col_double(),
  conf_int_n_divergences_lower = readr::col_double(),
  conf_int_n_divergences_upper = readr::col_double(),
  warning = readr::col_logical(),
  exact_ocs = readr::col_logical(),
  enumeration_omitted_mass = readr::col_double(),
  tie = readr::col_double(),
  mcse_tie = readr::col_double(),
  conf_int_tie_lower = readr::col_double(),
  conf_int_tie_upper = readr::col_double(),
  frequentist_power_at_equivalent_tie = readr::col_double(),
  frequentist_power_at_equivalent_tie_lower = readr::col_double(),
  frequentist_power_at_equivalent_tie_upper = readr::col_double(),
  frequentist_test = readr::col_character(),
  frequentist_reference_calibration = readr::col_character(),
  nominal_frequentist_power_separate = readr::col_double(),
  nominal_frequentist_power_separate_lower = readr::col_double(),
  nominal_frequentist_power_separate_upper = readr::col_double(),
  nominal_frequentist_power_pooling = readr::col_double(),
  nominal_frequentist_power_pooling_lower = readr::col_double(),
  nominal_frequentist_power_pooling_upper = readr::col_double()
)


bayesian_col_types <- readr::cols(
  case_study = readr::col_character(),
  method = readr::col_character(),
  target_sample_size_control = readr::col_double(),
  target_sample_size_treatment = readr::col_double(),
  parameters = readr::col_character(),
  source_denominator_change_factor = readr::col_double(),
  target_to_source_std_ratio = readr::col_double(),
  design_prior_type = readr::col_character(),
  prior_proba_success = readr::col_double(),
  prior_proba_no_benefit = readr::col_double(),
  prior_proba_benefit = readr::col_double(),
  prepost_proba_FP = readr::col_double(),
  prepost_proba_TP = readr::col_double(),
  average_tie = readr::col_double(),
  average_power = readr::col_double(),
  upper_bound_proba_FP = readr::col_double(),
  source_treatment_effect_estimate = readr::col_double(),
  source_standard_error = readr::col_double(),
  endpoint = readr::col_character(),
  summary_measure_likelihood = readr::col_character(),
  source_sample_size_control = readr::col_double(),
  source_sample_size_treatment = readr::col_double(),
  equivalent_source_sample_size_per_arm = readr::col_double()
)


sweet_spot_col_types <- readr::cols(
  method = readr::col_character(),
  parameters = readr::col_character(),
  target_sample_size_control = readr::col_double(),
  target_sample_size_treatment = readr::col_double(),
  control_drift = readr::col_double(),
  source_denominator = readr::col_double(),
  source_denominator_change_factor = readr::col_double(),
  case_study = readr::col_character(),
  parallelization = readr::col_logical(),
  sampling_approximation = readr::col_logical(),
  source_treatment_effect_estimate = readr::col_double(),
  target_to_source_std_ratio = readr::col_double(),
  dropout_probability = readr::col_double(),
  event_time_distribution = readr::col_character(),
  treatment_delay = readr::col_double(),
  theta_0 = readr::col_double(),
  null_space = readr::col_character(),
  summary_measure_likelihood = readr::col_character(),
  source_standard_error = readr::col_double(),
  source_sample_size_control = readr::col_double(),
  source_sample_size_treatment = readr::col_double(),
  equivalent_source_sample_size_per_arm = readr::col_double(),
  endpoint = readr::col_character(),
  source_control_rate = readr::col_double(),
  source_treatment_rate = readr::col_double(),
  sweet_spot_lower = readr::col_double(),
  sweet_spot_upper = readr::col_double(),
  sweet_spot_width = readr::col_double(),
  drift_range_upper = readr::col_double(),
  drift_range_lower = readr::col_double(),
  metric = readr::col_character()
)


# The time-to-event design axes. They were added to the scenario grid after some
# results had already been written, and a results file that predates them
# describes a single design, so keys built from one leave them out - see the
# intersect() calls in the analysis.
time_to_event_design_columns <- c("dropout_probability", "event_time_distribution",
                                  "treatment_delay")

# The columns holding the target arm sizes. Together they identify the target
# sample size of a scenario.
target_size_columns <- c("target_sample_size_control", "target_sample_size_treatment")

#' Replace the per-arm target sample size of older results by the arm sizes
#'
#' @description Results written before the target arms could differ carry a
#'   single `target_sample_size_per_arm` column. Both arms had that size, so it
#'   fills `target_sample_size_control` and `target_sample_size_treatment`.
#'
#' @param df A results data frame.
#' @return `df` with the two arm-size columns and without the old one.
#' @noRd
upgrade_target_size_columns <- function(df) {
  if (!"target_sample_size_per_arm" %in% names(df)) {
    return(df)
  }
  for (column in target_size_columns) {
    if (!column %in% names(df)) {
      df[[column]] <- df$target_sample_size_per_arm
    }
  }
  df$target_sample_size_per_arm <- NULL
  return(df)
}

# The columns that, together with the drift, identify one simulated scenario.
# They are the key the analysis joins a method's results to the separate
# analysis on, so a design axis missing here would make that join fan out.
scenario_columns <- c("target_sample_size_control", "target_sample_size_treatment", "control_drift", "source_denominator", "source_denominator_change_factor",
                      "case_study", "sampling_approximation", "source_treatment_effect_estimate",
                      "target_to_source_std_ratio", "dropout_probability", "event_time_distribution",
                      "treatment_delay", "theta_0", "null_space",
                      "summary_measure_likelihood", "source_standard_error", "source_sample_size_control",
                      "source_sample_size_treatment", "equivalent_source_sample_size_per_arm", "endpoint",
                      "source_control_rate", "source_treatment_rate")

unique_scenario_columns <- c("case_study", "target_sample_size_control", "target_sample_size_treatment", "source_denominator_change_factor",
                      "target_to_source_std_ratio", "dropout_probability", "event_time_distribution",
                      "treatment_delay")


expected_colnames_scenario <- c(
  "method",
  "parameters",
  "target_sample_size_control", "target_sample_size_treatment",
  "drift",
  "treatment_drift",
  "control_drift",
  "source_denominator",
  "source_denominator_change_factor",
  "case_study",
  "parallelization",
  "sampling_approximation",
  "source_treatment_effect_estimate",
  "target_treatment_effect",
  "target_to_source_std_ratio",
  "dropout_probability",
  "event_time_distribution",
  "treatment_delay",
  "theta_0",
  "null_space"
)

expected_colnames_results <- c(
  "success_proba",
  "mcse_success_proba",
  "conf_int_success_proba_lower",
  "conf_int_success_proba_upper",
  "coverage",
  "conf_int_coverage_lower",
  "conf_int_coverage_upper",
  "mse",
  "conf_int_mse_lower",
  "conf_int_mse_upper",
  "bias",
  "conf_int_bias_lower",
  "conf_int_bias_upper",
  "posterior_mean",
  "conf_int_posterior_mean_lower",
  "conf_int_posterior_mean_upper",
  "posterior_median",
  "conf_int_posterior_median_lower",
  "conf_int_posterior_median_upper",
  "precision",
  "conf_int_precision_lower",
  "conf_int_precision_upper",
  "interval_score",
  "conf_int_interval_score_lower",
  "conf_int_interval_score_upper",
  "credible_interval_lower",
  "credible_interval_upper",
  "posterior_parameters",
  "ess_moment",
  "conf_int_ess_moment_lower",
  "conf_int_ess_moment_upper",
  "ess_precision",
  "conf_int_ess_precision_lower",
  "conf_int_ess_precision_upper",
  "ess_elir",
  "conf_int_ess_elir_lower",
  "conf_int_ess_elir_upper",
  "rhat",
  "conf_int_rhat_lower",
  "conf_int_rhat_upper",
  "mcmc_ess",
  "conf_int_mcmc_ess_lower",
  "conf_int_mcmc_ess_upper",
  "n_divergences",
  "conf_int_n_divergences_lower",
  "conf_int_n_divergences_upper",
  "warning",
  "exact_ocs",
  "enumeration_omitted_mass"
)

expected_colnames_source <- c(
  "source_treatment_effect_estimate",
  "source_standard_error",
  "source_sample_size_control",
  "source_sample_size_treatment",
  "equivalent_source_sample_size_per_arm",
  "summary_measure_likelihood",
  "endpoint",
  "source_control_rate",
  "source_treatment_rate"
)


#' Column types expected across the simulation results
#'
#' @description Maps the column names shared by the scenario, source and
#'   results frames to the type each is expected to hold. Only columns whose
#'   type is load-bearing downstream are listed; `check_colnames()` ignores
#'   columns absent from the spec.
#'
#' @keywords internal
expected_coltypes <- c(
  # Scenario
  method = "character",
  parameters = "list|character",
  target_sample_size_control = "numeric",
  target_sample_size_treatment = "numeric",
  drift = "numeric",
  treatment_drift = "numeric",
  control_drift = "numeric",
  source_denominator = "numeric",
  source_denominator_change_factor = "numeric",
  case_study = "character",
  sampling_approximation = "logical",
  source_treatment_effect_estimate = "numeric",
  target_treatment_effect = "numeric",
  target_to_source_std_ratio = "numeric",
  theta_0 = "numeric",
  null_space = "character",
  # Source data
  source_standard_error = "numeric",
  source_sample_size_control = "numeric",
  source_sample_size_treatment = "numeric",
  equivalent_source_sample_size_per_arm = "numeric",
  summary_measure_likelihood = "character",
  endpoint = "character",
  # Results
  success_proba = "numeric",
  mcse_success_proba = "numeric",
  conf_int_success_proba_lower = "numeric",
  conf_int_success_proba_upper = "numeric",
  coverage = "numeric",
  mse = "numeric",
  bias = "numeric",
  posterior_mean = "numeric",
  posterior_median = "numeric",
  precision = "numeric",
  interval_score = "numeric",
  credible_interval_lower = "numeric",
  credible_interval_upper = "numeric",
  posterior_parameters = "list|character",
  ess_moment = "numeric",
  ess_precision = "numeric",
  ess_elir = "numeric",
  rhat = "numeric",
  mcmc_ess = "numeric",
  n_divergences = "numeric",
  exact_ocs = "logical",
  enumeration_omitted_mass = "numeric",
  # Derived by the analysis layer
  tie = "numeric",
  conf_int_tie_lower = "numeric",
  conf_int_tie_upper = "numeric",
  frequentist_power_at_equivalent_tie = "numeric",
  frequentist_power_at_equivalent_tie_lower = "numeric",
  frequentist_power_at_equivalent_tie_upper = "numeric",
  nominal_frequentist_power_separate = "numeric",
  nominal_frequentist_power_separate_lower = "numeric",
  nominal_frequentist_power_separate_upper = "numeric",
  nominal_frequentist_power_pooling = "numeric",
  nominal_frequentist_power_pooling_lower = "numeric",
  nominal_frequentist_power_pooling_upper = "numeric"
)


#' Columns every consumer of a results frame relies on
#'
#' @description The identity columns that the analysis, plot and table layers
#'   read by name from whichever results frame they are handed. `drift` is
#'   deliberately absent: the Bayesian operating-characteristics frame does not
#'   carry it, and these columns must hold for both.
#'
#' @keywords internal
required_colnames_consumer <- c(
  "case_study",
  "method",
  "parameters",
  "target_sample_size_control", "target_sample_size_treatment",
  "source_denominator_change_factor"
)


#' Look up the expected types for a set of columns
#'
#' @param expected_colnames The column names to look up.
#'
#' @return The subset of [expected_coltypes] describing those columns. Columns
#'   the spec does not describe are omitted rather than reported, so that a
#'   frame carrying extra bookkeeping columns still passes.
#'
#' @keywords internal
default_coltypes <- function(expected_colnames) {
  return(expected_coltypes[intersect(names(expected_coltypes), expected_colnames)])
}


#' Predicates used to check a column against a declared type
#'
#' @description The vocabulary `check_colnames()` accepts in its `types`
#'   argument. `numeric` deliberately accepts integer columns, since sample
#'   sizes and replicate counts are read back as either.
#'
#' @keywords internal
column_type_predicates <- list(
  numeric = is.numeric,
  integer = is.integer,
  character = is.character,
  logical = is.logical,
  factor = is.factor,
  list = is.list
)


#' Describe the type a column actually holds
#'
#' @param column A dataframe column.
#'
#' @return A single string naming the column's type, for use in error messages.
#'
#' @keywords internal
describe_column_type <- function(column) {
  return(paste(class(column), collapse = "/"))
}


#' Check columns in a dataframe
#'
#' This function checks if a dataframe contains the expected columns and only the expected columns.
#' When `types` is supplied it additionally checks that those columns hold the
#' declared type, reporting every mismatch at once.
#'
#' @param df A dataframe to check.
#' @param expected_colnames A character vector of column names the dataframe should contain.
#' @param types An optional named character vector mapping column names to the
#'   type they are expected to hold, using the vocabulary of
#'   [column_type_predicates]. Names must be among `expected_colnames`; columns
#'   the spec does not mention are not type checked.
#'
#' @return No return value, called for side effects.
#'
#' @keywords internal
check_colnames <- function(df,
                           expected_colnames,
                           types = default_coltypes(expected_colnames)) {
  missing_cols <- setdiff(expected_colnames, colnames(df))
  if (length(missing_cols) > 0) {
    error_message <- paste(
      "The dataframe is missing the following columns:",
      paste(missing_cols, collapse = ", ")
    )
    stop(error_message)
  }

  missing_cols <- setdiff(colnames(df), expected_colnames)
  if (length(missing_cols) > 0) {
    error_message <- paste(
      "The dataframe has unexpected columns:",
      paste(missing_cols, collapse = ", ")
    )
    stop(error_message)
  }

  return(check_coltypes(df, types, expected_colnames))
}


#' Check that a dataframe carries the columns a consumer relies on
#'
#' @description Unlike [check_colnames()], this tolerates additional columns.
#'   Use it at the entry point of the analysis, plot and table layers, which
#'   receive frames enriched with derived columns but read a known subset of
#'   them by name.
#'
#' @param df A dataframe to check.
#' @param required_colnames A character vector of column names that must be
#'   present.
#' @param types An optional named character vector of expected types; defaults
#'   to the entries of [expected_coltypes] describing the required columns.
#' @param context An optional string naming the caller, included in the error
#'   so that a schema drift points at the consumer that noticed it.
#'
#' @return No return value, called for side effects.
#'
#' @keywords internal
check_required_colnames <- function(df,
                                    required_colnames,
                                    types = default_coltypes(required_colnames),
                                    context = NULL) {
  prefix <- if (is.null(context)) "" else paste0(context, ": ")

  missing_cols <- setdiff(required_colnames, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste0(
      prefix, "the dataframe is missing the following required columns: ",
      paste(missing_cols, collapse = ", ")
    ))
  }

  return(tryCatch(
    check_coltypes(df, types, required_colnames),
    error = function(e) stop(paste0(prefix, conditionMessage(e)), call. = FALSE)
  ))
}


#' Check the types of a dataframe's columns
#'
#' @param df A dataframe whose columns have already been checked for presence.
#' @param types A named character vector mapping column names to expected types.
#' @param expected_colnames The column names `types` is allowed to mention.
#'
#' @return No return value, called for side effects.
#'
#' @keywords internal
check_coltypes <- function(df, types, expected_colnames = names(types)) {
  if (length(types) == 0) {
    return(invisible(NULL))
  }

  if (is.null(names(types)) || any(names(types) == "")) {
    stop("The type specification must be a named character vector.")
  }

  alternatives <- strsplit(types, "|", fixed = TRUE)
  unknown_types <- setdiff(unique(unlist(alternatives)), names(column_type_predicates))
  if (length(unknown_types) > 0) {
    stop(paste(
      "The type specification names types that cannot be checked:",
      paste(unknown_types, collapse = ", ")
    ))
  }

  unknown_cols <- setdiff(names(types), expected_colnames)
  if (length(unknown_cols) > 0) {
    stop(paste(
      "The type specification names columns that are not expected:",
      paste(unknown_cols, collapse = ", ")
    ))
  }

  mismatches <- character(0)
  for (column_name in names(types)) {
    column <- df[[column_name]]

    # A column holding nothing but NA carries no type information: R types it
    # as logical regardless of what the column means, so checking it would
    # reject frames that are merely empty for this scenario.
    if (is.atomic(column) && length(column) > 0 && all(is.na(column))) {
      next
    }

    accepted <- alternatives[[column_name]]
    if (!any(vapply(accepted, function(type) column_type_predicates[[type]](column), logical(1)))) {
      mismatches <- c(mismatches, paste0(
        column_name, " (expected ", types[[column_name]],
        ", found ", describe_column_type(column), ")"
      ))
    }
  }

  if (length(mismatches) > 0) {
    stop(paste(
      "The dataframe has columns of the wrong type:",
      paste(mismatches, collapse = "; ")
    ))
  }
  return(invisible(NULL))
}

#' Remove columns from a dataframe
#'
#' This function removes specified columns from a dataframe.
#'
#' @param df A dataframe to modify.
#' @param to_remove A character vector of column names to remove.
#'
#' @return A dataframe with specified columns removed.
#'
#' @keywords internal
remove_columns_from_df <- function(df, to_remove) {
  for (key in to_remove) {
    if (key %in% colnames(df)) {
      df <- df[, -which(colnames(df) == key)]
    }
  }
  return(df)
}

#' Concatenate all simulation results to a global dataframe
#'
#' This function reads all CSV files with names starting with "results" from the specified
#' directory, concatenates them, and saves the result to a new CSV file.
#'
#' @param results_dir A character string specifying the path of the results folder.
#' @param ocs_filename A character string specifying the filename of the results dataframe.
#'
#' @return No return value, called for side effects.
#'
#' @importFrom dplyr bind_rows
#'
#' @keywords internal
concatenate_simulation_results <- function(results_dir, ocs_filename) {
  # Initialize an empty data frame to store the concatenated results
  df_results <- data.frame()

  if (grepl("frequentist", ocs_filename)) {
    ocs_type <- "frequentist"
  } else if (grepl("bayesian", ocs_filename)) {
    ocs_type <- "bayesian"
  } else {
    stop("Wrong filename defined in output config.")
  }
  folder_dir <- file.path(results_dir, ocs_type)

  # List all files in the directory tree
  files <- list.files(folder_dir, recursive = TRUE, full.names = TRUE)

  for (file_path in files) {
    if (grepl("\\.csv$", file_path)) {
      if (startsWith(basename(file_path), "results")) {
        # Read the CSV file
        df <- upgrade_target_size_columns(read.csv(file_path, stringsAsFactors = FALSE))

        # Concatenate to the global data frame
        df_results <- dplyr::bind_rows(df_results, df)
      }
    }
  }

  # Save the global data frame to a new CSV file
  return(readr::write_csv(df_results, file.path(results_dir, ocs_filename)))
}


#' Format simulation output as a printable table
#'
#' @param output A simulation output list (as returned by e.g.
#'   `Model$estimate_frequentist_operating_characteristics()`).
#'
#' @return A data frame with one row per metric, suitable for printing.
#'
#' @keywords internal
format_simulation_output_table <- function(output) {
  # Extract the data
  data <- list(
    "Success Probability" = c(
      output$success_proba,
      output$conf_int_success_proba_lower,
      output$conf_int_success_proba_upper
    ),
    "Coverage" = c(
      output$coverage,
      output$conf_int_coverage_lower,
      output$conf_int_coverage_upper
    ),
    "MSE" = c(
      output$mse,
      output$conf_int_mse_lower,
      output$conf_int_mse_upper
    ),
    "Bias" = c(
      output$bias,
      output$conf_int_bias_lower,
      output$conf_int_bias_upper
    ),
    "Posterior Mean" = c(
      output$posterior_mean,
      output$conf_int_posterior_mean_lower,
      output$conf_int_posterior_mean_upper
    ),
    "Posterior Median" = c(
      output$posterior_median,
      output$conf_int_posterior_median_lower,
      output$conf_int_posterior_median_upper
    ),
    "Precision" = c(
      output$precision,
      output$conf_int_precision_lower,
      output$conf_int_precision_upper
    ),
    "Interval Score" = c(
      output$interval_score,
      output$conf_int_interval_score_lower,
      output$conf_int_interval_score_upper
    ),
    "Credible Interval" = c(
      NA,
      output$credible_interval_lower,
      output$credible_interval_upper
    ),
    "ESS Moment" = c(
      output$ess_moment,
      output$conf_int_ess_moment_lower,
      output$conf_int_ess_moment_upper
    ),
    "ESS Precision" = c(
      output$ess_precision,
      output$conf_int_ess_precision_lower,
      output$conf_int_ess_precision_upper
    ),
    "ESS ELIR" = c(
      output$ess_elir,
      output$conf_int_ess_elir_lower,
      output$conf_int_ess_elir_upper
    )
  )

  # Create a data frame
  df <- data.frame(
    Metric = names(data),
    Value = sapply(data, `[`, 1),
    CI_Lower = sapply(data, `[`, 2),
    CI_Upper = sapply(data, `[`, 3)
  )

  # Format the values
  df$Value <- format(df$Value, digits = 4, scientific = FALSE)
  df$CI_Lower <- format(df$CI_Lower, digits = 4, scientific = FALSE)
  df$CI_Upper <- format(df$CI_Upper, digits = 4, scientific = FALSE)

  # Create the confidence interval column
  df$`95% CI` <- paste0("[", df$CI_Lower, ", ", df$CI_Upper, "]")

  # Select and rename columns
  df <- df[, c("Metric", "Value", "95% CI")]

  # Print the table
  return(print(df, row.names = FALSE))
}

format_parameters_to_json <- function(json_parameters, escape = FALSE) {
  json_parameters <- jsonlite::toJSON(json_parameters, pretty = TRUE)
  json_parameters <- gsub("\"", "\'", json_parameters)
  json_parameters <- as.character(json_parameters)
  if (escape == TRUE) {
    json_parameters <- paste0("[\n  ", json_parameters, "\n]")
  }
  return(json_parameters)
}

load_data <- function(results_row, type, reload_data_objects = FALSE, case_studies_config_dir = NULL) {
  if (!(type %in% c("target", "source"))){
    stop("type must either be 'target' or 'source'.")
  }

  if (reload_data_objects) {
    if (is.null(case_studies_config_dir)) {
      case_studies_config_dir <- paste0(system.file("conf/case_studies", package = "BExTE"), "/")
    }
    case_study_config <- yaml::yaml.load_file(file.path(case_studies_config_dir, paste0(results_row$case_study, ".yml")))

    source_data <- SourceData$new(
      case_study_config = case_study_config,
      source_denominator = results_row$source_denominator
    )
    if (type == "target") {
      # The time-to-event design axes each describe a different trial, so they
      # are read back like the drift. Results written before those axes existed
      # do not carry them, and such a run only ever had the default design.
      design_value <- function(name, default) {
        if (name %in% names(results_row) && !is.na(results_row[[name]][1])) {
          return(results_row[[name]][1])
        } else {
          return(default)
        }
      }

      data <- TargetDataFactory$new()
      data <- data$create(
        source_data = source_data,
        case_study_config = case_study_config,
        target_sample_size_control = results_row$target_sample_size_control,
        target_sample_size_treatment = results_row$target_sample_size_treatment,
        treatment_drift = results_row$treatment_drift,
        control_drift = results_row$control_drift,
        summary_measure_likelihood = source_data$summary_measure_likelihood,
        target_to_source_std_ratio = results_row$target_to_source_std_ratio,
        dropout_probability = design_value("dropout_probability", 0),
        event_time_distribution = design_value("event_time_distribution", "exponential"),
        treatment_delay = design_value("treatment_delay", 0)
      )
    } else if (type == "source") {
      data <- source_data
    }
  } else {
    if (type == "target") {
      data <- list(
        sample_size_control = results_row$target_sample_size_control,
        sample_size_treatment = results_row$target_sample_size_treatment,
        treatment_effect = results_row$target_treatment_effect,
        standard_deviation = results_row$target_standard_deviation,
        summary_measure_likelihood = results_row$summary_measure_likelihood,
        treatment_rate = results_row$target_treatment_rate,
        control_rate = results_row$target_control_rate
      )
    } else if (type == "source") {
      data <- list(
        treatment_effect_estimate = results_row$source_treatment_effect_estimate,
        standard_error = results_row$source_standard_error,
        sample_size_control = results_row$source_sample_size_control,
        sample_size_treatment = results_row$source_sample_size_treatment,
        equivalent_source_sample_size_per_arm = results_row$equivalent_source_sample_size_per_arm,
        treatment_rate = results_row$source_treatment_rate,
        control_rate = results_row$source_control_rate
      )

      if (is.null(data$equivalent_source_sample_size_per_arm) ||
          anyNA(data$equivalent_source_sample_size_per_arm)) {
        data$equivalent_source_sample_size_per_arm <- 2 * data$sample_size_control * data$sample_size_treatment / (data$sample_size_control + data$sample_size_treatment)
      }
    }
  }
  return(data)
}

# Function to create boolean filter
create_boolean_filter <- function(df, conditions, exclude_key = NULL) {
  if (!is.null(exclude_key)) {
    conditions <- conditions[!names(conditions) %in% exclude_key]
  }

  # Convert all columns to numeric where possible
  df <- data.frame(lapply(df, function(x)
    as.numeric(x)))


  filter <- rep(TRUE, nrow(df))
  for (col in names(conditions)) {
    condition <- conditions[[col]]
    filter <- (df[[col]] == condition) & filter
  }
  return(filter)
}


check_confidence_intervals <- function(df, metrics) {
  # Loop through each metric and check the confidence intervals
  for (metric in metrics) {
    upper_col <- paste0("conf_int_", metric, "_upper")
    lower_col <- paste0("conf_int_", metric, "_lower")

    # Check if the metric and its CI bounds exist in the dataframe
    if (all(c(metric, lower_col, upper_col) %in% names(df))) {

      # Check if the lower bound is indeed lower than the upper bound
      invalid_bounds <- which(df[[lower_col]] > df[[upper_col]])

      if (length(invalid_bounds) > 0) {
        issue <- paste0("Lower bound of the CI is larger than the upper bound for ", metric)
        warning(issue)
      }

      # Check if the metric value lies within the confidence interval
      out_of_bounds <- which(df[[metric]] < df[[lower_col]] | df[[metric]] > df[[upper_col]])
      if (length(out_of_bounds) > 0) {
        issue <- paste0("Mean value outside the CI bounds for ", metric)
        warning(issue)
      }
    }
  }
  return(invisible(NULL))
}


# Used to extract parameters as nested list
extract_nested_parameter <- function(parameters){
  # Initialize an empty list to store the nested list
  nested_list <- list()

  # Iterate over the column names of the parameters to construct the nested list
  for (colname in colnames(parameters)) {

    # Split the column name by period (".") to find the hierarchy
    split_names <- strsplit(colname, "\\.")[[1]]

    if (length(split_names) == 2) {
      if (!split_names[1] %in% names(nested_list)) {
        nested_list[[split_names[1]]] <- list()
      }

      x <- parameters[[colname]]
      x_numeric <- suppressWarnings(as.numeric(x))

      # Replace NA values with the original input
      x_clean <- ifelse(is.na(x_numeric), x, x_numeric)

      # Assign values dynamically to the second level (e.g., family or std_dev)
      nested_list[[split_names[1]]][[split_names[2]]] <- x_clean

      # For 'initial_prior', assign the value directly
    } else {
      nested_list[colname] <- parameters[[colname]]
    }
  }
  return(nested_list)
}

# Function to log errors globally
global_error_handler <- function() {
  err <- geterrmessage()  # Get the error message
  futile.logger::flog.error("Global error occurred: %s", err)

  # Optionally log other debugging info such as the call stack
  return(futile.logger::flog.error("Call stack:\n%s", paste(deparse(sys.calls()), collapse = "\n")))

  # Log additional variables of interest (if needed)
  # For instance, if you're in a loop or a function with certain variables
  # futile.logger::flog.error("Variable state at error - x: %s, y: %s", x, y)  # Customize as needed
}

generate_log_filename <- function(base_name = "error_log.log", suffix_type = "timestamp") {
  if (!file.exists(base_name)) {
    return(base_name)  # Return the base name if no file exists
  }

  # If the file exists, create a new filename with a suffix
  if (suffix_type == "timestamp") {
    # Add a timestamp to the filename
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    new_name <- sub(".log$", paste0("_", timestamp, ".log"), base_name)
  } else if (suffix_type == "counter") {
    # Use a counter to create unique filenames
    counter <- 1
    repeat {
      new_name <- sub(".log$", paste0("_", counter, ".log"), base_name)
      if (!file.exists(new_name)) break
      counter <- counter + 1
    }
  }

  return(new_name)
}

save_state <- function(iteration, scenario, worker_id = NULL, env) {
  if (is.null(worker_id)){
    file = paste0("./logs/", env, "/checkpoints/checkpoint_iter_", iteration, ".RData")
  } else {
    file = paste0("./logs/", env, "/checkpoints/checkpoint_iter_", worker_id, "_iter_", iteration, ".RData")
  }
  # Extract the directory path from the file path
  dir_path <- dirname(file)

  # Check if the directory exists, and if not, create it
  if (!file.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }

  return(save(scenario, file = file))
}


is_approx_equal <- function(x, y, tolerance = .Machine$double.eps^0.5) {
  return(abs(x - y) < tolerance)
}

# Function to compare rows, ignoring NAs
compare_ignore_na <- function(row, ref_row) {
  # Only compare non-NA elements in the reference row
  non_na_indices <- !is.na(ref_row)
  return(all(row[non_na_indices] == ref_row[non_na_indices]))
}

#' Split a total target sample size between the arms
#'
#' @param total Total target sample size.
#' @param allocation_ratio Ratio of the treatment to the control arm size.
#' @return A list with whole-number `control` and `treatment`, as doubles,
#'   since the generation cache hashes them. For a ratio of 1,
#'   each arm gets `floor(total / 2)`, as balanced designs always have.
#' @noRd
split_target_sample_size <- function(total, allocation_ratio = 1) {
  assert_single_number(total)
  assert_single_number(allocation_ratio)
  if (!(allocation_ratio > 0) || !is.finite(allocation_ratio)) {
    stop("The target allocation ratio must be a positive number.", call. = FALSE)
  }
  if (allocation_ratio == 1) {
    control <- floor(total / 2)
    treatment <- control
  } else {
    control <- round(total / (1 + allocation_ratio))
    treatment <- round(total) - control
  }
  if (control < 2 || treatment < 2) {
    stop(sprintf(
      "A total target sample size of %s with allocation ratio %s leaves fewer than 2 patients in an arm.",
      total, allocation_ratio
    ), call. = FALSE)
  }
  return(list(control = control, treatment = treatment))
}

#' Label of the target arm sizes for titles and file names
#'
#' @param n_control,n_treatment Arm sizes.
#' @param style "title" for figure and table titles, "latex" for LaTeX
#'   tables, or "file" for file names.
#' @return A character vector.
#' @noRd
target_size_label <- function(n_control, n_treatment, style = c("title", "latex", "file")) {
  style <- match.arg(style)
  return(switch(style,
    title = return(sprintf("n_C = %s, n_T = %s", n_control, n_treatment)),
    latex = return(sprintf("$n_C = %s$, $n_T = %s$", n_control, n_treatment)),
    file = return(sprintf("nC=%s_nT=%s", n_control, n_treatment))
  ))
}

#' Unwrap the length-1 list columns that rbind() over lists leaves behind
#'
#' `do.call(rbind, list_of_lists)` returns a matrix of mode list, so
#' `data.frame()` on it gives a frame whose every column is a list of
#' length-1 values rather than a vector of scalars. Columns whose cells are
#' all scalars are unwrapped into the vector they stand for; a column that
#' genuinely holds several values per row, or one wrapped in I(), is left
#' alone.
#'
#' @param df A dataframe.
#'
#' @return `df`, with its scalar list columns replaced by plain vectors.
#' @noRd
unwrap_scalar_list_columns <- function(df) {
  df[] <- lapply(df, function(column) {
    if (is.list(column) && !inherits(column, "AsIs") && all(lengths(column) == 1L)) {
      return(unlist(column, use.names = FALSE))
    } else {
      return(column)
    }
  })

  return(df)
}

## Copies of two small BExTE validators, which BExTE does not export.
check_probability_value <- function(p) {
  if (length(p) == 0) {
    return(NA)
  }

  if (is.na(p)) {
    return(NA)
  }

  if (p < 0) {
    warning("Invalid value for the probability, should be between 0 and 1. Rounding the value to 0.")
    futile.logger::flog.warn("Invalid value for the probability, should be between 0 and 1. Rounding the value to 0.")
    p <- 0
  }

  if (p > 1) {
    warning("Invalid value for the probability, should be between 0 and 1. Rounding the value to 1.")
    futile.logger::flog.warn("Invalid value for the probability, should be between 0 and 1. Rounding the value to 1.")
    p <- 1
  }
  return(p)
}

#' Check that a value is a single number
#'
#' A like-for-like replacement for assertions::assert_number() in code that
#' builds an object per scenario or per result row. It accepts exactly what
#' assert_number() accepts - one numeric value, NA and Inf included - but
#' costs about 0.4us against 190us, because it does not inspect its own call
#' or dispatch over a list of validators. The argument is only deparsed when
#' the check fails, so the message still names the offending value.
#'
#' @param x Value to check.
#' @param arg_name Name to report; defaults to the deparsed argument.
#'
#' @return `x`, invisibly. Called for the error it raises.
#' @noRd
assert_single_number <- function(x, arg_name = deparse(substitute(x))) {
  if (!is.numeric(x) || length(x) != 1L) {
    detail <- if (!is.numeric(x)) {
      paste0("class is ", class(x)[1], ", not numeric")
    } else {
      paste0("length is ", length(x), ", not 1")
    }
    stop(sprintf("'%s' is not a number! (%s)", arg_name, detail), call. = FALSE)
  }

  return(invisible(x))
}
