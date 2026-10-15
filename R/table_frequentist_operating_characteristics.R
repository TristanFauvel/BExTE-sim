#' Generate a metric vs drift table
#'
#' @param metric The metric to be plotted
#' @param results_metrics_df A dataframe containing the results metrics
#' @param theta_0 The null hypothesis value
#' @param case_study The case study being analyzed
#' @param method The method being used
#' @param category The category of analysis ('parameters' or 'target_sample_size')
#' @param control_drift Boolean indicating whether to control for drift
#' @param target_sample_size_control,target_sample_size_treatment The target
#'   arm sizes
#' @param parameters_combinations The combinations of parameters
#' @param xvars A list containing x-axis variable information
#' @param source_denominator_change_factor The source denominator change factor of the scenario to keep.
#' @param target_to_source_std_ratio Ratio between the target and source studies sampling standard deviation.
#' @param wide_table Whether to pivot the table so that each parameter
#'   combination gets its own column.
#'
#' @return Generates and saves tables in HTML, PDF, and LaTeX formats
#'
#' @importFrom dplyr filter mutate select arrange rename group_by ungroup
#' @importFrom rlang sym
#' @keywords internal
table_metric_vs_drift <- function(metric,
                                  results_metrics_df,
                                  theta_0,
                                  case_study,
                                  method,
                                  category,
                                  control_drift,
                                  target_sample_size_control,
                                  target_sample_size_treatment,
                                  parameters_combinations,
                                  xvars,
                                  source_denominator_change_factor = 1,
                                  target_to_source_std_ratio = 1,
                                  wide_table = TRUE) {
  check_table_packages()
  # Filter results by case study and target sample size per arm
  results_metrics_df <- results_metrics_df %>%
    dplyr::filter(case_study == !!case_study) %>%
    dplyr::filter(method == !!method)


  if (control_drift) {
    results_metrics_df <- results_metrics_df %>% dplyr::filter(treatment_drift == 0)
    xvar <- xvars$control_drift
  } else {
    results_metrics_df <- results_metrics_df %>% dplyr::filter(control_drift == 0)
    xvar <- xvars$drift
  }

  results_df <- results_metrics_df %>% dplyr::filter(method == !!method)
  if (category == "target_sample_size") {
    results_df$target_sample_size <- target_size_label(
      results_df$target_sample_size_control, results_df$target_sample_size_treatment, "latex"
    )
  }

  categories <- unique(results_df[[category]])

  # Filter on the scenarios of main interest
  if (category == "parameters") {
    results_df <- results_df %>%
      dplyr::filter(
        target_sample_size_control == !!target_sample_size_control,
        target_sample_size_treatment == !!target_sample_size_treatment,
        matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor),
        matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio),
      )

    # Make sure that the data are sorted according to the parameters values
    params <- get_parameters(results_df[, "parameters", drop = FALSE])
    column_names <- colnames(params)
    sorted_index <- do.call(order, lapply(params[column_names], as.factor))
    # Reorder the data :
    results_df[results_df$method == method, ] <- results_df[sorted_index, ]
  } else if (category == "source_denominator") {
    results_df <- results_df %>%
      dplyr::filter(
        target_sample_size_control == !!target_sample_size_control,
        target_sample_size_treatment == !!target_sample_size_treatment,
        matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio),
      )
  } else if (category == "target_to_source_std_ratio") {
    results_df <- results_df %>%
      dplyr::filter(
        target_sample_size_control == !!target_sample_size_control,
        target_sample_size_treatment == !!target_sample_size_treatment,
        matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor)
      )
  } else if (category == "target_sample_size") {
    results_df <- results_df %>%
      dplyr::filter(
        matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio),
        matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor)
      )
  } else {
    stop("Category is not supported")
  }

  if (nrow(results_df) == 0){
    return()
  }

  results_df$rows <- seq(1, nrow(results_df))

  parameters_labels <- NA
  if (category == "target_sample_size") {
    results_df$label <- results_df$target_sample_size
    labels_title <- "Target sample size"
  } else if (category == "parameters") {
    parameters_labels <- unname(sapply(seq_len(nrow(results_df)), function(i) {
      return(process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE))
    }))

    labels_title <- "Parameters"
  } else if (category == "source_denominator") {
    results_df$label <- results_df$source_denominator
    labels_title <- "Source denominator change factor"
  } else if (category == "target_to_source_std_ratio") {
    results_df$label <- results_df$target_to_source_std_ratio
    labels_title <- latex2exp::TeX("$\\sigma_T/\\sigma_S$") # Ratio between target and source standard deviation
  } else {
    stop("Not implemented for this category")
  }

  results_df <- cbind(results_df, get_parameters(results_df[, "parameters", drop = FALSE]))

  if (category != "parameters") {
    # Filter rows that match parameters_combinations
    results_df <- filter_to_parameter_combinations(results_df, parameters_combinations)
    param_values <- unlist(unique(results_df[, category]))
  } else {
    param_values <- unlist(get_parameters(unique(results_df[, category, drop = FALSE])))
  }

  # Do not make the plot if all values for the parameters of interest are NA
  if (all(sapply(param_values, is.na))) {
    return()
  }

  if (nrow(results_df) == 0) {
    warning("Dataframe is empty")
    return()
  }

  xvar_name <- rlang::sym(xvar$name)
  xvar_label <- rlang::sym(xvar$label)

  if (metric %in% names(frequentist_metrics)) {
    selected_metric <- frequentist_metrics[[metric]]
  } else if (metric %in% names(inference_metrics)) {
    selected_metric <- inference_metrics[[metric]]
  } else if (metric %in% names(bayesian_metrics)) {
    selected_metric <- bayesian_metrics[[metric]]
  } else {
    stop("Metric is not supported")
  }

  if (is.null(selected_metric$name)) {
    stop("Metric name is NULL")
  }
  selected_metric_name <- rlang::sym(selected_metric$name)
  selected_metric_uncertainty_lower <- rlang::sym(paste0(selected_metric$metric_uncertainty, "_lower"))
  selected_metric_uncertainty_upper <- rlang::sym(paste0(selected_metric$metric_uncertainty, "_upper"))
  selected_label <- selected_metric$label

  labels <- list()
  for (i in seq_along(categories)) {
    combination <- categories[i]
    if (category == "source_denominator") {
      if (is.na(combination)) {
        next
      } else {
        results_df[results_df[, "source_denominator"] == combination, "label"] <- combination
      }
    }
  }


  if (!(category == "parameters")){
    parameters_label <- make_labels_from_parameters(parameters_combinations, method)
    parameters_str <- convert_params_to_str(methods_dict[[method]], parameters_combinations)

    if (parameters_label == "") {
      parameters_label_title <- ""
    } else {
      parameters_label_title <- paste0(", ", parameters_label)
    }

    results_df$labels <- parameters_labels
  }


  if (category == "target_sample_size") {
    title <- sprintf(
      "%s, %s, %s %s",
      frequentist_metrics[[metric]]$label,
      stringr::str_to_title(case_study),
      methods_labels[[method]]$label,
      parameters_label_title
    )

    filename <- append_parameters_str(paste0(
      case_study,
      "_",
      method,
      "_",
      selected_metric_name,
      "_vs_",
      xvar_name,
      "_cat_target_sample_size_parameters"
    ), parameters_str)
  } else if (category == "source_denominator") {
    title <- sprintf(
      "%s, %s, %s %s, %s",
      frequentist_metrics[[metric]]$label,
      stringr::str_to_title(case_study),
      methods_labels[[method]]$label,
      parameters_label_title,
      target_size_label(target_sample_size_control, target_sample_size_treatment, "latex")
    )
    filename <- append_parameters_str(paste0(
      case_study,
      "_",
      method,
      "_",
      selected_metric_name,
      "_vs_",
      xvar_name,
      "_cat_source_denominator_",
      target_size_label(target_sample_size_control, target_sample_size_treatment, "file")
    ), parameters_str)

    source_denominator_change_factor <- NA

  } else if (category == "parameters") {
    title <- sprintf(
      "%s, %s, %s, %s",
      frequentist_metrics[[metric]]$label,
      stringr::str_to_title(case_study),
      methods_labels[[method]]$label,
      target_size_label(target_sample_size_control, target_sample_size_treatment, "latex")
    )
    filename <- paste0(
      case_study,
      "_",
      method,
      "_",
      selected_metric_name,
      "_vs_",
      xvar_name,
      "_cat_parameters_",
      target_size_label(target_sample_size_control, target_sample_size_treatment, "file")
    )
  } else if (category == "target_to_source_std_ratio") {
    title <- sprintf(
      "%s, %s, %s %s, %s",
      frequentist_metrics[[metric]]$label,
      stringr::str_to_title(case_study),
      methods_labels[[method]]$label,
      parameters_label_title,
      target_size_label(target_sample_size_control, target_sample_size_treatment, "latex")
    )
    filename <- append_parameters_str(paste0(
      case_study,
      "_",
      method,
      "_",
      selected_metric_name,
      "_vs_",
      xvar_name,
      "_cat_target_to_source_std_ratio",
      "_",
      target_size_label(target_sample_size_control, target_sample_size_treatment, "file")
    ), parameters_str)

    target_to_source_std_ratio <- NA
  } else {
    stop("Not implemented for this category")
  }

  title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor, as_latex = TRUE)
  filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


  if (category == "target_sample_size") {
    results_df$label <- results_df$target_sample_size
    labels_title <- "Target sample size"
  } else if (category == "parameters") {
    results_df$label <- format_results_df_parameters(results_df, include_method_name = FALSE)
    labels_title <- "Parameters"
  } else if (category == "source_denominator") {
    results_df$label <- results_df$source_denominator
    labels_title <- "Source denominator change factor"
  } else if (category == "target_to_source_std_ratio") {
    results_df$label <- results_df$target_to_source_std_ratio
    labels_title <- expression(Tex("$\\sigma_T/\\sigma_S")) # Ratio between target and source standard deviation
  } else {
    stop("Not implemented for this category")
  }

  data_table <- results_df %>%
    dplyr::select(
      label,
      !!xvar_name,
      !!selected_metric_name,
      !!selected_metric_uncertainty_lower,
      !!selected_metric_uncertainty_upper
    ) %>%
    dplyr::mutate(
      !!xvar_label := format_num(!!xvar_name),
      !!selected_metric$label := format_num(!!selected_metric_name),
      error_low = format_num(!!selected_metric_uncertainty_lower),
      error_upper = format_num(!!selected_metric_uncertainty_upper),
      !!selected_metric$uncertainty_label := paste0("[", error_low, ", ", error_upper, "]")
    ) %>%
    dplyr::select(
      label,
      !!xvar_label,
      !!selected_metric$label,
      !!selected_metric$uncertainty_label
    ) %>%
    dplyr::rename(
      "$\\delta$" = "Drift in treatment effect",
    )


  # Merge Mean and CI columns
  data_table <- data_table %>%
    unite(!!selected_metric$label, selected_metric$label, selected_metric$uncertainty_label, sep = " ")
  if (wide_table && !(all(is.na(data_table$label) | data_table$label == ""))){
    # Pivot the table so that Parameters become columns and Sample Size per Arm is in rows
    data_table_formatted <- data_table %>%
      tidyr::pivot_wider(
        names_from = label,  # The column to pivot (label)
        values_from = !!rlang::sym(selected_metric$label)  # The values for the new columns
      )
    longtable = FALSE
  } else {
    data_table_formatted <- data_table
    longtable = TRUE
  }


  directory <- file.path(tables_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  file_path <- file.path(directory, filename)


  return(export_table(data_table = data_table_formatted, ncollapses = 1, title = title, file_path = file_path))
}

#' Plot metric vs parameters
#'
#' @param results_metrics_df A dataframe containing the results metrics
#' @param metric The metric to be plotted
#' @param case_study The case study being analyzed
#' @param method The method being used
#' @param target_sample_size_control,target_sample_size_treatment The target
#'   arm sizes
#' @param theta_0 The null hypothesis value
#' @param source_denominator_change_factor The source denominator change factor of the scenario to keep.
#' @param target_to_source_std_ratio Ratio between the target and source studies sampling standard deviation.
#'
#' @return Generates and saves tables in HTML, PDF, and LaTeX formats
#'
#' @importFrom dplyr filter mutate select arrange rename group_by ungroup
#' @importFrom rlang sym
#' @keywords internal
table_metric_vs_parameters <- function(results_metrics_df,
                                       metric,
                                       case_study = "belimumab",
                                       method = "RMP",
                                       target_sample_size_control,
                                       target_sample_size_treatment,
                                       theta_0 = NULL,
                                       source_denominator_change_factor = 1,
                                       target_to_source_std_ratio = 1) {
  check_table_packages()
  results_df <- results_metrics_df %>%
    dplyr::filter(
      case_study == !!case_study,
      target_sample_size_control == !!target_sample_size_control,
      target_sample_size_treatment == !!target_sample_size_treatment,
      method == !!method,
      matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor),
      matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio)
    )

  if (nrow(results_df) == 0) {
    warning("Dataframe is empty")
    return()
  }

  case_study_config <- yaml::read_yaml(file.path(case_studies_config_dir, paste0(case_study, ".yml")))

  source_treatment_effect_estimate <- unique(results_df$source_treatment_effect_estimate)[1]
  mandatory_drift_values <- important_drift_values(source_treatment_effect_estimate, case_study_config)

  if ("drift" %in% colnames(results_df)) {
    # Find the closest values in df$target_treatment_effect to important_target_treatment_effects
    closest_values <- sapply(mandatory_drift_values, function(x) {
      return(results_df$drift[which.min(abs(results_df$drift - x))])
    })

    # Filter the dataframe to keep only the rows with the closest values
    results_df <- results_df %>%
      dplyr::filter(drift %in% closest_values) %>%
      dplyr::arrange(drift)

    if (length(unique(results_df$drift)) != 3) {
      stop("Only three main treatment effect values should be selected")
    }
  }

  parameters_df <- get_parameters(results_df[, "parameters", drop = FALSE])

  methods_parameters <- methods_dict[[method]]

  # Effects label
  treatment_effects_names <- c("No effect",
                               "Partially consistent effect",
                               "Consistent effect")

  # Extract the metric of interest
  if (metric %in% names(frequentist_metrics)) {
    selected_metric <- frequentist_metrics[[metric]]
  } else if (metric %in% names(inference_metrics)) {
    selected_metric <- inference_metrics[[metric]]
  } else if (metric %in% names(bayesian_metrics)) {
    selected_metric <- bayesian_metrics[[metric]]
  } else {
    stop("Metric is not supported")
  }

  if (is.null(selected_metric$name)) {
    stop("Metric name is NULL")
  }

  selected_metric_name <- rlang::sym(selected_metric$name)
  selected_metric_uncertainty_lower <- rlang::sym(paste0(selected_metric$metric_uncertainty, "_lower"))
  selected_metric_uncertainty_upper <- rlang::sym(paste0(selected_metric$metric_uncertainty, "_upper"))
  selected_label <- selected_metric$label

  for (i in seq_len(ncol(parameters_df))) {
    # Select the other parameters
    other_parameters <- unique(parameters_df[,-i])

    if (!has_other_parameters(other_parameters)){
      n_params_loop = 1
    } else {
      n_params_loop = nrow(other_parameters)
    }
    for (j in seq(n_params_loop)){
      if (!has_other_parameters(other_parameters)){
        other_params_label <- ""
        other_params_str <- ""
        parameters_subdf <- parameters_df
        results_subdf <- results_df
      } else {
        other_params_label <- make_labels_from_parameters(parameter = other_parameters[j,], method = method)
        other_params_label <- paste0(other_params_label, ", ")
        other_params_str <- convert_params_to_str(methods_dict[[method]], other_parameters[j,])

        # Filter the dataframe on the value of these other parameters
        matching_filter <- apply(parameters_df[,-i], 1, function(row) all(row == other_parameters[j,]))
        parameters_subdf <- parameters_df[matching_filter,]

        results_subdf <- results_df[matching_filter,]
      }

      if (length(unique(parameters_subdf[, i])) < 2) {
        next
      } else {
        parameter_values <- as.numeric(parameters_subdf[, i])
      }

      results_subdf$parameter_values <- parameter_values

      data_table <- results_subdf %>%
        dplyr::mutate(
          means = format_num(!!rlang::sym(selected_metric$name)),
          error_low = format_num(!!rlang::sym(
            paste0(selected_metric$metric_uncertainty, "_lower")
          )),
          error_upper = format_num(!!rlang::sym(
            paste0(selected_metric$metric_uncertainty, "_upper")
          )),
          CI = paste0("[", error_low, ", ", error_upper, "]")
        ) %>%
        dplyr::select(parameter_values, means, CI) %>%
        dplyr::arrange(parameter_values)

      data_table <- data_table %>%
        # dplyr::group_by(parameter_values) %>%
        dplyr::rename(
          "Parameter" = parameter_values,
          !!selected_metric$label := means,
          !!selected_metric$uncertainty_label := CI
        )

      data_table <- data_table %>%
        unite(!!selected_metric$label, selected_metric$label, selected_metric$uncertainty_label, sep = " ")


      directory <- file.path(tables_dir, case_study)
      if (!dir.exists(directory)) {
        dir.create(directory,
                   showWarnings = FALSE,
                   recursive = TRUE)
      }

      filename <- paste0(
        case_study,
        "_",
        method,
        "_metrics_vs_parameters_",
        target_size_label(target_sample_size_control, target_sample_size_treatment, "file")
      )

      if (other_params_str != ""){
        filename <- paste0(filename, "_", other_params_str)
      }

      title <- sprintf(
        "%s, %s, %s%s",
        stringr::str_to_title(case_study),
        methods_labels[[method]]$label,
        other_params_label,
        target_size_label(target_sample_size_control, target_sample_size_treatment, "latex")
      )
      title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)
      filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

      file_path <- file.path(directory, filename)

      export_table(data_table = data_table, ncollapses = 1, title = title, file_path = file_path)

    }
  }
  return(invisible(NULL))
}

#' Plot metric vs sample size
#'
#' @param metric The metric to be plotted
#' @param results_metrics_df A dataframe containing the results metrics
#' @param case_study The case study being analyzed
#' @param method The method being used
#' @param theta_0 The null hypothesis value
#' @param source_denominator_change_factor The source denominator change factor of the scenario to keep.
#' @param target_to_source_std_ratio Ratio between the target and source studies sampling standard deviation.
#' @param wide_table Whether to pivot the table so that each parameter
#'   combination gets its own column.
#'
#' @return Generates and saves tables in HTML, PDF, and LaTeX formats
#'
#' @importFrom dplyr filter mutate select arrange rename group_by ungroup
#' @importFrom rlang sym
#' @keywords internal
table_metric_vs_sample_size <- function(metric,
                                        results_metrics_df,
                                        case_study,
                                        method,
                                        theta_0 = 0,
                                        source_denominator_change_factor = 1,
                                        target_to_source_std_ratio = 1,
                                        wide_table = TRUE) {
  check_table_packages()



  results_df <- results_metrics_df %>%
    dplyr::filter(
      case_study == !!case_study,
      method == !!method,
      matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor),
      matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio)
    )

  results_df <- cbind(results_df, get_parameters(results_df[, "parameters", drop = FALSE]))

  case_study_config <- yaml::read_yaml(file.path(case_studies_config_dir, paste0(case_study, ".yml")))

  source_treatment_effect_estimate <- unique(results_df$source_treatment_effect_estimate)[1]
  mandatory_drift_values <- important_drift_values(source_treatment_effect_estimate, case_study_config)

  if ("drift" %in% colnames(results_df)) {
    # Find the closest values in df$target_treatment_effect to important_target_treatment_effects
    closest_values <- sapply(mandatory_drift_values, function(x) {
      return(results_df$drift[which.min(abs(results_df$drift - x))])
    })

    # Filter the dataframe to keep only the rows with the closest values
    results_df <- results_df %>%
      dplyr::filter(drift %in% closest_values) %>%
      dplyr::arrange(drift)

    if (length(unique(results_df$drift)) != 3) {
      stop("Only three main treatment effect values should be selected")
    }
  }

  # Effects label
  effects <- c("No effect",
               "Partially consistent effect",
               "Consistent effect")

  # Factorize target treatment effect
  results_df$target_treatment_effect <- factor(
    results_df$target_treatment_effect,
    levels = unique(results_df$target_treatment_effect),
    labels = effects
  )


  # Extract the metric of interest
  if (metric %in% names(frequentist_metrics)) {
    selected_metric <- frequentist_metrics[[metric]]
  } else if (metric %in% names(inference_metrics)) {
    selected_metric <- inference_metrics[[metric]]
  } else if (metric %in% names(bayesian_metrics)) {
    selected_metric <- bayesian_metrics[[metric]]
  } else {
    stop("Metric is not supported")
  }

  if (is.null(selected_metric$name)) {
    stop("Metric name is NULL")
  }

  results_df$Parameters <- format_results_df_parameters(results_df, include_method_name = FALSE)

  for (treatment_effect in effects){
    # Filter on the treatment effect
    results_df_filtered <- results_df %>%
      dplyr::filter(target_treatment_effect %in% treatment_effect)

    data_table <- results_df_filtered %>%
      dplyr::mutate(
        means = format_num(!!rlang::sym(selected_metric$name)),
        error_low = format_num(!!rlang::sym(
          paste0(selected_metric$metric_uncertainty, "_lower")
        )),
        error_upper = format_num(!!rlang::sym(
          paste0(selected_metric$metric_uncertainty, "_upper")
        )),
        CI = paste0("[", error_low, ", ", error_upper, "]"),
        combined_mean_CI = paste0(means, " ", CI)  # Combine mean and CI into one column
      ) %>%
      dplyr::arrange(target_sample_size_control + target_sample_size_treatment,
                     target_sample_size_control) %>%
      dplyr::select(Parameters, target_sample_size_control, target_sample_size_treatment, combined_mean_CI)

    data_table <- data_table %>%
      dplyr::rename(
        "$n_C$" = target_sample_size_control,
        "$n_T$" = target_sample_size_treatment,
        !!selected_metric$label := combined_mean_CI  # Rename the combined column
      )

    if (wide_table && !(all(is.na(data_table$Parameters) | data_table$Parameters == ""))){
      # Pivot the table so that Parameters become columns and Sample Size per Arm is in rows
      data_table_formatted <- data_table %>%
        tidyr::pivot_wider(
          names_from = Parameters,  # The column to pivot (Parameters)
          values_from = !!rlang::sym(selected_metric$label)  # The values for the new columns
        )
      longtable = FALSE
    } else {
      data_table_formatted <- data_table
      longtable = TRUE
    }


    directory <- file.path(tables_dir, case_study)
    if (!dir.exists(directory)) {
      dir.create(directory, showWarnings = FALSE, recursive = TRUE)
    }
    filename <- paste0(
      case_study,
      "_",
      method,
      "_",
      frequentist_metrics[[metric]]$name,
      "_vs_target_sample_size_",
      gsub(" ", "_", treatment_effect)
    )

    title <-  latex2exp::TeX(
      sprintf(
        "%s, %s, %s, %s",
        frequentist_metrics[[metric]]$label,
        stringr::str_to_title(case_study),
        methods_labels[[method]]$label,
        treatment_effect
      )
    )

    title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor, as_latex = TRUE)
    filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


    file_path <- file.path(directory, filename)

    export_table(data_table = data_table_formatted, ncollapses = 1, title = title, file_path = file_path, longtable = longtable)
  }
  return(invisible(NULL))
}


#' Generate tables for methods operating characteristics
#'
#' @param results_metrics_df A dataframe containing the results metrics
#' @param metrics A list of metrics to be analyzed
#'
#' @return Generates and saves multiple tables for different metrics and methods
#'
#' @importFrom dplyr filter
#' @export
tables_metric_vs_scenario <- function(results_metrics_df, metrics) {
  check_table_packages()
  case_studies <- unique(results_metrics_df$case_study)

  # Get the list of methods
  methods <- unique(results_metrics_df$method)

  for (case_study in case_studies) {
    case_study_config <- yaml::yaml.load_file(system.file(
      paste0("conf/case_studies/", case_study, ".yml"),
      package = "BExTE"
    ))

    filtered_results_metrics_df <- results_metrics_df %>%
      dplyr::filter(case_study == !!case_study)

    if (nrow(filtered_results_metrics_df) == 0) {
      warning("Dataframe is empty")
      next
    }

    theta_0 <- case_study_config$theta_0

    target_sample_sizes <- unique(filtered_results_metrics_df[, target_size_columns, drop = FALSE])

    for (method in methods) {
      if (sum(filtered_results_metrics_df$method == method) == 0) {
        next
      }
      parameters_combinations <- unique(get_parameters(filtered_results_metrics_df[filtered_results_metrics_df$method == method, "parameters", drop = FALSE]))


      for (metric in names(metrics)) {
        for (size_index in seq_len(nrow(target_sample_sizes))) {

          # table_metric_vs_drift(
          #   metric = metric,
          #   results_metrics_df = filtered_results_metrics_df,
          #   case_study = case_study,
          #   method = method,
          #   theta_0 = theta_0,
          #   category = "target_to_source_std_ratio",
          #   control_drift = FALSE,
          #   target_sample_size_control = target_sample_sizes$target_sample_size_control[size_index],
          #   target_sample_size_treatment = target_sample_sizes$target_sample_size_treatment[size_index],
          #   parameters_combinations = NULL,
          #   xvars = xvars
          # )
        }

        for (source_denominator_change_factor in unique(filtered_results_metrics_df$source_denominator_change_factor)) {
          for (target_to_source_std_ratio in unique(filtered_results_metrics_df$target_to_source_std_ratio)) {
            for (size_index in seq_len(nrow(target_sample_sizes))) {
              for (i in seq_len(nrow(parameters_combinations))) {
                table_metric_vs_drift(
                  metric = metric,
                  results_metrics_df = filtered_results_metrics_df,
                  case_study = case_study,
                  method = method,
                  theta_0 = theta_0,
                  category = "parameters",
                  control_drift = FALSE,
                  target_sample_size_control = target_sample_sizes$target_sample_size_control[size_index],
                  target_sample_size_treatment = target_sample_sizes$target_sample_size_treatment[size_index],
                  parameters_combinations = NULL,
                  xvars = xvars
                )
              }
            }

            table_metric_vs_sample_size(
              metric = metric,
              results_metrics_df = filtered_results_metrics_df,
              case_study = case_study,
              method = method,
              theta_0 = theta_0,
              wide_table = TRUE
            )
          }
        }
      }
    }
  }
  return(invisible(NULL))
}
