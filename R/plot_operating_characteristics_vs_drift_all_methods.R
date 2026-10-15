plot_element_metric_vs_xvar_methods <- function(plt,
                                        data,
                                        xvar_name,
                                        selected_metric_name,
                                        selected_metric_uncertainty_lower,
                                        selected_metric_uncertainty_upper,
                                        cap_size,
                                        markersize,
                                        labels_title,
                                        join_points,
                                        color = NULL,
                                        is_baseline = FALSE,
                                        dodging = TRUE,
                                        parameters_labels = NULL) {
  if (nrow(data) == 0) {
    stop("The dataframe is empty.")
  }


  parameters_labels <- setNames(object = parameters_labels, nm = data$parameters_labels)

  if (dodging){
    if (!is.null(selected_metric_uncertainty_lower)){
      plt <- plt +
        ggplot2::geom_errorbar(
          data = data,
          ggplot2::aes(
            x = !!xvar_name,
            y = !!selected_metric_name,
            ymin = !!selected_metric_uncertainty_lower,
            ymax = !!selected_metric_uncertainty_upper,
            color = factor(parameters_labels, levels = unique(data$parameters_labels)),
            # Used to keep the parameters ordered
            group = factor(parameters_labels, levels = unique(data$parameters_labels)),
            shape = factor(parameters_labels, levels = unique(data$parameters_labels))
          ),
          position = ggplot2::position_dodge(width = cap_size),
          width = cap_size
        )
    }
    plt <- plt +
      ggplot2::geom_point(
        data = data,
        ggplot2::aes(
          x = !!xvar_name,
          y = !!selected_metric_name,
          color = factor(parameters_labels, levels = unique(data$parameters_labels)),
          # Used to keep the parameters ordered
          group = factor(parameters_labels, levels = unique(data$parameters_labels)),
          shape = factor(parameters_labels, levels = unique(data$parameters_labels))
        ),
        position = ggplot2::position_dodge(width = cap_size),
        size = markersize,
      )
    } else {
      if (!is.null(selected_metric_uncertainty_lower)){
        plt <- plt +
          ggplot2::geom_errorbar(
            data = data,
            ggplot2::aes(
              x = !!xvar_name,
              y = !!selected_metric_name,
              ymin = !!selected_metric_uncertainty_lower,
              ymax = !!selected_metric_uncertainty_upper,
              color = factor(parameters_labels, levels = unique(data$parameters_labels)),
              # Used to keep the parameters ordered
              group = factor(parameters_labels, levels = unique(data$parameters_labels)),
              shape = factor(parameters_labels, levels = unique(data$parameters_labels))
            ),
            width = cap_size
          )
      }
      plt <- plt +
        ggplot2::geom_point(
          data = data,
          ggplot2::aes(
            x = !!xvar_name,
            y = !!selected_metric_name,
            color = factor(parameters_labels, levels = unique(data$parameters_labels)),
            # Used to keep the parameters ordered
            group = factor(parameters_labels, levels = unique(data$parameters_labels)),
            shape = factor(parameters_labels, levels = unique(data$parameters_labels))
          ),
          size = markersize
        )
  }

  plt <- plt + ggplot2::labs(color = labels_title) + ggplot2::scale_color_viridis_d(labels = parameters_labels, guide = "legend") + ggplot2::scale_fill_viridis_d(labels = parameters_labels, guide = "legend") +
    ggplot2::scale_shape_manual(values = seq_along(parameters_labels), labels = parameters_labels, guide = "legend") +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Methods"), shape = ggplot2::guide_legend(title = "Methods"))


  if (join_points == TRUE) {
    if (dodging){
      plt <- plt + ggplot2::geom_line(
        #  Add lines between points
        data = data,
        ggplot2::aes(
          x = !!xvar_name,
          y = !!selected_metric_name,
          color = factor(parameters_labels, levels = unique(data$parameters_labels)),
          # Used to keep the parameters ordered
          group = factor(parameters_labels, levels = unique(data$parameters_labels)),
        ),
        position = ggplot2::position_dodge(width = cap_size)
      )
    } else {
      plt <- plt + ggplot2::geom_line(
        #  Add lines between points
        data = data,
        ggplot2::aes(
          x = !!xvar_name,
          y = !!selected_metric_name,
          color = factor(parameters_labels, levels = unique(data$parameters_labels)),
          group = factor(parameters_labels, levels = unique(data$parameters_labels)),
        )
      )
    }
  }

  return(plt)
}


#' Function to generate a metric vs drift plot
#'
#' @param metric The metric to plot
#' @param results_metrics_df The dataframe containing the results and metrics
#' @param theta_0 The true treatment effect
#' @param case_study The case study name
#' @param target_sample_size_control,target_sample_size_treatment The target
#'   arm sizes.
#' @param xvars A list of x-variables for control drift and treatment drift
#' @param target_to_source_std_ratio Ratio between the target and source studies sampling standard deviation.
#' @param join_points Whether to join the point with a line or not.
#' @param source_denominator_change_factor The source denominator change factor of the scenario to keep.
#'
#' @return None
#'
#' @keywords internal
plot_metric_vs_drift_methods <- function(metric,
                                 results_metrics_df,
                                 theta_0,
                                 case_study,
                                 target_sample_size_control,
                                 target_sample_size_treatment,
                                 xvars,
                                 target_to_source_std_ratio = 1,
                                 join_points = FALSE,
                                 source_denominator_change_factor = 1) {
  check_plot_packages()

  if (!(metric %in% colnames(results_metrics_df))) {
    warning("Metric not in input dataframe.")
    return()
  }

  # Remove non-zero control drift
  results_df <- results_metrics_df %>% dplyr::filter(control_drift == 0)

  results_df <- results_df %>%
    dplyr::arrange(method, drift)

  xvar <- xvars$drift

  results_df <- filter_to_important_parameter_values(results_df)

  color_map <- viridis::viridis(n = 256)

  results_df$rows <- seq(1, nrow(results_df))

  # Process the row of a results dataframe to create a Method + Parameters label
  parameters_labels <- sapply(seq_len(nrow(results_df)), function(i) {
    return(process_method_parameters_label(results_df[i, ], methods_labels))
  })

  results_df$parameters_labels <- lapply(seq_len(nrow(results_df)), function(i) {
    return(process_method_parameters_label(results_df[i, ], methods_labels))
  })

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

  # Calculate the range of the x-axis
  x_range <- range(results_df[, xvar$name])
  x_width <- diff(x_range)
  # Set a relative cap size (necessary, otherwise the cap size will depend on the x-axis width)

  cap_size <- x_width * relative_error_cap_width # Adjust relative_error_cap_width to control the relative cap size

  labels_title <- "Methods"
  plt <- ggplot2::ggplot()
  plt <- plot_element_metric_vs_xvar_methods(
    plt,
    results_df,
    xvar_name,
    selected_metric_name,
    selected_metric_uncertainty_lower,
    selected_metric_uncertainty_upper,
    cap_size,
    markersize,
    labels_title,
    join_points,
    color,
    dodging = FALSE,
    parameters_labels =  parameters_labels
  )

  if(metric %in% c("tie", "success_proba")){
    plt <- plt +
      ggplot2::geom_hline(yintercept = analysis_config$nominal_tie, color = "black", linetype = "dashed") + # Add horizontal line at y = 0.05
      ggplot2::scale_y_continuous(breaks = nominal_tie_breaks(analysis_config$nominal_tie))  # Add y-tick at y = 0.05
  }

  case_study_config <- yaml::read_yaml(file.path(case_studies_config_dir, paste0(case_study, ".yml")))
  source_treatment_effect_estimate <- unique(results_df$source_treatment_effect_estimate)[1]
  # Not sorted: important_drift_values() returns the no-effect, partially
  # consistent and consistent drifts in that order, which is the order of the
  # labels below. Sorting reversed them for a negative source effect
  # (teriflunomide, mepolizumab), putting "No effect" on the consistent drift.
  mandatory_drift_values <- important_drift_values(source_treatment_effect_estimate, case_study_config)

  mandatory_drift_values <- data.frame(
    drift_value = mandatory_drift_values,
    effect_label = c("No effect", "Partially consistent effect", "Consistent effect")
  )

  plt <- plt +
    ggplot2::geom_vline(data = mandatory_drift_values,
               ggplot2::aes(xintercept = drift_value, linetype = effect_label), color = "black") +
    ggplot2::scale_linetype_manual(values = c("No effect" = "dotted",
                                   "Partially consistent effect" = "dashed",
                                   "Consistent effect" = "solid")) +
    ggplot2::labs(linetype = "Treatment effect")

  title <- sprintf(
          "%s, %s",
          stringr::str_to_title(case_study),
          target_size_label(target_sample_size_control, target_sample_size_treatment, "latex")
        )

  figure_name <- paste0(
    case_study,
    "_",
    selected_metric_name,
    "_vs_",
    xvar_name,
    "_all_methods_target_",
    target_size_label(target_sample_size_control, target_sample_size_treatment, "file")
  )

  title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)
  figure_name <- format_filename(filename = figure_name, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


  # Add colors, labels and theme
  plt <- plt +
    ggplot2::labs(x = xvar_label, y = selected_label, title = title)

  # Save plot
  directory <- paste0(figures_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  plot.size <- set_size(426)
  fig_width_in <- plot.size[1]
  fig_height_in <- plot.size[2]*1.15

  file_path <- file.path(directory, figure_name)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
    return()
  }
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf")
  return(export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png"))


}

#' Function to plot metric vs sample size
#'
#' @description This function plots a metric against the sample size for a given method, case study, and parameters combinations.
#'
#' @param metric The metric to plot.
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param case_study The case study to plot the metric for.
#' @param theta_0 Boundary of the null hypothesis space.
#' @param target_to_source_std_ratio Ratio between the target and source studies sampling standard deviation.
#' @param source_denominator_change_factor The source denominator change factor of the scenario to keep.
#' @param dodging Currently unused.
#' @param join_points Whether to join the point with a line or not.
#' @param target_treatment_effect Which of the three main treatment effects to keep:
#'   "No effect", "Partially consistent effect" or "Consistent effect".
#'
#' @return None
#'
#' @keywords internal
plot_metric_vs_sample_size_methods <- function(metric,
                                       results_metrics_df,
                                       case_study,
                                       theta_0 = 0,
                                       target_to_source_std_ratio = 1,
                                       source_denominator_change_factor = 1,
                                       dodging = FALSE,
                                       join_points = FALSE,
                                       target_treatment_effect = target_treatment_effect) {
  check_plot_packages()


  if (!(metric %in% colnames(results_metrics_df))) {
    warning("Metric not in input dataframe.")
    return()
  }
  # Remove non-zero control drift
  results_df <- results_metrics_df %>% dplyr::filter(control_drift == 0)

  # Filter results_metrics_df on the selected case study and method
  results_df <- results_metrics_df %>%
    dplyr::filter(
      case_study == !!case_study,
      matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor),
      matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio),
    )

  if (nrow(results_df) == 0) {
    warning("Dataframe is empty")
    return()
  }



  case_study_config <- yaml::read_yaml(file.path(case_studies_config_dir, paste0(case_study, ".yml")))

  # Filter on important drift value
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


    closest_values <- sort(closest_values)
    results_df <- results_df %>%
      mutate(effect_name = case_when(
        drift == closest_values[1] ~ "No effect",
        drift == closest_values[2] ~ "Partially consistent effect",
        drift == closest_values[3] ~ "Consistent effect"
      ))
    # Filter on the target_treatmnent_effect
    results_df <- results_df %>%
      dplyr::filter(effect_name == !! target_treatment_effect)
  }

  if (nrow(results_df) == 0) {
    warning("Dataframe is empty")
    return()
  }


  results_df <- filter_to_important_parameter_values(results_df)


  # Process the row of a results dataframe to create a Method + Parameters label
  parameters_labels <- sapply(seq_len(nrow(results_df)), function(i) {
    return(process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE))
  })

  results_df$parameters_labels <- lapply(seq_len(nrow(results_df)), function(i) {
    return(process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE))
  })

  ## The plotmath labels cannot be read back for the numbers they carry, and
  ## method_parameter_colors() needs those to place a value on its ramp.
  results_df$parameters_labels_not_latex <- lapply(seq_len(nrow(results_df)), function(i) {
    return(process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE, as_latex = FALSE))
  })

  results_df$methods_parameters_labels <- lapply(seq_len(nrow(results_df)), function(i) {
    return(process_method_parameters_label(results_df[i, ], methods_labels, method_name = TRUE, as_latex = FALSE))
  })

  parameters_labels <- setNames(object = parameters_labels, nm = results_df$parameters_labels)

  methods <- format_results_df_methods(results_df)
  results_df$method <- methods

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

  # The x-axis is the total target sample size.
  results_df$target_sample_size_total <- results_df$target_sample_size_control +
    results_df$target_sample_size_treatment

  # Calculate the range of the x-axis
  x_range <- range(results_df$target_sample_size_total)
  x_width <- diff(x_range)
  # Set a relative cap size (necessary, otherwise the cap size will depend on the x-axis width)
  cap_size <- x_width * relative_error_cap_width # Adjust relative_error_cap_width to control the relative cap size

  # Calculate minimum difference between x-values
  x_values <- sort(unique(results_df$target_sample_size_total))
  if (length(x_values) > 1) {
    x_diffs <- diff(x_values)
    min_x_diff <- min(x_diffs)
  } else {
    min_x_diff <- 1
  }

  # Determine number of groups
  n_groups <- length(unique(results_df$parameters_labels))

  # Set fractions for cap width and dodge width
  cap_width_fraction <- 0.5
  dodge_width_fraction <- 0.8

  # Calculate cap_width and dodge_width
  cap_width <- min_x_diff * cap_width_fraction
  total_dodge_width <- min_x_diff * dodge_width_fraction
  dodge_width <- total_dodge_width / (n_groups - 1)

  results_df <- results_df %>%
    group_by(parameters_labels) %>%
    mutate(
      jitter_offset = jitter(0, amount = 4)[1],  # Generate a single jitter offset per group
      jittered_x = target_sample_size_total + jitter_offset  # Add the same offset to all rows in the group
    ) %>%
    ungroup()

  # Split data by 'method'
  results_df_split <- split(results_df, results_df$method)

  # Unique parameters and methods
  unique_parameters <- unique(results_df$parameters_labels)
  unique_methods <- sort(unique(results_df$method))

  # Shapes for methods, keyed by method rather than handed out by position
  shapes <- method_shape_map(unique_methods)


  # Initialize the plot
  plt <- ggplot2::ggplot(mapping = ggplot2::aes(x = target_sample_size_total, y = !!selected_metric_name))

  # Iterate over each method and add layers
  plt <- plt + purrr::imap(results_df_split, function(df_subset, method_name) {
    df_subset <- df_subset[!is.na(df_subset[[selected_metric_name]]), ]


    # Get unique parameters for this method
    method_parameters <- unique(df_subset$parameters_labels)

    if (length(method_parameters) == 0){
      return()
    }

    # Shades of this method's hue, ordered by the parameter value
    method_colors <- method_parameter_color_map(df_subset, method_name)

    # Get the shape for this method
    method_shape <- shapes[method_name]

    # Create the plot layers for this method
    layers <- list(
      # Add points for this method
      ggplot2::geom_point(
        data = df_subset,
        ggplot2::aes(x = jittered_x, color = factor(parameters_labels, levels = method_parameters)),
        shape = method_shape,
        size = 2
      ),
      ggplot2::geom_line(
        data = df_subset,
        ggplot2::aes(x = jittered_x, color = factor(parameters_labels, levels = method_parameters))
      ),
      # Add vertical error bars
      ggplot2::geom_errorbar(
        data = df_subset,
        ggplot2::aes(
          x = jittered_x,
          ymin = !!ggplot2::sym(paste0("conf_int_", selected_metric_name, "_lower")),
          ymax = !!ggplot2::sym(paste0("conf_int_", selected_metric_name, "_upper")),
          color = factor(parameters_labels, levels = method_parameters)
        ),
        width = cap_width
      ),
      # Define the color scale for parameters within this method
      ggplot2::scale_color_manual(
        values = method_colors,
        labels = method_parameters,
        name = method_name,
        drop = TRUE,
        guide = ggplot2::guide_legend(override.aes = list(shape = method_shape), ncol = 1)
      )
    )


    # Reset the color scale for the next method
    layers <- c(layers, list(ggnewscale::new_scale_color()))

    return(layers)
  })

  if(metric %in% c("tie", "success_proba")){
    plt <- plt +
      ggplot2::geom_hline(yintercept = analysis_config$nominal_tie, color = "black", linetype = "dashed") + # Add horizontal line at y = 0.05
      ggplot2::scale_y_continuous(breaks = nominal_tie_breaks(analysis_config$nominal_tie))  # Add y-tick at y = 0.05
  }

  plot_title <- sprintf(
    "%s, %s",
    stringr::str_to_title(case_study),
    target_treatment_effect
  )
  plot_title <- format_title(title = plot_title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

  # Add shape scale for methods and labels
  plt <- plt +
    ggplot2::scale_shape_manual(
      values = shapes,
      name = "Methods",
      guide = ggplot2::guide_legend(order = 1, nrow = 2)
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "Target sample size (total)",
      y = selected_label
    )

  # Apply theme settings
  plt <- plt + ggplot2::theme_bw() + ggplot2::theme(
    axis.text = ggplot2::element_text(family = font, size = text_size),
    axis.text.y = ggplot2::element_text(family = font, size = small_text_size),
    axis.text.x = ggplot2::element_text(family = font, size = small_text_size),
    axis.title = ggplot2::element_text(family = font, size = text_size),
    plot.title = ggplot2::element_text(family = font, size = text_size),
    legend.text = ggplot2::element_text(family = font, size = small_text_size),
    legend.title = ggplot2::element_text(family = font, size = text_size),
    legend.key.size = grid::unit(0.1, "cm"),  # Key size
    legend.spacing.x = grid::unit(0.03, "cm"),
    legend.spacing.y = grid::unit(0.01, "cm"),  # Narrow vertical spacing
    legend.position = "bottom",
    legend.direction = "vertical",
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.box = "horizontal",             # Stack legends vertically
    legend.box.just = "center"             # Align legends to the left
  )

  plot.size <- set_size(textwidth)
  fig_width_in <- plot.size[1]*1.35
  fig_height_in <- plot.size[2]*1.5

  # plt <- plt +
  #   ggplot2::labs(title = title, x = "Target sample size (total)",
  #                 y = selected_label)

  directory <- paste0(figures_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  figure_name <- paste0(
    case_study,
    "_", target_treatment_effect,
    "_",
    selected_metric_name,
    "_vs_target_sample_size_all_methods"
  )

  figure_name <- format_filename(filename = figure_name, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)


  # plot.size <- set_size(426)
  # fig_width_in <- plot.size[1]
  # fig_height_in <- plot.size[2]

  file_path <- file.path(directory, figure_name)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
    return()
  }
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf", adjust_theme = FALSE)
  return(export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png", adjust_theme = FALSE))

}


#' Plot methods operating characteristics
#'
#' @description This function generates plots for the operating characteristics of different methods.
#'
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param metrics A list of metrics to be plotted.
#'
#' @return None
#'
#' @examples NA
plot_metric_vs_scenario_methods <- function(results_metrics_df, metrics) {
  check_plot_packages()
  check_required_colnames(
    results_metrics_df,
    required_colnames_consumer,
    context = "plot_metric_vs_scenario_methods"
  )

  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)

  treatment_effects_labels <- c("No effect", "Partially consistent effect", "Consistent effect")


  for (case_study in case_studies) {
    case_study_config <- yaml::yaml.load_file(system.file(
      paste0("conf/case_studies/", case_study, ".yml"),
      package = "BExTE"
    ))

    results_metrics_df1 <- results_metrics_df %>%
      dplyr::filter(case_study == !!case_study)

    if (nrow(results_metrics_df1) == 0) {
      warning("Dataframe is empty")
      next
    }

    theta_0 <- case_study_config$theta_0

    target_sample_sizes <- unique(results_metrics_df1[, target_size_columns, drop = FALSE])

    for (metric in names(metrics)) {
      if (!(metric %in% c("tie"))){
        for (source_denominator_change_factor in unique(results_metrics_df1$source_denominator_change_factor)) {
          results_metrics_df2 <- results_metrics_df1 %>%
            dplyr::filter(matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor))

          for (target_to_source_std_ratio in unique(results_metrics_df2$target_to_source_std_ratio)) {
            results_metrics_df3 <- results_metrics_df2 %>%
              dplyr::filter(
                matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio),
              )

            for (size_index in seq_len(nrow(target_sample_sizes))) {
              target_sample_size_control <- target_sample_sizes$target_sample_size_control[size_index]
              target_sample_size_treatment <- target_sample_sizes$target_sample_size_treatment[size_index]
              results_metrics_df4 <- results_metrics_df3 %>%
                dplyr::filter(
                  target_sample_size_control == !!target_sample_size_control,
                  target_sample_size_treatment == !!target_sample_size_treatment
                )


              if (length(unique(results_metrics_df4$method))<2){
                next
              }
              plot_metric_vs_drift_methods(
                metric = metric,
                results_metrics_df = results_metrics_df4,
                case_study = case_study,
                theta_0 = theta_0,
                source_denominator_change_factor = source_denominator_change_factor,
                target_to_source_std_ratio = target_to_source_std_ratio,
                target_sample_size_control = target_sample_size_control,
                target_sample_size_treatment = target_sample_size_treatment,
                xvars = xvars,
                join_points= FALSE
              )
          }

          for (target_treatment_effect in treatment_effects_labels){
            plot_metric_vs_sample_size_methods(
              metric = metric,
              results_metrics_df = results_metrics_df3,
              case_study = case_study,
              theta_0 = theta_0,
              target_to_source_std_ratio = target_to_source_std_ratio,
              target_treatment_effect = target_treatment_effect,
              join_points = TRUE,
              dodging = TRUE
            )
            }
          }
        }
      }
    }
  }
  return(invisible(NULL))
}
