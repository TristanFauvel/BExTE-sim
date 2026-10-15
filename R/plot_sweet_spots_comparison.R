#' Create a forest plot
#'
#' @description This function creates a forest plot using the provided data.
#'
#' @param data The data for creating the forest plot.
#' @param title The title of the forest plot.
#' @param x_metric_name Name of the metric on the x-axis
#' @param sweet_spot_lower Lower limit of the metric on the x-axis
#' @param sweet_spot_upper Upper limit of the metric on the x-axis
#' @param x_metric_label Label of the metric on the x-axis
#'
#' @return A ggplot2 object representing the forest plot.
#' @keywords internal
forest_subplot_sweet_spot <- function(data,
                                      title,
                                      sweet_spot_lower,
                                      sweet_spot_upper,
                                      x_range,
                                      x_metric_label,
                                      x_metric_name,
                                      methods_labels,
                                      legend = FALSE,
                                      case_study,
                                      sort_by = "values") {

  # Unnest the sweet spot bounds
  data <- data %>%
    unnest(cols = c(!!rlang::sym(sweet_spot_lower), !!rlang::sym(sweet_spot_upper)))


  tryCatch({
    # Convert to numeric
    data[[sweet_spot_lower]] <- as.numeric(unlist(data[[sweet_spot_lower]]))
    data[[sweet_spot_upper]] <- as.numeric(unlist(data[[sweet_spot_upper]]))
  }, error = function(e) {
    message("Error occurred: ", e$message)
    message("Attempting to fix by unnesting columns...")

    # Attempt unnesting and retry
    data <- data %>%
      unnest(cols = c(!!rlang::sym(sweet_spot_lower), !!rlang::sym(sweet_spot_upper)))
    return(invisible(data))

  })

  data <- data %>%
    unnest(cols = c(!!rlang::sym(sweet_spot_lower), !!rlang::sym(sweet_spot_upper)))

  data <- data %>%
    mutate(
      label = purrr::map_chr(seq_len(nrow(.)), function(i) {
        return(as.character(process_method_parameters_label(data[i, ], methods_labels)))
      })
    )

  # Group by label (method/parameter combination)
  data_grouped <- data %>%
    group_by(label)

  # Assign row numbers based on sorting criteria
  if (sort_by == "values") {
    # Calculate total width per group
    widths <- data_grouped %>%
      summarize(
        total_width = sum(!!rlang::sym(sweet_spot_upper) - !!rlang::sym(sweet_spot_lower), na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      arrange(total_width) %>%
      mutate(row = row_number())
  } else if (sort_by == "methods_parameters") {
    # Assign row numbers based on methods/parameters
    widths <- data_grouped %>%
      summarize() %>%
      arrange(label) %>%
      mutate(row = row_number())
  } else {
    # Default assignment if sort_by is not recognized
    widths <- data_grouped %>%
      summarize() %>%
      mutate(row = row_number())
  }

  # Join back to data to get 'row' assigned to each label
  data <- data %>%
    left_join(widths %>% select(label, row), by = "label")

  # Parse labels into expressions for plotting
  labels <- parse(text = widths$label)

  # Build the plot
  plt <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      y = row,
      xmin = !!rlang::sym(sweet_spot_lower),
      xmax = !!rlang::sym(sweet_spot_upper)
    )
  ) +
    ggplot2::geom_errorbarh(height = 0.2, linewidth = 1) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(size = 12),
      text = ggplot2::element_text(size = 14),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::scale_y_continuous(breaks = widths$row, labels = labels)

  case_study_config <- yaml::read_yaml(file.path(case_studies_config_dir, paste0(case_study, ".yml")))
  source_treatment_effect_estimate <- unique(data$source_treatment_effect_estimate)[1]
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

 return(plt)
}


#' Generate a forest plot
#'
#' @description This function generates a forest plot based on the provided results dataframe, selected case study, and selected target arm sizes.
#'
#' @param results_freq_df The results dataframe.
#' @param x_metric Metric on the x-axis
#'
#' @return None
forest_plot_sweet_spot <- function(results_freq_df, x_metric) {
  selected_case_study <- unique(results_freq_df$case_study)
  selected_target_sizes <- unique(results_freq_df[, target_size_columns, drop = FALSE])


  directory <- file.path(figures_dir, selected_case_study)

  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }


  if (x_metric %in% names(frequentist_metrics)) {
    metric <- frequentist_metrics[[x_metric]]

    if (metric$label == "TIE"){
      return()
    }

    metric$label = paste0("Sweet spot relative to ", metric$label)
  } else if (x_metric %in% names(inference_metrics)) {
    metric <- inference_metrics[[x_metric]]
  } else if (x_metric %in% names(bayesian_metrics)) {
    metric <- bayesian_metrics[[x_metric]]
  } else {
    if (x_metric == "success_proba_smaller_than_nominal_TIE"){
      metric = frequentist_metrics[["tie"]]
      metric$name <- "success_proba_smaller_than_nominal_TIE"
      metric$label <- "Drift values for which Pr(Study success) < nominal TIE"
    } else if (x_metric == "power_larger_than_nominal"){
      metric = frequentist_metrics[["power"]]
      metric$name <- "power_larger_than_nominal"
      metric$label <- "Drift values for which Pr(Study success) > nominal power"
    } else {
      stop("Metric is not supported")
    }
  }

  x_metric_name <- rlang::sym(paste0("sweet_spot_", metric$name))
  filename <- paste0(
    selected_case_study,
    "_",
    x_metric_name,
    "_forest_plot_target_",
    target_size_label(
      selected_target_sizes$target_sample_size_control,
      selected_target_sizes$target_sample_size_treatment,
      "file"
    )
  )

  case_study <- unique(results_freq_df$case_study)
  target_to_source_std_ratio <- unique(results_freq_df$target_to_source_std_ratio)
  source_denominator_change_factor <- unique(results_freq_df$source_denominator_change_factor)

  filename <- format_filename(
    filename = filename,
    case_study = case_study,
    target_to_source_std_ratio = target_to_source_std_ratio,
    source_denominator_change_factor = source_denominator_change_factor
  )

  file_path <- file.path(directory, filename)
  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
    return()
  }

  id_cols <- c("method", "parameters", target_size_columns, "control_drift", "source_denominator", "source_denominator_change_factor", "case_study", "sampling_approximation", "source_treatment_effect_estimate",  "target_to_source_std_ratio", "theta_0", "null_space", "summary_measure_likelihood", "source_standard_error", "source_sample_size_control", "source_sample_size_treatment", "equivalent_source_sample_size_per_arm", "endpoint", "source_control_rate", "source_treatment_rate")

  results_freq_df <- results_freq_df %>%
    pivot_wider(id_cols = id_cols, names_from = metric, values_from = c("sweet_spot_lower", "sweet_spot_upper", "sweet_spot_width", "drift_range_upper", "drift_range_lower"))

  if (!(paste0("sweet_spot_width_", metric$name) %in% colnames(results_freq_df))) {
    warning("Metric not in input dataframe.")
    return()
  }

  sweet_spot_lower <- rlang::sym(paste0("sweet_spot_lower_", metric$name))
  sweet_spot_upper <- rlang::sym(paste0("sweet_spot_upper_", metric$name))
  x_metric_label <- metric$label

  df <- results_freq_df %>%
    dplyr::mutate(
      lower_bound = !!sweet_spot_lower,
      lower_bound = !!sweet_spot_upper
    )

  title <- sprintf(
    "%s, %s, %s",
    x_metric_label,
    stringr::str_to_title(unique(df$case_study)),
    target_size_label(
      unique(df$target_sample_size_control),
      unique(df$target_sample_size_treatment),
      "latex"
    )
  )

  title <- format_title(title = title, case_study = unique(df$case_study), target_to_source_std_ratio = unique(df$target_to_source_std_ratio), source_denominator_change_factor = unique(df$source_denominator_change_factor))

  x_range <- c(unique(df[,paste0("drift_range_lower_", metric$name)]), unique(df[,paste0("drift_range_upper_", metric$name)]))

  plt <- forest_subplot_sweet_spot(
    df,
    title,
    sweet_spot_lower,
    sweet_spot_upper,
    x_range = x_range,
    x_metric_label = x_metric_label,
    x_metric_name = x_metric_name,
    methods_labels = methods_labels,
    case_study = selected_case_study,
    legend = FALSE
  )

  # Note that adding the Tikz export of the export_plots function does not work.
  if (plots_to_latex == TRUE) {
    tikzDevice::tikz(file = paste0(file_path, ".tex"), width = 5.92)
  }

  if (plots_to_latex == TRUE) {
    dev.off()
  }

  plot_size <- set_size(textwidth)
  fig_width_in <- plot_size[1]
  fig_height_in <- plot_size[2]*1.2

  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf")
  return(export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png"))
}


#' Plot methods comparison
#'
#' @description This function plots the comparison of methods based on the provided results metrics dataframe.
#'
#' @param sweet_spots_df The results metrics dataframe.
#' @param metrics Metrics for which to create forest plots.
#'
#' @return None
forest_plot_sweet_spots_comparison <- function(sweet_spots_df, metrics) {
  check_plot_packages()

  case_studies <- unique(sweet_spots_df$case_study)


  metrics <- unique(sweet_spots_df$metric)
  for (metric in metrics) {
    for (case_study in case_studies) {
      sweet_spots_df1 <- sweet_spots_df[sweet_spots_df$case_study == case_study,]

      target_sample_sizes <- unique(sweet_spots_df1[, target_size_columns, drop = FALSE])
      for (size_index in seq_len(nrow(target_sample_sizes))) {
        sweet_spots_df2 <- dplyr::semi_join(
          sweet_spots_df1, target_sample_sizes[size_index, , drop = FALSE], by = target_size_columns
        )

        target_to_source_std_ratios <- unique(sweet_spots_df2$target_to_source_std_ratio)

        for (target_to_source_std_ratio in target_to_source_std_ratios) {
          sweet_spots_df3 <- sweet_spots_df2 %>%
            dplyr::filter(
              matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio)
            )

          source_denominator_change_factors <- unique(sweet_spots_df3$source_denominator_change_factor)

          for (source_denominator_change_factor in source_denominator_change_factors){
            sweet_spots_df4 <- sweet_spots_df3 %>%
              dplyr::filter(
                matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor)
              )
            forest_plot_sweet_spot(sweet_spots_df4, metric)
          }
        }
      }
    }
  }
  return(invisible(NULL))
}


plot_sweet_spot_width_metric_vs_sample_size <- function(metric,
                                                 sweet_spot_df,
                                                 case_study,
                                                 theta_0 = 0,
                                                 target_to_source_std_ratio = 1,
                                                 source_denominator_change_factor = 1,
                                                 dodging = FALSE,
                                                 join_points = FALSE) {

    metric_name <- metric$name
    sweet_spot_df <- sweet_spot_df %>% dplyr::filter(metric == metric_name)

    # Remove unused columns
    sweet_spot_df <- sweet_spot_df[, !colnames(sweet_spot_df) %in% c("sweet_spot_lower", "sweet_spot_upper", "sweet_spot_width")]

    # Remove non-zero control drift
    results_df <- sweet_spot_df %>% dplyr::filter(control_drift == 0)

    # Filter sweet_spot_df on the selected case study and method
    results_df <- sweet_spot_df %>%
      dplyr::filter(
        case_study == !!case_study,
        matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor),
        matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio),
      )

    if (nrow(results_df) == 0) {
      warning("Dataframe is empty")
      return()
    }


    # Unnest the sweet spot bounds
    results_df <- results_df %>%
      unnest(cols = "sweet_spot_total_width")

    tryCatch({
      # Convert to numeric
      results_df[["sweet_spot_total_width"]] <- as.numeric(unlist(results_df[["sweet_spot_total_width"]]))
      results_df[["sweet_spot_total_width"]] <- as.numeric(unlist(results_df[["sweet_spot_total_width"]]))
    }, error = function(e) {
      message("Error occurred: ", e$message)
      message("Attempting to fix by unnesting columns...")

      # Attempt unnesting and retry
      results_df <- results_df %>%
        unnest(cols = "sweet_spot_total_width")
      return(invisible(results_df))
    })

    # Keep only unique rows
    results_df <- unique(results_df)

    case_study_config <- yaml::read_yaml(file.path(case_studies_config_dir, paste0(case_study, ".yml")))

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

    dodge_width_fraction <- 0.8

    total_dodge_width <- min_x_diff * dodge_width_fraction
    dodge_width <- total_dodge_width / (n_groups - 1)

    selected_metric_name <- "sweet_spot_total_width"

    results_df_split <- split(results_df, results_df$method)
    unique_methods <- sort(unique(results_df$method))
    shapes <- method_shape_map(unique_methods)


    title <- sprintf(
      "%s",
      stringr::str_to_title(case_study)
    )

    title <- format_title(title = title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

    # Initialize the plot
    plt <- ggplot2::ggplot(mapping = ggplot2::aes_string(x = "target_sample_size_total", y = "sweet_spot_total_width"))

    plt <- plt + purrr::imap(results_df_split, function(df_subset, method_name) {
      df_subset <- df_subset[!is.na(df_subset[["sweet_spot_total_width"]]), ]

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
          ggplot2::aes(color = factor(parameters_labels, levels = method_parameters)),
          shape = method_shape,
          size = 2
          #position = position_jitter(width = 4, height = 0)
        ),
        ggplot2::geom_line(
          data = df_subset,
          ggplot2::aes(color = factor(parameters_labels, levels = method_parameters))
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

    # Add shape scale for methods and labels
    plt <- plt +
      ggplot2::scale_shape_manual(
        values = shapes,
        name = "Methods",
        guide = ggplot2::guide_legend(order = 1, nrow = 2, byrow = TRUE)
      ) +
      ggplot2::labs(
        title =  title,
        x = "Target sample size (total)",
        y = metric$label
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
      legend.spacing.x = grid::unit(0.01, "cm"),
      legend.spacing.y = grid::unit(0.01, "cm"),  # Narrow vertical spacing
      legend.position = "bottom",
      legend.direction = "vertical",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.box = "horizontal"             # Stack legends vertically
    )

    plot.size <- set_size(textwidth)
    fig_width_in <- plot.size[1]*1.35
    fig_height_in <- plot.size[2]*1.5

    directory <- paste0(figures_dir, case_study)
    if (!dir.exists(directory)) {
      dir.create(directory, showWarnings = FALSE, recursive = TRUE)
    }

    figure_name <- paste0(
      case_study,
      "_sweet_spot_width_",
      metric_name,
      "_vs_target_sample_size_all_methods"
    )

    figure_name <- format_filename(filename = figure_name, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

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
#' @param sweet_spot_df The data frame containing the results and metrics.
#' @param sweet_spots_metrics A list of metrics to be plotted.
#'
#' @return None
#'
#' @examples NA
plot_sweet_spot_width_vs_scenario <- function(sweet_spot_df, sweet_spots_metrics) {
  check_plot_packages()
  # Get the list of case studies
  case_studies <- unique(sweet_spot_df$case_study)

  for (case_study in case_studies) {
    case_study_config <- yaml::yaml.load_file(system.file(
      paste0("conf/case_studies/", case_study, ".yml"),
      package = "BExTE"
    ))

    sweet_spot_df1 <- sweet_spot_df %>%
      dplyr::filter(case_study == !!case_study)

    if (nrow(sweet_spot_df1) == 0) {
      warning("Dataframe is empty")
      next
    }

    theta_0 <- case_study_config$theta_0

    for (metric in sweet_spots_metrics) {
      for (source_denominator_change_factor in unique(sweet_spot_df1$source_denominator_change_factor)) {
        sweet_spot_df2 <- sweet_spot_df1 %>%
          dplyr::filter(matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor))

        for (target_to_source_std_ratio in unique(sweet_spot_df2$target_to_source_std_ratio)) {
          sweet_spot_df3 <- sweet_spot_df2 %>%
            dplyr::filter(
              matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio),
            )

          plot_sweet_spot_width_metric_vs_sample_size(
            metric = metric,
            sweet_spot_df = sweet_spot_df3,
            case_study = case_study,
            theta_0 = theta_0,
            target_to_source_std_ratio = target_to_source_std_ratio,
            join_points = TRUE,
            dodging = TRUE
          )
        }
      }
    }
  }
  return(invisible(NULL))
}
