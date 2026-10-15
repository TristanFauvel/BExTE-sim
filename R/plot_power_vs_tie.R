#' Function to generate a power vs tie plot
#'
#' @param results_metrics_df The dataframe containing the results and metrics
#' @param case_study The case study name
#' @param target_sample_size_control,target_sample_size_treatment The target
#'   arm sizes.
#' @param treatment_effect The treatment effect type ("consistent", "no_effect", "partially_consistent")
#' @param power_difference Logical indicating whether to calculate power difference
#' @param metrics The metrics to include in the plot
#' @param source_denominator_change_factor The source denominator change factor of the scenario to keep.
#' @param target_to_source_std_ratio Ratio between the target and source studies sampling standard deviation.
#'
#' @return None
#'
#' @keywords internal
power_vs_tie <- function(results_metrics_df,
                         case_study,
                         target_sample_size_control,
                         target_sample_size_treatment,
                         treatment_effect,
                         power_difference,
                         metrics,
                         source_denominator_change_factor,
                         target_to_source_std_ratio) {
  check_plot_packages()

  directory <- paste0(figures_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  results_df <- results_metrics_df %>%
    dplyr::filter(
      target_sample_size_control == !!target_sample_size_control,
      target_sample_size_treatment == !!target_sample_size_treatment,
      case_study == !!case_study,
      matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor),
      matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio)
    )

  labels <- format_results_df_parameters(results_df)
  methods <- format_results_df_methods(results_df)

  results_df$label <- labels
  results_df$method <- methods

  theta_0 <- unique(results_df$theta_0)
  results_df_tie <- results_df %>% dplyr::filter(target_treatment_effect == theta_0)

  results_df <- filter_to_treatment_effect(results_df, treatment_effect)

  if (power_difference) {
    results_df <- add_power_difference_columns(results_df)
    power_label <- "power_difference"
  } else {
    power_label <- "power"
  }

  filename <- paste0(
    case_study,
    "_",
    power_label,
    "_",
    target_size_label(target_sample_size_control, target_sample_size_treatment, "file"),
    "_treatment_effect=",
    treatment_effect
  )

  filename <- format_filename(filename = filename, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

  file_path <- file.path(directory, filename)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE){
    return()
  }


  if (power_difference == FALSE) {
    ymin <- 0
    ymax <- 1
  } else if (power_difference == TRUE) {
    y_range <- sort(range(results_df[, "success_proba"]))
    ymin <- y_range[1]
    ymax <- y_range[2]
  }
  # xmin <- 0
  # xmax <- 1
  x_range <- sort(range(results_df[, "tie"]))
  xmin <- x_range[1]
  xmax <- x_range[2]

  x_width <- xmax - xmin
  y_width <- ymax - ymin

  y_cap_size <- x_width * relative_error_cap_width / 2 # Adjust relative_error_cap_width to control the relative cap size
  x_cap_size <- y_width * relative_error_cap_width / 2 # Adjust relative_error_cap_width to control the relative cap size

  cap_size <- min(x_cap_size, y_cap_size)
  x_cap_size <- cap_size
  y_cap_size <- cap_size

  treatment_effects_labels <- list(
    no_effect = "No effect",
    partially_consistent = "Partially consistent effect",
    consistent = "Consistent effect"
  )

  plot_title <- sprintf(
      "%s, %s, %s",
      stringr::str_to_title(case_study),
      target_size_label(target_sample_size_control, target_sample_size_treatment, "latex"),
      treatment_effects_labels[[treatment_effect]]
    )

  plot_title <- format_title(title = plot_title, case_study = case_study, target_to_source_std_ratio = target_to_source_std_ratio, source_denominator_change_factor = source_denominator_change_factor)

  plt <- ggplot2::ggplot(results_df, ggplot2::aes(x = tie)) +
    # coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
    ggplot2::geom_point(ggplot2::aes(
      y = success_proba,
      color = as.factor(method),
      shape = as.factor(method)
    )) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = conf_int_success_proba_lower,
        ymax = conf_int_success_proba_upper,
        color =  as.factor(method)
      ),
      width = y_cap_size
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = conf_int_tie_lower,
        xmax = conf_int_tie_upper,
        y = success_proba,
        color = as.factor(method)
      ),
      height = x_cap_size
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "TIE",
      y = ifelse(
        power_difference,
        "Power difference compared to a frequentist test",
        "Power"
      ),
      color = "Methods",
      # Title for the color legend
      shape = "Methods" # Title for the shape legend
    )

  if (power_difference == FALSE) {
    # Add the nominal frequentist power and power at equivalent TIE
    plt <- plt +
      ggplot2::geom_point(
      ggplot2::aes(
        y = !!ggplot2::sym("nominal_frequentist_power_separate"),
        color = "Nominal frequentist power",
        shape = "Nominal frequentist power"
      )
    ) +
      ggplot2::geom_point(
        ggplot2::aes(
          y = !!ggplot2::sym("frequentist_power_at_equivalent_tie"),
          color = "Frequentist power at equivalent TIE",
          shape = "Frequentist power at equivalent TIE"
        )
      )

    plt <- plt +
      ggplot2::geom_vline(xintercept = analysis_config$nominal_tie, color = "black", linetype = "dashed") + # Add horizontal line at x = 0.025
      ggplot2::scale_x_continuous(breaks = nominal_tie_breaks(analysis_config$nominal_tie))  # Add x-tick at x = 0.025
  }

  plt <- plt +
    ggplot2::scale_color_viridis_d() # Apply viridis color scale

  # Save plot

  plot.size <- set_size(textwidth)
  fig_width_in <- plot.size[1]
  fig_height_in <- plot.size[2]

  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf")
  return(export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png"))
}

#' Plot frequentist methods operating characteristics
#'
#' @description This function generates plots for the operating characteristics of frequentist methods.
#'
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param metrics A list of metrics to be plotted.
#'
#' @return None
#'
#' @examples NA
power_vs_tie_plots <- function(results_metrics_df, metrics) {
  check_plot_packages()
  check_required_colnames(
    results_metrics_df,
    required_colnames_consumer,
    context = "power_vs_tie_plots"
  )

  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)
  treatment_effects <- c("partially_consistent", "consistent")

  for (case_study in case_studies) {
    results_metrics_df1 <- results_metrics_df[results_metrics_df$case_study == case_study,]

    # Get the list of target sample sizes
    target_sample_sizes <- unique(results_metrics_df1[results_metrics_df1$case_study == case_study, target_size_columns, drop = FALSE])

    for (size_index in seq_len(nrow(target_sample_sizes))) {
      target_sample_size <- target_sample_sizes[size_index, , drop = FALSE]
      results_metrics_df2 <- dplyr::semi_join(results_metrics_df1, target_sample_size, by = target_size_columns)

      for (treatment_effect in treatment_effects) {

        target_to_source_std_ratios <- unique(results_metrics_df2$target_to_source_std_ratio)

        for (target_to_source_std_ratio in target_to_source_std_ratios) {
          results_metrics_df3 <- results_metrics_df2 %>%
            dplyr::filter(
              matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio)
            )

          source_denominator_change_factors <- unique(results_metrics_df3$source_denominator_change_factor)

          for (source_denominator_change_factor in source_denominator_change_factors){
            results_metrics_df4 <- results_metrics_df3 %>%
              dplyr::filter(
                matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor)
              )

            power_vs_tie(
              results_metrics_df4,
              case_study,
              target_sample_size$target_sample_size_control,
              target_sample_size$target_sample_size_treatment,
              treatment_effect = treatment_effect,
              power_difference = FALSE,
              metrics = metrics,
              source_denominator_change_factor = source_denominator_change_factor,
              target_to_source_std_ratio = target_to_source_std_ratio
            )

            power_vs_tie(
              results_metrics_df4,
              case_study,
              target_sample_size$target_sample_size_control,
              target_sample_size$target_sample_size_treatment,
              treatment_effect = treatment_effect,
              power_difference = TRUE,
              metrics = metrics,
              source_denominator_change_factor = source_denominator_change_factor,
              target_to_source_std_ratio = target_to_source_std_ratio
            )
          }
        }
      }
    }
  }
  return(invisible(NULL))
}
