plot_metric_vs_ess <- function(results_metrics_df,
                         case_study,
                         target_sample_size_control,
                         target_sample_size_treatment,
                         treatment_effect,
                         ess_method,
                         metric,
                         source_denominator_change_factor = 1,
                         target_to_source_std_ratio = 1) {

  if (!(metric$name %in% c("mse", "bias", "coverage", "precision", "interval_score"))){
   return()
  }

  directory <- paste0(figures_dir, case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  filename <- paste0(
    case_study,
    "_",
    metric$name,
    "_vs_",
    ess_method$name,
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

  if (treatment_effect == "consistent") {
    results_df <- results_df %>%
      dplyr::filter(target_treatment_effect == source_treatment_effect_estimate)
  } else if (treatment_effect == "no_effect") {
    theta_0 <- unique(results_df$theta_0)
    results_df <- results_df %>%
      dplyr::filter(target_treatment_effect == theta_0)
  } else if (treatment_effect == "partially_consistent") {
    results_df <- results_df %>%
      dplyr::filter(abs(
        target_treatment_effect - source_treatment_effect_estimate / 2
      ) < 1e-4)
  } else {
    stop('treatment_effect must be "consistent", "no_effect" or "partially_consistent"')
  }

  y_range <- sort(range(results_df[, metric$name]))
  ymin <- y_range[1]
  ymax <- y_range[2]

  x_range <- sort(range(results_df[,  ess_method$name]))
  xmin <- x_range[1]
  xmax <- x_range[2]

  x_width <- xmax - xmin
  y_width <- ymax - ymin

  y_cap_size <- x_width * relative_error_cap_width / 2 # Adjust relative_error_cap_width to control the relative cap size
  x_cap_size <- y_width * relative_error_cap_width / 2 # Adjust relative_error_cap_width to control the relative cap size

  # cap_size <- min(x_cap_size, y_cap_size)
  # x_cap_size <- cap_size
  # y_cap_size <- cap_size

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


  metric_name <- rlang::sym(metric$name)
  metric_uncertainty_lower <- rlang::sym(paste0(metric$metric_uncertainty, "_lower"))
  metric_uncertainty_upper <- rlang::sym(paste0(metric$metric_uncertainty, "_upper"))
  metric_label <- metric$label

  ess_name <- rlang::sym(ess_method$name)
  ess_uncertainty_lower <- rlang::sym(paste0(ess_method$metric_uncertainty, "_lower"))
  ess_uncertainty_upper <- rlang::sym(paste0(ess_method$metric_uncertainty, "_upper"))
  ess_label <- ess_method$label


  plt <- ggplot2::ggplot(results_df, ggplot2::aes(x = !!ess_name))

  if (metric$name == "mse" && ess_method$name == "ess_moment") {
    # The ESS counts patients per arm of a balanced study, so the target is
    # drawn at its equivalent sample size per arm.
    plt <- plt + ggplot2::geom_vline(
      xintercept = effective_arm_size(target_sample_size_control, target_sample_size_treatment),
      color = "black",
      linetype = "dashed"
    )
  }

  plt <- plt +
    ggplot2::geom_point(ggplot2::aes(
      y = !!metric_name,
      color = as.factor(method),
      shape = as.factor(method)
    )) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = !!metric_uncertainty_lower,
        ymax = !!metric_uncertainty_upper,
        color =  as.factor(method)
      ),
      width = y_cap_size
    ) +
    ggplot2::labs(
      title = plot_title,
      x = ess_label,
      y = metric_label,
      color = "Methods",
      # Title for the color legend
      shape = "Methods" # Title for the shape legend
    )

  # MSE-versus-ESS panels show uncertainty only in the MSE estimate.
  if (metric$name != "mse") {
    plt <- plt + ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = !!ess_uncertainty_lower,
        xmax = !!ess_uncertainty_upper,
        y = !!metric_name,
        color = as.factor(method)
      ),
      height = x_cap_size
    )
  }

  ## Keyed by method, not by position among the methods this figure happens to
  ## contain: viridis_d and the default shape scale both hand out their values
  ## in level order, so a case study missing a method recoloured every method
  ## after it. Here a method carries the same colour and symbol everywhere.
  plt <- plt +
    ggplot2::scale_color_manual(values = method_hue_map(results_df$method), name = "Methods") +
    ggplot2::scale_shape_manual(values = method_shape_map(results_df$method), name = "Methods")

  if (metric$name == "mse") {
    plt <- plt +
      ggplot2::guides(
        color = ggplot2::guide_legend(ncol = 5, byrow = TRUE),
        shape = ggplot2::guide_legend(ncol = 5, byrow = TRUE)
      )
  }

  # Save plot
  plot.size <- set_size(textwidth)
  fig_width_in <- plot.size[1]
  fig_height_in <- plot.size[2]

  if (metric$name == "mse") {
    fig_width_in <- fig_width_in * 1.35
    fig_height_in <- fig_height_in + 1.7
  }

  legend_theme <- if (metric$name == "mse") {
    ggplot2::theme(legend.position = "bottom", legend.title.position = "top")
  } else {
    NULL
  }

  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf",
               theme_extra = legend_theme)
  return(export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png",
               theme_extra = legend_theme))
}



#' Plot methods operating characteristics
#'
#' @description This function generates plots for the operating characteristics of different methods.
#'
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param frequentist_metrics List of frequentist metrics, as defined by
#'   sourcing `metrics_config.R`.
#' @param inference_metrics List of inference metrics, as defined by sourcing
#'   `metrics_config.R`.
#'
#' @return None
#'
#' @examples NA
plot_metrics_vs_ess <- function(results_metrics_df, frequentist_metrics, inference_metrics) {
  check_plot_packages()
  check_required_colnames(
    results_metrics_df,
    required_colnames_consumer,
    context = "plot_metrics_vs_ess"
  )

  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)

  ess_methods <- c("ess_elir", "ess_moment", "ess_precision")

  treatment_effects <-  treatment_effects <- c("no_effect", "partially_consistent", "consistent")

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


    target_sample_sizes <- unique(filtered_results_metrics_df[, target_size_columns, drop = FALSE])

    for (metric in frequentist_metrics) {
      for (size_index in seq_len(nrow(target_sample_sizes))) {
        for (ess_method in  ess_methods){
          for (treatment_effect in treatment_effects){
            plot_metric_vs_ess(
              results_metrics_df = filtered_results_metrics_df,
              case_study = case_study,
              ess_method = inference_metrics[[ess_method]],
              target_sample_size_control = target_sample_sizes$target_sample_size_control[size_index],
              target_sample_size_treatment = target_sample_sizes$target_sample_size_treatment[size_index],
              treatment_effect = treatment_effect,
              metric = metric,
            )
          }
        }
      }
    }
  }
  return(invisible(NULL))
}
