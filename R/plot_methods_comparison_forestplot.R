# Sort the rows of each method by the values of its parameters, keeping the
# rows of each method in the positions they had.
sort_within_methods <- function(data) {
  for (method in unique(data$method)) {
    data_subset <- data[data$method == method, ]
    params <- get_parameters(data_subset[, "parameters", drop = FALSE])
    sorted_index <- do.call(order, lapply(params[colnames(params)], as.factor))
    data[data$method == method, ] <- data_subset[sorted_index, ]
  }
  data
}

#' Forest-plot colours for a BExTE-app colour scheme
#'
#' @description The forest plots are the only charts BExTE-app renders as static
#' images (they are multi-panel gtables, so they cannot go through
#' bexte_plotly() like every other chart) and therefore cannot inherit the
#' page's stylesheet. These helpers take one of the app's bexte_palette() schemes and
#' recolour a subplot for it. `palette = NULL` is the publication figure,
#' unchanged.
#'
#' @param palette A colour scheme from bexte_palette(), or NULL.
#'
#' @return A ggplot2 theme, or NULL when there is no palette.
#' @keywords internal
forest_palette_theme <- function(palette) {
  if (is.null(palette)) {
    return(NULL)
  }
  return(ggplot2::theme(
    text = ggplot2::element_text(colour = palette$ink),
    axis.text = ggplot2::element_text(colour = palette$ink),
    axis.title = ggplot2::element_text(colour = palette$ink),
    plot.title = ggplot2::element_text(colour = palette$ink),
    legend.text = ggplot2::element_text(colour = palette$ink),
    legend.title = ggplot2::element_text(colour = palette$ink),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_line(colour = palette$ink_muted),
    panel.background = ggplot2::element_rect(fill = palette$surface, colour = NA),
    panel.border = ggplot2::element_rect(fill = NA, colour = palette$hairline),
    panel.grid.major = ggplot2::element_line(colour = palette$hairline, linetype = "solid",
                                    linewidth = 0.15),
    strip.background = ggplot2::element_rect(fill = palette$hairline, colour = NA),
    strip.text = ggplot2::element_text(colour = palette$ink),
    plot.background = ggplot2::element_blank()
  ))
}

#' The colours of the forest plots' reference lines, defaulting to the literal
#' values the publication figures use so that `palette = NULL` output is
#' untouched.
#'
#' @param palette A colour scheme from bexte_palette(), or NULL.
#'
#' @return A list with `target`, `source` and `ink` colours.
#' @keywords internal
forest_reference_colours <- function(palette) {
  if (is.null(palette)) {
    return(list(target = "black", source = "red", ink = "black"))
  }
  return(list(target = palette$ref_target, source = palette$ref_source, ink = palette$ink))
}

#' A point/pointrange layer carrying the palette's ink. Without a palette the
#' colour argument is omitted altogether, leaving ggplot2's default.
#'
#' @param geom The ggplot2 layer constructor to call.
#' @param palette A colour scheme from bexte_palette(), or NULL.
#' @param ... Further arguments for `geom`.
#'
#' @return A ggplot2 layer.
#' @keywords internal
forest_ink_layer <- function(geom, palette, ...) {
  args <- list(...)
  if (!is.null(palette)) {
    args$colour <- palette$ink
  }
  return(do.call(geom, args))
}

#' Express a metric as a ratio to the separate analysis
#'
#' @description Divides each row's metric by the value the `separate` (no
#'   borrowing) method takes in the same scenario. The scenario key is the full
#'   set of coordinates that identify a simulated configuration, so the ratio
#'   never collapses across drift or sample size.
#'
#'   A zero or missing denominator makes the ratio undefined. Those rows are
#'   dropped with a warning rather than reported as `Inf` or as some finite
#'   stand-in.
#'
#' @param df A frequentist results frame containing `separate` method rows.
#' @param metric_columns Character vector of columns to divide - typically the
#'   metric and its two confidence bounds.
#'
#' @return `df` with `metric_columns` divided by the separate analysis's value,
#'   in the original column order.
#'
#' @keywords internal
scale_by_separate <- function(df, metric_columns) {
  key <- c(
    "case_study", target_size_columns, "drift",
    "target_to_source_std_ratio", "source_denominator_change_factor"
  )
  missing_columns <- setdiff(c(key, metric_columns), names(df))
  if (length(missing_columns) > 0) {
    stop(
      "scale_by_separate() needs these columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  baseline <- df[df$method == "separate", c(key, metric_columns[1]), drop = FALSE]
  if (nrow(baseline) == 0) {
    stop("No `separate` method rows to use as the denominator.")
  }
  names(baseline)[names(baseline) == metric_columns[1]] <- ".separate_value"
  baseline <- unique(baseline)

  original_columns <- names(df)
  out <- merge(df, baseline, by = key, all.x = TRUE, sort = FALSE)

  usable <- !is.na(out$.separate_value) & out$.separate_value != 0
  if (any(!usable)) {
    warning(
      sum(!usable),
      " row(s) dropped: the separate analysis's value is zero or missing, ",
      "so the ratio is undefined."
    )
    out <- out[usable, , drop = FALSE]
  }

  for (column in metric_columns) {
    out[[column]] <- out[[column]] / out$.separate_value
  }

  return(out[, original_columns, drop = FALSE])
}


## Half the source study's equivalent per-arm sample size: the value the
## moment-based ESS colour scale is centred on. White therefore means
## borrowing half of what the source study is worth, blue less and red more,
## which is fixed per case study rather than relative to whatever range a
## particular figure happens to span.
forest_ess_midpoint <- function(data) {
  source_per_arm <- data$equivalent_source_sample_size_per_arm
  if (is.null(source_per_arm) || all(is.na(source_per_arm))) {
    source_per_arm <- 2 * data$source_sample_size_control * data$source_sample_size_treatment /
      (data$source_sample_size_control + data$source_sample_size_treatment)
  }
  return(unique(stats::median(source_per_arm, na.rm = TRUE)) / 2)
}

## The moment-based ESS each row achieved, or NULL when the run did not
## record it. Returning NULL rather than erroring keeps the forest plot
## usable on a results frame that predates the column.
forest_ess_values <- function(data) {
  if (is.null(data$ess_moment) || all(!is.finite(data$ess_moment))) {
    return(NULL)
  }
  return(data$ess_moment)
}

## The point-and-interval layer. The marker is filled by moment-based ESS and
## outlined thinly, so a near-white fill still reads as a point; `fill` is
## used rather than `colour` because the ESS panels already spend their
## colour scale on the target and source reference lines.
forest_ess_layer <- function(data, palette, ess_limits) {
  if (is.null(forest_ess_values(data)) || !is.null(palette)) {
    return(forest_ink_layer(ggplot2::geom_pointrange, palette, size = 0.001, linewidth = 1))
  }
  return(ggplot2::geom_pointrange(
    ggplot2::aes(fill = ess_moment),
    ## No explicit colour: the outline stays at ggplot2's default so the
    ## dark-theme palette can still set it, which is the contract
    ## forest_ink_layer() exists to honour.
    shape = 21, stroke = 0.2, size = 0.61, linewidth = 1.2
  ))
}

## The shared diverging fill scale. Limits are passed in so all three panels
## of a figure share one scale and one colourbar.
forest_ess_scale <- function(data, ess_limits, ess_midpoint) {
  if (is.null(forest_ess_values(data))) {
    return(NULL)
  }
  return(ggplot2::scale_fill_gradient2(
    name = "Moment-based ESS",
    low = "#3B4CC0", mid = "white", high = "#E8000B",
    midpoint = if (is.null(ess_midpoint)) forest_ess_midpoint(data) else ess_midpoint,
    limits = ess_limits,
    oob = scales::squish,
    guide = ggplot2::guide_colourbar(
      title.position = "top", title.hjust = 0.5,
      ## A border, so the pale middle of the scale still reads as part
      ## of the bar rather than as a gap in it.
      frame.colour = "grey20", frame.linewidth = 0.3,
      ticks.colour = "grey20",
      barwidth = grid::unit(2.2, "in"), barheight = grid::unit(0.12, "in")
    )
  ))
}

## Scientific notation on the x axis: a single power of ten factored out of
## the break labels and named in the axis title, so the ticks read 1.6, 2.0,
## 2.4 under "MSE (10^-2)" instead of 0.016, 0.020, 0.024.
forest_scientific_x <- function(values, label) {
  finite <- values[is.finite(values) & values != 0]
  if (length(finite) == 0) {
    return(list(scale = NULL, label = label))
  }
  exponent <- floor(log10(max(abs(finite))))
  ## Leave ordinary-sized numbers alone; a factor of 1 helps nobody.
  if (exponent >= -1 && exponent <= 2) {
    return(list(scale = NULL, label = label))
  }
  factor <- 10^exponent
  return(list(
    scale = ggplot2::scale_x_continuous(
      labels = function(x) formatC(x / factor, format = "fg", digits = 2)
    ),
    label = latex2exp::TeX(sprintf("%s ($10^{%d}$)", label, exponent))
  ))
}

## The shared legend row: the moment-based ESS colourbar, and on the ESS
## panels the target and source reference-line keys beside it. Taken from the
## panel that carries both guides, so the three panels can then be drawn
## without any legend of their own.
forest_legend_grob <- function(plt) {
  close_device <- use_font_capable_device()
  on.exit(close_device(), add = TRUE)

  with_legend <- plt + ggplot2::theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal"
  )
  built <- ggplot2::ggplotGrob(with_legend)
  boxes <- which(vapply(built$grobs, function(g) grepl("guide-box", g$name), logical(1)))
  for (index in boxes) {
    if (!inherits(built$grobs[[index]], "zeroGrob")) {
      return(built$grobs[[index]])
    }
  }
  return(NULL)
}

## The width a forest panel spends on everything that is not the plotting
## region: y-axis labels, ticks, axis title and margins. The first panel
## carries the method labels and so spends far more than the other two.
##
## Laying the three panels out on fixed fractions of the page gives them
## unequal plotting regions, because each fraction has to cover its own
## labels first. Handing each panel its own overhead plus an equal share of
## what is left makes the three grids the same width, so a distance along
## the x axis means the same thing in each.
forest_nonpanel_width <- function(grob) {
  panel_columns <- grob$layout$l[grepl("^panel", grob$layout$name)]
  if (length(panel_columns) == 0) {
    return(grid::unit(0, "cm"))
  }
  spanned <- seq(min(panel_columns), max(panel_columns))
  others <- setdiff(seq_along(grob$widths), spanned)
  if (length(others) == 0) {
    return(grid::unit(0, "cm"))
  }
  return(grid::convertWidth(sum(grob$widths[others]), "cm"))
}

## Each panel's own overhead, plus an equal share of the remaining width,
## resolved to absolute inches.
##
## The widths have to be absolute: `unit(x, "cm") + unit(1, "null")` is a
## unit arithmetic object, and arrangeGrob cannot distribute leftover
## space across those - the null parts get nothing and every panel
## collapses to its labels.
forest_panel_widths <- function(grobs, available_in) {
  overheads <- vapply(grobs, function(grob) {
    return(as.numeric(grid::convertWidth(forest_nonpanel_width(grob), "in")))
  }, numeric(1))
  ## A floor wide enough for the text that sits over and under a panel:
  ## its title ("Consistent (no drift)") and its axis label ("Relative
  ## Power (vs. separate analysis)"). Both are drawn at the panel's own
  ## width, so a panel narrower than they are lets them spill into the
  ## neighbouring panel. Scaled by the font, since that is what sets how
  ## wide the text runs; the page is then widened to hold the result.
  ## small_text_size is one of the style globals conf/plots_config.R
  ## supplies; fall back to its usual value rather than failing when the
  ## caller has not sourced them.
  font_size <- if (exists("small_text_size", inherits = TRUE)) {
    get("small_text_size")
  } else {
    12
  }
  panel_floor_in <- 2.2 * font_size / 12
  panel_share <- max((available_in - sum(overheads)) / length(grobs), panel_floor_in)
  return(grid::unit(overheads + panel_share, "in"))
}

#' Create a forest plot
#'
#' @description This function creates a forest plot using the provided data.
#'
#' @param data The data for creating the forest plot.
#' @param title The title of the forest plot.
#' @param ylabel A logical value indicating whether to display the y-axis label.
#' @param x_metric_name Name of the metric on the x-axis
#' @param x_metric_uncertainty_lower Lower limit of the metric on the x-axis
#' @param x_metric_uncertainty_upper Upper limit of the metric on the x-axis
#' @param x_metric_label Label of the metric on the x-axis
#' @param methods_labels Display labels of the methods, as defined by sourcing
#'   `methods_plots_config.R`.
#' @param legend Whether to draw a legend for the reference lines.
#' @param sort_by `"methods_parameters"` to order each method's rows by its
#'   parameter values, `"value"` to order the rows by the metric, or FALSE to
#'   keep the input order.
#' @param palette A colour scheme from bexte_palette() for BExTE-app's
#'   dark mode, or NULL for the publication figure.
#' @param reference_line x position of a dotted vertical reference line, or NULL for none.
#' @param ess_limits Limits of the moment-based ESS fill scale, so that several
#'   panels can share one scale, or NULL for the range of `data`.
#' @param ess_midpoint Midpoint of the ESS fill scale, or NULL for
#'   `forest_ess_midpoint(data)`.
#' @param nominal_tie_line x position of a dashed line marking the nominal type
#'   I error rate, or NULL for none.
#'
#' @return A ggplot object representing the forest plot.
#' @keywords internal
forest_subplot <- function(data,
                           title,
                           ylabel,
                           x_metric_name,
                           x_metric_uncertainty_lower,
                           x_metric_uncertainty_upper,
                           x_metric_label,
                           methods_labels,
                           legend = FALSE,
                           sort_by = FALSE,
                           palette = NULL,
                           reference_line = NULL,
                           ess_limits = NULL,
                           ess_midpoint = NULL,
                           nominal_tie_line = NULL) {
  refs <- forest_reference_colours(palette)
  ## A ratio plot needs its own reference at 1; the metric-name-driven
  ## reference lines further down do not cover it.
  extra_reference <- reference_line
  if (sort_by == "methods_parameters"){
    # Make sure that the data are sorted according to the parameters values
    data <- sort_within_methods(data)
  } else if (sort_by == "value"){
    x_metric_name_str <- as.character(x_metric_name)
    data <- data[order(data[[x_metric_name_str]]), ]
  }


  # Process the row of a results dataframe to create a Method + Parameters label
  labels <- sapply(seq_len(nrow(data)), function(i) {
    return(process_method_parameters_label(data[i, ], methods_labels))
  })

  data$rows <- seq(1, nrow(data))

  ## A single power of ten factored out of the tick labels and named in
  ## the axis title, so the ticks stay short.
  scientific_x <- forest_scientific_x(data[[as.character(x_metric_name)]], x_metric_label)

  plt <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      y = rows,
      x = !!x_metric_name,
      xmin = !!x_metric_uncertainty_lower,
      xmax = !!x_metric_uncertainty_upper
    )
  ) +
    forest_ess_layer(data, palette, ess_limits) +
    ggplot2::labs(title = title, x = scientific_x$label) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      text = ggplot2::element_text(family = font, size = small_text_size),
      axis.text.y = ggplot2::element_text(size = small_text_size),
      strip.text.x = ggplot2::element_text(size = small_text_size),
      legend.key = ggplot2::element_blank(),
      # strip.background = element_blank(),
      panel.grid.major = ggplot2::element_line(color = "gray80", linetype = "solid",  size = 0.15),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = small_text_size),
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.background = ggplot2::element_blank()
    ) +
    forest_palette_theme(palette) +
    ggplot2::scale_y_continuous(breaks = data$rows, labels = labels) + # LaTeX labels
    forest_ess_scale(data, ess_limits, ess_midpoint) +
    scientific_x$scale
    # scale_x_continuous(
    #   breaks = pretty(data[[x_metric_name]], n = 5), # Set 10 evenly spaced breaks
    #   labels = function(x) format(x, nsmall = 2)     # Optional: Format to 2 decimal places
    # )
  if (!ylabel) {
    plt <- plt + ggplot2::theme(axis.text.y = ggplot2::element_blank(), # Removes y-axis tick labels
                       axis.ticks.y = ggplot2::element_blank()) # Removes y-axis ticks
  }
  if (x_metric_name == "posterior_mean" ||
      x_metric_name == "posterior_median" || x_metric_name == "point_estimate_mean") {
    theta_0 <- unique(data$theta_0)
    if (legend) {
      plt <- plt + ggplot2::geom_vline(ggplot2::aes(xintercept = theta_0, color = "theta_0"), linetype = "dotted") +
        ggplot2::scale_color_manual(
          name = "",
          values = c("theta_0" = refs$ink),
          labels = expression(theta[0])
        )
    } else {
      plt <- plt + ggplot2::geom_vline(xintercept = theta_0, linetype = "dotted", color = refs$ink)
    }
  } else if (x_metric_name == "ess_moment" ||
             x_metric_name == "ess_precision" ||
             x_metric_name == "ess_elir") {
    if (is.null(data$equivalent_source_sample_size_per_arm) ||
        anyNA(data$equivalent_source_sample_size_per_arm)) {
      equivalent_source_sample_size_per_arm <- 2 * data$source_sample_size_control * data$source_sample_size_treatment / (data$source_sample_size_control + data$source_sample_size_treatment)
    } else {
      equivalent_source_sample_size_per_arm <- data$equivalent_source_sample_size_per_arm
    }
    # The ESS counts patients per arm of a balanced study, so the target is
    # drawn at its equivalent sample size per arm.
    target_equivalent_size <- unique(effective_arm_size(
      data$target_sample_size_control, data$target_sample_size_treatment
    ))
    source_sample_size_per_arm <- unique(round(equivalent_source_sample_size_per_arm))

    if (legend) {
      plt <- plt +
        ggplot2::geom_vline(
          ggplot2::aes(xintercept = target_equivalent_size, color = "Target equivalent sample size per arm"),
          linetype = "dotted"
        ) +
        ggplot2::geom_vline(
          ggplot2::aes(xintercept = source_sample_size_per_arm, color = "Source sample size per arm"),
          linetype = "dotted"
        ) +
        ggplot2::scale_color_manual(
          name = "",
          values = c(
            "Target equivalent sample size per arm" = refs$target,
            "Source sample size per arm" = refs$source
          ),
          labels = c(
            "Target equivalent sample size per arm" = expression(n[T]), # nolint: T_and_F_symbol_linter.
            "Source sample size per arm" = expression(N[S]/2)
          )
        )
    } else {
      plt <- plt +
        ggplot2::geom_vline(
          xintercept = target_equivalent_size,
          linetype = "dotted",
          "color" = refs$target
        ) +
        ggplot2::geom_vline(
          xintercept = source_sample_size_per_arm,
          linetype = "dotted",
          "color" = refs$source
        )
    }
  }
  # } else if (x_metric_name == "success_proba" |
  #            x_metric_name == "average_tie" |
  #            x_metric_name == "average_power" |
  #            x_metric_name == "prior_proba_no_benefit" |
  #            x_metric_name == "prior_proba_benefit" |
  #            x_metric_name == "prepost_proba_FP" |
  #            x_metric_name == "prepost_proba_TP" |
  #            x_metric_name == "upper_bound_proba_FP" |
  #            x_metric_name == "prior_proba_success") {
  #   # # Make sure the x range are within (0,1)
  #   # plot_build <- ggplot_build(p)
  #   # default_x_limits <- plot_build$layout$panel_params[[1]]$x.range
  #   #
  #   # # Default x-axis limits
  #   # xl <- default_x_limits[1]
  #   # xu <- default_x_limits[2]
  #   #
  #   # # Compute new limits
  #   # new_xl <- max(xl, 0)
  #   # new_xu <- min(xu, 1)
  #
  #
  #   # plt <- plt + scale_x_continuous(limits = c(new_xl, new_xu))
  # }
  if (!is.null(extra_reference)) {
    plt <- plt + ggplot2::geom_vline(
      xintercept = extra_reference,
      linetype = "dotted",
      color = refs$ink
    )
  }
  ## The nominal type I error rate, on the panel where the probability of
  ## declaring success is the type I error rate. Dashed, to read apart
  ## from the dotted reference lines above.
  if (!is.null(nominal_tie_line)) {
    plt <- plt + ggplot2::geom_vline(
      xintercept = nominal_tie_line,
      linetype = "dashed",
      color = refs$ink
    )
  }
  return(plt)
}


#' Create a forest plot
#'
#' @description This function creates a forest plot using the provided data.
#'
#' @param data The data for creating the forest plot.
#' @param title The title of the forest plot.
#' @param ylabel A logical value indicating whether to display the y-axis label.
#' @param x_metric_name Name of the metric on the x-axis
#' @param x_metric_label Label of the metric on the x-axis
#' @param palette A colour scheme from bexte_palette() for BExTE-app's
#'   dark mode, or NULL for the publication figure.
#'
#' @return A ggplot2::ggplot( object representing the forest plot.
#' @keywords internal
forest_subplot_no_uncertainty <- function(data,
                                          title,
                                          ylabel,
                                          x_metric_name,
                                          x_metric_label,
                                          methods_labels,
                                          legend = FALSE,
                                          sort_by = FALSE,
                                          palette = NULL) {

  if (sort_by == "methods_parameters"){
    # Make sure that the data are sorted according to the parameters values
    data <- sort_within_methods(data)
  } else if (sort_by == "value"){
    x_metric_name_str <- as.character(x_metric_name)
    data <- data[order(data[[x_metric_name_str]]), ]
  }

  labels <- sapply(seq_len(nrow(data)), function(i) {
    return(process_method_parameters_label(data[i, ], methods_labels))
  })

  data$rows <- seq(1, nrow(data))

  plt <- ggplot2::ggplot(data, ggplot2::aes(y = rows, x = !!x_metric_name)) +
    forest_ink_layer(ggplot2::geom_point, palette, size = 0.001) + # Scatter plot
    ggplot2::labs(title = title, x = x_metric_label) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      text = ggplot2::element_text(family = font, size = text_size),
      axis.text.y = ggplot2::element_text(size = text_size),
      strip.text.x = ggplot2::element_text(size = text_size),
      legend.key = ggplot2::element_blank(),
      # strip.background = element_blank(),
      panel.grid.major = ggplot2::element_line(color = "gray80", linetype = "solid",  size = 0.15),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = text_size),
      plot.title = ggplot2::element_text(hjust = 0.5)
      # Center the title
      # plot.background = element_blank()
    ) +
    forest_palette_theme(palette) +
    ggplot2::scale_y_continuous(breaks = data$rows, labels = labels) # LaTeX labels

  if (!ylabel) {
    plt <- plt + ggplot2::theme(axis.text.y = ggplot2::element_blank(), # Removes y-axis tick labels
                       axis.ticks.y = ggplot2::element_blank()) # Removes y-axis ticks
  }

  return(plt)
}

forest_combined_plot <- function(data,
                                 x_metric_name,
                                 x_metric_uncertainty_lower,
                                 x_metric_uncertainty_upper,
                                 x_metric_label,
                                 methods_labels,
                                 legend = TRUE,
                                 sorting_by = "value",
                                 palette = NULL) {

  # Sort data based on the specified sorting method
  if (sorting_by == "methods_parameters") {
    data <- sort_within_methods(data)
  } else if (sorting_by == "value") {
    x_metric_name_str <- as.character(x_metric_name)
    data <- data[order(data[[x_metric_name_str]]), ]
  }

  # Generate labels for the y-axis
  labels <- sapply(seq_len(nrow(data)), function(i) {
    return(process_method_parameters_label(data[i, ], methods_labels))
  })
  data$rows <- seq(1, nrow(data))

  # Create the plot
  plt <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      y = rows,
      x = !!x_metric_name,
      xmin = !!x_metric_uncertainty_lower,
      xmax = !!x_metric_uncertainty_upper,
      color = effect_type
    )
  ) +
    ggplot2::geom_pointrange(size = 0.001, linewidth = 1) +
    ggplot2::labs(
      title = "Add title",
      x = x_metric_label,
      color = ""
    ) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      text = ggplot2::element_text(family = font, size = small_text_size),
      axis.text.y = ggplot2::element_text(size = small_text_size),
      strip.text.x = ggplot2::element_text(size = small_text_size),
      legend.key = ggplot2::element_blank(),
      # strip.background = element_blank(),
      panel.grid.major = ggplot2::element_line(color = "gray80", linetype = "solid",  size = 0.15),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = small_text_size),
      plot.title = ggplot2::element_text(hjust = 0.5)
      # plot.background = element_blank()
    ) +
    forest_palette_theme(palette) +
    ggplot2::scale_y_continuous(breaks = data$rows, labels = labels) +
    ggplot2::scale_color_manual(values = c(
      "No effect" = "blue",
      "Partially consistent" = "orange",
      "Consistent (no drift)" = "green"
    ))

  # Add vertical lines for specific cases
  if (x_metric_name == "posterior_mean" || x_metric_name == "posterior_median" || x_metric_name == "point_estimate_mean") {
    theta_0 <- unique(data$theta_0)
    plt <- plt + ggplot2::geom_vline(xintercept = theta_0, linetype = "dotted", color = "black")
  } else if (x_metric_name %in% c("ess_moment", "ess_precision", "ess_elir")) {
    target_equivalent_size <- unique(effective_arm_size(
      data$target_sample_size_control, data$target_sample_size_treatment
    ))
    equivalent_source_sample_size_per_arm <- ifelse(
      is.null(data$equivalent_source_sample_size_per_arm) | anyNA(data$equivalent_source_sample_size_per_arm),
      2 * data$source_sample_size_control * data$source_sample_size_treatment / (data$source_sample_size_control + data$source_sample_size_treatment),
      data$equivalent_source_sample_size_per_arm
    )
    source_sample_size_per_arm <- unique(round(equivalent_source_sample_size_per_arm))

    plt <- plt +
      ggplot2::geom_vline(xintercept = target_equivalent_size, linetype = "dotted", color = "black") +
      ggplot2::geom_vline(xintercept = source_sample_size_per_arm, linetype = "dotted", color = "red")
  }

  return(plt)
}


#' Generate a forest plot
#'
#' @description This function generates a forest plot based on the provided results dataframe, selected case study, and selected target arm sizes.
#'
#' @param results_freq_df The results dataframe.
#' @param x_metric Metric on the x-axis
#' @param panels If TRUE, draw one panel per treatment effect (no effect,
#'   partially consistent, consistent); otherwise draw the three in a
#'   single combined plot.
#' @param palette A colour scheme from bexte_palette() for BExTE-app's
#'   dark mode, or NULL for the publication figure.
#' @param relative_to_separate When TRUE, divide the metric and its confidence
#'   bounds by the separate analysis's value in the same scenario, label the
#'   axis accordingly, and draw a reference line at 1. Used by supplementary
#'   figures S9, S14, S17, S27, S32 and S37.
#' @param filename_suffix Appended to the generated filename, to keep apart
#'   plots of the same slice under different time-to-event designs (figures
#'   S45-S47), which would otherwise share one file.
#' @param subtitle Optional line drawn above the panels.
#'
#' @return None
#' @export
forest_plot <- function(results_freq_df, x_metric, panels = TRUE, palette = NULL,
                        relative_to_separate = FALSE, filename_suffix = "",
                        subtitle = NULL) {
  check_plot_packages()
  close_device <- use_font_capable_device()
  on.exit(close_device(), add = TRUE)

  selected_case_study <- unique(results_freq_df$case_study)
  selected_target_sizes <- unique(results_freq_df[, target_size_columns, drop = FALSE])

  if (x_metric %in% names(frequentist_metrics)) {
    metric <- frequentist_metrics[[x_metric]]
  } else if (x_metric %in% names(inference_metrics)) {
    metric <- inference_metrics[[x_metric]]
  } else if (x_metric %in% names(bayesian_metrics)) {
    metric <- bayesian_metrics[[x_metric]]
  } else {
    stop("Metric is not supported")
  }

  if (!(metric$name %in% colnames(results_freq_df))) {
    warning("Metric not in input dataframe.")
    return()
  }

  x_metric_name <- rlang::sym(metric$name)
  x_metric_uncertainty_lower <- rlang::sym(paste0(metric$metric_uncertainty, "_lower"))
  x_metric_uncertainty_upper <- rlang::sym(paste0(metric$metric_uncertainty, "_upper"))
  x_metric_label <- metric$label

  if (relative_to_separate) {
    results_freq_df <- scale_by_separate(
      results_freq_df,
      c(metric$name, as.character(x_metric_uncertainty_lower),
        as.character(x_metric_uncertainty_upper))
    )
    x_metric_label <- paste0(x_metric_label, " relative to a separate analysis")
  }

  directory <- file.path(figures_dir, selected_case_study)

  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  filename <- paste0(
    selected_case_study,
    "_",
    if (relative_to_separate) "relative_" else "",
    x_metric,
    "_forest_plot_target_",
    target_size_label(
      selected_target_sizes$target_sample_size_control,
      selected_target_sizes$target_sample_size_treatment,
      "file"
    ),
    filename_suffix
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

  if (!is.null(palette)) {
    filename <- paste0(filename, "_dark")
  }

  file_path <- file.path(directory, filename)

  if (file.exists(paste0(file_path, ".png")) &&
      (!is.null(palette) || file.exists(paste0(file_path, ".pdf"))) &&
      remake_figures == FALSE){
    return(invisible(paste0(file_path, ".png")))
  }

  if (x_metric %in% names(bayesian_metrics)) {
    df <- results_freq_df %>%
      dplyr::mutate(means = !!x_metric_name)
  } else {
    df <- results_freq_df %>%
      dplyr::mutate(
        means = !!x_metric_name,
        error_low = !!x_metric_uncertainty_lower,
        error_upper = !!x_metric_uncertainty_upper
      )
  }

  round_decimal <- 4
  # Calculate important target treatment effects
  important_target_treatment_effects <- round(
    df %>%
      pull(source_treatment_effect_estimate) %>%
      unique() %>%
      first() %>%
      `*`(c(0, 1 / 2, 1)),
    round_decimal
  )

  # The following is to handle problems with the way R handles floating point precision.
  # Find the closest values in df$target_treatment_effect to important_target_treatment_effects
  closest_values <- sapply(important_target_treatment_effects, function(x) {
    return(df$target_treatment_effect[which.min(abs(df$target_treatment_effect - x))])
  })

  # Filter the dataframe based on approximate equality
  df <- df %>%
    rowwise() %>%
    filter(any(sapply(closest_values, function(val) is_approx_equal(target_treatment_effect, val)))) %>%
    ungroup() %>%
    ## By distance from no effect, not by the signed value. The three
    ## scenarios are the source effect times 0, 1/2 and 1, and the
    ## labels below are assigned in whatever order this leaves. Sorting
    ## the signed value puts the largest benefit first whenever the
    ## source effect is negative - teriflunomide and mepolizumab, where
    ## benefit is a negative log rate ratio - which labelled the full
    ## effect "No treatment effect" and the null "Consistent", swapping
    ## two of the three panels.
    dplyr::arrange(abs(target_treatment_effect))

  # Custom function to find the closest value
  closest_value <- function(x, values) {
    return(values[which.min(abs(values - x))])
  }

  # Assign closest values to target_treatment_effect
  df <- df %>%
    mutate(target_treatment_effect= sapply(target_treatment_effect, closest_value, values = closest_values))


  # # Filter the dataframe to keep only the rows with the closest values
  # df <- df %>%
  #   dplyr::filter(target_treatment_effect %in% closest_values) %>%
  #   dplyr::arrange(target_treatment_effect)

  if (length(unique(df$target_treatment_effect)) != 3) {
    stop("Only three main treatment effect values should be selected")
  }

  # Factorize target treatment effect
  df$target_treatment_effect <- factor(
    df$target_treatment_effect,
    levels = unique(df$target_treatment_effect),
    labels = c(
      "No treatment effect",
      "Partially consistent treatment effect",
      "Consistent treatment effect"
    )
  )

  no_effect_data <- df[df[, "target_treatment_effect"] == "No treatment effect", ]
  partially_consistent_effect_data <- df[df[, "target_treatment_effect"] == "Partially consistent treatment effect", ]
  consistent_effect_data <- df[df[, "target_treatment_effect"] == "Consistent treatment effect", ]

  if (nrow(no_effect_data) != nrow(partially_consistent_effect_data) || nrow(consistent_effect_data) != nrow(partially_consistent_effect_data) || nrow(no_effect_data) != nrow(consistent_effect_data)){
    warning("Inconsistency in the number of rows in the partially consistent, consistent, and no-effect data. Skipping this plot.")
    return()
  }

  # Ensure that the methods and parameters are correctly ordered
  no_effect_data <- no_effect_data %>%
    arrange(method, parameters)

  partially_consistent_effect_data <- partially_consistent_effect_data %>%
    arrange(method, parameters)

  consistent_effect_data <- consistent_effect_data %>%
    arrange(method, parameters)

  # Sort the dataframes by the value of interest
  x_metric_name_str <- as.character(x_metric_name)
  new_order <- order(no_effect_data[[x_metric_name_str]])
  no_effect_data <- no_effect_data[new_order, ]
  partially_consistent_effect_data <- partially_consistent_effect_data[new_order, ]
  consistent_effect_data <- consistent_effect_data[new_order, ]


  ## The probability of declaring success is the type I error rate when
  ## there is no effect and the power when there is one, so the three
  ## panels do not share an axis label.
  if (x_metric_name == "success_proba"){
    x_metric_label_no_effect <- "TIE"
    x_metric_label_partially_consistent <- "Power"
    x_metric_label_consistent <- "Power"
  } else {
    x_metric_label_no_effect <- x_metric_label
    x_metric_label_partially_consistent <- x_metric_label
    x_metric_label_consistent <- x_metric_label
  }

  ## A relative plot's axis is a ratio to the separate analysis, not the
  ## metric itself, and saying so is the difference between reading 2 as
  ## "twice the power of no borrowing" and as a power of 2. The break
  ## keeps the longer label inside a narrow panel.
  ## Under no effect the probability of declaring success is the type I
  ## error rate, so that panel - and only that panel - can carry the
  ## nominal rate as a reference. A relative plot's axis is a ratio, on
  ## which the nominal rate has no meaning, so it is left off there.
  nominal_tie <- if (x_metric_name == "success_proba" && !relative_to_separate &&
                     exists("analysis_config", inherits = TRUE)) {
    get("analysis_config")$nominal_tie
  } else {
    NULL
  }

  if (relative_to_separate) {
    as_relative <- function(label) {
      return(paste0("Relative ", label, "\n(vs. separate analysis)"))
    }
    x_metric_label_no_effect <- as_relative(x_metric_label_no_effect)
    x_metric_label_partially_consistent <-
      as_relative(x_metric_label_partially_consistent)
    x_metric_label_consistent <- as_relative(x_metric_label_consistent)
  }


  if (panels == TRUE){
    # Create forest plots for each effect
    ## One scale across the three panels, so a colour means the same thing
    ## in each and a single colourbar can describe all of them.
    ess_all <- c(no_effect_data$ess_moment,
                 partially_consistent_effect_data$ess_moment,
                 consistent_effect_data$ess_moment)
    ess_all <- ess_all[is.finite(ess_all)]
    ess_limits <- if (length(ess_all) > 0) range(ess_all) else NULL
    ess_midpoint <- forest_ess_midpoint(no_effect_data)

    plot_no_effect <- forest_subplot(
      no_effect_data,
      "No effect",
      ylabel = TRUE,
      x_metric_name,
      x_metric_uncertainty_lower,
      x_metric_uncertainty_upper,
      x_metric_label_no_effect,
      methods_labels = methods_labels,
      legend = FALSE,
      sort_by = FALSE,
      reference_line = if (relative_to_separate) 1 else NULL,
      palette = palette,
      ess_limits = ess_limits,
      ess_midpoint = ess_midpoint,
      nominal_tie_line = nominal_tie
    )
    plot_partially_consistent_effect <- forest_subplot(
      partially_consistent_effect_data,
      "Partially consistent",
      ylabel = FALSE,
      x_metric_name,
      x_metric_uncertainty_lower,
      x_metric_uncertainty_upper,
      x_metric_label_partially_consistent,
      methods_labels = methods_labels,
      legend = FALSE,
      sort_by = FALSE,
      reference_line = if (relative_to_separate) 1 else NULL,
      palette = palette,
      ess_limits = ess_limits,
      ess_midpoint = ess_midpoint
    )
    plot_consistent_effect <- forest_subplot(
      consistent_effect_data,
      "Consistent (no drift)",
      ylabel = FALSE,
      x_metric_name,
      x_metric_uncertainty_lower,
      x_metric_uncertainty_upper,
      x_metric_label_consistent,
      methods_labels = methods_labels,
      legend = TRUE,
      sort_by = FALSE,
      reference_line = if (relative_to_separate) 1 else NULL,
      palette = palette,
      ess_limits = ess_limits,
      ess_midpoint = ess_midpoint
    )

    ## Widths first: the page is then sized to hold them, rather than the
    ## panels being squeezed into a fixed page. At 12pt the method labels
    ## are wide enough that a fixed width overflows - titles collide and
    ## the last panel runs off the edge.
    panel_widths <- NULL

    ## One legend row for the whole figure, taken from the panel that
    ## carries every guide, before the panels are stripped of theirs.
    legend_row <- forest_legend_grob(plot_consistent_effect)
    drop_legend <- ggplot2::theme(legend.position = "none")

    # Convert to grobs
    grob_no_effect <- ggplot2::ggplotGrob(plot_no_effect + drop_legend)
    grob_partially_consistent_effect <- ggplot2::ggplotGrob(plot_partially_consistent_effect + drop_legend)
    grob_consistent_effect <- ggplot2::ggplotGrob(plot_consistent_effect + drop_legend)

    panel_widths <- forest_panel_widths(
      list(grob_no_effect, grob_partially_consistent_effect, grob_consistent_effect),
      set_size(textwidth)[1]
    )

    # Note that adding the Tikz export to the export_plots function does not work.
    if (plots_to_latex == TRUE) {
      # grid.arrange() (unlike arrangeGrob()) also draws to the current
      # graphics device as a side effect - that's exactly what's wanted here,
      # to capture the drawing into the open tikz() device below.
      grid::grid.newpage()
      tikzDevice::tikz(file = paste0(file_path, ".tex"), width = 5.92)
      plt <- gridExtra::grid.arrange(
        grob_no_effect,
        grob_partially_consistent_effect,
        grob_consistent_effect,
        ncol = 3,
        widths = panel_widths
      )
      dev.off()
    } else {
      # arrangeGrob() only builds the combined grob, with no drawing side
      # effect - unlike grid.arrange(), it doesn't require (or open) a
      # graphics device, which matters when this runs headless (e.g. from
      # BExTE-app), where the implicitly-opened default device can't
      # render the plot theme's font and grid.arrange() would error.
      panels <- gridExtra::arrangeGrob(
        grob_no_effect,
        grob_partially_consistent_effect,
        grob_consistent_effect,
        ncol = 3,
        widths = panel_widths
      )
      plt <- if (is.null(legend_row)) {
        panels
      } else {
        gridExtra::arrangeGrob(panels, legend_row, ncol = 1,
                               heights = grid::unit.c(
                                 grid::unit(1, "null"),
                                 grid::grobHeight(legend_row)
                               ))
      }
    }
  } else {
    # Combine data from all three cases into one data frame
    combined_data <- rbind(
      cbind(no_effect_data, effect_type = "No effect"),
      cbind(partially_consistent_effect_data, effect_type = "Partially consistent"),
      cbind(consistent_effect_data, effect_type = "Consistent (no drift)")
    )

    plt <- forest_combined_plot(
      combined_data,
      x_metric_name = "posterior_mean",
      x_metric_uncertainty_lower = "lower_bound",
      x_metric_uncertainty_upper = "upper_bound",
      x_metric_label = "Posterior Mean",
      methods_labels = methods_labels,
      legend = TRUE,
      palette = palette
    )
  }
  plot_size <- set_size(textwidth)
  ## Wide enough for what the panels were actually given, never narrower
  ## than the text width.
  fig_width_in <- if (exists("panel_widths", inherits = FALSE) && !is.null(panel_widths)) {
    max(plot_size[1], sum(as.numeric(grid::convertWidth(panel_widths, "in"))))
  } else {
    plot_size[1]
  }
  ## A fixed height crams however many methods were simulated into the
  ## same space, which is why the y-axis labels overlapped once the runs
  ## carried all eleven methods and their parameter grids. Give every row
  ## its own strip of the page instead, and never less than the old
  ## height.
  ## As tight as the labels allow: a row need only be as tall as the text
  ## in it, which is small_text_size points, plus a quarter of that for
  ## the descenders and subscripts the method labels are full of. Tying it
  ## to the font rather than a fixed number keeps the spacing minimal if
  ## the text size is ever changed.
  label_height_in <- (small_text_size) / 72
  fig_height_in <- max(
    plot_size[2],
    ## 1.05, not 1.25: a row is the height of its own text plus a hair, so
    ## the labels sit directly under one another with no gap. Below about
    ## 1.0 they would start to touch.
    1.05 * label_height_in * nrow(no_effect_data) + 1.6
  )

  if (!is.null(subtitle)) {
    if (ggplot2::is.ggplot(plt)) {
      plt <- plt + ggplot2::labs(subtitle = subtitle)
    } else {
      subtitle_grob <- grid::textGrob(
        subtitle, x = grid::unit(0.5, "lines"), hjust = 0,
        gp = grid::gpar(fontfamily = font, fontsize = small_text_size)
      )
      plt <- gridExtra::arrangeGrob(
        subtitle_grob, plt, ncol = 1,
        heights = grid::unit.c(
          grid::grobHeight(subtitle_grob) + grid::unit(0.5, "lines"),
          grid::unit(1, "null")
        )
      )
      fig_height_in <- fig_height_in + 0.3
    }
  }

  if (is.null(palette)) {
    export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf", forest_plot = TRUE)
  }
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png", forest_plot = TRUE,
               bg = palette$surface)
  return(invisible(paste0(file_path, ".png")))
}


#' Generate a forest plot for bayesian metrics
#'
#' @description This function generates a forest plot based on the provided results dataframe, selected case study, and selected target arm sizes.
#'
#' @param results_bayes_df The results dataframe.
#' @param x_metric Metric on the x-axis
#' @param palette A colour scheme from bexte_palette() for BExTE-app's
#'   dark mode, or NULL for the publication figure.
#'
#' @return None
#' @export
forest_plot_bayesian <- function(results_bayes_df, x_metric, palette = NULL) {
  check_plot_packages()
  close_device <- use_font_capable_device()
  on.exit(close_device(), add = TRUE)

  selected_case_study <- unique(results_bayes_df$case_study)
  selected_target_sizes <- unique(results_bayes_df[, target_size_columns, drop = FALSE])


  directory <- file.path(figures_dir, selected_case_study)
  if (!dir.exists(directory)) {
    dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  }

  if (x_metric %in% names(bayesian_metrics)) {
    metric <- bayesian_metrics[[x_metric]]
  } else {
    stop("Metric is not supported")
  }

  if (!(x_metric %in% colnames(results_bayes_df))) {
    warning("Metric not in input dataframe.")
    return()
  }

  x_metric_name <- rlang::sym(metric$name)
  x_metric_label <- metric$label


  filename <- paste0(
    selected_case_study,
    "_",
    x_metric,
    "_forest_plot_target_",
    target_size_label(
      selected_target_sizes$target_sample_size_control,
      selected_target_sizes$target_sample_size_treatment,
      "file"
    )
  )

  if (!is.null(palette)) {
    filename <- paste0(filename, "_dark")
  }

  file_path <- file.path(directory, filename)

  if (file.exists(paste0(file_path, ".png")) &&
      (!is.null(palette) || file.exists(paste0(file_path, ".pdf"))) &&
      remake_figures == FALSE){
    return(invisible(paste0(file_path, ".png")))
  }

  round_decimal <- 4

  df <- results_bayes_df %>%
    dplyr::mutate(means = !!x_metric_name)

  ui_design_prior_data <- df[df[, "design_prior_type"] == "ui_design_prior", ]
  analysis_design_prior_data <- df[df[, "design_prior_type"] == "analysis_prior", ]
  source_posterior_design_prior_data <- df[df[, "design_prior_type"] == "source_posterior", ]


  # Sort the dataframes by the value of interest
  x_metric_name_str <- as.character(x_metric_name)
  new_order <- order(ui_design_prior_data[[x_metric_name_str]])
  ui_design_prior_data <- ui_design_prior_data[new_order, ]
  analysis_design_prior_data <- analysis_design_prior_data[new_order, ]
  source_posterior_design_prior_data <- source_posterior_design_prior_data[new_order, ]

  # Create forest plots for each effect
  plot_ui_design_prior <- forest_subplot_no_uncertainty(
    ui_design_prior_data,
    "UI design prior",
    ylabel = TRUE,
    x_metric_name,
    x_metric_label,
    methods_labels = methods_labels,
    legend = FALSE,
    palette = palette
  )
  plot_analysis_design_prior <- forest_subplot_no_uncertainty(
    analysis_design_prior_data,
    "Analysis prior",
    ylabel = FALSE,
    x_metric_name,
    x_metric_label,
    methods_labels = methods_labels,
    legend = FALSE,
    palette = palette
  )
  plot_source_posterior_design_prior <- forest_subplot_no_uncertainty(
    source_posterior_design_prior_data,
    "Source posterior as design prior",
    ylabel = FALSE,
    x_metric_name,
    x_metric_label,
    methods_labels = methods_labels,
    legend = FALSE,
    palette = palette
  )

  # Convert to grobs
  grob_ui_design_prior <- ggplot2::ggplotGrob(plot_ui_design_prior)
  grob_analysis_design_prior <- ggplot2::ggplotGrob(plot_analysis_design_prior)
  grob_source_posterior_design_prior <- ggplot2::ggplotGrob(plot_source_posterior_design_prior)

  # arrangeGrob() (unlike grid.arrange()) only builds the combined grob, with
  # no drawing side effect, so it doesn't require (or open) a graphics
  # device - which matters when this runs headless (e.g. from
  # BExTE-app), where the implicitly-opened default device can't render the plot
  # theme's font and grid.arrange() would error.
  plt <- gridExtra::arrangeGrob(
    grob_ui_design_prior,
    grob_analysis_design_prior,
    grob_source_posterior_design_prior,
    ncol = 3,
    widths = c(0.46, 0.23, 0.31)
  )

  plot.size <- set_size(textwidth)
  fig_width_in <- plot.size[1]
  fig_height_in <- plot.size[2]

  if (is.null(palette)) {
    export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf", forest_plot = TRUE)
  }
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png", forest_plot = TRUE,
               bg = palette$surface)
  return(invisible(paste0(file_path, ".png")))
}


#' Plot methods comparison
#'
#' @description This function plots the comparison of methods based on the provided results metrics dataframe.
#'
#' @param results_metrics_df The results metrics dataframe.
#' @param metrics Metrics for which to create forest plots.
#'
#' @return None
forest_plot_methods_comparison <- function(results_metrics_df, metrics) {
  check_plot_packages()
  check_required_colnames(
    results_metrics_df,
    required_colnames_consumer,
    context = "forest_plot_methods_comparison"
  )

  case_studies <- unique(results_metrics_df$case_study)

  results_metrics_df <- results_metrics_df[(results_metrics_df[, "source_denominator_change_factor"] == 1 &
                                              !is.na(results_metrics_df[, "source_denominator_change_factor"])), ]

  for (metric in names(metrics)) {
    if (metric %in% c('tie', 'power')){
      next
    }

    for (case_study in case_studies) {
      results_metrics_df1 <- results_metrics_df[results_metrics_df$case_study == case_study,]

      target_sample_sizes <- unique(results_metrics_df1[, target_size_columns, drop = FALSE])
      target_sample_sizes <- target_sample_sizes[stats::complete.cases(target_sample_sizes), , drop = FALSE]
      for (size_index in seq_len(nrow(target_sample_sizes))) {
        results_metrics_df2 <- dplyr::semi_join(
          results_metrics_df1, target_sample_sizes[size_index, , drop = FALSE], by = target_size_columns
        )

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

            forest_plot(results_metrics_df4, metric)
          }
        }
      }
    }
  }
  return(invisible(NULL))
}

#' Plot methods comparison
#'
#' @description This function plots the comparison of methods based on the provided results metrics dataframe.
#'
#' @param results_metrics_df The results metrics dataframe.
#' @param metrics Metrics for which to create forest plots.
#'
#' @return None
forest_plot_methods_comparison_bayesian <- function(results_metrics_df, metrics) {
  check_plot_packages()
  case_studies <- unique(results_metrics_df$case_study)

  results_metrics_df <- results_metrics_df[(results_metrics_df[, "source_denominator_change_factor"] == 1 &
                                              !is.na(results_metrics_df[, "source_denominator_change_factor"])), ]

  for (metric in names(metrics)) {
    for (case_study in case_studies) {
      results_metrics_df1 <- results_metrics_df[results_metrics_df$case_study == case_study,]

      target_sample_sizes <- unique(results_metrics_df1[, target_size_columns, drop = FALSE])
      target_sample_sizes <- target_sample_sizes[stats::complete.cases(target_sample_sizes), , drop = FALSE]
      for (size_index in seq_len(nrow(target_sample_sizes))) {
        results_metrics_df2 <- dplyr::semi_join(
          results_metrics_df1, target_sample_sizes[size_index, , drop = FALSE], by = target_size_columns
        )

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
            forest_plot_bayesian(results_metrics_df4, metric)
          }
        }
      }
    }
  }
  return(invisible(NULL))
}
