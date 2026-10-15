## Read title/label text through ggplot's nested title grobs.
vs_tie_legend_text <- function(grob) {
  if (!is.null(grob$label)) {
    label <- grob$label
    if (is.expression(label)) label <- label[[1]]
    if (is.language(label)) return(paste(deparse(label), collapse = ""))
    return(paste(as.character(label), collapse = ""))
  }
  children <- c(grob$grobs, as.list(grob$children))
  return(unlist(lapply(children, vs_tie_legend_text), use.names = FALSE))
}

vs_tie_legend_title <- function(legend) {
  title <- which(legend$layout$name == "title")
  if (!length(title)) return(NA_character_)
  text <- vs_tie_legend_text(legend$grobs[[title[1]]])
  if (!length(text)) return(NA_character_)
  return(trimws(text[1]))
}

vs_tie_key_labels <- function(method, parameters) {
  if (length(parameters) != 1) return(parameters)
  if (vs_tie_key_is_bare(parameters)) return(method)
  parameter <- parameters[[1]]
  if (is.expression(parameter)) parameter <- parameter[[1]]
  return(as.expression(list(bquote(.(method) ~ .(parameter)))))
}

## Does this method contribute a single key with nothing to say about it?
## Such a key is better labelled with the method name than headed by it
## and left blank. The labels are plotmath expressions, so emptiness is
## checked on the deparsed form rather than with nzchar().
vs_tie_key_is_bare <- function(method_parameters) {
  if (length(method_parameters) != 1) {
    return(FALSE)
  }
  text <- paste(deparse(method_parameters[[1]]), collapse = "")
  return(grepl('^(expression\\()?\\s*""\\s*\\)?$', text))
}

## Where each method's legend sits in the versus-type-I-error legend block.
##
## Entries are the labels the legends carry, and the grid is read as drawn:
## a label repeated across neighbouring cells spans them. Methods present in
## the data but absent from this table are appended in rows of their own
## rather than dropped, so adding a method to a figure cannot silently lose
## its legend.
PAPER_VS_TIE_LEGEND_LAYOUT <- rbind(
  c("Pooling",  "RMP", "Conditional PP"),
  c("Separate", "RMP", "NPP"),
  c("EBPP",     "RMP", "p-PP"),
  c("Com. PP",  "Com. PP", "Empirical RMP"),
  c("Com. prior", "Com. prior", "NPP (KL)"),
  c("TtP (diff.)", "TtP (eq.)", "TtP (eq.)")
)

## Lay the method legends out as a grid of groups rather than a single row
## or a single column.
##
## Each method contributes its own guide, and ggplot places multiple guides
## either all in a row - which runs off the page once there are eight of
## them at 12pt - or all in a column, which is complete but far taller than
## the panel. The guide box is a gtable whose children are the individual
## legends, so they can be pulled out and arranged in columns.
vs_tie_legend_grid <- function(plt, labels, available_in = Inf) {
  close_device <- use_font_capable_device()
  on.exit(close_device(), add = TRUE)

  built <- ggplot2::ggplotGrob(
    plt + ggplot2::theme(
      legend.position = "bottom", legend.box = "vertical",
      legend.title.position = "top", legend.margin = ggplot2::margin(2, 4, 2, 4)
    )
  )
  boxes <- which(vapply(built$grobs, function(g) grepl("guide-box", g$name), logical(1)))

  legends <- list()
  for (index in boxes) {
    box <- built$grobs[[index]]
    if (!inherits(box, "gtable")) {
      next
    }
    ## The zeroGrob placeholders are not legends.
    legends <- c(legends, box$grobs[vapply(box$grobs, function(g) {
      return(inherits(g, "gtable"))
    }, logical(1))])
  }
  if (length(legends) == 0) {
    return(NULL)
  }
  ## Match actual titles (or inline method labels), never guide order.
  labels <- vapply(legends, function(g) {
    title <- vs_tie_legend_title(g)
    if (!is.na(title)) return(title)
    text <- vs_tie_legend_text(g)
    matches <- labels[vapply(labels, function(label) {
      return(any(vapply(text, function(value) {
        return(identical(value, label) || startsWith(value, paste0('"', label, '" ~')))
      }, logical(1))))
    }, logical(1))]
    if (length(matches) != 1) stop("Cannot identify a vs-TIE legend")
    return(matches)
  }, character(1))
  ## Place each legend where PAPER_VS_TIE_LEGEND_LAYOUT asks for it. A
  ## label repeated across cells spans them, which is how the wide
  ## commensurate-prior legend runs along the bottom.
  layout <- PAPER_VS_TIE_LEGEND_LAYOUT
  placed <- intersect(as.vector(layout), labels)
  ## Anything the table does not mention gets a row to itself.
  leftover <- setdiff(labels, placed)
  for (label in leftover) {
    layout <- rbind(layout, rep(label, ncol(layout)))
  }

  layout_matrix <- matrix(
    match(layout, labels), nrow = nrow(layout), ncol = ncol(layout)
  )
  ## Cells naming a method this figure does not draw stay empty.
  keep_rows <- apply(layout_matrix, 1, function(r) any(!is.na(r)))
  layout_matrix <- layout_matrix[keep_rows, , drop = FALSE]
  if (nrow(layout_matrix) == 0) {
    return(NULL)
  }

  keep_cols <- apply(layout_matrix, 2, function(x) any(!is.na(x)))
  layout_matrix <- layout_matrix[, keep_cols, drop = FALSE]
  used <- sort(unique(as.vector(layout_matrix)[!is.na(as.vector(layout_matrix))]))

  ## Satisfy single-cell sizes first, then only add the space a spanning
  ## legend still needs. RMP's height and Com. PP's width count once.
  dimensions <- function(axis, measure, convert) {
    sizes <- numeric(dim(layout_matrix)[axis])
    spans <- lapply(used, function(i) {
      return(unique(which(layout_matrix == i, arr.ind = TRUE)[, axis]))
    })
    for (k in order(lengths(spans))) {
      cells <- spans[[k]]
      needed <- convert(measure(legends[[used[k]]]), "in", valueOnly = TRUE)
      deficit <- max(0, needed - sum(sizes[cells]))
      sizes[cells] <- sizes[cells] + deficit / length(cells)
    }
    return(sizes)
  }
  heights <- dimensions(1, grid::grobHeight, grid::convertHeight)
  widths <- dimensions(2, grid::grobWidth, grid::convertWidth)
  if (is.finite(available_in) && sum(widths) < available_in) {
    widths <- widths + (available_in - sum(widths)) / length(widths)
  }
  layout_matrix[] <- match(layout_matrix, used)
  aligned <- lapply(legends[used], function(g) {
    return(grid::grobTree(g, vp = grid::viewport(
      x = 0, y = 1, just = c("left", "top"),
      width = grid::grobWidth(g), height = grid::grobHeight(g)
    )))
  })
  return(gridExtra::arrangeGrob(
    grobs = aligned, layout_matrix = layout_matrix,
    heights = grid::unit(heights, "in"), widths = grid::unit(widths, "in")
  ))
}

## Shared by frequentist and Bayesian metric plots. Grow the export when
## necessary; narrowing grid cells does not shrink the text inside them.
vs_tie_compose <- function(plt, labels, width, height) {
  close_device <- use_font_capable_device()
  on.exit(close_device(), add = TRUE)

  legend <- vs_tie_legend_grid(plt, labels, available_in = width)
  if (is.null(legend)) return(list(plot = plt, width = width, height = height))
  legend_height <- grid::convertHeight(grid::grobHeight(legend), "in", valueOnly = TRUE)
  width <- max(width, grid::convertWidth(grid::grobWidth(legend), "in", valueOnly = TRUE))
  panel <- ggplot2::ggplotGrob(plt + ggplot2::theme(legend.position = "none"))
  return(list(
    plot = gridExtra::arrangeGrob(panel, legend, ncol = 1,
      heights = grid::unit.c(grid::unit(1, "null"), grid::grobHeight(legend))),
    width = width, height = max(height, height * 0.7 + legend_height)
  ))
}

#' Plot an operating characteristic against the type I error rate
#'
#' @description Draws one point per method and parameter combination, placing
#'   the operating characteristic against the type I error rate that method
#'   incurs, so that methods can be compared at the error rate they actually
#'   spend rather than at their nominal one.
#'
#' @param results_metrics_df The dataframe containing the results and metrics.
#' @param case_study The case study name.
#' @param target_sample_size_control,target_sample_size_treatment The target
#'   arm sizes.
#' @param treatment_effect The treatment effect scenario ("consistent",
#'   "no_effect" or "partially_consistent").
#' @param operating_characteristic The metric to plot, as an entry of
#'   `frequentist_metrics` or `inference_metrics`.
#' @param source_denominator_change_factor The source denominator change factor.
#' @param target_to_source_std_ratio The target to source standard deviation ratio.
#' @param show_tie_error_bars Whether to draw the Monte Carlo interval on the
#'   type I error axis.
#'
#' @return None
#'
#' @keywords internal
operating_characteristic_vs_tie <- function(
    results_metrics_df,
    case_study,
    target_sample_size_control,
    target_sample_size_treatment,
    treatment_effect,
    operating_characteristic,
    source_denominator_change_factor,
    target_to_source_std_ratio,
    show_tie_error_bars = FALSE
) {
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

  theta_0 <- unique(results_df$theta_0)
  results_df_tie <- results_df %>% dplyr::filter(target_treatment_effect == theta_0)

  if (treatment_effect == "consistent") {
    results_df <- results_df %>%
      dplyr::filter(target_treatment_effect == source_treatment_effect_estimate)
  } else if (treatment_effect == "no_effect") {
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

  # Process the row of a results dataframe to create a Method + Parameters label
  results_df$parameters_labels <- lapply(seq_len(nrow(results_df)), function(i) {
    return(process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE))
  })

  results_df$parameters_labels_not_latex <- lapply(seq_len(nrow(results_df)), function(i) {
    return(process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE, as_latex = FALSE))
  })

  results_df$methods_parameters_labels <- lapply(seq_len(nrow(results_df)), function(i) {
    return(process_method_parameters_label(results_df[i, ], methods_labels, method_name = TRUE, as_latex = FALSE))
  })

  parameters_labels <- setNames(object = results_df$parameters_labels, nm = results_df$parameters_labels)
  methods <- format_results_df_methods(results_df)
  results_df$method <- methods
  methods <- unique(methods)

  filename <- paste0(
    case_study,
    "_",
    operating_characteristic$name,
    "_vs_tie_",
    target_size_label(target_sample_size_control, target_sample_size_treatment, "file"),
    "_treatment_effect=",
    treatment_effect
  )

  filename <- format_filename(
    filename = filename,
    case_study = case_study,
    target_to_source_std_ratio = target_to_source_std_ratio,
    source_denominator_change_factor = source_denominator_change_factor
  )
  file_path <- file.path(directory, filename)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE) {
    return()
  }

  y_range <- sort(range(results_df[, operating_characteristic$name]))
  ymin <- y_range[1]
  ymax <- y_range[2]

  x_range <- sort(range(results_df[, "tie"]))
  xmin <- x_range[1]
  xmax <- x_range[2]

  x_width <- xmax - xmin
  y_width <- ymax - ymin

  # geom_errorbar() measures its cap width on the TIE (x) axis.
  y_cap_size <- x_width * relative_error_cap_width / 2
  x_cap_size <- y_width * relative_error_cap_width / 2

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

  plot_title <- format_title(
    title = plot_title,
    case_study = case_study,
    target_to_source_std_ratio = target_to_source_std_ratio,
    source_denominator_change_factor = source_denominator_change_factor
  )


  # Split data by 'method'
  results_df_split <- split(results_df, results_df$method)

  # Unique parameters and methods
  unique_parameters <- unique(results_df$parameters_labels)
  unique_methods <- sort(unique(results_df$method))

  # Shapes for methods, keyed by method rather than handed out by position
  shapes <- method_shape_map(unique_methods)

  # Initialize the plot
  plt <- ggplot2::ggplot(mapping = ggplot2::aes(x = tie, y = !!ggplot2::sym(operating_characteristic$name)))

  # Iterate over each method and add layers
  plt <- plt + purrr::imap(results_df_split, function(df_subset, method_name) {
    df_subset <- df_subset[!is.na(df_subset[[operating_characteristic$name]]), ]

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
        ggplot2::aes(
          x = tie,
          y = !!ggplot2::sym(operating_characteristic$name),
          color = factor(parameters_labels, levels = method_parameters)
        ),
        shape = method_shape,
        size = 2
      ),
      # Add vertical error bars
      ggplot2::geom_errorbar(
        data = df_subset,
        ggplot2::aes(
          x = tie,
          ymin = !!ggplot2::sym(paste0("conf_int_", operating_characteristic$name, "_lower")),
          ymax = !!ggplot2::sym(paste0("conf_int_", operating_characteristic$name, "_upper")),
          color = factor(parameters_labels, levels = method_parameters)
        ),
        width = y_cap_size
      ),
      ## The horizontal bars carry the Monte Carlo error on the type I
      ## error rate. With one point per method-parameter combination they
      ## cross each other and obscure the points they belong to, so they
      ## are off unless asked for.
      if (show_tie_error_bars) {
        ggplot2::geom_errorbarh(
          data = df_subset,
          ggplot2::aes(
            y = !!ggplot2::sym(operating_characteristic$name),
            xmin = conf_int_tie_lower,
            xmax = conf_int_tie_upper,
            color = factor(parameters_labels, levels = method_parameters)
          ),
          height = x_cap_size
        )
      } else {
        NULL
      },
      # Define the color scale for parameters within this method
      ## A method with one key needs no separate heading: the method name
      ## goes on the key itself, so it reads "Pooling" rather than a
      ## "Pooling" title above a lone unlabelled symbol. Where that key
      ## also carries parameters they follow the name.
      ggplot2::scale_color_manual(
        values = method_colors,
        labels = vs_tie_key_labels(method_name, method_parameters),
        name = if (length(method_parameters) == 1) NULL else method_name,
        drop = TRUE,
        ## Up to three entries per row: a method with nine parameter values
        ## (RMP's weights) becomes a 3x3 block rather than a column nine
        ## tall, which is what made the legend taller than the panel.
        guide = ggplot2::guide_legend(
          override.aes = list(shape = method_shape),
          ncol = min(3, length(method_parameters)),
          byrow = TRUE
        )
      )
    )

    # Reset the color scale for the next method
    layers <- c(layers, list(ggnewscale::new_scale_color()))

    return(layers)
  })

  # Add vertical line
  plt <- plt +
    ggplot2::geom_vline(
      xintercept = analysis_config$nominal_tie,
      color = "black",
      linetype = "dashed"
    ) +
    ggplot2::scale_x_continuous(
      breaks = nominal_tie_breaks(analysis_config$nominal_tie)
    )

  # A metric measured against a reference - a difference whose zero means "no
  # better and no worse" - names it, and gets a horizontal line there.
  if (!is.null(operating_characteristic$reference_line)) {
    plt <- plt +
      ggplot2::geom_hline(
        yintercept = operating_characteristic$reference_line,
        color = "black",
        linetype = "dotted"
      )
  }

  # Add shape scale for methods
  plt <- plt +
    ggplot2::scale_shape_manual(
      values = shapes,
      name = "Methods",
      guide = ggplot2::guide_legend(order = 1, nrow = 2)
    ) +
    ggplot2::labs(
      title = plot_title,
      x = "TIE",
      y = operating_characteristic$label
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
    legend.key.size = grid::unit(0.1, "cm"),
    legend.spacing.x = grid::unit(0.03, "cm"),
    legend.spacing.y = grid::unit(0.01, "cm"),
    legend.position = "bottom",
    legend.direction = "vertical",
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    ## The eight method groups laid side by side make a legend row wider
    ## than any page - at 12pt it was clipped at both ends however wide
    ## the figure was made. Stacking the groups wraps the legend instead
    ## of running it off the edge.
    legend.box = "vertical",
    legend.box.just = "left"
  )

  plot.size <- set_size(textwidth)
  ## The 1.35 was sized for the old 6pt legend text. The legend sits under
  ## the panel and is the widest thing on the page, so its width tracks
  ## the font: at 12pt it needs about half again as much room, and
  ## without this it is clipped at both ends.
  legend_scale <- max(1, small_text_size / 8)
  fig_width_in <- plot.size[1] * 1.35 * legend_scale
  ## The stacked legend is tall - one row per method group - so the
  ## page has to grow with it, or the panel is squeezed to a sliver.
  ## The legend is a compact grid rather than one group per row, so the
  ## page only needs the font's worth of extra room, not the large
  ## allowance a stacked legend wanted.
  fig_height_in <- plot.size[2] * 1.5 * legend_scale

  composed <- vs_tie_compose(plt, names(results_df_split), fig_width_in, fig_height_in)
  plt <- composed$plot
  fig_width_in <- composed$width
  fig_height_in <- composed$height * 0.85

  # Export the plots
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf", adjust_theme = FALSE)
  return(export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png", adjust_theme = FALSE))
}


#' Plot methods' operating characteristics
#'
#' @description OCs vs TIE
#'
#' @param results_metrics_df The data frame containing the results and metrics.
#' @param metrics A list of metrics to be plotted.
#'
#' @return None
#'
#' @examples NA
operating_characteristics_vs_tie_plots <- function(results_metrics_df, metrics) {
  check_plot_packages()
  check_required_colnames(
    results_metrics_df,
    required_colnames_consumer,
    context = "operating_characteristics_vs_tie_plots"
  )

  # Get the list of case studies
  case_studies <- unique(results_metrics_df$case_study)

  treatment_effects <- c("partially_consistent", "consistent" ) #, "no_effect"


  for (case_study in case_studies) {
    results_metrics_df1 <- results_metrics_df[results_metrics_df$case_study == case_study,]

    # Get the list of target sample sizes
    target_sample_sizes <- unique(results_metrics_df1[, target_size_columns, drop = FALSE])

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

          for (source_denominator_change_factor in source_denominator_change_factors) {
            results_metrics_df4 <- results_metrics_df3 %>%
              dplyr::filter(
                matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor)
              )

            for (operating_characteristic in metrics){
              operating_characteristic_vs_tie(
                results_metrics_df4,
                case_study,
                target_sample_size$target_sample_size_control,
                target_sample_size$target_sample_size_treatment,
                treatment_effect = treatment_effect,
                operating_characteristic = operating_characteristic,
                source_denominator_change_factor = source_denominator_change_factor,
                target_to_source_std_ratio = target_to_source_std_ratio
              )
            }
          }
        }
      }
    }
  }
  return(invisible(NULL))
}

bayesian_operating_characteristic_vs_tie <- function(results_metrics_df,
                                                     case_study,
                                                     target_sample_size_control,
                                                     target_sample_size_treatment,
                                                     operating_characteristic,
                                                     source_denominator_change_factor,
                                                     target_to_source_std_ratio,
                                                     tie_type,
                                                     design_prior_type,
                                                     design_prior_type_tie) {

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

  if (tie_type == "average_tie") {
    results_df_tie <- results_metrics_df %>%
      dplyr::filter(
        design_prior_type == !!design_prior_type_tie
      )

    results_df <- results_df %>% select(-average_tie)
    results_df <- results_df %>%
      dplyr::filter(
        design_prior_type == !!design_prior_type
      )

    results_df$average_tie <- results_df_tie$average_tie

    if (design_prior_type_tie == "ui_design_prior") {
      design_prior_tie_label <- "UI design prior"
    } else if (design_prior_type_tie == "analysis_prior") {
      design_prior_tie_label <- "Analysis prior"
    } else if (design_prior_type_tie == "source_posterior") {
      design_prior_tie_label <- "Source posterior as design prior"
    } else {
      stop("Design prior type not supported.")
    }
  }

  # Process the rows to create parameter labels
  results_df$parameters_labels <- lapply(seq_len(nrow(results_df)), function(i) {
    return(process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE))
  })

  results_df$parameters_labels_not_latex <- lapply(seq_len(nrow(results_df)), function(i) {
    return(process_method_parameters_label(results_df[i, ], methods_labels, method_name = FALSE, as_latex = FALSE))
  })

  results_df$methods_parameters_labels <- lapply(seq_len(nrow(results_df)), function(i) {
    return(process_method_parameters_label(results_df[i, ], methods_labels, method_name = TRUE, as_latex = FALSE))
  })

  parameters_labels <- setNames(object = results_df$parameters_labels, nm = results_df$parameters_labels)
  methods <- format_results_df_methods(results_df)
  results_df$method <- methods
  methods <- unique(methods)

  if (tie_type == "average_tie") {
    filename <- paste0(
      case_study,
      "_", design_prior_type, "_",
      operating_characteristic$name,
      "_vs_", design_prior_type_tie, "_", tie_type, "_",
      target_size_label(target_sample_size_control, target_sample_size_treatment, "file")
    )
  } else {
    filename <- paste0(
      case_study,
      "_",
      operating_characteristic$name,
      "_vs_", tie_type, "_",
      target_size_label(target_sample_size_control, target_sample_size_treatment, "file")
    )
  }

  filename <- format_filename(
    filename = filename,
    case_study = case_study,
    target_to_source_std_ratio = target_to_source_std_ratio,
    source_denominator_change_factor = source_denominator_change_factor
  )
  file_path <- file.path(directory, filename)

  if (file.exists(paste0(file_path, ".pdf")) && file.exists(paste0(file_path, ".png")) && remake_figures == FALSE) {
    return()
  }

  if (!(operating_characteristic$name %in% colnames(results_df))) {
    return()
  }

  y_range <- sort(range(results_df[, operating_characteristic$name]))
  ymin <- y_range[1]
  ymax <- y_range[2]

  x_range <- sort(range(results_df[, tie_type]))
  xmin <- x_range[1]
  xmax <- x_range[2]

  x_width <- xmax - xmin
  y_width <- ymax - ymin

  y_cap_size <- x_width * relative_error_cap_width / 2
  x_cap_size <- y_width * relative_error_cap_width / 2
  cap_size <- min(x_cap_size, y_cap_size)
  x_cap_size <- cap_size
  y_cap_size <- cap_size

  if (design_prior_type == "ui_design_prior") {
    design_prior_label <- "UI design prior"
  } else if (design_prior_type == "analysis_prior") {
    design_prior_label <- "Analysis prior"
  } else if (design_prior_type == "source_posterior") {
    design_prior_label <- "Source posterior as design prior"
  } else {
    stop("Design prior type not supported.")
  }

  plot_title <- sprintf(
    "%s, %s, %s",
    stringr::str_to_title(case_study),
    target_size_label(target_sample_size_control, target_sample_size_treatment, "latex"),
    design_prior_label
  )

  plot_title <- format_title(
    title = plot_title,
    case_study = case_study,
    target_to_source_std_ratio = target_to_source_std_ratio,
    source_denominator_change_factor = source_denominator_change_factor
  )

  # Split data by 'method'
  results_df_split <- split(results_df, results_df$method)

  # Unique parameters and methods
  unique_parameters <- unique(results_df$parameters_labels)
  unique_methods <- sort(unique(results_df$method))

  # Shapes for methods, keyed by method rather than handed out by position
  shapes <- method_shape_map(unique_methods)

  # Determine x-axis variable and label
  if (tie_type == "tie") {
    x_var <- "tie"
    x_label <- "TIE"
  } else {
    x_var <- "average_tie"
    x_label <- paste0("Average TIE with ", design_prior_tie_label)
  }

  if (all(is.na(results_df[[x_var]])) || all(is.na(results_df[[operating_characteristic$name]]))){
    return()
  }

  # Initialize the plot
  plt <- ggplot2::ggplot(mapping = ggplot2::aes_string(x = x_var, y = operating_characteristic$name))

  # Iterate over each method and add layers
  plt <- plt + purrr::imap(results_df_split, function(df_subset, method_name) {
    df_subset <- df_subset[!is.na(df_subset[[operating_characteristic$name]]), ]

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
      ),
      # Define the color scale for parameters within this method
      ## A method with one key needs no separate heading: the method name
      ## goes on the key itself, so it reads "Pooling" rather than a
      ## "Pooling" title above a lone unlabelled symbol. Where that key
      ## also carries parameters they follow the name.
      ggplot2::scale_color_manual(
        values = method_colors,
        labels = vs_tie_key_labels(method_name, method_parameters),
        name = if (length(method_parameters) == 1) NULL else method_name,
        drop = TRUE,
        ## Up to three entries per row: a method with nine parameter values
        ## (RMP's weights) becomes a 3x3 block rather than a column nine
        ## tall, which is what made the legend taller than the panel.
        guide = ggplot2::guide_legend(
          override.aes = list(shape = method_shape),
          ncol = min(3, length(method_parameters)),
          byrow = TRUE
        )
      )
    )


    # Add horizontal error bars if tie_type is "tie"
    if (tie_type == "tie") {
      layers <- c(layers, list(
        ggplot2::geom_errorbarh(
          data = df_subset,
          ggplot2::aes(
            xmin = conf_int_tie_lower,
            xmax = conf_int_tie_upper,
            y = !!ggplot2::sym(operating_characteristic$name),
            color = factor(parameters_labels, levels = method_parameters)
          ),
          height = x_cap_size
        )
      ))
    }

    # Reset the color scale for the next method
    layers <- c(layers, list(ggnewscale::new_scale_color()))

    return(layers)
  })

  # Add shape scale for methods and labels
  plt <- plt +
    ggplot2::scale_shape_manual(
      values = shapes,
      name = "Methods",
      guide = ggplot2::guide_legend(order = 1, nrow = 2)
    ) +
    ggplot2::labs(
      title = plot_title,
      x = x_label,
      y = operating_characteristic$label
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
  ## The 1.35 was sized for the old 6pt legend text. The legend sits under
  ## the panel and is the widest thing on the page, so its width tracks
  ## the font: at 12pt it needs about half again as much room, and
  ## without this it is clipped at both ends.
  legend_scale <- max(1, small_text_size / 8)
  fig_width_in <- plot.size[1] * 1.35 * legend_scale
  ## The stacked legend is tall - one row per method group - so the
  ## page has to grow with it, or the panel is squeezed to a sliver.
  ## The legend is a compact grid rather than one group per row, so the
  ## page only needs the font's worth of extra room, not the large
  ## allowance a stacked legend wanted.
  fig_height_in <- plot.size[2] * 1.5 * legend_scale


  composed <- vs_tie_compose(plt, names(results_df_split), fig_width_in, fig_height_in)
  plt <- composed$plot
  fig_width_in <- composed$width
  fig_height_in <- composed$height * 0.85

  # Export the plots
  export_plots(plt, file_path, fig_width_in, fig_height_in, type = "pdf", adjust_theme = FALSE)
  return(export_plots(plt, file_path, fig_width_in, fig_height_in, type = "png", adjust_theme = FALSE))

}


#' Plot methods operating characteristics
#'
#' @description OCs vs TIE.
#'
#'
#' @param results_bayes_df The Bayesian results data frame.
#' @param results_freq_df The frequentist results data frame, from which the type I error rates are read.
#' @param bayesian_metrics List of Bayesian metrics, as defined by sourcing
#'   `metrics_config.R`.
#' @return None
#'
#' @examples NA
bayesian_operating_characteristics_vs_tie_plots <- function(results_bayes_df, results_freq_df, bayesian_metrics) {
  check_plot_packages()
 common_columns <- intersect(names(results_bayes_df), names(results_freq_df))

  # Select the results that correspond to TIE computation.
  results_freq_df_tie <- results_freq_df[results_freq_df$target_treatment_effect == results_freq_df$theta_0, c(
    common_columns,
    c(
      "tie",
      "conf_int_tie_lower",
      "conf_int_tie_upper"
    )
  )]

  common_columns <- intersect(names(results_freq_df_tie),common_columns)

  # Now merge results_bayes_df with the filtered results_freq_selected on the common columns specified in 'by'
  results_bayes_df_merged <- results_bayes_df %>%
  left_join(results_freq_df_tie, by =  common_columns)

  # Get the list of case studies
  case_studies <- unique(results_bayes_df_merged$case_study)

  for (case_study in case_studies) {
    results_metrics_df1 <- results_bayes_df_merged[results_bayes_df_merged$case_study == case_study,]

    # Get the list of target sample sizes
    target_sample_sizes <- unique(results_metrics_df1[, target_size_columns, drop = FALSE])

    for (size_index in seq_len(nrow(target_sample_sizes))) {
      target_sample_size <- target_sample_sizes[size_index, , drop = FALSE]
      results_metrics_df2 <- dplyr::semi_join(results_metrics_df1, target_sample_size, by = target_size_columns)

      target_to_source_std_ratios <- unique(results_metrics_df2$target_to_source_std_ratio)

      for (target_to_source_std_ratio in target_to_source_std_ratios) {
        results_metrics_df3 <- results_metrics_df2 %>%
          dplyr::filter(
            matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio)
          )

        source_denominator_change_factors <- unique(results_metrics_df3$source_denominator_change_factor)

        for (source_denominator_change_factor in source_denominator_change_factors) {
          results_metrics_df4 <- results_metrics_df3 %>%
            dplyr::filter(
              matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor)
            )

          design_prior_types <- unique(results_metrics_df4$design_prior_type)

          #design_prior_types <- "analysis_prior" # FIXME
          for (operating_characteristic in bayesian_metrics){
            for (design_prior_type in design_prior_types){
              bayesian_operating_characteristic_vs_tie(
                results_metrics_df4,
                case_study,
                target_sample_size$target_sample_size_control,
                target_sample_size$target_sample_size_treatment,
                operating_characteristic = operating_characteristic,
                source_denominator_change_factor = source_denominator_change_factor,
                target_to_source_std_ratio = target_to_source_std_ratio,
                tie_type = "tie",
                design_prior_type = design_prior_type
              )

              for (design_prior_type_tie in design_prior_types){
                bayesian_operating_characteristic_vs_tie(
                  results_metrics_df4,
                  case_study,
                  target_sample_size$target_sample_size_control,
                  target_sample_size$target_sample_size_treatment,
                  operating_characteristic = operating_characteristic,
                  source_denominator_change_factor = source_denominator_change_factor,
                  target_to_source_std_ratio = target_to_source_std_ratio,
                  tie_type = "average_tie",
                  design_prior_type = design_prior_type,
                  design_prior_type_tie = design_prior_type_tie
                )
              }
            }
          }
        }
      }
    }
  }
  return(invisible(NULL))
}
