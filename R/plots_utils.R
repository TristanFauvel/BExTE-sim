# Function to format numbers in scientific notation with 2 significant digits
format_num <- function(x, digits = 2, scientific = FALSE) {
  return(format(x, digits = digits, scientific = scientific))
}

#' Open a graphics device that can render the plot font
#'
#' @description Composing a figure measures text: `grid::convertWidth()` and
#'   friends do it directly, and `ggplot2::ggplotGrob()` does it again while it
#'   assembles the guides. That measurement runs on whatever graphics device is
#'   current. Setting `options(device = )` is not enough, because it only
#'   decides what gets opened when no device is open at all; a device the
#'   caller already had open is used as it is, and a plain `pdf()` device
#'   cannot load the Computer Modern CID font the plot themes ask for. A stray
#'   device left behind by earlier code - another test file, an interactive
#'   session, a previous figure - then breaks a figure that is itself fine.
#'
#'   Opening a cairo device here makes the measurement independent of what the
#'   caller left behind. It is a no-op when cairo is unavailable, or when the
#'   current device is already a cairo one, so nesting these costs nothing.
#'
#' @return A function that closes the device this opened and makes the
#'   caller's previous device current again. Call it from `on.exit()`.
#'
#' @keywords internal
use_font_capable_device <- function() {
  no_op <- function() invisible(NULL)
  if (!isTRUE(capabilities("cairo")) ||
      identical(names(grDevices::dev.cur()), "cairo_pdf")) {
    return(no_op)
  }

  previous <- grDevices::dev.cur()
  path <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(filename = path)
  opened <- grDevices::dev.cur()

  return(function() {
    if (opened %in% grDevices::dev.list()) {
      grDevices::dev.off(opened)
    }
    if (previous > 1 && previous %in% grDevices::dev.list()) {
      grDevices::dev.set(previous)
    }
    unlink(path)
    return(invisible(NULL))
  })
}


#' Are there parameters left to condition a plot or table on?
#'
#' @description The metric-versus-parameters plots and tables loop over one of
#'   a method's parameters at a time and condition on the rest, taken as
#'   `parameters[, -i]`. Whether any rest remain is not the same question as
#'   whether that frame has rows. A method carrying a single parameter -
#'   separate and pooling carry only `initial_prior` - loses its only column,
#'   and a data frame with no columns keeps every one of its rows, so `nrow()`
#'   still reports something to loop over. The loop then asks for a label for a
#'   row that holds nothing, which fails with `invalid subscript type 'list'`.
#'
#'   A remainder that has collapsed to a bare vector also counts as nothing,
#'   which is what these call sites have always done.
#'
#' @param other_parameters The remaining parameters, as returned by
#'   `parameters[, -i]`.
#'
#' @return `TRUE` when there is at least one parameter left to condition on.
#'
#' @keywords internal
has_other_parameters <- function(other_parameters) {
  return(!is.null(nrow(other_parameters)) &&
    nrow(other_parameters) > 0 &&
    ncol(other_parameters) > 0)
}


#' Function to set figure dimensions
#'
#' @param width The desired width of the figure (in pts).
#' @param fraction The fraction of the width to use (default is 1).
#' @param aspect_ratio The desired aspect ratio of the figure (default is golden ratio).
#' @return A vector containing the width and height of the figure.
#' @keywords internal
set_size <- function(width,
                     fraction = 1,
                     aspect_ratio = NULL) {
  # Width of figure (in pts)
  fig_width_pt <- width * fraction

  # Convert from pt to inches
  inches_per_pt <- 1 / 72.27

  # Golden ratio to set aesthetic figure height
  # https://disq.us/p/2940ij3
  golden_ratio <- (5 ^ 0.5 - 1) / 2

  if (is.null(aspect_ratio)) {
    aspect_ratio <- golden_ratio
  }

  # Figure width in inches
  fig_width_in <- fig_width_pt * inches_per_pt
  # Figure height in inches
  fig_height_in <- fig_width_in * aspect_ratio

  fig_dim <- c(fig_width_in, fig_height_in)

  return(fig_dim)
}


#' Function to convert parameters in dataframe to string
#'
#' @param method The method name.
#' @param parameters The parameters dataframe.
#' @return The string representation of the parameters.
#' @keywords internal
convert_params_to_str <- function(method, parameters) {
  # Round floating point parameters to 2 decimal places
  parameters <- purrr::map_if(parameters, is.numeric, round, 2)

  if (is.null(names(parameters))) {
    stop("parameters must be a dataframe")
  }
  # Convert parameters to string representation
  labels <- paste(Filter(Negate(is.null), sapply(names(parameters), function(key) {
    if (length(method[[key]]$range) > 1) {
      return(paste0(method[[key]]$parameter_label, "=", parameters[[key]]))
    }
    return(invisible(NULL))
  })), collapse = "_")
  return(labels)
}

#' Append a parameter string to a figure or table filename
#'
#' @description Joins the two with an underscore, so that a value ending the
#'   stem ("sample_size=123") stays separate from the first parameter
#'   ("gamma=0.5"). A method without a varied parameter has an empty string,
#'   and gets no trailing underscore.
#'
#' @param stem The filename so far.
#' @param parameters_str The output of [convert_params_to_str()].
#' @return The filename with the parameters appended.
#' @keywords internal
append_parameters_str <- function(stem, parameters_str) {
  if (length(parameters_str) == 0 || !nzchar(parameters_str)) {
    return(stem)
  }
  return(paste0(stem, "_", parameters_str))
}

#' The methods whose parameters nest under a heterogeneity prior
#'
#' @description The commensurate power prior and the plain commensurate prior
#' both take their `heterogeneity_prior` as a nested list of a family and its
#' hyperparameters, rather than as a scalar. That shape is what the figure and
#' table code has to special-case, so it is what this names - the two methods
#' differ in whether they also carry a power parameter, which is irrelevant
#' everywhere this predicate is used.
#'
#' @keywords internal
COMMENSURATE_METHODS <- c("commensurate_power_prior", "commensurate_prior")

#' Whether a method's parameters nest under a heterogeneity prior
#'
#' @param method Method name.
#' @return A single logical.
#' @keywords internal
is_commensurate_method <- function(method) {
  return(method %in% COMMENSURATE_METHODS)
}

#' Function to make labels from parameters. Return a label formatted in Tex, for example "$\\xi_\\gamma$ = 0.5, $\\sigma_\\gamma$ = 0.1"
#'
#' @param parameter The parameter dataframe.
#' @param method The method name.
#' @return The labels generated from the parameters.
#' @keywords internal
make_labels_from_parameters <- function(parameter, method) {
  # The heterogeneity-prior labels read their values from here, before the
  # rounding below: rounded to 2 decimal places, alpha = 1/1000 prints as 0.
  unrounded <- parameter
  colnames(unrounded) <- gsub("heterogeneity_prior\\.", "", colnames(unrounded))
  prior_value <- function(name) {
    return(format(signif(as.numeric(unrounded[[name]]), 2), scientific = FALSE, drop0trailing = TRUE))
  }

  # Round floating point parameters to 2 decimal places
  parameter[sapply(parameter, is.numeric)] <- lapply(parameter[sapply(parameter, is.numeric)], round, 2)

  # Select parameters with non-singleton ranges or display keyword = TRUE
  selection <- list()

  colnames(parameter) <- gsub("heterogeneity_prior\\.", "", colnames(parameter))
  if (is_commensurate_method(method)){
    if (is.null(parameter$family)){
      if (!is.null(parameter$alpha) && !is.na(parameter$alpha)){
        parameter$family = "inverse_gamma"
      } else if (!is.null(parameter$std_dev) && !is.na(parameter$std_dev)){
        parameter$family = "half_normal"
      } else if (!is.null(parameter$location) && !is.na(parameter$location)){
        parameter$family = "cauchy"
      }
    }
    # Each label names the quantity its prior is placed on, as the Stan
    # model and the quadrature in R/vectorised_commensurate_prior.R do: the
    # half-normal on tau, the inverse gamma on tau^2, the Cauchy on log(tau).
    if (parameter$family == "half_normal"){
      label <- paste0("$\\tau \\sim HN(", prior_value("std_dev"), ")$")
    } else if (parameter$family == "inverse_gamma"){
      label <- paste0("$\\tau^2 \\sim IG(\\alpha = ", prior_value("alpha"), ", \\beta = ", prior_value("beta"), ")$")
    } else if (parameter$family == "cauchy"){
      label <- paste0("$\\log \\tau \\sim Cauchy(", prior_value("location"), ", ", prior_value("scale"), ")$")
    } else {
      stop("Heterogeneity prior family not implemented.")
    }
  } else {
    for (parameter_name in names(parameter)) {
      if (is.null(methods_dict[[method]][[parameter_name]])){
        stop("Parameter not listed in the method configuration.")
      }

      if (!is.null(methods_dict[[method]][[parameter_name]]$display)) {
        selection <- c(selection, methods_dict[[method]][[parameter_name]]$display)
      } else {
        selection <- c(selection, length(methods_dict[[method]][[parameter_name]]$range) > 1)
      }
    }

    selection <- unlist(selection)

    label <- ""
    if (!is.null(selection)) {
      i <- 0
      parameters_to_display <- parameter[, selection, drop = FALSE]
      for (parameter_name in colnames(parameters_to_display)) {
        # Combine column names and values and convert selected parameters to label

        new_label <- sprintf(
          "%s = %s",
          methods_dict[[method]][[parameter_name]]$parameter_notation,
          as.character(round(as.numeric(parameter[[parameter_name]]),2))
        )

        if (i > 0) {
          label <- paste(label, new_label, sep = ", ")
        } else {
          label <- paste0(label, new_label)
        }

        i <- i + 1
      }
    }
  }

  return(label)
}



process_method_parameters_label <- function(row,
                                            methods_labels,
                                            method_name = TRUE,
                                            parameters_colname = "parameters",
                                            as_latex = TRUE){
  # Process the row of a results dataframe to create a Method + Parameters label
  method <- unlist(row["method"])

  if (is.null(methods_labels[[method]])) {
    stop(
      paste0(
        "Method label is not defined. Add a label in methods_config.R for method ",
        method
      )
    )
  }

  parameters_df <- get_parameters(row[parameters_colname])

  # Return a label formatted in Tex, for example "$\\xi_\\gamma$ = 0.5, $\\sigma_\\gamma$ = 0.1"
  param_label <- make_labels_from_parameters(parameter = parameters_df, method = method)


  if (method_name == TRUE) {
    # Add the method name to the label
    if (param_label == "") {
      full_label <- methods_labels[[method]]$label
    } else {
      full_label <- paste(methods_labels[[method]]$label, param_label, sep = ", ")
    }
  } else {
    if (param_label == "") {
      full_label <- ""
    } else {
      full_label <- param_label
    }
  }

  if (as_latex){
    full_label <- latex2exp::TeX(full_label)
  }

  if (is.null(full_label)) {
    stop("Label is null.")
  }

  return(full_label)
}


#' Function to return a dataframe of parameters from json strings
#'
#' @description This function returns, for a dataframe containing parameters (in json strings) for a given method, an R dataframe.
#'
#' @param parameters_df The dataframe containing parameters in JSON strings.
#' @return An R dataframe with parsed parameters.
#' @export
get_parameters <- function(parameters_df) {
  if ("posterior_parameters" %in% colnames(parameters_df)) {
    # Process the posterior parameters
    parameters_df$parameters <- parameters_df$posterior_parameters
    parameters_df <- parameters_df[, !(names(parameters_df) %in% c("posterior_parameters")), drop = FALSE]
  }

  if (is.null(parameters_df$parameters)){
    stop("The input parameters dataframe does not contain parameters.")
  }

  if (all(unique(parameters_df$parameters) == "[]")) {
    return(parameters_df)
  }

  N <- nrow(parameters_df)

  # Check that the input is a dataframe
  if (!is.data.frame(parameters_df)) {
    if (length(parameters_df) == 1) {
      parameters_df <- gsub("'", "\"", parameters_df)
      json_parameters_df <- jsonlite::fromJSON(parameters_df, simplifyDataFrame = TRUE)
      if (length(names(json_parameters_df)) == 1 &&
          (json_parameters_df) == "parameters") {
        json_parameters_df <- unlist(unlist(json_parameters_df$parameters))
      } else {
        json_parameters_df <- unlist(json_parameters_df)
      }
      return(json_parameters_df)
    } else {
      stop("Input is not a dataframe")
    }
  } else {
    # Replace single quotes with double quotes
    parameters_df$parameters <- gsub("'", "\"", parameters_df$parameters)
  }

  results_list <- list()

  if (nrow(parameters_df) > 0) {
    # A results frame repeats each parameter string across every scenario it
    # was run on (hundreds of rows each), so each distinct string is parsed
    # once and its values shared by the rows that carry it.
    distinct_parameters <- unique(parameters_df$parameters)
    parsed <- lapply(distinct_parameters, function(parameters) {
      json_parameters_df <- data.frame(jsonlite::fromJSON(parameters, simplifyDataFrame = TRUE))

      if (is.null(names(json_parameters_df))) {
        return(NA)
      } else if (length(names(json_parameters_df)) == 1 &&
                 names(json_parameters_df) == "parameters") {
        return(unlist(unlist(json_parameters_df$parameters)))
      } else {
        return(unlist(json_parameters_df))
      }
    })

    results_list <- parsed[match(parameters_df$parameters, distinct_parameters)]
  }

  # Combine the list of data frames into a single data frame, removing any NULL entries
  results_list <- results_list[!sapply(results_list, is.null)]

  # Get the number of columns for each element in the list
  num_cols <- lengths(results_list)

  # Check if all elements have the same number of columns
  if (length(results_list) != 0 && length(unique(num_cols)) != 1) {

    # Step 1: Extract all unique column names across elements in the list
    all_columns <- unique(unlist(lapply(results_list, names)))

    # Step 2: Ensure each element has the same columns by adding missing ones as NA
    results_list_aligned <- lapply(results_list, function(x) {
      missing_cols <- setdiff(all_columns, names(x))
      x[missing_cols] <- NA
      x <- x[all_columns]  # Reorder to keep column order consistent
      return(x)
    })

    # Step 3: Bind rows without error
    results_list <- results_list_aligned

    # Now, results_df will have consistent columns

    # stop(
    #   "Error: The elements in the list have different numbers of columns. The function should only be applied to a dataframe containing parameters for a single method."
    # )
  }

  results_df <- do.call(rbind, results_list)

  if (!is.data.frame(results_df)) {
    results_df <- data.frame(results_df)
  }

  if (nrow(results_df) != N) {
    stop("Output does not have the same number of rows as the input df.")
  }

  return(results_df)
}

# Define a function to convert strings to numeric or boolean if possible
convert_if_possible <- function(x) {
  output <- vector("list", length(x)) # Initialize output as a list
  for (i in seq_along(x)) {
    if (is.character(x[i])) {
      # Try to convert to numeric
      numeric_value <- suppressWarnings(as.numeric(x[i]))
      if (!is.na(numeric_value)) {
        output[[i]] <- numeric_value
        next
      }

      # Try to convert to boolean
      if (tolower(x[i]) == "true") {
        output[[i]] <- TRUE
      } else if (tolower(x[i]) == "false") {
        output[[i]] <- FALSE
      } else {
        output[[i]] <- x[i] # Keep as character if no conversion occurs
      }
    }
  }
  # Return the original value if no conversion was possible
  return(unlist(output))
}

export_plots <- function(plt,
                         file_path,
                         fig_width_in,
                         fig_height_in,
                         type = "pdf", forest_plot = FALSE, adjust_theme = TRUE,
                         bg = NULL, theme_extra = NULL) {

  if (!(forest_plot) && adjust_theme == TRUE){
    plt <- plt + ggplot2::theme_bw() + ggplot2::theme(
      axis.text = ggplot2::element_text(family = font, size = text_size),
      axis.text.y = ggplot2::element_text(family = font, size = small_text_size),
      axis.text.x = ggplot2::element_text(family = font, size = small_text_size),
      axis.title = ggplot2::element_text(family = font, size = text_size),
      plot.title = ggplot2::element_text(family = font, size = text_size),
      legend.text = ggplot2::element_text(family = font, size = small_text_size),
      legend.title = ggplot2::element_text(family = font, size = text_size),
      legend.key = ggplot2::element_blank(),
      # strip.background = element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  }

  # Apply caller styling after theme_bw(), which resets legend placement.
  if (!is.null(theme_extra)) {
    plt <- plt + theme_extra
  }

  if (type == "pdf") {
    return(ggplot2::ggsave(
      filename = paste0(file_path, ".pdf"),
      plot = plt,
      device = cairo_pdf,
      width = fig_width_in,
      height = fig_height_in,
      dpi = dpi,
      bg = bg
    ))
  } else {
    ## PNGs are the copies that get dropped into slides, documents and
    ## issue threads, where a transparent background renders as whatever is
    ## behind it - black, in a dark viewer. ggsave() falls back to the
    ## theme's background when bg is NULL, and the forest plots skip the
    ## theme_bw() adjustment above, so they came out transparent. White
    ## unless the caller asked for something specific.
    return(ggplot2::ggsave(
      filename = paste0(file_path, ".png"),
      plot = plt,
      width = fig_width_in,
      height = fig_height_in,
      dpi = dpi,
      bg = if (is.null(bg)) "white" else bg
    ))
  }

}

format_title <- function(title, case_study, target_to_source_std_ratio = NA, source_denominator_change_factor = NA, as_latex = FALSE){
  ## A ratio of 1 is the default - the target and source standard
  ## deviations agree - so naming it lengthens every title without
  ## saying anything, the same way a denominator change factor of 1 does.
  if (case_study %in% c("botox", "dapagliflozin") &&
        !is.na(target_to_source_std_ratio) &&
        target_to_source_std_ratio != 1){
    title <- paste0(title, ", $\\sigma_T/\\sigma_S = $", target_to_source_std_ratio)
  }

  ## A change factor of 1 is the default - the source denominator used as
  ## it stands - so saying so in the title adds length without adding
  ## information. Only a factor that actually changes something is named.
  if (!is.na(source_denominator_change_factor) &&
        source_denominator_change_factor != 1){
    if (!(case_study %in% c("botox", "dapagliflozin", "aprepitant"))){
      title = paste0(title, ", Source denominator change factor = ",source_denominator_change_factor)
    }
  }

  if (as_latex == FALSE){
    title <- latex2exp::TeX(title)
  }

  return(title)
}

format_filename <- function(filename, case_study, target_to_source_std_ratio = NA, source_denominator_change_factor = NA){
  if (case_study %in% c("botox", "dapagliflozin")){
    filename <- paste0(filename,
                       "_target_to_source_std_ratio=",
                       target_to_source_std_ratio
    )
  }

  if (!is.na(source_denominator_change_factor)){
    if (!(case_study %in% c("botox", "dapagliflozin", "aprepitant"))){
      filename <- paste0(filename, "_source_denominator_change_factor=",  source_denominator_change_factor)
    }
  }

  # Convert to lowercase and replace spaces with underscores
  filename <- tolower(filename)
  filename <- gsub(" ", "_", filename)

  return(filename)
}

format_results_df_parameters <- function(results_df, include_method_name = TRUE){
  labels <- c()
  for (i in seq_len(nrow(results_df))) {
    method <- unlist(results_df[i, "method"])

    parameters_df <- get_parameters(results_df[i, "parameters", drop = FALSE])

    param_label <- make_labels_from_parameters(parameter = parameters_df, method = method)

    if (include_method_name){
      if (param_label == "") {
        full_label <- methods_labels[[method]]$label
      } else {
        full_label <- paste(methods_labels[[method]]$label, param_label, sep = ", ")
      }
    } else {
      full_label <- param_label
    }
    # Append the label to the labels vector
    labels <- c(labels, full_label)
  }
  return(labels)
}

format_results_df_methods <- function(results_df){
  methods <- c()
  for (i in seq_len(nrow(results_df))) {
    method <- unlist(results_df[i, "method"])
    methods <- c(methods, methods_labels[[method]]$label)
  }
  return(methods)
}


convert_CI_columns <- function(table_data_df) {
  # Get the names of all columns that start with "conf_int"
  conf_int_cols <- grep("^conf_int", colnames(table_data_df), value = TRUE)

  # Extract unique parameter names (e.g., "power_parameter_mean")
  unique_params <- unique(sub("conf_int_(lower|upper)_", "", conf_int_cols))

  # Dynamically combine parameter values with their confidence intervals
  for (param in unique_params) {
    # Find the lower and upper confidence interval columns
    lower_col <- paste0("conf_int_lower_", param)
    upper_col <- paste0("conf_int_upper_", param)

    # Check if the parameter column exists
    if (param %in% colnames(table_data_df)) {
      # Combine the value and confidence interval into a single string
      table_data_df[[param]] <- paste0(
        table_data_df[[param]], " [", table_data_df[[lower_col]], ", ", table_data_df[[upper_col]], "]"
      )
    } else {
      # If the parameter column doesn't exist, create a new one with just the CI
      table_data_df[[param]] <- paste0(
        "[", table_data_df[[lower_col]], ", ", table_data_df[[upper_col]], "]"
      )
    }
  }
  # Drop the original confidence interval columns
  cols_to_keep <- !grepl("^conf_int", colnames(table_data_df))
  table_data_df <- table_data_df[, cols_to_keep, drop = FALSE]
  return(table_data_df)
}


## ---------------------------------------------------------------------------
## Method plotting style
##
## A method has to look the same in every figure it appears in, otherwise a
## reader cannot carry a symbol from one panel to the next. The shape carries
## the method; within a method, the parameter values are shades of the method's
## hue. Both come from the methods_style table rather than from the contents of
## the figure being drawn.
## ---------------------------------------------------------------------------

#' Resolve a configuration object the plot code keeps in the global environment
#'
#' @description conf/methods_plots_config.R and the run's methods configuration
#'   assign with `<<-`, so the plot code reads them as free variables the way it
#'   already reads `font`, `dpi` and `textwidth`. Returning NULL rather than
#'   erroring lets the style helpers be called with an explicit table in tests,
#'   where no configuration has been loaded.
#'
#' @param name The object to look up.
#' @return The object, or NULL when no configuration has been loaded.
#' @keywords internal
style_config_object <- function(name) {
  return(mget(name, envir = environment(), ifnotfound = list(NULL), inherits = TRUE)[[1]])
}

style_methods_labels <- function() style_config_object("methods_labels")

style_methods_dict <- function() style_config_object("methods_dict")

#' Resolve a method to its configuration key
#'
#' @description The plot code overwrites `results_df$method` with the display
#'   label before it builds the scales, so a method reaches the style helpers
#'   as either `"conditional_power_prior"` or `"Conditional PP"`. Both have to
#'   land on the same entry.
#'
#' @param method A method key or display label.
#' @param labels The method label table.
#' @return The method key, or the input unchanged when no label matches.
#' @keywords internal
method_key <- function(method, labels = style_methods_labels()) {
  method <- as.character(method)[1]

  if (is.null(labels) || method %in% names(labels)) {
    return(method)
  }

  display <- vapply(labels, function(entry) {
    if (is.null(entry$label)) return(NA_character_) else return(as.character(entry$label)[1])
  }, character(1))

  hit <- names(labels)[!is.na(display) & display == method]

  if (length(hit) == 1) return(hit) else return(method)
}

#' Fixed shape and hue for a method
#'
#' @param method A method key or display label.
#' @param styles The style table, from conf/methods_plots_config.R.
#' @param labels The method label table, used to resolve a display label.
#' @return A list with `shape` and `hue`.
#' @keywords internal
method_style <- function(method,
                         styles = methods_style,
                         labels = style_methods_labels()) {
  style <- styles[[method_key(method, labels = labels)]]

  if (is.null(style)) {
    stop(sprintf(
      "No plot style defined for method '%s'. Add an entry to methods_style in conf/methods_plots_config.R.",
      as.character(method)[1]
    ))
  }

  return(style)
}

#' Shapes for a set of methods, keyed by method
#'
#' @description A named vector, so `scale_shape_manual()` matches by name
#'   instead of by position. The codes used to be handed out as
#'   `shape_codes[1:length(unique_methods)]` over the methods present in one
#'   figure, which meant a case study missing a method shifted the shape of
#'   every method after it.
#'
#' @param methods Method keys or display labels.
#' @param styles The style table.
#' @param labels The method label table.
#' @return A named numeric vector of plotting characters.
#' @keywords internal
method_shape_map <- function(methods,
                             styles = methods_style,
                             labels = style_methods_labels()) {
  methods <- unique(as.character(methods))

  return(stats::setNames(
    vapply(
      methods,
      function(method) as.numeric(method_style(method, styles = styles, labels = labels)$shape),
      numeric(1)
    ),
    methods
  ))
}

#' Blend a colour towards another
#'
#' @param color The colour to start from.
#' @param towards The colour to move towards.
#' @param weight How far to move, between 0 and 1.
#' @return A hex colour string.
#' @keywords internal
style_blend <- function(color, towards, weight) {
  mixed <- (1 - weight) * grDevices::col2rgb(color)[, 1] +
    weight * grDevices::col2rgb(towards)[, 1]

  return(grDevices::rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255))
}

#' A light-to-dark ramp through a method's hue
#'
#' @param hue The method's base colour, which sits at the middle of the ramp.
#' @param light_weight How far towards white the light end sits.
#' @param dark_weight How far towards black the dark end sits.
#' @return A function mapping positions in \eqn{[0, 1]} to hex colours.
#' @keywords internal
method_hue_ramp <- function(hue, light_weight = 0.75, dark_weight = 0.45) {
  ramp <- grDevices::colorRamp(c(
    style_blend(hue, "#FFFFFF", light_weight),
    hue,
    style_blend(hue, "#000000", dark_weight)
  ))

  return(function(positions) {
    channels <- ramp(pmin(pmax(positions, 0), 1))
    return(grDevices::rgb(channels[, 1], channels[, 2], channels[, 3], maxColorValue = 255))
  })
}

#' Numbers a parameter label assigns to its parameters
#'
#' @description Only numbers introduced by "=" count. A notation such as the
#'   initial prior's carries a digit of its own (pi nought) that has nothing to
#'   do with the value the parameter took.
#'
#' @param label A parameter label, as plain text.
#' @return The numbers in the order they appear, possibly none.
#' @keywords internal
style_label_numbers <- function(label) {
  matches <- regmatches(label, gregexpr("=\\s*-?[0-9]+\\.?[0-9]*", label))[[1]]

  if (length(matches) == 0) {
    return(numeric(0))
  }

  return(suppressWarnings(as.numeric(sub("^=\\s*", "", matches))))
}

#' The numeric parameter ranges a method's configuration declares
#'
#' @description Mirrors the selection `make_labels_from_parameters()` makes, so
#'   the ranges line up with the numbers that end up in the labels: a parameter
#'   is shown when it says so through `display`, or when its range holds more
#'   than one value. Non-numeric ranges are skipped - a categorical parameter
#'   has no ramp to sit on.
#'
#' @param key The method key.
#' @param dict The run's methods configuration.
#' @return A list of numeric vectors, in configuration order.
#' @keywords internal
style_numeric_ranges <- function(key, dict) {
  spec <- dict[[key]]

  if (is.null(spec)) {
    return(list())
  }

  ranges <- list()

  for (name in names(spec)) {
    parameter <- spec[[name]]
    values <- suppressWarnings(as.numeric(unlist(parameter$range)))

    if (length(values) == 0 || anyNA(values)) {
      next
    }

    displayed <- if (!is.null(parameter$display)) {
      isTRUE(parameter$display)
    } else {
      length(parameter$range) > 1
    }

    if (displayed) {
      ranges[[name]] <- values
    }
  }

  return(ranges)
}

#' A tuple of parameter values as a lookup key
#'
#' @description Rounded to the two decimals `make_labels_from_parameters()`
#'   rounds to, so a configured value and the number parsed back out of its
#'   label agree despite floating point.
#'
#' @param values A numeric vector.
#' @return A single string.
#' @keywords internal
style_tuple_key <- function(values) {
  return(paste(format(round(values, 2), nsmall = 2, trim = TRUE), collapse = "|"))
}

#' Where each parameter label sits on its method's ramp
#'
#' @description Anchored on the range the configuration declares rather than on
#'   the values one figure happens to show. A figure plotting w in \{0.1, 0.5,
#'   0.9\} and one plotting all nine weights put w = 0.5 at the same place, and
#'   so give it the same shade; ranking the values present would not.
#'
#'   Labels carrying no recoverable number - categorical parameters - fall back
#'   to sorted order, which is at least reproducible between runs.
#'
#' @param key The method key.
#' @param parameter_labels The labels to place, as plain text.
#' @param dict The run's methods configuration.
#' @return A numeric vector of positions in \eqn{[0, 1]}.
#' @keywords internal
method_parameter_positions <- function(key, parameter_labels, dict) {
  values <- lapply(parameter_labels, style_label_numbers)
  width <- unique(lengths(values))

  if (length(width) == 1 && width > 0) {
    ranges <- style_numeric_ranges(key, dict)

    if (length(ranges) == width) {
      if (width == 1) {
        span <- range(ranges[[1]])

        if (is.finite(diff(span)) && diff(span) > 0) {
          observed <- vapply(values, function(value) value[1], numeric(1))
          return(pmin(pmax((observed - span[1]) / diff(span), 0), 1))
        }
      } else {
        ## More than one parameter varies, so there is no single number to
        ## normalise. Walk the configured grid in order instead and take each
        ## combination's place in it.
        grid <- expand.grid(ranges, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        grid <- grid[do.call(order, as.list(grid)), , drop = FALSE]

        grid_keys <- apply(grid, 1, function(row) style_tuple_key(as.numeric(row)))
        label_keys <- vapply(values, style_tuple_key, character(1))
        place <- match(label_keys, grid_keys)

        if (!anyNA(place) && length(grid_keys) > 1) {
          return((place - 1) / (length(grid_keys) - 1))
        }
      }
    }
  }

  ranks <- rank(parameter_labels, ties.method = "first")

  return((ranks - 1) / max(length(ranks) - 1, 1))
}

#' Colours for the parameter values of one method
#'
#' @description Shades of the method's hue, light for low parameter values and
#'   dark for high ones. Returned named by the labels asked for, so that
#'   `scale_color_manual()` is a lookup rather than a position - which is what
#'   keeps a parameter value the same colour from one figure to the next. These
#'   colours used to come from an unseeded `sample()` over a 50-colour palette,
#'   so they differed between two runs of the same figure.
#'
#' @param method A method key or display label.
#' @param parameter_labels The labels to colour, as plain text.
#' @param styles The style table.
#' @param dict The run's methods configuration.
#' @param labels The method label table.
#' @return A named character vector of hex colours.
#' @keywords internal
method_parameter_colors <- function(method,
                                    parameter_labels,
                                    styles = methods_style,
                                    dict = style_methods_dict(),
                                    labels = style_methods_labels()) {
  hue <- method_style(method, styles = styles, labels = labels)$hue
  parameter_labels <- unique(as.character(parameter_labels))

  if (length(parameter_labels) == 0) {
    return(stats::setNames(character(0), character(0)))
  }

  ## One key needs no ramp: Pooling and Separate take no parameters, and their
  ## single entry should be the colour the table names rather than a washed-out
  ## end of a ramp.
  if (length(parameter_labels) == 1) {
    return(stats::setNames(hue, parameter_labels))
  }

  positions <- method_parameter_positions(
    method_key(method, labels = labels),
    parameter_labels,
    dict = dict
  )

  return(stats::setNames(method_hue_ramp(hue)(positions), parameter_labels))
}

#' Base hues for a set of methods, keyed by method
#'
#' @description The companion to `method_shape_map()`, for the figures that
#'   colour by method rather than by parameter value.
#'
#' @param methods Method keys or display labels.
#' @param styles The style table.
#' @param labels The method label table.
#' @return A named character vector of hex colours.
#' @keywords internal
method_hue_map <- function(methods,
                           styles = methods_style,
                           labels = style_methods_labels()) {
  methods <- unique(as.character(methods))

  return(stats::setNames(
    vapply(
      methods,
      function(method) method_style(method, styles = styles, labels = labels)$hue,
      character(1)
    ),
    methods
  ))
}

#' Colours for the parameter values present in one method's rows
#'
#' @description The plot code carries two spellings of each parameter label: the
#'   plotmath expression the legend is drawn with, and the plain text it was
#'   built from. Only the plain text still holds the numbers that place a value
#'   on its method's ramp, while the scale has to be keyed by the expression the
#'   layers map to, so the two are paired up row by row here.
#'
#' @param df_subset One method's rows, carrying `parameters_labels` and
#'   `parameters_labels_not_latex`.
#' @param method The method key or display label.
#' @param ... Passed on to `method_parameter_colors()`.
#' @return A character vector of hex colours, named by the plotmath labels.
#' @keywords internal
method_parameter_color_map <- function(df_subset, method, ...) {
  display <- as.character(df_subset$parameters_labels)
  plain <- as.character(df_subset$parameters_labels_not_latex)

  if (length(plain) != length(display)) {
    stop("parameters_labels_not_latex must accompany parameters_labels.")
  }

  keep <- !duplicated(display)
  display <- display[keep]
  plain <- plain[keep]

  colors <- method_parameter_colors(method, plain, ...)

  ## Empty parameter labels are valid for Pooling, Separate and EBPP.
  ## Character indexing never matches an empty name; match() does.
  return(stats::setNames(unname(colors[match(plain, names(colors))]), display))
}

#' Axis breaks that always show the nominal type-I error
#'
#' `pretty()` picks breaks that ignore the nominal TIE, so forcing the
#' nominal value in alongside them leaves two labels sitting on top of each
#' other. `check.overlap` resolves that collision the wrong way round: it
#' draws labels leftmost-first and discards the nominal tick, which is the
#' one the plot is read against. Drop the neighbouring pretty breaks
#' instead, so the nominal value is the only label in its neighbourhood.
#'
#' @param nominal_tie The nominal type-I error rate, or `NULL` when the
#'   scale has no nominal value to mark.
#' @return A function of the scale limits returning a sorted break vector.
#' @keywords internal
nominal_tie_breaks <- function(nominal_tie) {
  force(nominal_tie)

  return(function(limits) {
    breaks <- pretty(limits)

    if (length(nominal_tie) != 1 || !is.finite(nominal_tie)) {
      return(breaks)
    }

    # Half a step is roughly the width of a break label at these font sizes.
    spacing <- if (length(breaks) > 1) min(diff(breaks)) else diff(range(limits))
    clear <- abs(breaks - nominal_tie) >= 0.5 * spacing

    return(sort(unique(c(breaks[clear], nominal_tie))))
  })
}

#' Whether a pair of bound columns carries Monte Carlo uncertainty
#'
#' The power baselines are estimated one of two ways. A closed-form power has
#' no Monte Carlo error, and the compute_freq_power() family reports that as a
#' degenerate interval, `c(power, power)`. Drawing those bounds would put a
#' zero-height error bar - a bare cap tick - through every marker, so a plot
#' asks this first and adds the layer only when the bounds say something.
#'
#' Bound columns missing from `data` count as no uncertainty: results written
#' before the bounds were recorded still have to plot.
#'
#' @param data The dataframe backing the layer.
#' @param lower,upper Names of the bound columns, or `NULL`.
#'
#' @return `TRUE` when both columns are present and some row spans a non-zero
#'   width, `FALSE` otherwise.
#' @keywords internal
has_monte_carlo_uncertainty <- function(data, lower, upper) {
  if (is.null(lower) || is.null(upper)) {
    return(FALSE)
  }

  lower <- as.character(lower)
  upper <- as.character(upper)

  if (!all(c(lower, upper) %in% colnames(data))) {
    return(FALSE)
  }

  width <- data[[upper]] - data[[lower]]

  return(any(is.finite(width) & width > 0))
}

#' Keep the rows of one treatment effect scenario
#'
#' @param results_df A data frame with `target_treatment_effect` and
#'   `source_treatment_effect_estimate`.
#' @param treatment_effect "consistent" (target equals the source estimate) or
#'   "partially_consistent" (target is half the source estimate).
#'
#' @return The matching rows of `results_df`.
#' @keywords internal
filter_to_treatment_effect <- function(results_df, treatment_effect) {
  if (treatment_effect == "consistent") {
    results_df %>%
      dplyr::filter(target_treatment_effect == source_treatment_effect_estimate)
  } else if (treatment_effect == "no_effect") {
    stop("The treatment effect should not be non-positive")
  } else if (treatment_effect == "partially_consistent") {
    results_df %>%
      dplyr::filter(abs(
        target_treatment_effect - source_treatment_effect_estimate / 2
      ) < 1e-4)
  } else {
    stop('treatment_effect must be "consistent", "no_effect" or "partially_consistent"')
  }
}

#' Keep the parameter combinations a figure shows
#'
#' Sorts the rows by method, then keeps for each method the rows whose
#' parameters take the important values in `methods_dict` (the full range when
#' none are set). Commensurate methods keep only the half-normal prior with
#' standard deviation 1 and the inverse-gamma prior with alpha 0.3333.
#'
#' @param results_df A data frame with `method` and `parameters` columns.
#'
#' @return The kept rows of `results_df`, sorted by method.
#' @keywords internal
filter_to_important_parameter_values <- function(results_df) {
  boolean_filter <- c()
  results_df <- results_df %>% arrange(method)
  for (method in unique(results_df$method)) {
    parameters <- get_parameters(results_df[results_df$method == method, "parameters", drop = FALSE])
    filter_per_method <- rep(TRUE, nrow(parameters))
    for (key in names(methods_dict[[method]])) {
      if (is_commensurate_method(method)) {
        filter <- (parameters$heterogeneity_prior.family == "half_normal" & parameters$heterogeneity_prior.std_dev == 1) |
          (parameters$heterogeneity_prior.family == "inverse_gamma" & parameters$heterogeneity_prior.alpha == 0.3333)
      } else {
        important_parameters_values <- methods_dict[[method]][[key]][["important_values"]]
        if (is.null(important_parameters_values)) {
          important_parameters_values <- methods_dict[[method]][[key]][["range"]]
        }
        filter <- parameters[, key] %in% important_parameters_values
      }
      filter_per_method <- filter & filter_per_method
    }
    boolean_filter <- c(boolean_filter, filter_per_method)
  }
  results_df[boolean_filter, ]
}

#' Call a function once for each scenario of a results table
#'
#' Loops over the case studies, the methods, optionally the target sample
#' sizes, the target to source standard deviation ratios and the source
#' denominator change factors present in `results_metrics_df`. For each
#' combination it calls `fn` with the rows of that scenario. Rows where a ratio
#' or a factor is `NA` belong to every scenario.
#'
#' @param results_metrics_df A data frame of results.
#' @param fn A function of two arguments: the rows of the scenario, and a list
#'   with `case_study`, `method`, `target_sample_size_control`,
#'   `target_sample_size_treatment` (both `NULL` unless `by_target_size`),
#'   `target_to_source_std_ratio` and `source_denominator_change_factor`.
#' @param by_target_size Loop over the target sample sizes as well.
#'
#' @return `NULL`, invisibly.
#' @keywords internal
for_each_scenario <- function(results_metrics_df, fn, by_target_size = FALSE) {
  for (case_study in unique(results_metrics_df$case_study)) {
    df1 <- results_metrics_df[results_metrics_df$case_study == case_study, ]
    for (method in unique(results_metrics_df$method)) {
      df2 <- df1[df1$method == method, ]

      if (by_target_size) {
        target_sample_sizes <- unique(df2[, target_size_columns, drop = FALSE])
        sizes <- lapply(seq_len(nrow(target_sample_sizes)), function(i) as.list(target_sample_sizes[i, , drop = FALSE]))
      } else {
        sizes <- list(list(target_sample_size_control = NULL, target_sample_size_treatment = NULL))
      }

      for (size in sizes) {
        df3 <- df2
        if (by_target_size) {
          df3 <- df2[
            df2$target_sample_size_control == size$target_sample_size_control &
              df2$target_sample_size_treatment == size$target_sample_size_treatment,
          ]
        }

        for (target_to_source_std_ratio in unique(df3$target_to_source_std_ratio)) {
          df4 <- df3 %>%
            dplyr::filter(
              matches_or_na(target_to_source_std_ratio, !!target_to_source_std_ratio)
            )

          for (source_denominator_change_factor in unique(df4$source_denominator_change_factor)) {
            df5 <- df4 %>%
              dplyr::filter(
                matches_or_na(source_denominator_change_factor, !!source_denominator_change_factor)
              )

            fn(df5, list(
              case_study = case_study,
              method = method,
              target_sample_size_control = size$target_sample_size_control,
              target_sample_size_treatment = size$target_sample_size_treatment,
              target_to_source_std_ratio = target_to_source_std_ratio,
              source_denominator_change_factor = source_denominator_change_factor
            ))
          }
        }
      }
    }
  }
  invisible(NULL)
}

#' Packages that the plotting code needs
#'
#' They are in Suggests, so that users who only fit priors and compute
#' operating characteristics do not have to install them.
#' @noRd
PLOT_PACKAGES <- c("ggplot2", "ggnewscale", "gridExtra", "latex2exp",
                   "purrr", "stringr", "viridis")

#' Check that the packages for plotting are installed
#'
#' Every plotting entry point calls this first.
#'
#' @return `NULL`, invisibly. Raises an error if a package is missing and the
#'   user does not install it.
#' @noRd
check_plot_packages <- function() {
  return(rlang::check_installed(PLOT_PACKAGES, reason = "to draw plots."))
}
