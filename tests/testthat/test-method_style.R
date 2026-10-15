## A method has to look the same in every figure it appears in, otherwise the
## reader cannot carry a symbol from one panel to the next. Both halves of the
## identity used to be derived from whatever happened to be in the frame:
## shapes were `shape_codes[1:length(unique_methods)]`, so a case study missing
## one method shifted every later method's shape, and the colours came from an
## unseeded `sample()`, so they differed between two runs of the same figure.
## The style table replaces both with a lookup keyed by the method.

methods_plots_config <- function() {
  env <- new.env()
  ## The file assigns with <<-, which needs an enclosing binding to find.
  eval(parse(system.file("conf/methods_plots_config.R", package = "BExTEsim")), envir = env)
  env
}

## A stand-in for the run's methods_dict, holding the two shapes the real
## configs take: one varying parameter (RMP's weight) and two (NPP's mean and
## standard deviation).
style_test_dict <- function() {
  list(
    RMP = list(
      prior_weight = list(
        range = seq(0, 1, length.out = 11),
        parameter_label = "w",
        parameter_notation = "$w$",
        type = "continuous"
      ),
      initial_prior = list(
        range = list("noninformative"),
        parameter_label = "pi_0",
        parameter_notation = "$\\pi_0$",
        type = "categorical"
      )
    ),
    NPP = list(
      power_parameter_mean = list(
        range = list(0.25, 0.5),
        parameter_label = "xi_gamma",
        parameter_notation = "$\\xi_\\gamma$",
        type = "continuous",
        display = TRUE
      ),
      power_parameter_std = list(
        range = list(0.1, 0.2, 0.4),
        parameter_label = "sigma_gamma",
        parameter_notation = "$\\sigma_\\gamma$",
        type = "continuous",
        display = TRUE
      )
    )
  )
}

test_that("the style table covers every labelled method", {
  ## shape_codes held ten entries for eleven methods, so a run with all of
  ## them silently handed the eleventh an NA shape.
  config <- methods_plots_config()
  styles <- get("methods_style", envir = config)
  labels <- get("methods_labels", envir = config)

  expect_setequal(names(styles), names(labels))
})

test_that("no two methods share a shape or a hue", {
  config <- methods_plots_config()
  styles <- get("methods_style", envir = config)

  shapes <- vapply(styles, function(style) style$shape, numeric(1))
  hues <- vapply(styles, function(style) style$hue, character(1))

  expect_equal(anyDuplicated(shapes), 0)
  expect_equal(anyDuplicated(hues), 0)
})

test_that("a method can be looked up by key or by display label", {
  ## The plot code overwrites results_df$method with the display label before
  ## it builds the scales, so both spellings have to resolve.
  styles <- list(RMP = list(shape = 16, hue = "#0072B2"))

  by_key <- method_style("RMP", styles = styles)
  by_label <- method_style("RMP", styles = styles)

  expect_equal(by_key$shape, 16)
  expect_equal(by_key, by_label)

  config <- methods_plots_config()
  real <- get("methods_style", envir = config)
  labels <- get("methods_labels", envir = config)

  expect_equal(
    method_style("conditional_power_prior", styles = real, labels = labels),
    method_style("Conditional PP", styles = real, labels = labels)
  )
})

test_that("an unknown method is an error rather than a silent NA", {
  expect_error(
    method_style("not_a_method", styles = list(RMP = list(shape = 16, hue = "#0072B2"))),
    "not_a_method"
  )
})

test_that("a method keeps its shape whichever other methods share the figure", {
  styles <- list(
    RMP = list(shape = 16, hue = "#0072B2"),
    NPP = list(shape = 17, hue = "#D55E00"),
    pooling = list(shape = 3, hue = "#000000")
  )

  full <- method_shape_map(c("RMP", "NPP", "pooling"), styles = styles)
  ## The case study that has no NPP results used to shift pooling onto NPP's
  ## shape, because the codes were handed out by position.
  without_npp <- method_shape_map(c("RMP", "pooling"), styles = styles)

  expect_equal(full[["pooling"]], 3)
  expect_equal(without_npp[["pooling"]], full[["pooling"]])
  expect_equal(without_npp[["RMP"]], full[["RMP"]])
})

test_that("a parameter value keeps its colour when the figure shows fewer values", {
  dict <- style_test_dict()
  styles <- list(RMP = list(shape = 16, hue = "#0072B2"))

  nine <- method_parameter_colors(
    "RMP",
    sprintf("$w$ = %s", c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)),
    styles = styles, dict = dict
  )
  three <- method_parameter_colors(
    "RMP",
    sprintf("$w$ = %s", c(0.1, 0.5, 0.9)),
    styles = styles, dict = dict
  )

  expect_equal(three[["$w$ = 0.5"]], nine[["$w$ = 0.5"]])
  expect_equal(three[["$w$ = 0.1"]], nine[["$w$ = 0.1"]])
  expect_equal(three[["$w$ = 0.9"]], nine[["$w$ = 0.9"]])
})

test_that("the shades run from light to dark with the parameter value", {
  dict <- style_test_dict()
  styles <- list(RMP = list(shape = 16, hue = "#0072B2"))

  colors <- method_parameter_colors(
    "RMP",
    sprintf("$w$ = %s", c(0.1, 0.5, 0.9)),
    styles = styles, dict = dict
  )

  luminance <- function(hex) {
    rgb <- grDevices::col2rgb(hex)[, 1]
    sum(rgb * c(0.2126, 0.7152, 0.0722))
  }

  expect_gt(luminance(colors[["$w$ = 0.1"]]), luminance(colors[["$w$ = 0.5"]]))
  expect_gt(luminance(colors[["$w$ = 0.5"]]), luminance(colors[["$w$ = 0.9"]]))
})

test_that("a method with two varying parameters is ordered on both", {
  ## NPP displays a prior mean and a prior standard deviation, so the shade
  ## has to key on the pair rather than on a single number.
  dict <- style_test_dict()
  styles <- list(NPP = list(shape = 17, hue = "#D55E00"))

  all_six <- method_parameter_colors(
    "NPP",
    c(
      "$\\xi_\\gamma$ = 0.25, $\\sigma_\\gamma$ = 0.1",
      "$\\xi_\\gamma$ = 0.25, $\\sigma_\\gamma$ = 0.2",
      "$\\xi_\\gamma$ = 0.25, $\\sigma_\\gamma$ = 0.4",
      "$\\xi_\\gamma$ = 0.5, $\\sigma_\\gamma$ = 0.1",
      "$\\xi_\\gamma$ = 0.5, $\\sigma_\\gamma$ = 0.2",
      "$\\xi_\\gamma$ = 0.5, $\\sigma_\\gamma$ = 0.4"
    ),
    styles = styles, dict = dict
  )
  subset <- method_parameter_colors(
    "NPP",
    c(
      "$\\xi_\\gamma$ = 0.25, $\\sigma_\\gamma$ = 0.2",
      "$\\xi_\\gamma$ = 0.5, $\\sigma_\\gamma$ = 0.4"
    ),
    styles = styles, dict = dict
  )

  expect_length(unique(all_six), 6)
  for (label in names(subset)) {
    expect_equal(subset[[label]], all_six[[label]], info = label)
  }
})

test_that("a method with a single parameter value gets its base hue", {
  ## Pooling and Separate take no parameters; their one key should be the
  ## colour the style table names, not a washed-out end of a ramp.
  styles <- list(pooling = list(shape = 3, hue = "#000000"))

  colors <- method_parameter_colors("pooling", "", styles = styles, dict = list())

  expect_equal(unname(colors), "#000000")
})

test_that("labels with no recoverable number still get stable colours", {
  ## Categorical parameters fall back to sorted labels; the guarantee there is
  ## reproducibility between runs, which the unseeded sample() never gave.
  styles <- list(RMP = list(shape = 16, hue = "#0072B2"))
  labels <- c("$\\pi_0$ = vague", "$\\pi_0$ = unit information")

  first <- method_parameter_colors("RMP", labels, styles = styles, dict = list())
  second <- method_parameter_colors("RMP", rev(labels), styles = styles, dict = list())

  expect_equal(first[labels], second[labels])
  expect_length(unique(first), 2)
})

test_that("colours are returned named by the labels they were asked for", {
  ## scale_color_manual() matches values to levels by name, which is what
  ## makes the scale a lookup instead of a position.
  dict <- style_test_dict()
  styles <- list(RMP = list(shape = 16, hue = "#0072B2"))
  labels <- sprintf("$w$ = %s", c(0.2, 0.8))

  colors <- method_parameter_colors("RMP", labels, styles = styles, dict = dict)

  expect_setequal(names(colors), labels)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", colors)))
})

test_that("parameterless methods retain visible colors in plot scales", {
  styles <- list(
    pooling = list(shape = 3, hue = "#000000"),
    separate = list(shape = 4, hue = "#7F7F7F"),
    EB_PP = list(shape = 8, hue = "#E69F00")
  )
  df <- data.frame(parameters_labels_not_latex = c("", ""))
  df$parameters_labels <- list(expression(""), expression(""))
  for (method in names(styles)) {
    colors <- method_parameter_color_map(df, method, styles = styles, dict = list())
    expect_equal(unname(colors), styles[[method]]$hue)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = 1, y = 1,
      colour = factor(parameters_labels))) +
      ggplot2::geom_point(shape = styles[[method]]$shape) +
      ggplot2::scale_colour_manual(values = colors)
    expect_equal(unique(ggplot2::ggplot_build(p)$data[[1]]$colour), styles[[method]]$hue)
  }
})
