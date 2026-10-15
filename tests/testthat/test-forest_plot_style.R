## The forest plots colour each point by the moment-based ESS it achieved,
## on a diverging scale centred on half the source study's equivalent per-arm
## size: white means borrowing half of what the source is worth, blue less,
## red more. Centring on a fixed, interpretable value rather than on whatever
## range a given figure spans is what makes the colours comparable between
## figures of the same case study.

test_that("the colour scale is centred on half the source study's per-arm size", {
  data <- data.frame(equivalent_source_sample_size_per_arm = rep(234, 5))

  expect_equal(forest_ess_midpoint(data), 117)
})

test_that("the midpoint falls back to the source arm sizes when the column is absent", {
  ## Older results directories predate equivalent_source_sample_size_per_arm.
  data <- data.frame(
    equivalent_source_sample_size_per_arm = rep(NA_real_, 3),
    source_sample_size_control = rep(235, 3),
    source_sample_size_treatment = rep(233, 3)
  )

  ## Harmonic mean of the two arms, halved.
  expected <- (2 * 235 * 233 / (235 + 233)) / 2
  expect_equal(forest_ess_midpoint(data), expected)
})

test_that("a results frame without an ESS column still plots", {
  ## Returning NULL rather than erroring keeps the plot usable on a frame that
  ## never recorded ess_moment; the layer falls back to plain black points.
  expect_null(forest_ess_values(data.frame(x = 1:3)))
  expect_null(forest_ess_values(data.frame(ess_moment = c(NA, NaN))))
  expect_equal(forest_ess_values(data.frame(ess_moment = c(1, 2))), c(1, 2))
})

test_that("small values get a power of ten factored into the axis label", {
  ## MSE around 0.02 reads as 1.5, 2.0, 2.5 under "MSE (10^-2)" rather than
  ## as 0.015, 0.020, 0.025 crammed under the ticks.
  scientific <- forest_scientific_x(c(0.015, 0.020, 0.025), "MSE")

  expect_false(is.null(scientific$scale))
  expect_false(identical(scientific$label, "MSE"))
})

test_that("ordinary-sized values are left alone", {
  ## A factor of 1 helps nobody, and a probability axis reads fine as is.
  for (values in list(c(0.2, 0.5, 0.9), c(1, 5, 20), c(50, 100, 200))) {
    scientific <- forest_scientific_x(values, "Pr(Study success)")
    expect_null(scientific$scale)
    expect_equal(scientific$label, "Pr(Study success)")
  }
})

test_that("an all-zero or empty axis does not trip the factoring", {
  expect_null(forest_scientific_x(numeric(0), "MSE")$scale)
  expect_null(forest_scientific_x(c(0, 0), "MSE")$scale)
  expect_null(forest_scientific_x(c(NA, Inf), "MSE")$scale)
})

## The three panels of a forest plot must have the same plotting region, so a
## distance along the x axis means the same thing in each. Only the first
## carries the method labels, so giving the panels fixed fractions of the page
## makes their grids unequal - the labelled one loses whatever its labels take.

forest_test_grob <- function(y_labels) {
  data <- data.frame(rows = 1:5, x = 1:5, lo = (1:5) - 0.1, hi = (1:5) + 0.1)
  plt <- ggplot2::ggplot(data, ggplot2::aes(y = rows, x = x, xmin = lo, xmax = hi)) +
    ggplot2::geom_pointrange()
  if (!y_labels) {
    plt <- plt + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  }
  ggplot2::ggplotGrob(plt)
}

test_that("the three panels get the same plotting width despite unequal labels", {
  grobs <- list(forest_test_grob(TRUE), forest_test_grob(FALSE), forest_test_grob(FALSE))
  available <- 6

  widths <- as.numeric(forest_panel_widths(grobs, available))
  overheads <- vapply(grobs, function(g) {
    as.numeric(grid::convertWidth(forest_nonpanel_width(g), "in"))
  }, numeric(1))
  panels <- widths - overheads

  ## The labelled panel needs more total width precisely because its labels
  ## take more room.
  expect_gt(overheads[1], overheads[2])
  expect_gt(widths[1], widths[2])

  expect_equal(panels[1], panels[2], tolerance = 1e-6)
  expect_equal(panels[2], panels[3], tolerance = 1e-6)
  ## The widths fill the page when there is room for them. When the
  ## per-panel floor binds - 12pt titles need about 2.2in each - the total
  ## legitimately exceeds the nominal width and the page is widened to
  ## hold it, which is what stops the panels overprinting each other.
  expect_gte(sum(widths), min(available, sum(widths)))
  if (sum(widths) <= available + 1e-6) {
    expect_equal(sum(widths), available, tolerance = 1e-6)
  }
})

test_that("the widths are absolute, not unit arithmetic", {
  ## `unit(x, "cm") + unit(1, "null")` is a unit arithmetic object, and
  ## arrangeGrob cannot distribute leftover space across those - every panel
  ## collapses onto its labels. Absolute inches are what keeps the layout.
  widths <- forest_panel_widths(list(forest_test_grob(TRUE), forest_test_grob(FALSE)), 6)

  expect_true(grid::is.unit(widths))
  expect_false(anyNA(as.numeric(widths)))
  expect_length(widths, 2)
})

test_that("labels wider than the page still leave a usable grid", {
  ## A floor, so an extreme label set cannot squeeze the panels out entirely.
  widths <- as.numeric(forest_panel_widths(list(forest_test_grob(TRUE)), 0.1))

  expect_gt(widths[1], 0.1)
})

## Under no effect the probability of declaring success IS the type I error
## rate, so that panel - and only that panel - carries the nominal rate as a
## dashed reference. On the power panels it would mark nothing, and on a
## relative plot the axis is a ratio, on which a rate has no meaning.

forest_style_fixture <- function() {
  data <- readRDS(testthat::test_path("fixtures", "forest_plot_freq.rds"))
  data <- data[data$drift == data$drift[1], , drop = FALSE]
  data$rows <- seq_len(nrow(data))
  data
}

vline_layers <- function(plt) {
  Filter(function(layer) inherits(layer$geom, "GeomVline"), plt$layers)
}

test_that("the no-effect panel marks the nominal type I error rate", {
  skip_if_not(exists("methods_labels"), "plot globals not sourced")
  data <- forest_style_fixture()

  plt <- forest_subplot(
    data, "No effect", TRUE, rlang::sym("success_proba"),
    "conf_int_success_proba_lower", "conf_int_success_proba_upper",
    "TIE", methods_labels, nominal_tie_line = 0.025
  )

  lines <- vline_layers(plt)
  expect_length(lines, 1)
  expect_equal(unique(lines[[1]]$data$xintercept), 0.025)
  ## Dashed, so it reads apart from the dotted reference lines.
  expect_equal(lines[[1]]$aes_params$linetype, "dashed")
})

test_that("the power panels carry no nominal-rate line", {
  skip_if_not(exists("methods_labels"), "plot globals not sourced")
  data <- forest_style_fixture()

  plt <- forest_subplot(
    data, "Partially consistent", FALSE, rlang::sym("success_proba"),
    "conf_int_success_proba_lower", "conf_int_success_proba_upper",
    "Power", methods_labels
  )

  expect_length(vline_layers(plt), 0)
})
