baseline_data <- data.frame(
  drift = c(-0.2, 0, 0.2),
  nominal_frequentist_power_separate = c(0.3, 0.5, 0.7),
  nominal_frequentist_power_separate_lower = c(0.27, 0.47, 0.67),
  nominal_frequentist_power_separate_upper = c(0.33, 0.53, 0.73)
)

has_errorbar_layer <- function(plt) {
  any(vapply(plt$layers, function(layer) inherits(layer$geom, "GeomErrorbar"), logical(1)))
}

draw_baseline <- function(data, ...) {
  plot_baseline_success_proba_vs_xvar(
    plt = ggplot2::ggplot(),
    data = data,
    xvar_name = "drift",
    selected_metric_name = "nominal_frequentist_power_separate",
    cap_size = 0.05,
    markersize = 2,
    join_points = FALSE,
    label = "Nominal Pr(Success)",
    ...
  )
}

bound_columns <- list(
  selected_metric_uncertainty_lower = "nominal_frequentist_power_separate_lower",
  selected_metric_uncertainty_upper = "nominal_frequentist_power_separate_upper"
)

test_that("Monte Carlo bounds are drawn as error bars", {
  plt <- do.call(draw_baseline, c(list(data = baseline_data), bound_columns))

  expect_true(has_errorbar_layer(plt))
})

test_that("closed-form bounds are left off the plot", {
  # compute_freq_power() returns c(power, power) when it has a closed form, so
  # every bar would be zero-height: a bare cap tick through each marker.
  degenerate_data <- baseline_data
  degenerate_data$nominal_frequentist_power_separate_lower <-
    degenerate_data$nominal_frequentist_power_separate
  degenerate_data$nominal_frequentist_power_separate_upper <-
    degenerate_data$nominal_frequentist_power_separate

  plt <- do.call(draw_baseline, c(list(data = degenerate_data), bound_columns))

  expect_false(has_errorbar_layer(plt))
})

test_that("results predating the bound columns still plot", {
  legacy_data <- baseline_data[, c("drift", "nominal_frequentist_power_separate")]

  plt <- do.call(draw_baseline, c(list(data = legacy_data), bound_columns))

  expect_false(has_errorbar_layer(plt))
  expect_true(any(vapply(plt$layers, function(layer) inherits(layer$geom, "GeomPoint"), logical(1))))
})

test_that("all-missing bounds are treated as absent", {
  na_data <- baseline_data
  na_data$nominal_frequentist_power_separate_lower <- NA_real_
  na_data$nominal_frequentist_power_separate_upper <- NA_real_

  plt <- do.call(draw_baseline, c(list(data = na_data), bound_columns))

  expect_false(has_errorbar_layer(plt))
})

test_that("callers that pass no bounds get the previous points-only baseline", {
  plt <- draw_baseline(data = baseline_data)

  expect_false(has_errorbar_layer(plt))
})

test_that("the baseline still rejects an empty dataframe", {
  expect_error(
    do.call(draw_baseline, c(list(data = baseline_data[0, ]), bound_columns)),
    "The dataframe is empty"
  )
})
