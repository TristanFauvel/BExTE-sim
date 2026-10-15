test_that("sweet_spot_determination finds one bounded interval", {
  result <- sweet_spot_determination(
    x_values = 1:10,
    y_values = c(3, 5, 7, 9, 11, 8, 6, 4, 2, 0),
    reference_value = 5,
    larger_is_better = TRUE
  )

  expected <- data.frame(
    start = 2,
    end = 7.5,
    width = 5.5,
    total_width = 5.5
  )
  expect_equal(result, expected)
})

test_that("sweet_spot_determination returns multiple intervals", {
  expect_warning(
    result <- sweet_spot_determination(
      x_values = 1:10,
      y_values = c(3, 5, 7, 9, 11, 8, 6, 4, 6, 8),
      reference_value = 6,
      larger_is_better = TRUE
    ),
    "More than one sweet spot"
  )

  expected <- data.frame(
    start = c(2.5, 9),
    end = c(7, 10),
    width = c(4.5, 1),
    total_width = c(5.5, 5.5)
  )
  expect_equal(result, expected)
})

test_that("sweet_spot_determination handles values below the reference", {
  result <- sweet_spot_determination(
    x_values = 1:10,
    y_values = c(3, 5, 7, 9, 11, 8, 6, 4, 2, 0),
    reference_value = 12,
    larger_is_better = TRUE
  )

  expected <- data.frame(start = NA, end = NA, width = 0, total_width = 0)
  expect_equal(result, expected)
})

test_that("sweet_spot_determination handles empty values", {
  expect_warning(
    result <- sweet_spot_determination(
      x_values = 1:10,
      y_values = numeric(),
      reference_value = 5,
      larger_is_better = TRUE
    ),
    "Input vectors must have non-zero length"
  )

  expected <- data.frame(start = NA, end = NA, width = NA)
  expect_equal(result, expected)
})

test_that("sweet_spot_determination treats equality as favorable", {
  result <- sweet_spot_determination(
    x_values = 1:10,
    y_values = rep(5, 10),
    reference_value = 5,
    larger_is_better = TRUE
  )

  expected <- data.frame(start = 1, end = 10, width = 9, total_width = 9)
  expect_equal(result, expected)
})

## A time-to-event case study simulates several designs - dropout, event time
## distribution, treatment delay - over the same drift grid. sweet_spot() used
## to pool them into one curve per method, which put several points at every
## drift and then failed binding one results row per design to sweet spots
## found on the mixture ("arguments imply differing number of rows: 6, 4").

two_design_results <- function() {
  drift <- seq(-1, 1, by = 0.25)
  grid <- expand.grid(
    drift = drift,
    method = c("separate", "borrowing"),
    dropout_probability = c(0, 0.1),
    stringsAsFactors = FALSE
  )
  # The borrowing method beats the separate analysis where |drift| is below
  # a design-specific bound: 0.5 without dropout, 0.25 with it.
  bound <- ifelse(grid$dropout_probability == 0, 0.5, 0.25)
  mse <- ifelse(grid$method == "separate", 1, 0.5 + abs(grid$drift) - bound + 0.5)

  data.frame(
    grid,
    parameters = ifelse(grid$method == "separate", "[]", "{'weight': [0.5]}"),
    target_treatment_effect = grid$drift,
    case_study = "example",
    target_sample_size_control = 100,
    target_sample_size_treatment = 100,
    source_denominator_change_factor = 1,
    target_to_source_std_ratio = 1,
    control_drift = 0,
    event_time_distribution = "exponential",
    treatment_delay = 0,
    mse = mse,
    conf_int_mse_lower = mse,
    conf_int_mse_upper = mse,
    stringsAsFactors = FALSE
  )
}

test_that("sweet_spot finds one sweet spot per time-to-event design", {
  result <- sweet_spot(
    two_design_results(),
    metrics = list(list(name = "mse", larger_is_better = FALSE))
  )
  borrowing <- result[result$method == "borrowing", ]
  borrowing <- borrowing[order(borrowing$dropout_probability), ]

  expect_equal(nrow(borrowing), 2)
  expect_equal(borrowing$dropout_probability, c(0, 0.1))
  expect_equal(unlist(borrowing$sweet_spot_lower), c(-0.5, -0.25))
  expect_equal(unlist(borrowing$sweet_spot_upper), c(0.5, 0.25))
})

test_that("several sweet spots of one curve are kept on one row", {
  collapsed <- collapse_sweet_spots(data.frame(
    start = c(-1, 0.5), end = c(-0.5, 1), width = c(0.5, 0.5), total_width = c(1, 1)
  ))

  expect_equal(nrow(collapsed), 1)
  expect_equal(collapsed$start[[1]], c(-1, 0.5))
  expect_equal(collapsed$end[[1]], c(-0.5, 1))

  single <- data.frame(start = -1, end = 1, width = 2, total_width = 2)
  expect_identical(collapse_sweet_spots(single), single)
})
