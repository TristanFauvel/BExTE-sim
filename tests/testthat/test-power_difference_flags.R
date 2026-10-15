power_row <- function(success_proba, lower, upper,
                      comparator, comparator_lower, comparator_upper) {
  data.frame(
    success_proba = success_proba,
    conf_int_success_proba_lower = lower,
    conf_int_success_proba_upper = upper,
    frequentist_power_at_equivalent_tie = comparator,
    frequentist_power_at_equivalent_tie_lower = comparator_lower,
    frequentist_power_at_equivalent_tie_upper = comparator_upper
  )
}

test_that("flag_power_differences adds a difference interval that accounts for both estimates", {
  flagged <- flag_power_differences(
    power_row(0.80, 0.76, 0.84, 0.60, 0.55, 0.65)
  )

  expect_equal(flagged$power_difference, 0.20)
  expect_equal(
    flagged$power_difference_upper - flagged$power_difference_lower,
    2 * sqrt(0.04^2 + 0.05^2)
  )
})

test_that("flag_power_differences declares a gain only when the difference interval excludes zero", {
  clear_gain <- flag_power_differences(
    power_row(0.80, 0.76, 0.84, 0.60, 0.55, 0.65)
  )
  overlapping <- flag_power_differences(
    power_row(0.62, 0.56, 0.68, 0.60, 0.55, 0.65)
  )

  expect_true(clear_gain$power_gain)
  expect_false(clear_gain$power_loss)
  expect_false(overlapping$power_gain)
  expect_false(overlapping$power_loss)
})

test_that("flag_power_differences declares a loss when the difference interval lies below zero", {
  flagged <- flag_power_differences(
    power_row(0.60, 0.55, 0.65, 0.80, 0.76, 0.84)
  )

  expect_true(flagged$power_loss)
  expect_false(flagged$power_gain)
})

test_that("flag_power_differences applies the same threshold to gains and to losses", {
  # A gain and the mirror-image loss of identical magnitude must be flagged
  # identically; the two analyses previously used different decision rules.
  gain <- flag_power_differences(
    power_row(0.80, 0.76, 0.84, 0.60, 0.55, 0.65)
  )
  loss <- flag_power_differences(
    power_row(0.60, 0.55, 0.65, 0.80, 0.76, 0.84)
  )

  expect_equal(gain$power_gain, loss$power_loss)
  expect_equal(gain$power_difference, -loss$power_difference)
  expect_equal(
    gain$power_difference_upper - gain$power_difference_lower,
    loss$power_difference_upper - loss$power_difference_lower
  )
})

test_that("flag_power_differences flags nothing when the comparator interval is missing", {
  flagged <- flag_power_differences(
    power_row(0.80, 0.76, 0.84, 0.60, NA_real_, NA_real_)
  )

  expect_false(flagged$power_gain)
  expect_false(flagged$power_loss)
  expect_true(is.na(flagged$power_difference_lower))
})

test_that("flag_power_differences returns an empty result for an empty input", {
  empty <- flag_power_differences(
    power_row(0.8, 0.76, 0.84, 0.6, 0.55, 0.65)[0, ]
  )

  expect_equal(nrow(empty), 0)
  expect_true(all(
    c("power_difference", "power_difference_lower", "power_difference_upper",
      "power_gain", "power_loss") %in% names(empty)
  ))
})
