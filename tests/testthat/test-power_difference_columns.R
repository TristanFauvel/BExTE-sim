power_difference_input <- data.frame(
  success_proba = 0.80,
  conf_int_success_proba_lower = 0.76,
  conf_int_success_proba_upper = 0.84,
  frequentist_power_at_equivalent_tie = 0.60,
  frequentist_power_at_equivalent_tie_lower = 0.55,
  frequentist_power_at_equivalent_tie_upper = 0.65
)

test_that("add_power_difference_columns subtracts the comparator point estimate", {
  result <- add_power_difference_columns(power_difference_input)

  expect_equal(result$success_proba, 0.20)
})

test_that("add_power_difference_columns widens the interval by the comparator uncertainty", {
  result <- add_power_difference_columns(power_difference_input)

  # Subtracting only the comparator point estimate would leave the interval
  # exactly as wide as the success probability interval.
  expect_gt(
    result$conf_int_success_proba_upper - result$conf_int_success_proba_lower,
    0.84 - 0.76
  )
  expect_equal(
    result$conf_int_success_proba_upper - result$conf_int_success_proba_lower,
    2 * sqrt(0.04^2 + 0.05^2)
  )
})

test_that("add_power_difference_columns brackets the point difference", {
  result <- add_power_difference_columns(power_difference_input)

  expect_lt(result$conf_int_success_proba_lower, result$success_proba)
  expect_gt(result$conf_int_success_proba_upper, result$success_proba)
})

test_that("add_power_difference_columns accepts a correlation between the two estimates", {
  independent <- add_power_difference_columns(power_difference_input)
  correlated <- add_power_difference_columns(power_difference_input, correlation = 0.5)

  expect_equal(correlated$success_proba, independent$success_proba)
  expect_lt(
    correlated$conf_int_success_proba_upper - correlated$conf_int_success_proba_lower,
    independent$conf_int_success_proba_upper - independent$conf_int_success_proba_lower
  )
})

test_that("add_power_difference_columns handles several rows at once", {
  two_rows <- rbind(power_difference_input, power_difference_input)
  two_rows$frequentist_power_at_equivalent_tie[2] <- 0.30
  two_rows$frequentist_power_at_equivalent_tie_lower[2] <- 0.25
  two_rows$frequentist_power_at_equivalent_tie_upper[2] <- 0.35

  result <- add_power_difference_columns(two_rows)

  expect_equal(result$success_proba, c(0.20, 0.50))
  expect_equal(nrow(result), 2)
})

test_that("add_power_difference_columns propagates missing comparator bounds", {
  missing_bounds <- power_difference_input
  missing_bounds$frequentist_power_at_equivalent_tie_lower <- NA_real_

  result <- add_power_difference_columns(missing_bounds)

  expect_equal(result$success_proba, 0.20)
  expect_true(is.na(result$conf_int_success_proba_upper))
})
