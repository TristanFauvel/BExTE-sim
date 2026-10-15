## get_parameters() parses each distinct parameter string once and shares the
## values across the rows that repeat it; the rows must still come back in
## their own order, one per input row.

test_that("get_parameters returns one row per input row, in order", {
  parameters_df <- data.frame(parameters = c(
    "{'gamma': 0.5}", "{'gamma': 0.9}", "{'gamma': 0.5}",
    "{'gamma': 0.9}", "{'gamma': 0.1}", "{'gamma': 0.5}"
  ))

  parsed <- get_parameters(parameters_df)

  expect_equal(nrow(parsed), 6)
  expect_equal(as.numeric(parsed$gamma), c(0.5, 0.9, 0.5, 0.9, 0.1, 0.5))
})

test_that("get_parameters aligns rows whose parameter sets differ", {
  parameters_df <- data.frame(parameters = c(
    "{'a': 1, 'b': 2}", "{'a': 3}", "{'a': 1, 'b': 2}"
  ))

  parsed <- get_parameters(parameters_df)

  expect_equal(nrow(parsed), 3)
  expect_equal(as.numeric(parsed$a), c(1, 3, 1))
  expect_equal(as.numeric(parsed$b), c(2, NA, 2))
})
