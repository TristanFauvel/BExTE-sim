test_that("matches_or_na keeps matching values and missing values", {
  expect_equal(matches_or_na(c(1, 2, NA), 1), c(TRUE, FALSE, TRUE))
})

test_that("matches_or_na works inside dplyr::filter with an unquoted value", {
  df <- data.frame(factor = c(1, 2, NA), id = 1:3)
  factor <- 2
  out <- dplyr::filter(df, matches_or_na(factor, !!factor))
  expect_equal(out$id, c(2L, 3L))
})
