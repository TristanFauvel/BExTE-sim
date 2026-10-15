test_that("check_required_colnames accepts a frame carrying extra columns", {
  df <- data.frame(
    method = "conjugate",
    drift = 0.5,
    extra_bookkeeping = 1,
    stringsAsFactors = FALSE
  )

  expect_silent(check_required_colnames(df, c("method", "drift")))
})


test_that("check_required_colnames reports the columns a frame is missing", {
  df <- data.frame(method = "conjugate", stringsAsFactors = FALSE)

  expect_error(
    check_required_colnames(df, c("method", "drift")),
    "drift"
  )
})


test_that("check_required_colnames type checks the required columns", {
  df <- data.frame(method = 1, drift = 0.5, extra = "x", stringsAsFactors = FALSE)

  expect_error(
    check_required_colnames(df, c("method", "drift")),
    "method.*expected character.*found numeric"
  )
})


test_that("check_required_colnames does not type check the extra columns", {
  df <- data.frame(method = "conjugate", drift = 0.5, stringsAsFactors = FALSE)
  # `case_study` is in the shared spec as character, but is not required here.
  df$case_study <- 1

  expect_silent(check_required_colnames(df, c("method", "drift")))
})


test_that("check_required_colnames names the caller in its error", {
  df <- data.frame(method = "conjugate", stringsAsFactors = FALSE)

  expect_error(
    check_required_colnames(df, c("method", "drift"), context = "power_vs_tie"),
    "power_vs_tie"
  )
})
