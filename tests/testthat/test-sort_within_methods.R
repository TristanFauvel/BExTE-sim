test_that("sort_within_methods sorts each method by its parameters in place", {
  df <- data.frame(
    method = c("a", "b", "a", "b"),
    parameters = c("a2", "b2", "a1", "b1"),
    stringsAsFactors = FALSE
  )
  testthat::local_mocked_bindings(
    get_parameters = function(parameters_df) {
      data.frame(p = substring(parameters_df$parameters, 2))
    }
  )
  out <- sort_within_methods(df)
  expect_equal(out$method, df$method)
  expect_equal(out$parameters, c("a1", "b1", "a2", "b2"))
})
