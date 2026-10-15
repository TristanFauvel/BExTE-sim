test_that("filter_to_important_parameter_values keeps important values, sorted by method", {
  assign("methods_dict", list(
    pooling = list(foo = list(range = c(1, 2, 3), important_values = c(1, 3))),
    separate = list(bar = list(range = c(5, 6)))
  ), envir = globalenv())
  withr::defer(rm("methods_dict", envir = globalenv()))

  df <- data.frame(
    method = c("separate", "pooling", "pooling", "pooling", "separate"),
    parameters = c("s5", "p1", "p2", "p3", "s6"),
    stringsAsFactors = FALSE
  )
  values <- c(s5 = 5, s6 = 6, p1 = 1, p2 = 2, p3 = 3)
  testthat::local_mocked_bindings(
    get_parameters = function(parameters_df) {
      v <- unname(values[parameters_df$parameters])
      data.frame(foo = v, bar = v)
    },
    is_commensurate_method = function(method) FALSE
  )

  out <- filter_to_important_parameter_values(df)
  expect_equal(out$parameters, c("p1", "p3", "s5", "s6"))
  expect_equal(out$method, c("pooling", "pooling", "separate", "separate"))
})
