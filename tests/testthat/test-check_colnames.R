test_that("check_colnames accepts a dataframe whose columns match the type spec", {
  df <- data.frame(
    method = "conjugate",
    drift = 0.5,
    stringsAsFactors = FALSE
  )

  expect_silent(
    check_colnames(
      df,
      c("method", "drift"),
      types = c(method = "character", drift = "numeric")
    )
  )
})


test_that("check_colnames reports a column whose type differs from the spec", {
  df <- data.frame(
    method = "conjugate",
    drift = "0.5",
    stringsAsFactors = FALSE
  )

  expect_error(
    check_colnames(
      df,
      c("method", "drift"),
      types = c(method = "character", drift = "numeric")
    ),
    "drift.*expected numeric.*character"
  )
})


test_that("check_colnames reports every mismatched column, not only the first", {
  df <- data.frame(
    method = 1,
    drift = "0.5",
    stringsAsFactors = FALSE
  )

  error_message <- tryCatch(
    check_colnames(
      df,
      c("method", "drift"),
      types = c(method = "character", drift = "numeric")
    ),
    error = function(e) conditionMessage(e)
  )

  expect_match(error_message, "method")
  expect_match(error_message, "drift")
})


test_that("check_colnames accepts an integer column where the spec asks for numeric", {
  df <- data.frame(target_sample_size_control = 50L, target_sample_size_treatment = 50L)

  expect_silent(
    check_colnames(
      df,
      c("target_sample_size_control", "target_sample_size_treatment"),
      types = c(target_sample_size_control = "numeric", target_sample_size_treatment = "numeric")
    )
  )
})


test_that("check_colnames accepts a list column where the spec asks for list", {
  df <- data.frame(method = "conjugate", stringsAsFactors = FALSE)
  df$parameters <- list(list(alpha = 1))

  expect_silent(
    check_colnames(
      df,
      c("method", "parameters"),
      types = c(method = "character", parameters = "list")
    )
  )
})


test_that("check_colnames ignores columns the type spec does not mention", {
  df <- data.frame(
    method = "conjugate",
    drift = 0.5,
    stringsAsFactors = FALSE
  )

  expect_silent(
    check_colnames(df, c("method", "drift"), types = c(method = "character"))
  )
})


test_that("check_colnames rejects a type spec naming a type it cannot check", {
  df <- data.frame(drift = 0.5)

  expect_error(
    check_colnames(df, "drift", types = c(drift = "duble")),
    "duble"
  )
})


test_that("check_colnames rejects a type spec naming a column it does not expect", {
  df <- data.frame(drift = 0.5)

  expect_error(
    check_colnames(df, "drift", types = c(dirft = "numeric")),
    "dirft"
  )
})


test_that("check_colnames still reports missing and unexpected columns", {
  df <- data.frame(method = "conjugate", drift = 0.5, stringsAsFactors = FALSE)

  expect_error(check_colnames(df, c("method", "drift", "case_study")), "missing")
  expect_error(check_colnames(df, "method"), "unexpected")
})


test_that("default_coltypes returns only the spec entries for the requested columns", {
  types <- default_coltypes(c("method", "drift"))

  expect_setequal(names(types), c("method", "drift"))
  expect_identical(types[["method"]], "character")
  expect_identical(types[["drift"]], "numeric")
})


test_that("default_coltypes ignores requested columns the spec does not describe", {
  types <- default_coltypes(c("method", "warning"))

  expect_setequal(names(types), "method")
})


test_that("check_colnames type checks against the shared spec by default", {
  df <- data.frame(method = 1, drift = 0.5)

  expect_error(check_colnames(df, c("method", "drift")), "method")
})


test_that("check_colnames skips type checking when given an empty spec", {
  df <- data.frame(method = 1, drift = 0.5)

  expect_silent(check_colnames(df, c("method", "drift"), types = character(0)))
})


test_that("check_colnames accepts an all-NA column whatever type the spec declares", {
  df <- data.frame(method = "conjugate", source_denominator = NA)

  expect_silent(
    check_colnames(
      df,
      c("method", "source_denominator"),
      types = c(method = "character", source_denominator = "numeric")
    )
  )
})


test_that("check_colnames still type checks a column holding some non-NA values", {
  df <- data.frame(method = "conjugate", source_denominator = c(TRUE, NA))

  expect_error(
    check_colnames(
      df,
      c("method", "source_denominator"),
      types = c(method = "character", source_denominator = "numeric")
    ),
    "source_denominator"
  )
})


test_that("check_colnames accepts any of the alternatives in a type spec", {
  serialised <- data.frame(parameters = "{'a': 1}", stringsAsFactors = FALSE)
  in_memory <- data.frame(row = 1)
  in_memory$parameters <- list(list(a = 1))
  in_memory <- in_memory["parameters"]

  expect_silent(
    check_colnames(serialised, "parameters", types = c(parameters = "list|character"))
  )
  expect_silent(
    check_colnames(in_memory, "parameters", types = c(parameters = "list|character"))
  )
})


test_that("check_colnames rejects a column matching none of the alternatives", {
  df <- data.frame(parameters = 1)

  expect_error(
    check_colnames(df, "parameters", types = c(parameters = "list|character")),
    "expected list\\|character, found numeric"
  )
})


test_that("check_colnames rejects an unknown type among the alternatives", {
  df <- data.frame(parameters = 1)

  expect_error(
    check_colnames(df, "parameters", types = c(parameters = "list|charcter")),
    "charcter"
  )
})
