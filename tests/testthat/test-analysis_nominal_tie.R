# The type I error screens compare against the nominal rate the analysis was
# run at, which the caller supplies, rather than a hard-coded 0.025.

local_output_path <- function() {
  path <- tempfile("analysis_nominal_tie_")
  dir.create(path)
  path
}

read_written <- function(output_path, file_name) {
  utils::read.csv(file.path(output_path, file_name))
}

power_loss_frame <- function(conf_int_tie_lower) {
  data.frame(
    method = "conjugate",
    null_space = "left",
    target_treatment_effect = 1,
    theta_0 = 0,
    success_proba = 0.3,
    mcse_success_proba = 0.01,
    conf_int_success_proba_lower = 0.28,
    conf_int_success_proba_upper = 0.32,
    frequentist_power_at_equivalent_tie = 0.6,
    frequentist_power_at_equivalent_tie_lower = 0.58,
    frequentist_power_at_equivalent_tie_upper = 0.62,
    conf_int_tie_lower = conf_int_tie_lower,
    stringsAsFactors = FALSE
  )
}

test_that("an inflated type I error is judged against the nominal rate given", {
  output_path <- local_output_path()
  file_name <- "power_loss_inflated_tie_cases.csv"
  df <- power_loss_frame(conf_int_tie_lower = 0.03)

  analyze_power_loss_inflated_tie(df, output_path, nominal_tie = 0.025)
  expect_equal(nrow(read_written(output_path, file_name)), 1)

  analyze_power_loss_inflated_tie(df, output_path, nominal_tie = 0.05)
  expect_equal(nrow(read_written(output_path, file_name)), 0)
})

test_that("a non-inflated type I error is judged against the nominal rate given", {
  output_path <- local_output_path()
  file_name <- "noninflated_tie_cases.csv"
  df <- data.frame(
    target_treatment_effect = 0,
    theta_0 = 0,
    conf_int_tie_upper = 0.04
  )

  analyze_noninflated_tie(df, output_path, nominal_tie = 0.025)
  expect_equal(nrow(read_written(output_path, file_name)), 0)

  analyze_noninflated_tie(df, output_path, nominal_tie = 0.05)
  expect_equal(nrow(read_written(output_path, file_name)), 1)
})

test_that("the null boundary is found at a non-zero theta_0 despite rounding", {
  output_path <- local_output_path()
  # 0.1 + 0.2 is not exactly 0.3 in floating point, which is how a target
  # effect built as drift + source estimate lands on theta_0.
  df <- data.frame(
    target_treatment_effect = c(0.1 + 0.2, 0.5),
    theta_0 = 0.3,
    conf_int_tie_upper = 0.01
  )

  analyze_noninflated_tie(df, output_path, nominal_tie = 0.025)
  written <- read_written(output_path, "noninflated_tie_cases.csv")

  expect_equal(nrow(written), 1)
  expect_equal(written$target_treatment_effect, 0.3)
})

test_that("the nominal rate must be supplied and be a probability", {
  df <- power_loss_frame(conf_int_tie_lower = 0.03)

  expect_error(analyze_power_loss_inflated_tie(df, tempdir()), "nominal_tie")
  expect_error(
    analyze_noninflated_tie(df, tempdir(), nominal_tie = 2.5),
    "strictly between 0 and 1"
  )
})
