power_gains_frame <- function(...) {
  df <- data.frame(
    method = "conjugate",
    null_space = "left",
    target_treatment_effect = 1,
    theta_0 = 0,
    success_proba = 0.8,
    mcse_success_proba = 0.01,
    conf_int_success_proba_lower = 0.7,
    conf_int_success_proba_upper = 0.9,
    frequentist_power_at_equivalent_tie = 0.6,
    frequentist_power_at_equivalent_tie_lower = 0.5,
    frequentist_power_at_equivalent_tie_upper = 0.7,
    nominal_frequentist_power_separate = 0.5,
    conf_int_tie_lower = 0.03,
    stringsAsFactors = FALSE
  )
  overrides <- list(...)
  for (name in names(overrides)) df[[name]] <- overrides[[name]]
  df
}


test_that("analyze_power_gains rejects a frame missing a column it decides on", {
  df <- power_gains_frame()
  df$nominal_frequentist_power_separate <- NULL

  expect_error(
    analyze_power_gains(df, tempdir()),
    "nominal_frequentist_power_separate"
  )
})


test_that("analyze_power_gains rejects a decision column of the wrong type", {
  df <- power_gains_frame(success_proba = "0.8")

  expect_error(analyze_power_gains(df, tempdir()), "success_proba")
})


test_that("analyze_power_loss rejects a frame missing a column it decides on", {
  df <- power_gains_frame()
  df$frequentist_power_at_equivalent_tie <- NULL

  expect_error(
    analyze_power_loss(df, tempdir()),
    "frequentist_power_at_equivalent_tie"
  )
})


test_that("analyze_power_loss_inflated_tie rejects a frame missing its TIE column", {
  df <- power_gains_frame()
  df$conf_int_tie_lower <- NULL

  expect_error(
    analyze_power_loss_inflated_tie(df, tempdir(), nominal_tie = 0.025),
    "conf_int_tie_lower"
  )
})


test_that("the schema error names the analysis that rejected the frame", {
  df <- power_gains_frame()
  df$theta_0 <- NULL

  expect_error(analyze_power_gains(df, tempdir()), "analyze_power_gains")
})
