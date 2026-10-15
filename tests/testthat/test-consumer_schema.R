test_that("the shared consumer columns are the ones every results frame carries", {
  expect_setequal(
    required_colnames_consumer,
    c("case_study", "method", "parameters", "target_sample_size_control",
      "target_sample_size_treatment", "source_denominator_change_factor")
  )
})


test_that("the shared consumer columns exclude drift, which Bayesian frames lack", {
  expect_false("drift" %in% required_colnames_consumer)
})


test_that("forest_plot_methods_comparison rejects a frame missing an identity column", {
  df <- data.frame(
    method = "conjugate",
    parameters = "{}",
    target_sample_size_control = 50,
    target_sample_size_treatment = 50,
    source_denominator_change_factor = 1,
    stringsAsFactors = FALSE
  )

  expect_error(
    forest_plot_methods_comparison(df, metrics = list()),
    "forest_plot_methods_comparison.*case_study"
  )
})


test_that("power_vs_tie_tables rejects a frame missing an identity column", {
  df <- data.frame(
    case_study = "botox",
    parameters = "{}",
    target_sample_size_control = 50,
    target_sample_size_treatment = 50,
    source_denominator_change_factor = 1,
    stringsAsFactors = FALSE
  )

  expect_error(
    power_vs_tie_tables(df, metrics = list()),
    "power_vs_tie_tables.*method"
  )
})


test_that("bayesian_ocs_plots rejects an identity column of the wrong type", {
  df <- data.frame(
    case_study = 1,
    method = "conjugate",
    parameters = "{}",
    target_sample_size_control = 50,
    target_sample_size_treatment = 50,
    source_denominator_change_factor = 1,
    stringsAsFactors = FALSE
  )

  expect_error(
    bayesian_ocs_plots(df, metrics = list()),
    "bayesian_ocs_plots.*case_study"
  )
})
