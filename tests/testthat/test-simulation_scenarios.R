
# Test for the important_drift_values function
test_that("important_drift_values computes the correct important drift values", {
  source_treatment_effect <- 1
  case_study_config <- list(summary_measure_likelihood = "normal")

  expected_important_drift_values <- c(-1, -0.5, 0)

  expect_equal(
    important_drift_values(source_treatment_effect, case_study_config),
    expected_important_drift_values
  )
})

# Test for the compute_control_drift_range function
test_that("compute_control_drift_range computes the correct control drift range", {
  source_treatment_effect <- 1

  expected_control_drift_range <- c(0)

  expect_equal(
    compute_control_drift_range(source_treatment_effect),
    expected_control_drift_range
  )
})

# Test for the compute_source_denominator_range function
test_that("compute_source_denominator_range computes the correct source denominator range", {
  source_data <- list(
    endpoint = "binary",
    control_rate = 0.2,
    treatment_effect_estimate = 0.5
  )
  simulation_config <- list(denominator_change_factor = c(0.8, 1, 1.2))
  case_study_config <- list(summary_measure_likelihood = "normal")

  expected_source_denominator_range <- list(
    source_denominator = c(0.20, 0.25, 0.30),
    source_denominator_change_factor = c(0.8, 1, 1.2)
  )
  expect_equal(
    compute_source_denominator_range(source_data, simulation_config, case_study_config),
    expected_source_denominator_range
  )
})


# Building the scenario grid with data.frame(do.call(rbind, list_of_lists))
# yields a matrix of mode list, so every column arrives as a list of length-1
# values rather than a vector of scalars. Nothing downstream wants that, and
# the results-frame schema rejects it outright.

test_that("unwrap_scalar_list_columns unwraps columns whose cells are scalars", {
  df <- data.frame(id = 1:3)
  df$x <- list(1, 2, 3)
  df$y <- list("a", "b", "c")

  out <- unwrap_scalar_list_columns(df)

  expect_false(any(vapply(out, is.list, logical(1))))
  expect_identical(out$x, c(1, 2, 3))
  expect_identical(out$y, c("a", "b", "c"))
  expect_identical(out$id, 1:3)
  expect_identical(nrow(out), nrow(df))
})

test_that("unwrap_scalar_list_columns leaves everything else alone", {
  # A column that genuinely holds several values per row must survive.
  df <- data.frame(id = 1:2)
  df$many <- list(c(1, 2), c(3, 4, 5))

  out <- unwrap_scalar_list_columns(df)

  expect_true(is.list(out$many))
  expect_identical(out$many, df$many)
  expect_identical(out$id, 1:2)

  # So must a column deliberately wrapped in I(), which is how the grid keeps
  # each row's method parameters together.
  wrapped <- data.frame(id = 1:2)
  wrapped$parameters <- I(list("a = 1", "a = 2"))
  expect_true(inherits(unwrap_scalar_list_columns(wrapped)$parameters, "AsIs"))

  # Atomic columns are untouched.
  plain <- data.frame(a = c(1.5, 2.5), b = c("x", "y"), stringsAsFactors = FALSE)
  expect_identical(unwrap_scalar_list_columns(plain), plain)
})

test_that("simulation_scenarios returns a grid of scalars, not list columns", {
  config_dir <- testthat::test_path("fixtures", "aprepitant_mcmc_config_light", "")
  skip_if(config_dir == "/", "packaged configuration not available")

  scenarios_config <- yaml::read_yaml(paste0(config_dir, "scenarios_config.yml"))
  cases <- simulation_scenarios(
    config_dir = config_dir,
    scenarios_config = scenarios_config,
    case_studies_config_dir = paste0(system.file("conf/case_studies", package = "BExTE"), "/")
  )

  # `parameters` holds each row's method parameters and is deliberately
  # wrapped in I(); every other column describes one scenario value per row.
  list_columns <- names(cases)[vapply(cases, is.list, logical(1))]
  expect_identical(list_columns, "parameters")

  expect_true(is.numeric(cases$target_sample_size_control))
  expect_true(is.numeric(cases$target_sample_size_treatment))
  expect_true(is.numeric(cases$drift))
  expect_true(is.character(cases$case_study))
  expect_true(is.character(cases$null_space))
  expect_true(is.numeric(cases$theta_0))
})
