test_that("append_parameters_str separates a trailing value from the first parameter", {
  expect_identical(
    append_parameters_str("teriflunomide_rmp_mse_vs_drift_target_sample_size_per_arm=123", "w=0.5"),
    "teriflunomide_rmp_mse_vs_drift_target_sample_size_per_arm=123_w=0.5"
  )
})

test_that("append_parameters_str adds no trailing underscore without a varied parameter", {
  expect_identical(append_parameters_str("stem", ""), "stem")
  expect_identical(append_parameters_str("stem", character(0)), "stem")
})

test_that("filenames built from convert_params_to_str keep their parameters separate", {
  # One parameter varied on the grid, one fixed: only the varied one is named.
  method <- list(
    power_parameter = list(range = c(0, 0.5, 1), parameter_label = "gamma"),
    initial_prior = list(range = list("noninformative"), parameter_label = "pi_0")
  )
  parameters <- list(power_parameter = 0.5, initial_prior = "noninformative")

  expect_identical(
    append_parameters_str("sample_size_per_arm=123", convert_params_to_str(method, parameters)),
    "sample_size_per_arm=123_gamma=0.5"
  )

  # A method with nothing varied, such as the EBPP, contributes nothing.
  fixed_method <- list(
    initial_prior = list(range = list("noninformative"), parameter_label = "pi_0")
  )
  expect_identical(
    append_parameters_str(
      "sample_size_per_arm=123",
      convert_params_to_str(fixed_method, list(initial_prior = "noninformative"))
    ),
    "sample_size_per_arm=123"
  )
})
