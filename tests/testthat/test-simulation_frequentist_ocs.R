test_that("simulation_frequentist_ocs requires case studies and methods", {
  common_arguments <- list(
    env = "test",
    simulation_config = NULL,
    analysis_config = NULL,
    config_dir = NULL,
    case_studies_config_dir = NULL,
    logging_file_path = NULL
  )

  expect_error(
    do.call(
      simulation_frequentist_ocs,
      c(
        common_arguments,
        list(scenarios_config = list(
          case_studies = character(),
          methods = "separate"
        ))
      )
    ),
    "Select at least one case study and one method"
  )

  expect_error(
    do.call(
      simulation_frequentist_ocs,
      c(
        common_arguments,
        list(scenarios_config = list(
          case_studies = "example",
          methods = character()
        ))
      )
    ),
    "Select at least one case study and one method"
  )
})
