

test_that("the shipped simulation and mcmc configs satisfy their schemas", {
  simulation_config <- yaml::read_yaml(
    system.file("conf/simulation_config.yml", package = "BExTEsim")
  )
  expect_silent(
    validate_config(simulation_config, simulation_config_schema, "simulation_config.yml")
  )
})


test_that("every shipped config satisfies its schema", {
  conf_dir <- system.file("conf", package = "BExTEsim")

  scenarios <- list.files(conf_dir, pattern = "^scenarios_config\\.yml$",
                          recursive = TRUE, full.names = TRUE)
  expect_gt(length(scenarios), 0)
  for (file in scenarios) {
    expect_silent(
      validate_config(yaml::read_yaml(file), scenarios_config_schema, file)
    )
  }

  mcmc <- list.files(conf_dir, pattern = "^mcmc_config\\.yml$",
                     recursive = TRUE, full.names = TRUE)
  expect_gt(length(mcmc), 0)
  for (file in mcmc) {
    expect_silent(validate_config(yaml::read_yaml(file), mcmc_config_schema, file))
  }
})


test_that("read_config returns the parsed config when it is valid", {
  path <- system.file("conf/simulation_config.yml", package = "BExTEsim")

  config <- read_config(path, simulation_config_schema)

  expect_identical(config$seed, 42L)
  expect_identical(config$critical_value, 0.975)
})


test_that("read_config rejects an invalid config and names the file", {
  path <- file.path(tempdir(), "bad_simulation_config.yml")
  yaml::write_yaml(list(critical_value = 3), path)

  expect_error(
    read_config(path, simulation_config_schema),
    "bad_simulation_config.yml"
  )
})


test_that("check_decision_threshold accepts a critical value matching the interval level", {
  expect_silent(check_decision_threshold(
    list(critical_value = 0.975, confidence_level = 0.95),
    "simulation_config.yml"
  ))
})


test_that("check_decision_threshold rejects a critical value decoupled from the interval level", {
  expect_error(
    check_decision_threshold(
      list(critical_value = 0.99, confidence_level = 0.95),
      "simulation_config.yml"
    ),
    "simulation_config.yml is not valid: critical_value should be"
  )
})


test_that("read_simulation_config rejects a decoupled decision threshold and names the file", {
  config <- yaml::read_yaml(system.file("conf/simulation_config.yml", package = "BExTEsim"))
  config$critical_value <- 0.99
  path <- file.path(withr::local_tempdir(), "decoupled_simulation_config.yml")
  yaml::write_yaml(config, path)

  expect_error(read_simulation_config(path), "decoupled_simulation_config.yml")
})


test_that("run_simulation_env refuses a decoupled decision threshold before touching the disk", {
  withr::local_dir(withr::local_tempdir())

  expect_error(
    run_simulation_env(
      env = "test_env",
      config_dir = "./",
      case_studies_config_dir = "./",
      simulation_config = list(critical_value = 0.99, confidence_level = 0.95),
      analysis_config = list(),
      frequentist_metrics = list()
    ),
    "critical_value should be"
  )
  expect_false(dir.exists("logs"))
})
