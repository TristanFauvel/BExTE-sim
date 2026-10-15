## An environment's scenarios_config.yml selects which analysis steps run
## after the simulation. Before that key existed the four steps were
## hardcoded, so a run that needed only the power baselines still paid for the
## sweet spot and the Bayesian operating characteristics.

test_that("a config without the key runs every step", {
  expect_equal(resolve_analysis_steps(list(), TRUE), ANALYSIS_STEPS)
  expect_equal(resolve_analysis_steps(list(analysis_steps = NULL), TRUE), ANALYSIS_STEPS)
})

test_that("a config naming a subset runs only that subset, in the order given", {
  config <- list(analysis_steps = list(
    "frequentist_power_at_equivalent_tie",
    "frequentist_power_at_nominal_tie"
  ))

  expect_equal(
    resolve_analysis_steps(config, TRUE),
    c("frequentist_power_at_equivalent_tie", "frequentist_power_at_nominal_tie")
  )

  # yaml::read_yaml gives a list for a block sequence and a character vector
  # for a single scalar, and both have to work.
  expect_equal(resolve_analysis_steps(list(analysis_steps = "sweet_spot"), TRUE),
               "sweet_spot")
})

test_that("a misspelled step is rejected where the config is read", {
  # Silently skipping it would surface much later, as a figure with nothing
  # to plot, and the run would have to be repeated.
  expect_error(
    resolve_analysis_steps(list(analysis_steps = list("sweet_spots")), TRUE),
    "steps that do not exist: sweet_spots"
  )
  expect_error(
    resolve_analysis_steps(list(analysis_steps = list("bayesian_ocs", "typo")), TRUE),
    "typo"
  )
})

test_that("the power steps are dropped for a case study outside the package", {
  # They call load_data() without threading case_studies_config_dir through,
  # so they only resolve package-shipped case studies.
  expect_equal(resolve_analysis_steps(list(), FALSE), c("sweet_spot", "bayesian_ocs"))

  config <- list(analysis_steps = list("frequentist_power_at_nominal_tie"))
  expect_equal(resolve_analysis_steps(config, FALSE), character(0))
})

test_that("the schema accepts the key and rejects a malformed one", {
  base <- list(
    n_replicates = 10, ndrift = 3, denominator_change_factor = list(1),
    sample_size_factors = list(2), case_studies = list("botox"),
    methods = list("RMP")
  )

  expect_silent(validate_config(base, scenarios_config_schema, "test"))
  expect_silent(validate_config(
    utils::modifyList(base, list(analysis_steps = list("sweet_spot"))),
    scenarios_config_schema, "test"
  ))
  expect_error(
    validate_config(
      utils::modifyList(base, list(analysis_steps = list(1, 2))),
      scenarios_config_schema, "test"
    ),
    "analysis_steps should be"
  )
})
