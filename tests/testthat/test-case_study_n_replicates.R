## A run can give one case study fewer replicates than the rest: the paper
## simulates aprepitant, whose binomial likelihood is fitted by MCMC, at 1000
## replicates and the normal-likelihood case studies at 10000.

test_that("a case study the map names gets its own replicate count", {
  scenarios_config <- list(
    n_replicates = 10000,
    case_study_n_replicates = list(aprepitant = 1000)
  )

  expect_equal(case_study_n_replicates(scenarios_config, "aprepitant"), 1000)
  expect_equal(case_study_n_replicates(scenarios_config, "botox"), 10000)
})

test_that("without the map every case study gets the run-wide count", {
  expect_equal(
    case_study_n_replicates(list(n_replicates = 500), "aprepitant"),
    500
  )
})

test_that("the map is accepted by the scenarios config schema", {
  withr::with_tempdir({
    yaml::write_yaml(
      list(
        n_replicates = 10000,
        ndrift = 30,
        denominator_change_factor = 1,
        sample_size_factors = c(2, 4),
        case_studies = c("botox", "aprepitant"),
        methods = "separate",
        case_study_n_replicates = list(aprepitant = 1000)
      ),
      "scenarios_config.yml"
    )
    config <- read_config("scenarios_config.yml", scenarios_config_schema)
    expect_equal(case_study_n_replicates(config, "aprepitant"), 1000)
  })
})
