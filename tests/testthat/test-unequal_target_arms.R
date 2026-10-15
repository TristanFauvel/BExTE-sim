## The target study may have arms of different sizes. A summary-based method
## sees them only through the standard error and the equivalent sample size
## per arm, 2 / (1 / n_C + 1 / n_T); an arm-level one reads each arm.

case_study_config <- function(case_study, sampling_approximation = NULL) {
  config <- yaml::yaml.load_file(
    system.file("conf/case_studies", paste0(case_study, ".yml"), package = "BExTE")
  )
  if (!is.null(sampling_approximation)) {
    config$sampling_approximation <- sampling_approximation
  }
  config
}

unequal_target <- function(case_study, n_control = 40, n_treatment = 120,
                           sampling_approximation = NULL) {
  config <- case_study_config(case_study, sampling_approximation)
  TargetDataFactory$new()$create(
    source_data = ObservedSourceData$new(config),
    case_study_config = config,
    target_sample_size_control = n_control,
    target_sample_size_treatment = n_treatment,
    treatment_drift = 0,
    summary_measure_likelihood = config$summary_measure_likelihood
  )
}

test_that("the total target sample size is split by the allocation ratio", {
  expect_equal(split_target_sample_size(93, 1), list(control = 46, treatment = 46))
  expect_equal(split_target_sample_size(300, 2), list(control = 100, treatment = 200))
  expect_equal(split_target_sample_size(301, 0.5), list(control = 201, treatment = 100))
  expect_error(split_target_sample_size(5, 10), "fewer than 2 patients")
  expect_error(split_target_sample_size(100, 0), "positive")
})

test_that("the allocation ratio is read per case study, then for the run", {
  expect_equal(target_allocation_ratio(list(), "botox"), 1)
  config <- list(target_allocation_ratio = 2,
                 case_study_target_allocation_ratio = list(belimumab = list(1.5)))
  expect_equal(target_allocation_ratio(config, "botox"), 2)
  expect_equal(target_allocation_ratio(config, "belimumab"), 1.5)
  expect_error(target_allocation_ratio(list(target_allocation_ratio = -1), "botox"),
               "positive")
})

test_that("the analytical binomial power tests the two arms at their own sizes", {
  target <- unequal_target("aprepitant")
  expected <- pwr::pwr.2p2n.test(
    h = pwr::ES.h(target$treatment_rate, target$control_rate),
    n1 = 120, n2 = 40, sig.level = 0.025, alternative = "greater"
  )$power
  expect_equal(
    analytical_power(0.025, target, "z-test", theta_0 = 0, alternative = "greater"),
    expected
  )
})

test_that("the pooled binomial power sums over each arm's own responder counts", {
  target <- unequal_target("aprepitant", n_control = 20, n_treatment = 60)
  source <- ObservedSourceData$new(case_study_config("aprepitant"))
  power <- pooled_binomial_power(0.025, target, source, "greater")
  expect_true(power > 0 && power < 1)

  # Brute force over the target responder counts.
  source_treatment <- counts_from_rate(source$treatment_rate, source$sample_size_treatment)
  source_control <- counts_from_rate(source$control_rate, source$sample_size_control)
  grid <- expand.grid(treatment = 0:60, control = 0:20)
  rate_treatment <- (source_treatment + grid$treatment) / (60 + source$sample_size_treatment)
  rate_control <- (source_control + grid$control) / (20 + source$sample_size_control)
  statistic <- (2 * asin(sqrt(rate_treatment)) - 2 * asin(sqrt(rate_control))) /
    sqrt(1 / (60 + source$sample_size_treatment) + 1 / (20 + source$sample_size_control))
  probability <- stats::dbinom(grid$treatment, 60, target$treatment_rate) *
    stats::dbinom(grid$control, 20, target$control_rate)
  expect_equal(power, sum(probability[statistic > stats::qnorm(0.975)]))
})

test_that("results written with a single per-arm size load with both arm sizes", {
  old <- data.frame(case_study = "botox", target_sample_size_per_arm = 91)
  new <- upgrade_target_size_columns(old)
  expect_equal(new$target_sample_size_control, 91)
  expect_equal(new$target_sample_size_treatment, 91)
  expect_false("target_sample_size_per_arm" %in% names(new))
  expect_identical(upgrade_target_size_columns(new), new)
})

test_that("the scenario grid splits each total by the allocation ratio", {
  config_dir <- testthat::test_path("fixtures", "aprepitant_mcmc_config_light", "")
  skip_if(config_dir == "/", "packaged configuration not available")
  scenarios_config <- yaml::read_yaml(paste0(config_dir, "scenarios_config.yml"))
  case_studies_config_dir <- paste0(system.file("conf/case_studies", package = "BExTE"), "/")

  balanced <- simulation_scenarios(config_dir, scenarios_config, case_studies_config_dir)
  scenarios_config$target_allocation_ratio <- 2
  unequal <- simulation_scenarios(config_dir, scenarios_config, case_studies_config_dir)

  expect_false("target_sample_size_per_arm" %in% names(unequal))
  expect_equal(balanced$target_sample_size_control, balanced$target_sample_size_treatment)
  expect_true(all(abs(unequal$target_sample_size_treatment -
                        2 * unequal$target_sample_size_control) <= 2))
  expect_true(all(unequal$target_sample_size_control +
                    unequal$target_sample_size_treatment >=
                    balanced$target_sample_size_control + balanced$target_sample_size_treatment))
})
