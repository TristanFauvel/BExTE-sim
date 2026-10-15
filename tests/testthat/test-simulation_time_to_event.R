# The target trial follows the modified fixed-calendar design: patients enter
# uniformly over an accrual window, each may be followed for at most L, and the
# database is closed F after the end of recruitment.

accrual_period <- 3.41
final_follow_up <- 48 / 52
max_follow_up_time <- 96 / 52
weibull_shape <- 0.72
relapse_free_probability <- 0.442
source_control_rate <- 0.534
source_effect <- -0.393

simulate_arm <- function(n_subjects,
                         n_replicates,
                         parameter,
                         event_time_distribution = "exponential",
                         dropout_rate = 0,
                         accrual = accrual_period,
                         final = final_follow_up) {
  BExTE:::sample_time_to_event_arm(
    n_subjects = n_subjects,
    n_replicates = n_replicates,
    parameter = parameter,
    event_time_distribution = event_time_distribution,
    weibull_shape = weibull_shape,
    accrual_period = accrual,
    final_follow_up = final,
    max_follow_up_time = max_follow_up_time,
    dropout_rate = dropout_rate
  )
}


test_that("the time-to-event grid axes are only simulated for that endpoint", {
  scenarios <- list(
    control_drift_range = list(-0.405, 0.405),
    dropout_probability = list(0, 0.05),
    event_time_distribution = list("exponential", "weibull")
  )

  time_to_event <- list(endpoint = "time_to_event")
  continuous <- list(endpoint = "continuous")

  expect_equal(
    compute_control_drift_range(0, scenarios, time_to_event),
    c(-0.405, 0, 0.405)
  )
  expect_equal(compute_control_drift_range(0, scenarios, continuous), 0)
  # A configuration that does not mention control drift keeps a single scenario.
  expect_equal(compute_control_drift_range(0, list(), time_to_event), 0)

  ranges <- compute_time_to_event_ranges(scenarios, time_to_event)
  expect_equal(ranges$dropout_probability, c(0, 0.05))
  expect_equal(ranges$event_time_distribution, c("exponential", "weibull"))

  defaults <- compute_time_to_event_ranges(scenarios, continuous)
  expect_equal(defaults$dropout_probability, 0)
  expect_equal(defaults$event_time_distribution, "exponential")
})


test_that("the primary design is no heterogeneity, no dropout and exponential times", {
  expect_true(is_primary_time_to_event_design(0, 0, "exponential"))
  expect_false(is_primary_time_to_event_design(0.405, 0, "exponential"))
  expect_false(is_primary_time_to_event_design(0, 0.05, "exponential"))
  expect_false(is_primary_time_to_event_design(0, 0, "weibull"))
})


test_that("the sensitivity reference is read only for the time-to-event endpoint", {
  scenarios <- list(sensitivity_reference = list(
    sample_size_factor = 6, denominator_change_factor = 1
  ))
  time_to_event <- list(endpoint = "time_to_event")

  reference <- time_to_event_sensitivity_reference(scenarios, time_to_event)
  expect_equal(reference$sample_size_factor, 6)
  expect_equal(reference$denominator_change_factor, 1)

  expect_null(time_to_event_sensitivity_reference(scenarios, list(endpoint = "continuous")))
  # Without the key every design crosses the whole grid, as it did before the
  # key existed.
  expect_null(time_to_event_sensitivity_reference(list(), time_to_event))

  expect_error(
    time_to_event_sensitivity_reference(
      list(sensitivity_reference = list(sample_size_factor = 6)), time_to_event
    ),
    "denominator_change_factor"
  )
})


test_that("the sensitivity designs are simulated only at the reference point", {
  config_dir <- testthat::test_path("fixtures", "combined_teriflunomide", "")
  case_studies_config_dir <- paste0(system.file("conf/case_studies", package = "BExTE"), "/")
  skip_if(!nzchar(config_dir) || !file.exists(paste0(config_dir, "scenarios_config.yml")))

  scenarios_config <- read_config(
    paste0(config_dir, "scenarios_config.yml"), scenarios_config_schema
  )
  scenarios_config$case_studies <- "teriflunomide"
  scenarios_config$methods <- "separate"

  grid <- unwrap_scalar_list_columns(simulation_scenarios(
    config_dir = config_dir,
    scenarios_config = scenarios_config,
    case_studies_config_dir = case_studies_config_dir
  ))

  primary <- mapply(
    is_primary_time_to_event_design,
    grid$control_drift, grid$dropout_probability, grid$event_time_distribution
  )

  reference <- scenarios_config$sensitivity_reference
  expected_size <- floor((1483 / reference$sample_size_factor) / 2)

  # The whole grid runs under the primary design.
  expect_setequal(grid$source_denominator_change_factor[primary],
                  unlist(scenarios_config$denominator_change_factor))
  expect_length(unique(grid$target_sample_size_control[primary]),
                length(unlist(scenarios_config$sample_size_factors)))

  # The sensitivity designs run at one trial size and one source denominator,
  # but still across the whole drift range.
  expect_equal(unique(grid$target_sample_size_control[!primary]), expected_size)
  expect_equal(unique(grid$source_denominator_change_factor[!primary]),
               reference$denominator_change_factor)
  expect_setequal(grid$drift[!primary], grid$drift[primary])
})


test_that("a sensitivity reference outside the sample size factors is refused", {
  config_dir <- testthat::test_path("fixtures", "combined_teriflunomide", "")
  case_studies_config_dir <- paste0(system.file("conf/case_studies", package = "BExTE"), "/")
  skip_if(!nzchar(config_dir) || !file.exists(paste0(config_dir, "scenarios_config.yml")))

  scenarios_config <- read_config(
    paste0(config_dir, "scenarios_config.yml"), scenarios_config_schema
  )
  scenarios_config$case_studies <- "teriflunomide"
  scenarios_config$methods <- "separate"
  scenarios_config$sensitivity_reference$sample_size_factor <- 99

  expect_error(
    simulation_scenarios(
      config_dir = config_dir,
      scenarios_config = scenarios_config,
      case_studies_config_dir = case_studies_config_dir
    ),
    "never be simulated"
  )
})


test_that("one axis at a time keeps the primary design and five sensitivity designs", {
  config_dir <- testthat::test_path("fixtures", "combined_teriflunomide", "")
  case_studies_config_dir <- paste0(system.file("conf/case_studies", package = "BExTE"), "/")
  skip_if(!nzchar(config_dir) || !file.exists(paste0(config_dir, "scenarios_config.yml")))

  scenarios_config <- read_config(
    paste0(config_dir, "scenarios_config.yml"), scenarios_config_schema
  )
  scenarios_config$case_studies <- "teriflunomide"
  scenarios_config$methods <- "separate"

  designs <- function(config) {
    grid <- unwrap_scalar_list_columns(simulation_scenarios(
      config_dir = config_dir,
      scenarios_config = config,
      case_studies_config_dir = case_studies_config_dir
    ))
    unique(grid[, c("control_drift", "dropout_probability", "event_time_distribution")])
  }

  # Crossed, every combination of the three axes: 3 x 3 x 2.
  expect_equal(nrow(designs(scenarios_config)), 18)

  scenarios_config$sensitivity_one_at_a_time <- TRUE
  one_at_a_time <- designs(scenarios_config)
  expect_equal(nrow(one_at_a_time), 6)
  expect_true(all(time_to_event_axes_off_primary(
    one_at_a_time$control_drift,
    one_at_a_time$dropout_probability,
    one_at_a_time$event_time_distribution
  ) <= 1))
})


teriflunomide_design <- function() {
  config <- yaml::read_yaml(
    system.file("conf/case_studies/teriflunomide.yml", package = "BExTE")
  )
  config$target
}


test_that("the delay is a sensitivity axis, simulated one at a time", {
  expect_false(is_primary_time_to_event_design(0, 0, "exponential", 0.230769))
  expect_equal(time_to_event_axes_off_primary(0, 0, "exponential", 0.230769), 1L)

  config_dir <- testthat::test_path("fixtures", "combined", "")
  case_studies_config_dir <- paste0(system.file("conf/case_studies", package = "BExTE"), "/")
  scenarios_config <- read_config(
    paste0(config_dir, "scenarios_config.yml"), scenarios_config_schema
  )
  scenarios_config$case_studies <- c("teriflunomide", "botox")
  scenarios_config$methods <- "separate"
  grid <- unwrap_scalar_list_columns(simulation_scenarios(
    config_dir = config_dir,
    scenarios_config = scenarios_config,
    case_studies_config_dir = case_studies_config_dir
  ))

  delayed <- grid[grid$treatment_delay > 0, ]
  expect_setequal(unique(delayed$treatment_delay), c(0.230769, 0.461538))
  expect_equal(unique(delayed$case_study), "teriflunomide")
  expect_equal(unique(delayed$target_sample_size_control), 123)
  # Nothing else moves off the primary design alongside the delay.
  expect_true(all(delayed$control_drift == 0 & delayed$dropout_probability == 0 &
                    delayed$event_time_distribution == "exponential"))
  # Seven sensitivity designs besides the primary one.
  teriflunomide <- grid[grid$case_study == "teriflunomide" &
                          grid$target_sample_size_control == 123 &
                          grid$source_denominator_change_factor == 1, ]
  expect_equal(nrow(unique(teriflunomide[, c("control_drift", "dropout_probability",
                                             "event_time_distribution", "treatment_delay")])), 8)
})
