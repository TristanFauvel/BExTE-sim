## load_data() rebuilds a scenario's target data from its results row. It used
## to leave out the time-to-event design axes, so the power baselines of every
## sensitivity design - dropout, a Weibull event time distribution, a delayed
## treatment effect - were simulated on the primary design's trials.

teriflunomide_row <- function(...) {
  data.frame(
    case_study = "teriflunomide",
    source_denominator = NA_real_,
    target_sample_size_control = 123,
    target_sample_size_treatment = 123,
    treatment_drift = 0,
    control_drift = 0,
    target_to_source_std_ratio = 1,
    ...
  )
}

test_that("load_data rebuilds the time-to-event design the row describes", {
  target_data <- load_data(
    teriflunomide_row(dropout_probability = 0.1,
                      event_time_distribution = "weibull",
                      treatment_delay = 0.25),
    type = "target",
    reload_data_objects = TRUE
  )

  expect_equal(target_data$dropout_probability, 0.1)
  expect_equal(target_data$event_time_distribution, "weibull")
  expect_equal(target_data$treatment_delay, 0.25)
})

test_that("load_data falls back to the primary design for older results", {
  for (row in list(
    teriflunomide_row(),
    teriflunomide_row(dropout_probability = NA_real_,
                      event_time_distribution = NA_character_,
                      treatment_delay = NA_real_)
  )) {
    target_data <- load_data(row, type = "target", reload_data_objects = TRUE)

    expect_equal(target_data$dropout_probability, 0)
    expect_equal(target_data$event_time_distribution, "exponential")
    expect_equal(target_data$treatment_delay, 0)
  }
})
