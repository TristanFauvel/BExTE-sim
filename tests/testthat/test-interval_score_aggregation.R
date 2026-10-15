# The interval score answers a reviewer's objection that an interval's width
# means nothing without its coverage, so it is only useful if it actually
# reaches the results frame alongside `precision` and `coverage`. These tests
# pin the aggregation step: that the scenario-level value is the mean of the
# per-replicate scores, that it carries a Monte Carlo interval like every other
# aggregated metric, and that it is charged against the same truth `coverage`
# uses.

canned_results_model <- function(results) {
  R6::R6Class(
    "IntervalScoreCannedResultsTestModel",
    inherit = Model,
    public = list(
      simulation_for_given_treatment_effect = function(...) results
    )
  )$new()
}

three_replicates <- list(
  fit_success = rep("Success", 3),
  test_decisions = c(1, 1, 0),
  posterior_means = c(1, 2, 3),
  posterior_medians = c(1, 2, 3),
  # Rows: one interval that misses low, one that covers, one that misses high.
  credible_intervals = matrix(
    c(-3, -1, 1, 3, 5, 7),
    ncol = 2, byrow = TRUE
  ),
  posterior_parameters = NULL,
  ess_moments = c(10, 20, 30),
  ess_precisions = c(5, 10, 15),
  ess_elir = c(1, 2, 3),
  rhat = c(1.00, 1.01, 1.02),
  mcmc_ess = c(500, 600, 700),
  n_divergences = c(0, 1, 2)
)

run_ocs <- function(results, treatment_effect = 2) {
  canned_results_model(results)$estimate_frequentist_operating_characteristics(
    theta_0 = 0,
    target_data = list(treatment_effect = treatment_effect),
    n_replicates = length(results$fit_success),
    critical_value = 0.975,
    confidence_level = 0.95,
    null_space = "left",
    n_samples_quantiles_estimation = 100,
    case_study = "unit_test",
    method = "separate"
  )
}


test_that("the printable summary reports the score next to the precision", {
  result <- run_ocs(three_replicates)
  table <- format_simulation_output_table(result)

  expect_true("Interval Score" %in% table$Metric)
  expect_equal(
    as.numeric(table$Value[table$Metric == "Interval Score"]),
    result$interval_score,
    tolerance = 1e-3
  )
})
