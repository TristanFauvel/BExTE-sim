## Table S5: precision, coverage and interval score per method and scenario.

precision_ecp_rows <- function(width = 0) {
  effects <- c(0, 0.25, 0.5)
  data.frame(
    case_study = "belimumab", target_sample_size_control = 140, target_sample_size_treatment = 140,
    method = "separate", parameters = "[[{'initial_prior': ['noninformative']}]]",
    target_to_source_std_ratio = NA, source_denominator_change_factor = 1,
    source_treatment_effect_estimate = 0.5, target_treatment_effect = effects,
    precision = c(0.40, 0.41, 0.42), conf_int_precision_lower = c(0.40, 0.41, 0.42) - width,
    conf_int_precision_upper = c(0.40, 0.41, 0.42) + width,
    coverage = c(0.95, 0.94, 0.93), conf_int_coverage_lower = c(0.95, 0.94, 0.93) - width,
    conf_int_coverage_upper = c(0.95, 0.94, 0.93) + width,
    interval_score = c(1.10, 1.20, 1.30), conf_int_interval_score_lower = c(1.10, 1.20, 1.30) - width,
    conf_int_interval_score_upper = c(1.10, 1.20, 1.30) + width
  )
}

## The method labels come from the run's methods_dict, a global the plot
## configuration defines; the columns are what is tested here.
local_plain_labels <- function(env = parent.frame()) {
  local_mocked_bindings(format_results_df_parameters = function(df) df$method, .env = env)
}

## export_table() warns that a longtable cannot be scaled to the page, for
## every long table; that is not what is tested here.
write_table <- function(rows) {
  withCallingHandlers(
    table_precision_ecp(rows, getwd()),
    warning = function(w) {
      if (grepl("Longtable cannot be resized", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
}

test_that("the table reports the interval score beside precision and coverage", {
  local_plain_labels()
  withr::with_tempdir({
    path <- write_table(precision_ecp_rows())
    tex <- paste(readLines(paste0(path, ".tex")), collapse = "\n")
    expect_match(tex, "Interval score")
    expect_match(tex, "1.300")
    ## Exact operating characteristics have zero-width intervals, so only the
    ## estimate is shown.
    expect_false(grepl("95\\\\% CI", tex))
    expect_false(grepl("[1.300, 1.300]", tex, fixed = TRUE))
  })
})

test_that("Monte Carlo intervals are kept", {
  local_plain_labels()
  withr::with_tempdir({
    path <- write_table(precision_ecp_rows(width = 0.01))
    tex <- paste(readLines(paste0(path, ".tex")), collapse = "\n")
    expect_match(tex, "1.300 [1.290, 1.310]", fixed = TRUE)
  })
})
