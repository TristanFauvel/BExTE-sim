test_that("deterministic Bayesian OC analysis handles an empty result set", {
  empty_results <- data.frame(case_study = character())

  result <- BExTEsim:::compute_bayesian_ocs(
    results_freq_df = empty_results,
    env = "pipeline_tests",
    config_dir = testthat::test_path("fixtures", "pipeline_tests", "")
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

test_that("the Bayesian OCs run on a cluster when the analysis does", {
  expect_true("parallelization" %in% names(formals(compute_bayesian_ocs)))

  # deparse() wraps long lines, so whitespace is squeezed before matching.
  body_text <- gsub("\\s+", " ", paste(deparse(body(compute_bayesian_ocs)), collapse = " "))
  expect_match(body_text, "analysis_uses_cluster(parallelization, length(jobs))", fixed = TRUE)
  # A foreach body is evaluated outside the package namespace, so the job
  # function travels as a variable rather than being called by name.
  expect_match(body_text, "run_job <- bayesian_ocs_for_combination", fixed = TRUE)
})
