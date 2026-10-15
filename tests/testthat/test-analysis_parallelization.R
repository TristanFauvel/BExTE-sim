## The post-processing analysis used to run single-threaded whatever the
## environment asked for: simulation_analysis() called
## frequentist_power_at_equivalent_tie() with parallelization = FALSE
## hardcoded, and that function's workers loaded BExTE with a bare
## library(BExTEsim), which fails when the package is only load_all()ed from
## source (main.R and BExTE-app both do exactly that).

test_that("analysis_runs_in_parallel reads both forms of the config entry", {
  expect_false(analysis_runs_in_parallel(NULL))
  expect_false(analysis_runs_in_parallel(FALSE))
  expect_true(analysis_runs_in_parallel(TRUE))

  # The per-method form: the analysis loops over every row at once, so any
  # method running in parallel means the analysis does too.
  expect_true(analysis_runs_in_parallel(list("RMP", "NPP")))
  expect_true(analysis_runs_in_parallel("RMP"))
})

test_that("analysis_runs_in_parallel rejects a malformed config entry", {
  expect_error(analysis_runs_in_parallel(c(TRUE, FALSE)), "single TRUE or FALSE")
  expect_error(analysis_runs_in_parallel(NA), "single TRUE or FALSE")
  expect_error(analysis_runs_in_parallel(list()), "single TRUE or FALSE")
  expect_error(analysis_runs_in_parallel(list(1, 2)), "single TRUE or FALSE")
})

test_that("simulation_analysis takes parallelization from the scenarios config", {
  expect_true("parallelization" %in% names(formals(simulation_analysis)))
  body_text <- paste(deparse(body(simulation_analysis)), collapse = " ")
  # The equivalent-TIE step must be handed the resolved setting, not FALSE.
  expect_match(body_text, "parallelization = run_in_parallel", fixed = TRUE)
})

test_that("workers can call BExTE functions when it is only loaded from source", {
  skip_on_cran()
  skip_if_not_installed("parallel")

  cl <- parallel::makeCluster(1L)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  load_bexte_in_workers(cl, packages = c("dplyr"))

  # get_parallel_worker_count() is internal to BExTE, so it only resolves in
  # the worker if the package itself really loaded there.
  worker_sees_bexte <- parallel::clusterEvalQ(cl, {
    is.function(BExTEsim:::get_parallel_worker_count) && is.function(dplyr::bind_rows)
  })
  expect_true(all(unlist(worker_sees_bexte)))
})

test_that("a small run skips the cluster, which costs more than it saves", {
  # ~12s to stand the cluster up against ~0.07s a row: below the threshold
  # the sequential path finishes first.
  expect_false(analysis_uses_cluster(TRUE, n_rows = 108))
  expect_false(analysis_uses_cluster(TRUE, n_rows = ANALYSIS_PARALLEL_MIN_ROWS - 1L))
  expect_true(analysis_uses_cluster(TRUE, n_rows = ANALYSIS_PARALLEL_MIN_ROWS))
  expect_true(analysis_uses_cluster(TRUE, n_rows = 10000))
})

test_that("the cluster is never used when parallelization is off", {
  expect_false(analysis_uses_cluster(FALSE, n_rows = 1e6))
  expect_false(analysis_uses_cluster(NULL, n_rows = 1e6))
  expect_false(analysis_uses_cluster(NA, n_rows = 1e6))
})

test_that("target data is validated once per row, not once per sampled alpha", {
  # compute_power_with_tie_ci() evaluates power at 1000 sampled alphas per
  # row; asserting inside compute_freq_power() made argument checking ~87% of
  # the analysis.
  hot <- paste(deparse(body(BExTEsim:::compute_freq_power)), collapse = " ")
  expect_false(grepl("assert_target_data_numbers", hot, fixed = TRUE))
  expect_false(grepl("assert_number", hot, fixed = TRUE))

  caller <- paste(deparse(body(BExTEsim:::compute_power_with_tie_ci)), collapse = " ")
  expect_true(grepl("assert_target_data_numbers", caller, fixed = TRUE))
})

test_that("assert_target_data_numbers still rejects malformed target data", {
  good <- list(treatment_effect = 0.3, standard_deviation = 1.1, sample_size_control = 50, sample_size_treatment = 50)
  expect_silent(assert_target_data_numbers(good))

  expect_error(assert_target_data_numbers(utils::modifyList(good, list(treatment_effect = "a"))))
  expect_error(assert_target_data_numbers(utils::modifyList(good, list(standard_deviation = c(1, 2)))))
  expect_error(assert_target_data_numbers(list(treatment_effect = 0.3, standard_deviation = 1.1)))
})

test_that("the nominal-TIE step is handed the resolved setting, like the equivalent-TIE one", {
  # It was the only analysis step with no parallel branch at all. On the
  # paper's 18,480-row environment that cost 26 hours single-threaded, while
  # the equivalent-TIE step over the same rows took 30 minutes on a cluster.
  expect_true("parallelization" %in% names(formals(frequentist_power_at_nominal_tie)))

  body_text <- paste(deparse(body(simulation_analysis)), collapse = " ")
  hits <- gregexpr("parallelization = run_in_parallel", body_text, fixed = TRUE)[[1]]
  # Equivalent TIE, nominal TIE and the Bayesian OCs.
  expect_equal(sum(hits > 0), 3L)

  # The cluster is worth standing up for the distinct designs still to
  # compute, which are the work, not for the number of rows, which is 56 times
  # larger. deparse() wraps long lines, so whitespace is squeezed first.
  nominal <- gsub("\\s+", " ", paste(deparse(body(frequentist_power_at_nominal_tie)), collapse = " "))
  expect_true(grepl("analysis_uses_cluster(parallelization, length(to_compute))",
                    nominal, fixed = TRUE))
})

test_that("both branches of the nominal-TIE step compute a row the same way", {
  # The two branches cannot be allowed to drift apart, and a worker cannot
  # see a copy of the computation that lives inside the loop body, so both
  # call the same helper rather than holding a copy of it.
  nominal <- paste(deparse(body(frequentist_power_at_nominal_tie)), collapse = " ")

  # Written once, as a closure, and called from each branch: the definition
  # plus the two call sites.
  calls <- gregexpr("compute_design", nominal, fixed = TRUE)[[1]]
  expect_gte(sum(calls > 0), 3L)

  hits <- gregexpr("nominal_tie_power_row", nominal, fixed = TRUE)[[1]]
  expect_equal(sum(hits > 0), 1L)

  expect_false(grepl("compute_freq_power(", nominal, fixed = TRUE))
  expect_false(grepl("compute_freq_power_pooling(", nominal, fixed = TRUE))
})

test_that("nominal_tie_power_row returns the six baseline columns and the reference's calibration", {
  mock_load_data <- function(results_row, type, reload_data_objects = FALSE) {
    list(type = type, treatment_effect = 0.5, standard_deviation = 1,
         sample_size_control = 30,
         sample_size_treatment = 30)
  }

  with_mocked_bindings(
    {
      row_power <- nominal_tie_power_row(
        row = data.frame(case_study = "example", theta_0 = 0, null_space = "left"),
        nominal_tie = 0.025,
        frequentist_test = "t-test",
        simulation_config = list(),
        n_replicates = 100
      )

      expect_equal(names(row_power), c(
        "nominal_frequentist_power_separate",
        "nominal_frequentist_power_separate_lower",
        "nominal_frequentist_power_separate_upper",
        "nominal_frequentist_power_pooling",
        "nominal_frequentist_power_pooling_lower",
        "nominal_frequentist_power_pooling_upper",
        "frequentist_reference_calibration"
      ))
      expect_equal(unlist(row_power[1:6], use.names = FALSE),
                   c(0.8, 0.75, 0.85, 0.9, 0.87, 0.93))
      # No null scenario was given, so the test stays at its nominal level.
      expect_equal(row_power$frequentist_reference_calibration, "none")

      # dplyr::bind_rows() is how the parallel branch turns these into
      # columns, and a named bound would carry its name into the column.
      expect_null(names(row_power$nominal_frequentist_power_separate_lower))
    },
    load_data = mock_load_data,
    compute_freq_power = function(...) list(power = 0.8, conf_int_power = c(lower = 0.75, upper = 0.85)),
    compute_freq_power_pooling = function(...) list(power = 0.9, conf_int_power = c(lower = 0.87, upper = 0.93)),
    .package = "BExTEsim"
  )
})
