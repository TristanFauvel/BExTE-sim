## The power at equivalent TIE used to simulate the separate analysis's trials
## once per result row, although they depend only on the design: every method
## and parameter combination of a design reseeded and drew the same trials
## again. On the paper's time-to-event case study that was ~2.5s a row, 14,784
## rows. The trials are now simulated once per design, and shared with the
## nominal-TIE step, which reads the same ones.

test_that("the vectorised t-test p-values equal BSDA::tsum.test exactly", {
  set.seed(3)
  estimate <- rnorm(200, 0.1, 0.3)
  standard_deviation <- runif(200, 0.5, 2)
  n <- sample(20:200, 200, replace = TRUE)

  for (alternative in c("greater", "less", "two.sided")) {
    expected <- vapply(seq_along(estimate), function(r) {
      suppressWarnings(BSDA::tsum.test(
        mean.x = estimate[r], mu = 0.05, alternative = alternative,
        s.x = standard_deviation[r], n.x = n[r]
      ))$p.value
    }, numeric(1))

    expect_identical(
      one_sample_t_test_p_values(estimate, 0.05, standard_deviation, n, alternative),
      expected
    )
  }
})

simulated_design_results <- function() {
  data.frame(
    method = rep(c("separate", "pooling"), each = 2),
    parameters = "{}",
    control_drift = 0,
    source_denominator = NA_real_,
    source_denominator_change_factor = 1,
    case_study = "example",
    target_to_source_std_ratio = 1,
    target_sample_size_control = 30,
    target_sample_size_treatment = 30,
    theta_0 = 0,
    null_space = "left",
    sampling_approximation = TRUE,
    summary_measure_likelihood = "normal",
    source_sample_size_treatment = 50,
    source_sample_size_control = 50,
    endpoint = "recurrent_event",
    source_standard_error = 0.2,
    source_treatment_effect_estimate = 0.5,
    equivalent_source_sample_size_per_arm = 50,
    target_treatment_effect = c(0, 0.5, 0, 0.5),
    success_proba = c(0.025, 0.8, 0.04, 0.85),
    mcse_success_proba = c(0.005, 0.02, 0.006, 0.02),
    conf_int_success_proba_lower = c(0.015, 0.76, 0.028, 0.81),
    conf_int_success_proba_upper = c(0.035, 0.84, 0.052, 0.89)
  )
}

counting_load_data <- function(counter) {
  function(results_row, type, reload_data_objects = FALSE) {
    if (type == "source") {
      return(list(
        treatment_effect_estimate = results_row$source_treatment_effect_estimate,
        standard_error = results_row$source_standard_error,
        equivalent_source_sample_size_per_arm = results_row$equivalent_source_sample_size_per_arm
      ))
    }
    list(
      sample_size_control = 30,
      sample_size_treatment = 30,
      treatment_effect = results_row$target_treatment_effect,
      standard_deviation = 1,
      summary_measure_likelihood = "normal",
      endpoint = "recurrent_event",
      generate = function(n_replicates) {
        counter$draws <- counter$draws + 1
        data.frame(
          treatment_effect_estimate = rnorm(n_replicates, results_row$target_treatment_effect, 1 / sqrt(30)),
          standard_deviation = 1,
          sample_size_control = 30,
          sample_size_treatment = 30
        )
      }
    )
  }
}

test_that("the trials are simulated once per design, not once per row", {
  counter <- new.env()
  counter$draws <- 0
  cache <- new.env(parent = emptyenv())

  with_mocked_bindings(
    {
      results <- frequentist_power_at_equivalent_tie(
        results = simulated_design_results(),
        analysis_config = list(frequentist_test = "t-test"),
        simulation_config = list(seed = 7),
        parallelization = FALSE,
        n_replicates = 200,
        trial_cache = cache
      )
    },
    load_data = counting_load_data(counter),
    .package = "BExTEsim"
  )

  # Four rows, two designs (the two treatment effects), two methods each.
  expect_equal(counter$draws, 2)
  expect_length(ls(cache), 2)
  expect_false(anyNA(results$frequentist_power_at_equivalent_tie))
  expect_true(all(results$frequentist_power_at_equivalent_tie_lower <=
                    results$frequentist_power_at_equivalent_tie_upper))
})

test_that("rows sharing a design and a type I error get the same power", {
  counter <- new.env()
  counter$draws <- 0
  results <- simulated_design_results()
  results$success_proba <- c(0.03, 0.8, 0.03, 0.85)
  results$mcse_success_proba <- c(0, 0.02, 0, 0.02)
  results$conf_int_success_proba_lower <- c(0.03, 0.76, 0.03, 0.81)
  results$conf_int_success_proba_upper <- c(0.03, 0.84, 0.03, 0.89)

  with_mocked_bindings(
    {
      first <- frequentist_power_at_equivalent_tie(
        results = results,
        analysis_config = list(frequentist_test = "t-test"),
        simulation_config = list(seed = 7),
        n_replicates = 500
      )
      second <- frequentist_power_at_equivalent_tie(
        results = results,
        analysis_config = list(frequentist_test = "t-test"),
        simulation_config = list(seed = 7),
        n_replicates = 500
      )
    },
    load_data = counting_load_data(counter),
    .package = "BExTEsim"
  )

  # A degenerate type I error leaves the simulated rejection count as the
  # only input, and both methods read it off the same trials.
  power <- first$frequentist_power_at_equivalent_tie
  expect_equal(power[1], power[3], tolerance = 0.05)
  # Seeded once per call, so a repeat reproduces every column.
  expect_identical(first, second)
})

test_that("the nominal-TIE step reads the trials the equivalent-TIE step cached", {
  counter <- new.env()
  counter$draws <- 0
  cache <- new.env(parent = emptyenv())
  results <- simulated_design_results()

  handed <- new.env()
  handed$p_values <- list()
  handed$trials <- list()
  record_separate <- function(..., p_values = NULL) {
    handed$p_values[[length(handed$p_values) + 1]] <- p_values
    list(power = 0.5, conf_int_power = c(0.4, 0.6))
  }
  record_pooling <- function(..., trials = NULL) {
    handed$trials[[length(handed$trials) + 1]] <- trials
    list(power = 0.5, conf_int_power = c(0.4, 0.6))
  }

  with_mocked_bindings(
    {
      frequentist_power_at_equivalent_tie(
        results = results,
        analysis_config = list(frequentist_test = "t-test"),
        simulation_config = list(seed = 7),
        n_replicates = 200,
        trial_cache = cache
      )
      frequentist_power_at_nominal_tie(
        results = results,
        analysis_config = list(frequentist_test = "t-test", nominal_tie = 0.025),
        simulation_config = list(seed = 7),
        n_replicates = 200,
        trial_cache = cache
      )
    },
    load_data = counting_load_data(counter),
    compute_freq_power = record_separate,
    compute_freq_power_pooling = record_pooling,
    .package = "BExTEsim"
  )

  expect_equal(counter$draws, 2)
  expect_length(handed$p_values, 2)
  expect_true(all(lengths(handed$p_values) == 200))
  # The pooled power reads the very same trials.
  expect_length(handed$trials, 2)
  expect_true(all(vapply(handed$trials, nrow, integer(1)) == 200))
})

test_that("without a cache, each design's trials are generated once for both baselines", {
  counter <- new.env()
  counter$draws <- 0

  with_mocked_bindings(
    powers <- frequentist_power_at_nominal_tie(
      results = simulated_design_results(),
      analysis_config = list(frequentist_test = "t-test", nominal_tie = 0.025),
      simulation_config = list(seed = 7),
      n_replicates = 200
    ),
    load_data = counting_load_data(counter),
    .package = "BExTEsim"
  )

  # Two designs; the separate and the pooled power used to draw them each.
  expect_equal(counter$draws, 2)
  expect_false(anyNA(powers$nominal_frequentist_power_separate))
  expect_false(anyNA(powers$nominal_frequentist_power_pooling))
})

test_that("shared trials give the powers a fresh simulation gives", {
  counter <- new.env()
  counter$draws <- 0
  cache <- new.env(parent = emptyenv())
  config <- list(frequentist_test = "t-test", nominal_tie = 0.025)

  with_mocked_bindings(
    {
      fresh <- frequentist_power_at_nominal_tie(
        simulated_design_results(), config, list(seed = 7), n_replicates = 300
      )
      frequentist_power_at_equivalent_tie(
        simulated_design_results(), config, list(seed = 7), n_replicates = 300,
        trial_cache = cache
      )
      shared <- frequentist_power_at_nominal_tie(
        simulated_design_results(), config, list(seed = 7), n_replicates = 300,
        trial_cache = cache
      )
    },
    load_data = counting_load_data(counter),
    .package = "BExTEsim"
  )

  columns <- grep("^nominal_frequentist_power", names(fresh), value = TRUE)
  expect_identical(shared[columns], fresh[columns])
})

test_that("a cache built for other settings is not reused", {
  keys <- trial_cache_key("design", list(seed = 1), 1000)
  expect_false(keys == trial_cache_key("design", list(seed = 2), 1000))
  expect_false(keys == trial_cache_key("design", list(seed = 1), 10000))
  expect_false(keys == trial_cache_key("other", list(seed = 1), 1000))
})

test_that("the design cluster is sized by the trials simulated, not the designs", {
  # 24 designs of 2,000 trials took 6.7s sequentially and 28s on a cluster.
  expect_false(analysis_uses_cluster(TRUE, 24 * 2000, min_rows = ANALYSIS_PARALLEL_MIN_TRIALS))
  expect_true(analysis_uses_cluster(TRUE, 264 * 10000, min_rows = ANALYSIS_PARALLEL_MIN_TRIALS))

  # deparse() wraps long lines, so whitespace is squeezed before matching.
  body_text <- gsub("\\s+", " ", paste(deparse(body(equivalent_tie_design_p_values)), collapse = " "))
  expect_match(body_text, "length(to_simulate) * n_replicates", fixed = TRUE)
})

test_that("the workers are handed the simulation as a variable, not by name", {
  # A foreach body is evaluated outside the package namespace, so an internal
  # function called by name there is "not found" in every worker.
  # deparse() wraps long lines, so whitespace is squeezed before matching.
  body_text <- gsub("\\s+", " ", paste(deparse(body(equivalent_tie_design_p_values)), collapse = " "))
  expect_match(body_text, "simulate_design <- design_trials", fixed = TRUE)
  expect_false(grepl("%dopar% { design_trials", body_text, fixed = TRUE))
})

test_that("the shared analysis cluster starts only when a step asks for it", {
  shared <- new_analysis_cluster()
  # Stopping one that never started is a no-op, which is what a run whose
  # steps all stayed sequential does on its way out.
  expect_null(shared$stop())

  body_text <- gsub("\\s+", " ", paste(deparse(body(simulation_analysis)), collapse = " "))
  expect_match(body_text, "cluster <- new_analysis_cluster()", fixed = TRUE)
  hits <- gregexpr("cluster = cluster", body_text, fixed = TRUE)[[1]]
  # Equivalent TIE, nominal TIE and the Bayesian OCs.
  expect_equal(sum(hits > 0), 3L)
})
