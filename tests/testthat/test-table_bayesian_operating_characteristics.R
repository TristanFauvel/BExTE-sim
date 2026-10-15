## Characterisation of the Bayesian tables: the data frame, title and file path
## each function hands to export_table().

bayes_rows <- function() {
  grid <- expand.grid(
    case_study = c("belimumab", "mepolizumab"),
    target_sample_size_control = c(140, 100),
    design_prior_type = c("analysis_prior", "ui_design_prior", "source_posterior"),
    stringsAsFactors = FALSE
  )
  n <- nrow(grid)
  grid$target_sample_size_treatment <- grid$target_sample_size_control
  grid$method <- "separate"
  grid$parameters <- "[[{'initial_prior': ['noninformative']}]]"
  grid$source_denominator_change_factor <- c(1, NA)[(seq_len(n) %% 2) + 1]
  grid$target_to_source_std_ratio <- 1
  grid$source_treatment_effect_estimate <- 0.5
  grid$source_standard_error <- 0.1
  grid$endpoint <- "e"
  grid$summary_measure_likelihood <- "normal"
  grid$source_sample_size_control <- 100
  grid$source_sample_size_treatment <- 100
  grid$equivalent_source_sample_size_per_arm <- 50
  grid$prior_proba_benefit <- seq(0.1, 0.9, length.out = n)
  grid$coverage <- seq(0.5, 0.99, length.out = n)
  grid
}

with_bayes_globals <- function(env = parent.frame(), metrics = c("coverage", "prior_proba_benefit")) {
  labels <- c(coverage = "Coverage", prior_proba_benefit = "Prior probability of benefit")
  assign("bayesian_metrics", lapply(metrics, function(m) list(name = m, label = labels[[m]])),
         envir = globalenv())
  assign("tables_dir", withr::local_tempdir(.local_envir = env), envir = globalenv())
  withr::defer({ rm("bayesian_metrics", envir = globalenv()); rm("tables_dir", envir = globalenv()) }, envir = env)
}

record_exports <- function(expr_fn, env = parent.frame()) {
  calls <- list()
  testthat::local_mocked_bindings(
    check_table_packages = function() invisible(NULL),
    format_results_df_parameters = function(df) df$method,
    export_table = function(data_table, title, file_path, ...) {
      calls[[length(calls) + 1L]] <<- list(
        data = as.data.frame(data_table), title = title, file = basename(file_path)
      )
    },
    .env = env
  )
  expr_fn()
  calls
}

test_that("table_bayesian_metrics exports one table per metric", {
  with_bayes_globals()
  calls <- record_exports(function() {
    table_bayesian_metrics(bayes_rows(), "belimumab", 140, 140)
  })
  expect_length(calls, 2L)
  expect_equal(sort(names(calls[[1]]$data)),
               sort(c("Method", "Analysis prior", "UI design prior", "Source posterior")))
  expect_equal(calls[[1]]$title, "Coverage, Belimumab, $n_C = 140$, $n_T = 140$")
  expect_equal(calls[[1]]$file, "coverage_belimumab_nc=140_nt=140_source_denominator_change_factor=1")
  expect_equal(names(calls[[2]]$data), c("Method", "Prior probability of benefit"))
  expect_equal(calls[[2]]$title, "Prior probability of benefit, Belimumab, $n_C = 140$, $n_T = 140$")
  expect_equal(calls[[2]]$file, "prior_proba_benefit_belimumab_nc=140_nt=140_source_denominator_change_factor=1")
})

test_that("the sample-size comparison has one column per target size, smallest first", {
  with_bayes_globals()
  calls <- record_exports(function() {
    table_bayesian_metrics_sample_sizes_comparison(bayes_rows(), "belimumab")
  })
  expect_length(calls, 2L)
  expect_equal(names(calls[[1]]$data), c("Method", "$n_C = 100$, $n_T = 100$", "$n_C = 140$, $n_T = 140$"))
  expect_equal(calls[[1]]$title, "Coverage, Belimumab")
  expect_equal(calls[[1]]$file, "belimumab_coverage_analysis_prior_source_denominator_change_factor=1")
  expect_equal(calls[[2]]$file, "belimumab_prior_proba_benefit_analysis_prior_source_denominator_change_factor=1")
})

test_that("the across-case-study table exports only the prior probabilities", {
  with_bayes_globals()
  calls <- record_exports(function() {
    table_bayesian_metrics_across_case_studies(bayes_rows())
  })
  expect_length(calls, 1L)
  expect_equal(names(calls[[1]]$data), c("Method", "Belimumab", "Mepolizumab"))
  expect_equal(calls[[1]]$title, "Prior probability of benefit")
  expect_equal(calls[[1]]$file, "prior_proba_benefit_analysis_prior_source_denominator_change_factor=1")
})

test_that("a metric after a prior probability still has one column per design prior", {
  with_bayes_globals(metrics = c("prior_proba_benefit", "coverage"))
  calls <- record_exports(function() {
    table_bayesian_metrics(bayes_rows(), "belimumab", 140, 140)
  })
  expect_equal(names(calls[[1]]$data), c("Method", "Prior probability of benefit"))
  expect_equal(sort(names(calls[[2]]$data)),
               sort(c("Method", "Analysis prior", "UI design prior", "Source posterior")))
})
