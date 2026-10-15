## Characterisation of the scenario loops of the posterior and empirical-Bayes
## plot and table wrappers: which leaf calls they make, with which scenario.

scenario_rows <- function() {
  grid <- expand.grid(
    case_study = c("A", "B"),
    method = c("m1", "m2"),
    target_sample_size_control = c(10, 20),
    target_to_source_std_ratio = c(1, 2, NA),
    source_denominator_change_factor = c(1, NA),
    parameters = c("p1", "p2"),
    stringsAsFactors = FALSE
  )
  grid$target_sample_size_treatment <- grid$target_sample_size_control
  grid
}

record_leaf_calls <- function(wrapper, leaf_names, env = parent.frame()) {
  calls <- list()
  assign("xvars", list(), envir = globalenv())
  withr::defer(suppressWarnings(rm("xvars", envir = globalenv())), envir = env)
  recorder <- function(leaf) {
    function(results_metrics_df, method, case_study,
             target_sample_size_control = NULL, target_sample_size_treatment = NULL, ...) {
      args <- list(...)
      calls[[length(calls) + 1L]] <<- list(
        leaf = leaf,
        case_study = case_study,
        method = method,
        size = target_sample_size_control,
        ratio = args$target_to_source_std_ratio,
        factor = args$source_denominator_change_factor,
        parameters = if (!is.null(args$parameters_combinations)) args$parameters_combinations$parameters,
        n_rows = nrow(results_metrics_df)
      )
      invisible(NULL)
    }
  }
  bindings <- stats::setNames(lapply(leaf_names, recorder), leaf_names)
  bindings$check_plot_packages <- function() invisible(NULL)
  bindings$check_table_packages <- function() invisible(NULL)
  bindings$check_required_colnames <- function(...) invisible(NULL)
  do.call(testthat::local_mocked_bindings, c(bindings, list(.env = env)))
  wrapper(scenario_rows())
  calls
}

calls_key <- function(calls) {
  df <- do.call(rbind, lapply(calls, function(x) {
    data.frame(
      leaf = x$leaf, case_study = x$case_study, method = x$method,
      size = if (is.null(x$size)) NA else x$size,
      ratio = x$ratio, factor = x$factor,
      parameters = if (is.null(x$parameters)) NA_character_ else x$parameters,
      stringsAsFactors = FALSE
    )
  }))
  df[do.call(order, unname(as.list(df))), ]
}

test_that("the empirical-Bayes plot and table wrappers visit each scenario once", {
  for (case in list(
    list(fn = empirical_bayes_parameters_plots, leaf = "plot_empirical_bayes_hyperparameters_vs_drift"),
    list(fn = function(df) table_empirical_bayes_hyperparameters_vs_scenario(df, NULL),
         leaf = "table_empirical_bayes_hyperparameters_vs_drift")
  )) {
    calls <- record_leaf_calls(case$fn, case$leaf)
    key <- calls_key(calls)
    # 2 case studies x 2 methods x 3 ratios x 2 factors x 2 parameter combinations
    expect_equal(nrow(key), 2 * 2 * 3 * 2 * 2)
    expect_false(anyDuplicated(key) > 0)
    expect_setequal(key$parameters, c("p1", "p2"))
    expect_setequal(key$ratio, c(1, 2, NA))
  }
})

test_that("the posterior plot and table wrappers visit each scenario and size once", {
  for (case in list(
    list(fn = posterior_parameters_plots,
         leaves = c("plot_posterior_parameters_vs_drift", "plot_posterior_vs_prior_parameters")),
    list(fn = function(df) posterior_parameters_tables(df, NULL),
         leaves = "table_posterior_parameters_vs_drift")
  )) {
    calls <- record_leaf_calls(case$fn, case$leaves)
    key <- calls_key(calls)
    # 2 case studies x 2 methods x 2 sizes x 3 ratios x 2 factors, per leaf
    expect_equal(nrow(key), 2 * 2 * 2 * 3 * 2 * length(case$leaves))
    expect_false(anyDuplicated(key) > 0)
    expect_setequal(key$size, c(10, 20))
  }
})

test_that("the empirical-Bayes wrappers hand each leaf the rows of its scenario", {
  data <- scenario_rows()
  calls <- record_leaf_calls(empirical_bayes_parameters_plots, "plot_empirical_bayes_hyperparameters_vs_drift")
  for (x in calls) {
    keep <- data$case_study == x$case_study & data$method == x$method &
      (is.na(data$target_to_source_std_ratio) | data$target_to_source_std_ratio %in% x$ratio) &
      (is.na(data$source_denominator_change_factor) | data$source_denominator_change_factor %in% x$factor)
    expect_equal(x$n_rows, sum(keep))
  }
})
