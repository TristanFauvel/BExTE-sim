## Neither baseline the nominal-TIE step computes depends on the borrowing
## method, but the results frame holds one row per design *and*
## method-parameter combination. The paper's environment repeats each of its
## 330 designs 56 times, which was 32 CPU-hours spent to produce 330 pairs of
## numbers. They are now computed once per design and copied.

design_frame <- function(sample_sizes, methods) {
  frame <- expand.grid(
    target_sample_size_control = sample_sizes,
    method = methods,
    stringsAsFactors = FALSE
  )
  frame$target_sample_size_treatment <- frame$target_sample_size_control
  frame$case_study <- "example"
  frame$parameters <- paste0('{"w":', seq_len(nrow(frame)), "}")
  frame$theta_0 <- 0
  frame$null_space <- "left"
  frame$target_treatment_effect <- 0.5
  frame
}

## Counting mock: the power is read off the design so that a value landing on
## the wrong row is visible, and every call is tallied so that a cache which
## silently recomputes cannot pass.
counting_mocks <- function() {
  tally <- new.env(parent = emptyenv())
  tally$separate <- 0L
  tally$pooling <- 0L

  list(
    tally = tally,
    load_data = function(results_row, type, reload_data_objects = FALSE) {
      list(type = type, treatment_effect = 0.5, standard_deviation = 1,
           sample_size_control = results_row$target_sample_size_control,
           sample_size_treatment = results_row$target_sample_size_treatment)
    },
    separate = function(alpha, target_data, ...) {
      tally$separate <- tally$separate + 1L
      power <- target_data$sample_size_control / 1000
      list(power = power, conf_int_power = c(power - 0.01, power + 0.01))
    },
    pooling = function(alpha, target_data, ...) {
      tally$pooling <- tally$pooling + 1L
      power <- target_data$sample_size_control / 500
      list(power = power, conf_int_power = c(power - 0.02, power + 0.02))
    }
  )
}

run_cached <- function(results, mocks) {
  final <- NULL
  with_mocked_bindings(
    capture.output(
      final <- frequentist_power_at_nominal_tie(
        results = results,
        analysis_config = list(nominal_tie = 0.025, frequentist_test = "t-test"),
        simulation_config = list(),
        n_replicates = 100
      )
    ),
    load_data = mocks$load_data,
    compute_freq_power = mocks$separate,
    compute_freq_power_pooling = mocks$pooling,
    .package = "BExTEsim"
  )
  final
}

test_that("rows that differ only in the method share one computation", {
  mocks <- counting_mocks()
  results <- design_frame(c(100, 200), c("RMP", "separate", "NPP"))
  expect_equal(nrow(results), 6)

  final <- run_cached(results, mocks)

  # Two designs, not six rows.
  expect_equal(mocks$tally$separate, 2L)
  expect_equal(mocks$tally$pooling, 2L)

  # And every row still carries the value for its own design.
  expect_equal(final$nominal_frequentist_power_separate,
               results$target_sample_size_control / 1000)
  expect_equal(final$nominal_frequentist_power_pooling,
               results$target_sample_size_control / 500)
  expect_equal(final$nominal_frequentist_power_separate_lower,
               results$target_sample_size_control / 1000 - 0.01)
  expect_equal(final$nominal_frequentist_power_pooling_upper,
               results$target_sample_size_control / 500 + 0.02)
})

test_that("rows that differ in a design column do not share a computation", {
  mocks <- counting_mocks()
  results <- design_frame(c(100, 200), c("RMP", "separate"))
  # A design axis the borrowing method has nothing to do with.
  results$control_drift <- c(0, 0, 0.3, 0.3)

  final <- run_cached(results, mocks)

  expect_equal(mocks$tally$separate, 4L)
  expect_equal(final$nominal_frequentist_power_separate,
               results$target_sample_size_control / 1000)
})

test_that("the cached value is the one the row would have computed itself", {
  # The cache must be exact, not close: run the same frame one row at a time
  # and compare. set.seed() inside simulate_test_p_values() is what makes
  # this true of the real computation too.
  mocks <- counting_mocks()
  results <- design_frame(c(100, 200, 300), c("RMP", "separate"))

  cached <- run_cached(results, mocks)
  one_at_a_time <- do.call(rbind, lapply(seq_len(nrow(results)), function(i) {
    run_cached(results[i, , drop = FALSE], counting_mocks())
  }))

  columns <- grep("^nominal_frequentist_power", names(cached), value = TRUE)
  expect_equal(length(columns), 6)
  for (column in columns) {
    expect_equal(cached[[column]], one_at_a_time[[column]], info = column)
  }
})

test_that("a missing entry never collides with the string NA", {
  design <- data.frame(
    endpoint = c(NA, "NA"),
    case_study = c("example", "example"),
    stringsAsFactors = FALSE
  )
  key <- nominal_tie_design_key(design)
  expect_equal(length(unique(key)), 2)

  same <- nominal_tie_design_key(data.frame(endpoint = c(NA, NA), case_study = c("a", "a")))
  expect_equal(length(unique(same)), 1)
})

test_that("the design columns cover everything the computation reads", {
  # A column read by nominal_tie_power_row() but absent here would merge rows
  # that genuinely differ, and the wrong baseline would be copied onto them.
  read_by <- function(fn, symbol) {
    text <- paste(deparse(body(fn)), collapse = " ")
    matches <- gregexpr(paste0(symbol, "\\$[A-Za-z_.][A-Za-z0-9_.]*"), text)
    unique(sub(paste0(symbol, "\\$"), "", regmatches(text, matches)[[1]]))
  }

  used <- union(read_by(nominal_tie_power_row, "row"),
                read_by(load_data, "results_row"))
  expect_true(length(used) > 0)
  expect_setequal(intersect(used, NOMINAL_TIE_DESIGN_COLUMNS), used)
})
