## The four case tables share one layout and differ in file and title.

power_cases_df <- function() {
  data.frame(
    case_study = c("belimumab", "belimumab", NA),
    target_sample_size_control = c(200, 140, 100),
    target_sample_size_treatment = c(200, 140, 100),
    method = "separate", parameters = "[[{'initial_prior': ['noninformative']}]]",
    source_denominator_change_factor = 1, target_to_source_std_ratio = 1
  )
}

capture_table <- function(table_fn) {
  out <- NULL
  testthat::local_mocked_bindings(
    check_table_packages = function() invisible(NULL),
    format_results_df_parameters = function(df) df$method,
    export_table = function(data_table, title, file_path, ...) {
      out <<- list(data = data_table, title = title, path = file_path)
    }
  )
  assign("tables_dir", "tables", envir = globalenv())
  withr::defer(rm("tables_dir", envir = globalenv()))
  table_fn(power_cases_df())
  out
}

capture_table_path <- function(table_fn) capture_table(table_fn)$path

test_that("power-loss tables write to distinct files", {
  p_loss <- capture_table_path(table_power_loss)
  p_inflated <- capture_table_path(table_power_loss_cases_inflated_tie)
  expect_false(identical(p_loss, p_inflated))
})

test_that("each case table keeps its file stem and title", {
  expected <- list(
    list(fn = table_noninflated_tie, stem = "noninflated_tie_cases", title = "no inflation of the TIE"),
    list(fn = table_power_gain, stem = "power_gain_cases", title = "power gains"),
    list(fn = table_power_loss, stem = "power_loss_cases", title = "power losss"),
    list(fn = table_power_loss_cases_inflated_tie, stem = "power_loss_cases_inflated_tie", title = "higher than the nominal TIE")
  )
  for (e in expected) {
    out <- capture_table(e$fn)
    expect_equal(basename(out$path), e$stem)
    expect_match(out$title, e$title, fixed = TRUE)
  }
})

test_that("case tables sort by total size, rename columns and keep six columns", {
  out <- capture_table(table_power_gain)$data
  expect_equal(names(out), c("Method", "Case study", "$n_C$", "$n_T$",
                             "Source denominator change factor", "$\\sigma_T/\\sigma_S$"))
  expect_equal(out[["$n_C$"]], c(100, 140, 200))
})

test_that("only the inflated-TIE power-loss table drops rows without a case study", {
  expect_equal(nrow(capture_table(table_power_loss)$data), 3L)
  expect_equal(nrow(capture_table(table_power_loss_cases_inflated_tie)$data), 2L)
})
