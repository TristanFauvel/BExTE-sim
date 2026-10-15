## Figures S9, S14, S17, S27, S32 and S37 show each method's operating
## characteristic as a ratio to the separate (no-borrowing) analysis in the
## same scenario. scale_by_separate() is that division; it matches on the full
## scenario key rather than collapsing across drift, and it refuses to invent a
## finite value when the denominator is zero or missing.

relative_fixture <- function() {
  readRDS(testthat::test_path("fixtures", "forest_plot_freq.rds"))
}

test_that("scale_by_separate divides by the separate row in the same scenario", {
  df <- relative_fixture()

  out <- scale_by_separate(df, "success_proba")

  expect_identical(names(out), names(df))
  ## Every separate row becomes exactly 1.
  separate_rows <- out[out$method == "separate", ]
  expect_equal(separate_rows$success_proba, rep(1, nrow(separate_rows)))

  ## A pooling row equals its own scenario's ratio, not a cross-drift average.
  ## (Index 2, not 1: the fixture's most extreme drift has a separate-method
  ## success_proba of exactly 0, which scale_by_separate() correctly drops as
  ## an undefined ratio, so that row is unavailable to compare against.)
  one_drift <- df$drift[df$method == "pooling"][2]
  numerator <- df$success_proba[df$method == "pooling" & df$drift == one_drift]
  denominator <- df$success_proba[df$method == "separate" & df$drift == one_drift]
  expect_equal(
    out$success_proba[out$method == "pooling" & out$drift == one_drift],
    numerator / denominator
  )
})

test_that("scale_by_separate scales the confidence bounds by the same denominator", {
  df <- relative_fixture()

  out <- scale_by_separate(
    df,
    c("success_proba", "conf_int_success_proba_lower", "conf_int_success_proba_upper")
  )

  one_drift <- df$drift[df$method == "pooling"][2]
  denominator <- df$success_proba[df$method == "separate" & df$drift == one_drift]
  expect_equal(
    out$conf_int_success_proba_lower[out$method == "pooling" & out$drift == one_drift],
    df$conf_int_success_proba_lower[df$method == "pooling" & df$drift == one_drift] / denominator
  )
})

test_that("scale_by_separate drops rows whose denominator is zero, never returning Inf", {
  df <- relative_fixture()
  zero_drift <- df$drift[df$method == "separate"][1]
  df$success_proba[df$method == "separate" & df$drift == zero_drift] <- 0

  expect_warning(out <- scale_by_separate(df, "success_proba"), "zero or missing")

  expect_false(any(is.infinite(out$success_proba)))
  expect_equal(sum(out$drift == zero_drift), 0)
})

test_that("scale_by_separate errors when there is no separate analysis to divide by", {
  df <- relative_fixture()
  df <- df[df$method != "separate", ]

  expect_error(scale_by_separate(df, "success_proba"), "separate")
})

setup_relative_plot_globals <- function(figures_dir) {
  source(system.file("conf/plots_config.R", package = "BExTEsim"))
  source(system.file("conf/methods_plots_config.R", package = "BExTEsim"))
  source(system.file("conf/metrics_config.R", package = "BExTEsim"))
  methods_env <- new.env()
  source(system.file("conf/full/methods_config.R", package = "BExTEsim"), local = methods_env)
  assign("methods_dict", methods_env$methods_dict, envir = .GlobalEnv)
  assign("figures_dir", figures_dir, envir = .GlobalEnv)
  assign("remake_figures", TRUE, envir = .GlobalEnv)
  invisible(NULL)
}

test_that("forest_plot writes a distinctly named file when relative_to_separate is TRUE", {
  skip_on_cran()
  df <- relative_fixture()
  withr::with_tempdir({
    setup_relative_plot_globals(file.path(getwd(), ""))

    forest_plot(df, "success_proba", relative_to_separate = TRUE)

    written <- list.files(getwd(), pattern = "\\.png$", recursive = TRUE)
    expect_equal(written, paste0(
      "botox/botox_relative_success_proba_forest_plot_",
      "target_nc=234_nt=234_target_to_source_std_ratio=1.png"
    ))
  })
})

test_that("forest_plot's default output is unchanged by the new argument", {
  df <- relative_fixture()
  withr::with_tempdir({
    setup_relative_plot_globals(file.path(getwd(), ""))

    forest_plot(df, "success_proba")

    written <- list.files(getwd(), pattern = "\\.png$", recursive = TRUE)
    expect_length(written, 1)
    expect_false(grepl("_relative_success_proba_", written, fixed = TRUE))
  })
})
