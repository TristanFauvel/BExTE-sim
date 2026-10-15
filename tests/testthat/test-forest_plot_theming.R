## Forest plots can be drawn for a dark surface (BExTE-app's Analyze tab
## renders them as static PNGs). The default ink only works on white:
## geom_pointrange() had no colour, so the point estimates came out black and
## vanished against the dark surface. forest_plot()/forest_plot_bayesian()
## therefore take an optional palette; NULL keeps the publication figures
## exactly as they were.

## The plot_*() functions resolve their style constants and methods_dict as
## free variables out of .GlobalEnv (see inst/scripts/plots.R), so a test that
## drives them has to stand those up first.
setup_forest_plot_globals <- function(figures_dir) {
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

## BExTE-app's dark scheme (bexte_palette("dark") in its theme.R).
dark_palette <- function() {
  list(ink = "#E5E9EF", ink_muted = "#94A3B8", hairline = "#2A3138",
       accent = "#2DD4BF", surface = "#171C22",
       ref_target = "#E5E9EF", ref_source = "#FF8A8A")
}

forest_plot_fixture <- function() {
  readRDS(testthat::test_path("fixtures", "forest_plot_freq.rds"))
}

test_that("forest_subplot draws point ranges in the palette ink when given a palette", {
  data <- forest_plot_fixture()
  data <- data[data$drift == data$drift[1], , drop = FALSE]
  data$rows <- seq_len(nrow(data))
  withr::with_tempdir({
    setup_forest_plot_globals(file.path(getwd(), ""))
    pal <- dark_palette()

    plt <- forest_subplot(
      data, "No effect",
      ylabel = TRUE,
      x_metric_name = rlang::sym("success_proba"),
      x_metric_uncertainty_lower = "conf_int_success_proba_lower",
      x_metric_uncertainty_upper = "conf_int_success_proba_upper",
      x_metric_label = "Success probability",
      methods_labels = methods_labels,
      legend = FALSE,
      sort_by = FALSE,
      palette = pal
    )

    pointrange <- Filter(
      function(l) inherits(l$geom, "GeomPointrange"),
      plt$layers
    )
    expect_length(pointrange, 1)
    expect_identical(pointrange[[1]]$aes_params$colour, pal$ink)
  })
})

test_that("forest_subplot leaves point ranges at the ggplot2 default without a palette", {
  data <- forest_plot_fixture()
  data <- data[data$drift == data$drift[1], , drop = FALSE]
  data$rows <- seq_len(nrow(data))
  withr::with_tempdir({
    setup_forest_plot_globals(file.path(getwd(), ""))

    plt <- forest_subplot(
      data, "No effect",
      ylabel = TRUE,
      x_metric_name = rlang::sym("success_proba"),
      x_metric_uncertainty_lower = "conf_int_success_proba_lower",
      x_metric_uncertainty_upper = "conf_int_success_proba_upper",
      x_metric_label = "Success probability",
      methods_labels = methods_labels,
      legend = FALSE,
      sort_by = FALSE
    )

    pointrange <- Filter(
      function(l) inherits(l$geom, "GeomPointrange"),
      plt$layers
    )
    expect_null(pointrange[[1]]$aes_params$colour)
  })
})

test_that("forest_plot writes a separate dark PNG and leaves the light figures untouched", {
  data <- forest_plot_fixture()
  withr::with_tempdir({
    figures_dir <- file.path(getwd(), "figures", "")
    dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
    setup_forest_plot_globals(figures_dir)

    light_path <- suppressWarnings(forest_plot(data, "success_proba"))
    expect_true(file.exists(light_path))
    light_pdf <- sub("\\.png$", ".pdf", light_path)
    expect_true(file.exists(light_pdf))
    light_png_before <- readBin(light_path, "raw", file.size(light_path))
    light_pdf_mtime <- file.mtime(light_pdf)

    dark_path <- suppressWarnings(
      forest_plot(data, "success_proba", palette = dark_palette())
    )

    ## A distinct file, so the publication figure is never overwritten.
    expect_match(basename(dark_path), "_dark\\.png$")
    expect_true(file.exists(dark_path))
    expect_false(identical(normalizePath(dark_path), normalizePath(light_path)))

    ## The light PNG is byte-identical and no second PDF was written.
    expect_identical(readBin(light_path, "raw", file.size(light_path)), light_png_before)
    expect_identical(file.mtime(light_pdf), light_pdf_mtime)
    expect_false(file.exists(sub("\\.png$", ".pdf", dark_path)))

    ## The dark render is actually different ink, not a copy.
    expect_false(identical(
      readBin(dark_path, "raw", file.size(dark_path)),
      light_png_before
    ))
  })
})
