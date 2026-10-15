## Figure titles name the method in full - "Conditional Power Prior", not
## "Conditional PP". The abbreviations stay in legends and axis labels, where
## space is tight and the same name is repeated many times over.

methods_plots_config <- function() {
  env <- new.env()
  ## The file assigns with <<-, which needs an enclosing binding to find.
  eval(parse(system.file("conf/methods_plots_config.R", package = "BExTEsim")), envir = env)
  env
}

test_that("the full names spell out the abbreviations rather than repeating them", {
  config <- methods_plots_config()
  labels <- get("methods_labels", envir = config)

  expect_equal(labels$conditional_power_prior$full_name, "Conditional Power Prior")
  expect_equal(labels$EB_PP$full_name, "Empirical Bayes Power Prior")
  expect_equal(labels$NPP$full_name, "Normalized Power Prior")
  expect_equal(labels$RMP$full_name, "Robust Mixture Prior")

  ## The short forms are still there for the legends.
  expect_equal(labels$conditional_power_prior$label, "Conditional PP")
  expect_equal(labels$commensurate_power_prior$label, "Com. PP")
})

test_that("no figure title is built from the abbreviated label", {
  ## The titles were switched to full_name in one pass; this keeps a new
  ## title from quietly reintroducing the abbreviation.
  plot_sources <- list.files(
    testthat::test_path("..", "..", "R"), pattern = "^plot_.*\\.R$", full.names = TRUE
  )
  skip_if(length(plot_sources) == 0, "package sources not available")

  offenders <- Filter(
    function(path) any(grepl("methods_labels\\[\\[method\\]\\]\\$label",
                             readLines(path, warn = FALSE))),
    plot_sources
  )

  expect_equal(basename(offenders), character(0))
})

test_that("the commensurate priors are named in short form in titles", {
  ## Their titles also carry the heterogeneity prior, e.g. tau^2 ~ IG(alpha =
  ## 0.14, beta = 1), and the spelt-out name pushed the sample size off the
  ## edge of the figure.
  config <- methods_plots_config()
  labels <- get("methods_labels", envir = config)

  expect_equal(labels$commensurate_power_prior$full_name, "Com. PP")
  expect_equal(labels$commensurate_prior$full_name, "Com. prior")
})

test_that("the p-value based power prior is named in short form in titles", {
  ## The spelt-out name runs long in a title that already carries the case
  ## study, the parameters and the sample size, so this one keeps the
  ## abbreviated form.
  config <- methods_plots_config()
  labels <- get("methods_labels", envir = config)

  expect_equal(labels$p_value_based_PP$full_name, "p-value based PP")
  ## The legend form is unchanged.
  expect_equal(labels$p_value_based_PP$label, "p-PP")
})

test_that("a source denominator change factor of 1 is left out of the title", {
  ## 1 means the source denominator is used as it stands, which is the
  ## default - naming it lengthens every title without saying anything.
  plain <- format_title("Belimumab", "belimumab",
                        source_denominator_change_factor = 1, as_latex = TRUE)
  expect_false(grepl("denominator", plain, ignore.case = TRUE))

  ## A factor that does change something is still named.
  changed <- format_title("Belimumab", "belimumab",
                          source_denominator_change_factor = 0.5, as_latex = TRUE)
  expect_match(changed, "Source denominator change factor = 0.5", fixed = TRUE)
})

test_that("case studies that never show the factor are unaffected", {
  for (case_study in c("botox", "dapagliflozin", "aprepitant")) {
    title <- format_title(case_study, case_study,
                          source_denominator_change_factor = 0.5, as_latex = TRUE)
    expect_false(grepl("denominator", title, ignore.case = TRUE), info = case_study)
  }
})

test_that("a target-to-source standard deviation ratio of 1 is left out", {
  ## 1 means the two standard deviations agree, which is the default.
  for (case_study in c("botox", "dapagliflozin")) {
    plain <- format_title(case_study, case_study,
                          target_to_source_std_ratio = 1, as_latex = TRUE)
    expect_false(grepl("sigma", plain, fixed = TRUE), info = case_study)

    ## A ratio that does change something is still named.
    changed <- format_title(case_study, case_study,
                            target_to_source_std_ratio = 2, as_latex = TRUE)
    expect_match(changed, "sigma_T/\\sigma_S", fixed = TRUE)
  }
})
