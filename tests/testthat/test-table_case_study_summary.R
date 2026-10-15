## Tables S1 and S3 are derived from the case study YAMLs, not from simulation
## output, so they can be generated before any run finishes.

config_dir <- function() {
  paste0(system.file("conf/case_studies", package = "BExTE"), "/")
}

## kableExtra's "striped" option wraps every cell of shaded rows in
## \cellcolor{gray!10}{...}; strip that so row assertions below don't have to
## depend on which rows happen to be shaded.
strip_shading <- function(contents) {
  gsub("\\\\cellcolor\\{gray!10\\}\\{([^}]*)\\}", "\\1", contents)
}

test_that("table_drift_ranges_and_sample_sizes lists every case study with its ranges and nominal sample sizes", {
  withr::with_tempdir({
    path <- table_drift_ranges_and_sample_sizes(config_dir(), getwd())

    expect_true(file.exists(path))
    contents <- strip_shading(paste(readLines(path), collapse = " "))

    expect_match(contents, "\\textbf{$N_S/6$}", fixed = TRUE)

    ## Match each row as a unit, cells in column order. The sample sizes are
    ## the nominal floor(N_S / k) the paper prints, not the doubled per-arm
    ## sizes: Botox N_S/4 is 117, although two arms of 58 are simulated.
    expect_match(
      contents,
      "Botox & $[-0.365, 0.365]$ & -0.2 & $[-0.165, 0.565]$ & 468 & 234 & 117 & 78",
      fixed = TRUE
    )
    ## Pooled source trials: 1125 for belimumab and 1483 for teriflunomide,
    ## not the single trials (577 and 761) an earlier table printed.
    expect_match(contents, "Belimumab & [^&]* & -0.481 & [^&]* & 1125 & 562 & 281 & 187")
    ## Teriflunomide's null lies to the right, so its range must reach the
    ## drift of 0.393 at which the target effect is null.
    expect_match(
      contents,
      "Teriflunomide & $[-0.652, 0.652]$ & 0.393 & $[-1.05, 0.259]$ & 1483 & 741 & 370 & 247",
      fixed = TRUE
    )
    ## The binomial case study takes its source effect from the response
    ## counts, 184/293 - 154/280 = 0.078, and its range is bounded by the
    ## rates staying in [0, 1].
    expect_match(
      contents,
      "Aprepitant & $[-0.628, 0.372]$ & -0.078 & $[-0.55, 0.45]$ & 573 & 286 & 143 & 95",
      fixed = TRUE
    )
  })
})

test_that("table_case_study_summary reports the source and target sizes as the sums of their arms", {
  withr::with_tempdir({
    path <- table_case_study_summary(c("botox", "aprepitant"), config_dir(), getwd())

    expect_true(file.exists(path))
    contents <- strip_shading(paste(readLines(path), collapse = " "))

    ## Column headers.
    expect_match(contents, "\\textbf{Source N}", fixed = TRUE)
    expect_match(contents, "\\textbf{Target N}", fixed = TRUE)

    ## Match each row as a unit, with Source N and Target N pinned to their
    ## column position (Endpoint/Summary measure/effect/SE cells wildcarded,
    ## since they aren't what this test is about) - so a value landing in
    ## the wrong column doesn't slip past a bare substring match.
    ##
    ## Botox source: 235 control + 233 treatment = 468. Target: 130 + 126 = 256.
    expect_match(contents, "Botox & Placebo & [^&]* & [^&]* & 468 & [^&]* & [^&]* & 256")

    ## Aprepitant source: 280 control + 293 treatment = 573, the arms of the
    ## adult trial as Jin et al. (2021) report them. Target: 52 + 57 = 109, the
    ## 125 mg arm of Salman et al. (2019), not the 55 Jin et al. print.
    expect_match(contents, "Aprepitant & Ondansetron & [^&]* & [^&]* & 573 & [^&]* & [^&]* & 109")
  })
})

test_that("table_case_study_summary covers every shipped case study", {
  withr::with_tempdir({
    case_studies <- sub(
      "\\.yml$", "", basename(list.files(config_dir(), pattern = "\\.yml$"))
    )

    ## The manifest's TS3 asks for all of them in one call, so a single YAML
    ## missing a field the table reads takes the whole table down - which is
    ## how teriflunomide.yml's absent `control:` surfaced, as "arguments
    ## imply differing number of rows: 1, 0" from data.frame() being handed a
    ## NULL. Naming two case studies explicitly, as the tests above do, can
    ## only ever catch that for the two that happen to be named.
    path <- table_case_study_summary(case_studies, config_dir(), getwd())

    expect_true(file.exists(path))
    contents <- strip_shading(paste(readLines(path), collapse = " "))
    expect_match(contents, "Teriflunomide & Placebo")

    ## Endpoints are YAML identifiers, and the two that carry an underscore
    ## belong to case studies neither of the tests above names - so a bare
    ## "_" reached the .tex and aborted the compile for the whole table.
    expect_match(contents, "time\\_to\\_event", fixed = TRUE)
    expect_match(contents, "recurrent\\_event", fixed = TRUE)
  })
})
