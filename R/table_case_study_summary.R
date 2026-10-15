## Supplementary tables S1 and S3. Both are derived from the case study YAMLs
## rather than from simulation output.

## LaTeX reads a bare "_" outside math mode as a subscript and errors out. The
## YAML values these tables print are identifiers, so they carry underscores
## routinely; export_table() escapes only "%", and only in column names, so
## escaping the cell values is left to the caller.
escape_latex_underscores <- function(x) {
  return(gsub("_", "\\\\_", x))
}

## The case studies in the order table S1 lists them.
PAPER_TABLE_S1_CASE_STUDIES <- c(
  "belimumab", "botox", "dapagliflozin", "mepolizumab", "aprepitant",
  "teriflunomide"
)

#' Drift ranges and target sample sizes for each case study (table S1)
#'
#' @description The paper's table S1, which merges what used to be two
#'   tables: the drift and treatment-effect ranges, and the target-study
#'   sample sizes. The ranges come from [compute_drift_range()], the function
#'   the simulation builds its drift grid with, so the table describes the
#'   grid the current configs produce rather than whatever an older results
#'   directory happens to hold.
#'
#'   The sample-size columns are the nominal totals `floor(N_S / k)`, as the
#'   paper prints them. The simulation puts `floor(N_S / (2k))` patients in
#'   each arm, so the total it actually simulates can be one or two lower.
#'
#'   Every case study is listed whatever the selection, because the paper's
#'   table always has all six rows.
#'
#' @param case_studies_config_dir Directory holding the case study YAMLs.
#' @param tables_dir Directory to write into.
#'
#' @return The path of the written `.tex` file.
#'
#' @keywords internal
table_drift_ranges_and_sample_sizes <- function(case_studies_config_dir, tables_dir) {
  check_table_packages()
  dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

  format_range <- function(bounds) {
    return(paste0("$[", signif(bounds[1], 3), ", ", signif(bounds[2], 3), "]$"))
  }

  rows <- lapply(PAPER_TABLE_S1_CASE_STUDIES, function(case_study) {
    config <- yaml::read_yaml(file.path(case_studies_config_dir, paste0(case_study, ".yml")))

    source_effect <- if (config$summary_measure_likelihood == "binomial") {
      config$source$responses$treatment / config$source$treatment -
        config$source$responses$control / config$source$control
    } else {
      config$source$treatment_effect
    }
    ## Only the end points are used, and they do not depend on ndrift.
    drift <- range(compute_drift_range(list(ndrift = 2), config))
    ## control + treatment, not the `total:` field.
    source_total <- config$source$control + config$source$treatment

    return(data.frame(
      `Case study` = escape_latex_underscores(config$name),
      `Drift range` = format_range(drift),
      ## As text, so kable's digits = 3 prints -0.2 rather than -0.200.
      `Drift when $\\theta_T^{(true)} = \\theta_0$` =
        as.character(signif(config$theta_0 - source_effect, 3)),
      `Treatment effect range` = format_range(drift + source_effect),
      `$N_S$` = source_total,
      `$N_S/2$` = floor(source_total / 2),
      `$N_S/4$` = floor(source_total / 4),
      `$N_S/6$` = floor(source_total / 6),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  })

  table_data <- do.call(rbind, rows)
  file_path <- file.path(tables_dir, "drift_ranges_and_sample_sizes")

  export_table(
    data_table = table_data,
    title = "Drift and treatment effect ranges, and target study sample sizes considered for each case study.",
    file_path = file_path,
    ncollapses = 0,
    longtable = FALSE
  )

  return(paste0(file_path, ".tex"))
}

#' Summary of the clinical case studies (table S3)
#'
#' @param case_studies Character vector of case study names.
#' @param case_studies_config_dir Directory holding the case study YAMLs.
#' @param tables_dir Directory to write into.
#'
#' @return The path of the written `.tex` file.
#'
#' @keywords internal
table_case_study_summary <- function(case_studies, case_studies_config_dir, tables_dir) {
  check_table_packages()
  dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

  rows <- lapply(case_studies, function(case_study) {
    config <- yaml::read_yaml(file.path(case_studies_config_dir, paste0(case_study, ".yml")))
    return(data.frame(
      `Case study` = escape_latex_underscores(config$name),
      `Control` = escape_latex_underscores(config$control),
      ## Endpoints are YAML identifiers, so several carry an underscore
      ## (`time_to_event`, `recurrent_event`). export_table() escapes only
      ## "%", and only in column names, so an unescaped one reaches the .tex
      ## as a subscript outside math mode and aborts the compile - taking the
      ## whole table down, not just that row.
      `Endpoint` = escape_latex_underscores(config$endpoint),
      `Summary measure` = escape_latex_underscores(config$summary_measure_likelihood),
      ## control + treatment, not the `total:` field.
      `Source N` = config$source$control + config$source$treatment,
      `Source effect` = round(config$source$treatment_effect, 4),
      `Source SE` = round(config$source$standard_error, 4),
      `Target N` = config$target$control + config$target$treatment,
      `Target effect` = round(config$target$treatment_effect, 4),
      `Target SE` = round(config$target$standard_error, 4),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  })

  table_data <- do.call(rbind, rows)
  file_path <- file.path(tables_dir, "case_study_summary")

  export_table(
    data_table = table_data,
    title = "Summary of the clinical case studies used to construct the simulation-study design.",
    file_path = file_path,
    longtable = FALSE
  )

  return(paste0(file_path, ".tex"))
}
