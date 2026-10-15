## Table S5 combines the interval diagnostics - precision, coverage and the
## interval score, which charges an interval for its width and for excluding
## the true effect in one number - for the same three treatment-effect
## scenarios used by the Belimumab forest plots.
table_precision_ecp <- function(results_df, tables_dir) {
  check_table_packages()
  required <- c(
    "case_study", "target_sample_size_control", "target_sample_size_treatment",
    "method", "parameters",
    "target_to_source_std_ratio", "source_denominator_change_factor",
    "source_treatment_effect_estimate", "target_treatment_effect",
    "precision", "conf_int_precision_lower", "conf_int_precision_upper",
    "coverage", "conf_int_coverage_lower", "conf_int_coverage_upper",
    "interval_score", "conf_int_interval_score_lower", "conf_int_interval_score_upper"
  )
  missing <- setdiff(required, names(results_df))
  if (length(missing)) {
    stop("Precision/ECP table is missing columns: ", paste(missing, collapse = ", "))
  }
  if (!nrow(results_df) || length(unique(results_df$case_study)) != 1L ||
      nrow(unique(results_df[, target_size_columns, drop = FALSE])) != 1L) {
    stop("Precision/ECP table needs one nonempty case-study and sample-size slice.")
  }

  source_effect <- unique(results_df$source_treatment_effect_estimate)
  if (length(source_effect) != 1L || !is.finite(source_effect)) {
    stop("Precision/ECP table needs one finite source treatment effect.")
  }
  target_effects <- c(0, source_effect / 2, source_effect)
  observed <- unique(results_df$target_treatment_effect)
  selected <- vapply(target_effects, function(effect) {
    return(observed[which.min(abs(observed - effect))])
  }, numeric(1))
  if (length(unique(selected)) != 3L) {
    stop("Precision/ECP table could not identify three distinct treatment effects.")
  }

  rows <- results_df[results_df$target_treatment_effect %in% selected, , drop = FALSE]
  if (anyNA(rows[, c("precision", "coverage", "interval_score",
                     "conf_int_precision_lower", "conf_int_precision_upper",
                     "conf_int_coverage_lower", "conf_int_coverage_upper",
                     "conf_int_interval_score_lower", "conf_int_interval_score_upper")])) {
    stop("Precision/ECP table has missing estimates or confidence limits.")
  }
  scenario <- match(rows$target_treatment_effect, selected)
  effect_labels <- c("No treatment effect", "Partially consistent", "Consistent")
  table_data <- data.frame(
    Method = format_results_df_parameters(rows),
    `Treatment effect` = effect_labels[scenario],
    check.names = FALSE
  )
  ## Operating characteristics computed exactly have intervals of zero width,
  ## which would only repeat the estimate twice; such a column shows the
  ## estimate alone, which also keeps the table within the page.
  add_column <- function(name, estimate, lower, upper) {
    exact <- all(sprintf("%.3f", lower) == sprintf("%.3f", upper))
    if (exact) {
      table_data[[name]] <<- sprintf("%.3f", estimate)
      return(invisible(table_data[[name]]))
    } else {
      table_data[[paste(name, "(95% CI)")]] <<- sprintf("%.3f [%.3f, %.3f]", estimate, lower, upper)
      return(invisible(table_data[[paste(name, "(95% CI)")]]))
    }
  }
  add_column("Precision", rows$precision,
             rows$conf_int_precision_lower, rows$conf_int_precision_upper)
  add_column("ECP", rows$coverage,
             rows$conf_int_coverage_lower, rows$conf_int_coverage_upper)
  add_column("Interval score", rows$interval_score,
             rows$conf_int_interval_score_lower, rows$conf_int_interval_score_upper)
  table_data <- table_data[order(table_data$Method, scenario), , drop = FALSE]
  rownames(table_data) <- NULL

  case_study <- unique(rows$case_study)
  n_control <- unique(rows$target_sample_size_control)
  n_treatment <- unique(rows$target_sample_size_treatment)
  directory <- file.path(tables_dir, case_study)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  filename <- paste0(case_study, "_precision_ecp_table_",
                     target_size_label(n_control, n_treatment, "file"))
  filename <- format_filename(
    filename, case_study,
    target_to_source_std_ratio = unique(rows$target_to_source_std_ratio),
    source_denominator_change_factor = unique(rows$source_denominator_change_factor)
  )
  file_path <- file.path(directory, filename)
  title <- paste0(
    "Belimumab: mean half-width of the 95% credible interval (precision), ",
    "empirical coverage probability (ECP) and interval score, ",
    target_size_label(n_control, n_treatment, "latex")
  )
  export_table(table_data, ncollapses = 0, title = title, file_path = file_path)
  return(invisible(file_path))
}
