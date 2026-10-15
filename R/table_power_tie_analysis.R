# Rename the scenario columns to table headers: "Case study", "$n_C$", ...
title_case_column_names <- function(columns) {
  columns <- gsub("_", " ", columns)                       # Replace underscores with spaces
  columns <- sub("^(\\w)", "\\U\\1", columns, perl = TRUE) # Capitalize only the first letter of the string
  columns
}

# Export a table of scenarios in which a method gains or loses power or TIE.
export_scenario_cases_table <- function(cases_df, file_stem, title, drop_na_case_study = FALSE) {
  check_table_packages()
  cases_df <- cases_df[, intersect(c(unique_scenario_columns, "method", "parameters"), names(cases_df))]

  file_path <- file.path(tables_dir, file_stem)

  cases_df <- cases_df[order(cases_df$target_sample_size_control + cases_df$target_sample_size_treatment,
                             cases_df$target_sample_size_control), ]

  if (drop_na_case_study) {
    cases_df <- cases_df %>%
      filter(!is.na(case_study))
  }

  cases_df$method <- format_results_df_parameters(cases_df)
  cases_df <- cases_df[, !names(cases_df) %in% "parameters"]
  colnames(cases_df) <- title_case_column_names(colnames(cases_df))
  cases_df <- cases_df[, c("Method", "Case study", "Target sample size control", "Target sample size treatment", "Source denominator change factor", "Target to source std ratio")]
  names(cases_df)[names(cases_df) == "Target to source std ratio"] <- "$\\sigma_T/\\sigma_S$"
  names(cases_df)[names(cases_df) == "Target sample size control"] <- "$n_C$"
  names(cases_df)[names(cases_df) == "Target sample size treatment"] <- "$n_T$"

  export_table(data_table = cases_df, title = title, file_path = file_path)
}

table_noninflated_tie <- function(noninflated_cases_df){
  export_scenario_cases_table(noninflated_cases_df, "noninflated_tie_cases",
                              'Cases in which no inflation of the TIE was observed')
}

table_power_gain <- function(power_gain_cases_df){
  export_scenario_cases_table(power_gain_cases_df, "power_gain_cases",
                              'Cases in which power gains were observed at equivalent TIE compared to the frequentist method')
}

table_power_loss <- function(power_loss_cases_df){
  export_scenario_cases_table(power_loss_cases_df, "power_loss_cases",
                              'Cases in which power losss were observed at equivalent TIE compared to the frequentist method')
}

table_power_loss_cases_inflated_tie <- function(power_loss_cases_inflated_tie_df){
  export_scenario_cases_table(power_loss_cases_inflated_tie_df, "power_loss_cases_inflated_tie",
                              'Cases in which power losss were observed at equivalent TIE compared to the frequentist method, for a TIE higher than the nominal TIE',
                              drop_na_case_study = TRUE)
}
