# Rows whose target treatment effect lies in the alternative space.
alternative_space_rows <- function(results_freq_df) {
  return(subset(results_freq_df,
         (null_space == "left" & target_treatment_effect > theta_0) |
           (null_space == "right" & target_treatment_effect < theta_0)))
}

analyze_power_gains <- function(results_freq_df, output_path){
  check_required_colnames(
    results_freq_df,
    c(  "null_space",
    "target_treatment_effect",
    "theta_0",
    "success_proba",
    "frequentist_power_at_equivalent_tie",
      "method",
      "conf_int_success_proba_lower",
      "nominal_frequentist_power_separate"),
    context = "analyze_power_gains"
  )

  # Check whether there are cases where the success proba is higher in the alternative space than the frequentist power at equivalent TIE
  subdf <- alternative_space_rows(results_freq_df)

  # Determine whether the power of the method of interest is higher than the
  # power of a separate analysis at equivalent TIE. The decision rests on an
  # interval for the difference that recovers a variance from each of the two
  # intervals, rather than on standard errors reconstructed by assuming those
  # intervals were symmetric and normal.
  subdf <- flag_power_differences(subdf)

  power_gain_cases <- subset(subdf, power_gain)


  # Determine whether the power of the method of interest is higher than the power of the frequentist method at the nominal TIE (otherwise there is no point in borrowing)
  power_gain_cases <- subset(power_gain_cases, conf_int_success_proba_lower > nominal_frequentist_power_separate)


  # Detect anomalies where a separate analysis leads to power gains
  separate_analysis_gain <- subset(power_gain_cases, method == "separate")

  write.csv(file = paste0(output_path, "/", "power_gain_cases.csv"), power_gain_cases)#, power_gain_cases[,expected_colnames_scenario])

  return(write.csv(file = paste0(output_path, "/", "power_gain_cases_separate_analysis.csv"), separate_analysis_gain))

}

analyze_power_loss <- function(results_freq_df, output_path){
  check_required_colnames(
    results_freq_df,
    c(  "null_space",
    "target_treatment_effect",
    "theta_0",
    "success_proba",
    "frequentist_power_at_equivalent_tie"),
    context = "analyze_power_loss"
  )

  # Check whether there are cases where the success proba is lower in the alternative space than the frequentist power at equivalent TIE while TIE is inflated compared to the nominal TIE.
  subdf <- alternative_space_rows(results_freq_df)

  # The same difference interval decides losses as decides gains, so the two
  # analyses no longer apply different levels of stringency to the same
  # comparison.
  subdf <- flag_power_differences(subdf)
  power_loss_cases <- subset(subdf, power_loss)


  return(write.csv(file = paste0(output_path, "/","power_loss_cases.csv"), power_loss_cases[, intersect(expected_colnames_scenario, names(power_loss_cases))]))
}

# The results carry no record of the critical value they were simulated at, so
# the nominal type I error rate has to be passed in, from the same
# analysis_config$nominal_tie the figures use, rather than assumed.
check_nominal_tie <- function(nominal_tie) {
  if (!is.numeric(nominal_tie) || length(nominal_tie) != 1 ||
      !is.finite(nominal_tie) || nominal_tie <= 0 || nominal_tie >= 1) {
    stop("nominal_tie must be a single number strictly between 0 and 1.")
  }
  return(invisible(NULL))
}

analyze_power_loss_inflated_tie <- function(results_freq_df, output_path, nominal_tie){
  check_nominal_tie(nominal_tie)
  check_required_colnames(
    results_freq_df,
    c(  "null_space",
    "target_treatment_effect",
    "theta_0",
    "success_proba",
    "frequentist_power_at_equivalent_tie",
      "conf_int_tie_lower"),
    context = "analyze_power_loss_inflated_tie"
  )

  # Check whether there are cases where the success proba is lower in the alternative space than the frequentist power at equivalent TIE while TIE is inflated compared to the nominal TIE.
  subdf <- alternative_space_rows(results_freq_df)

  subdf <- subset(subdf, conf_int_tie_lower > nominal_tie)

  subdf <- flag_power_differences(subdf)
  power_loss_cases <- subset(subdf, power_loss)


  return(write.csv(file = paste0(output_path, "/", "power_loss_inflated_tie_cases.csv"), power_loss_cases[, intersect(expected_colnames_scenario, names(power_loss_cases))]))

}

analyze_noninflated_tie <- function(results_freq_df, output_path, nominal_tie){
  check_nominal_tie(nominal_tie)
  check_required_colnames(
    results_freq_df,
    c("conf_int_tie_upper",
      "target_treatment_effect",
      "theta_0"),
    context = "analyze_noninflated_tie"
  )

  # Check whether there are cases where the TIE is not inflated

  non_inflated_tie_cases <- subset(results_freq_df, conf_int_tie_upper < nominal_tie)

  # Scenarios whose target effect sits on the boundary of the null hypothesis
  # space. The target effect is drift + source estimate, so this is the drift
  # of -source estimate only when theta_0 is 0, and the sum is not exactly
  # theta_0 in floating point either way.
  at_null_boundary <- abs(
    non_inflated_tie_cases$target_treatment_effect - non_inflated_tie_cases$theta_0
  ) <= sqrt(.Machine$double.eps) * pmax(1, abs(non_inflated_tie_cases$theta_0))
  non_inflated_tie_cases <- non_inflated_tie_cases[at_null_boundary, , drop = FALSE]

  # Remove the cases where there was no borrowing
  # filter <- non_inflated_tie_cases$conf_int_ess_moment_upper >= 0 & non_inflated_tie_cases$conf_int_ess_moment_lower >= 0
  # non_inflated_tie_cases <- non_inflated_tie_cases[!filter,]


  return(write.csv(file = paste0(output_path, "/", "noninflated_tie_cases.csv"), non_inflated_tie_cases))

}
