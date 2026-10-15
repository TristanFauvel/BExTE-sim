## Names that R CMD check cannot resolve, because they are not R objects of
## the package:
## - column names in dplyr verbs and ggplot2::aes(), and the `.` of the
##   magrittr pipe;
## - the style and path objects (text_size, font, figures_dir, ...) that the
##   plot and table code sources from inst/conf into the global environment;
## - the loop variables of foreach (task, job).
utils::globalVariables(c(
  ".", "analysis_config", "bayesian_metrics", "case_studies_config_dir",
  "case_study", "categories", "CI", "CI_low", "CI_power", "CI_tie",
  "CI_upper", "color", "combined_mean_CI", "conf_int_lower",
  "conf_int_success_proba_lower", "conf_int_success_proba_upper",
  "conf_int_tie_lower", "conf_int_tie_upper", "conf_int_upper",
  "control_drift", "design_prior_type", "dpi", "drift", "drift_value",
  "dropout_probability", "effect", "effect_label", "effect_name",
  "effect_type", "empirical_bayes_hyperparameters", "error_low",
  "error_low_power", "error_low_tie", "error_upper", "error_upper_power",
  "error_upper_tie", "ess_moment", "event_time_distribution", "figures_dir",
  "font", "frequentist_metrics", "hyperparameter", "inference_metrics",
  "jitter_offset", "jittered_x", "job", "label", "markersize",
  "mcse_success_proba", "means", "Method", "method", "methods_dict",
  "methods_labels", "methods_style", "nominal_frequentist_power_separate",
  "null_space", "parameters", "Parameters", "plots_to_latex", "power_gain",
  "power_loss", "relative_error_cap_width", "remake_figures", "rows",
  "small_text_size", "source_denominator_change_factor",
  "source_treatment_effect_estimate", "success_proba", "tables_dir",
  "target_sample_size", "target_sample_size_control",
  "target_sample_size_total", "target_sample_size_treatment", "target_size",
  "target_to_source_std_ratio", "target_treatment_effect", "task",
  "text_size", "textwidth", "theta_0", "tie", "total_width",
  "treatment_delay", "treatment_drift", "x", "xvars", "y", "ymax", "ymin"
))
