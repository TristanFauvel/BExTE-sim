#' Recover the number of Monte Carlo replicates behind a binomial estimate
#'
#' @description Operating characteristics are stored as summaries (estimate,
#'   Monte Carlo standard error, exact confidence interval) rather than as the
#'   underlying counts. Propagating their uncertainty correctly requires the
#'   replicate count, which this function recovers from whichever summary is
#'   informative.
#'
#'   The Monte Carlo standard error identifies the replicate count exactly
#'   whenever it is available and non-zero. It is zero precisely when no
#'   replicate succeeded or all of them did, and in that case one bound of the
#'   exact interval is a closed-form function of the replicate count alone. As a
#'   last resort, when no standard error was stored for an interior estimate,
#'   the count is approximated from the width of the interval.
#'
#' @param estimate The estimated proportion.
#' @param mcse The Monte Carlo standard error of the estimate, or `NA`.
#' @param conf_int_lower Lower bound of the exact (Clopper-Pearson) interval.
#' @param conf_int_upper Upper bound of the exact (Clopper-Pearson) interval.
#' @param conf_level The confidence level the interval was computed at.
#'
#' @return The number of replicates, or `NA_real_` when no summary identifies it.
#'
#' @keywords internal
binomial_replicate_count <- function(estimate,
                                     mcse,
                                     conf_int_lower,
                                     conf_int_upper,
                                     conf_level = 0.95) {
  if (is.na(estimate)) {
    return(NA_real_)
  }

  if (!is.na(mcse) && mcse > 0 && estimate > 0 && estimate < 1) {
    return(round(estimate * (1 - estimate) / mcse^2))
  }

  if (is.na(conf_int_lower) || is.na(conf_int_upper)) {
    return(NA_real_)
  }

  tail_probability <- (1 - conf_level) / 2

  # With no observed success the exact interval is [0, 1 - tail^(1/n)], and with
  # no observed failure it is [tail^(1/n), 1]; either bound inverts to n.
  if (estimate <= 0) {
    if (conf_int_upper <= 0 || conf_int_upper >= 1) {
      return(NA_real_)
    }
    return(round(log(tail_probability) / log1p(-conf_int_upper)))
  }

  if (estimate >= 1) {
    if (conf_int_lower <= 0 || conf_int_lower >= 1) {
      return(NA_real_)
    }
    return(round(log(tail_probability) / log(conf_int_lower)))
  }

  interval_width <- conf_int_upper - conf_int_lower
  if (interval_width <= 0) {
    return(NA_real_)
  }

  # Approximate fallback for results stored without a Monte Carlo standard
  # error: the replicate count a Wald interval of this width would imply.
  quantile_value <- stats::qnorm(1 - tail_probability)
  return(round(estimate * (1 - estimate) * (2 * quantile_value / interval_width)^2))
}


#' Draw from the posterior of a binomial proportion
#'
#' @description Draws from the Jeffreys posterior `Beta(x + 1/2, n - x + 1/2)`
#'   implied by observing `estimate * n_replicates` successes in `n_replicates`
#'   replicates.
#'
#'   This replaces approximating an exact binomial interval by a symmetric
#'   normal and discarding the draws that fall outside `(0, 1)`. The posterior
#'   is supported on `(0, 1)` by construction, so no draw is ever discarded, and
#'   it retains the skewness that the normal approximation removes -- which
#'   matters most for the small type I errors these comparisons are made at.
#'
#' @param n_samples Number of draws to return.
#' @param estimate The estimated proportion.
#' @param n_replicates The number of replicates the estimate is based on.
#'
#' @return A numeric vector of `n_samples` draws in `(0, 1)`.
#'
#' @keywords internal
sample_binomial_proportion <- function(n_samples, estimate, n_replicates) {
  if (is.na(estimate) || is.na(n_replicates) || n_replicates <= 0) {
    return(rep(NA_real_, n_samples))
  }

  successes <- estimate * n_replicates
  return(stats::rbeta(
    n_samples,
    shape1 = successes + 0.5,
    shape2 = n_replicates - successes + 0.5
  ))
}


#' Confidence interval for a difference, by variance estimate recovery
#'
#' @description Combines two separate intervals into an interval for the
#'   difference of the two estimates, following Newcombe's MOVER (Method of
#'   Variance Estimates Recovery) approach.
#'
#'   Subtracting only the comparator's point estimate leaves an interval exactly
#'   as wide as the first estimate's, which understates the uncertainty of the
#'   difference. MOVER recovers a variance from each interval at the bound that
#'   matters for the corresponding bound of the difference, so the result keeps
#'   the asymmetry of exact binomial intervals instead of symmetrising them.
#'
#' @param estimate_1,lower_1,upper_1 Estimate and interval bounds of the first
#'   quantity.
#' @param estimate_2,lower_2,upper_2 Estimate and interval bounds of the second
#'   quantity, which is subtracted from the first.
#' @param correlation Correlation between the two estimates. Defaults to `0`,
#'   which assumes them independent; a positive value narrows the interval.
#'
#' @return A named numeric vector with the `lower` and `upper` bounds of the
#'   interval for `estimate_1 - estimate_2`, or `NA` bounds when any input is
#'   missing.
#'
#' @keywords internal
mover_difference_ci <- function(estimate_1,
                                lower_1,
                                upper_1,
                                estimate_2,
                                lower_2,
                                upper_2,
                                correlation = 0) {
  inputs <- c(estimate_1, lower_1, upper_1, estimate_2, lower_2, upper_2, correlation)
  if (anyNA(inputs)) {
    return(c(lower = NA_real_, upper = NA_real_))
  }

  difference <- estimate_1 - estimate_2

  # The lower bound of the difference pairs the lower arm of the first interval
  # with the upper arm of the second, and vice versa for the upper bound.
  lower_arm_1 <- estimate_1 - lower_1
  upper_arm_2 <- upper_2 - estimate_2
  upper_arm_1 <- upper_1 - estimate_1
  lower_arm_2 <- estimate_2 - lower_2

  lower_variance <- lower_arm_1^2 + upper_arm_2^2 -
    2 * correlation * lower_arm_1 * upper_arm_2
  upper_variance <- upper_arm_1^2 + lower_arm_2^2 -
    2 * correlation * upper_arm_1 * lower_arm_2

  return(c(
    lower = difference - sqrt(max(0, lower_variance)),
    upper = difference + sqrt(max(0, upper_variance))
  ))
}


#' Interval for the difference between success probability and comparator power
#'
#' @description Applies [mover_difference_ci()] row by row to the success
#'   probability and the frequentist power at equivalent type I error.
#'
#' @param results_df A data frame carrying `success_proba`, its interval bounds,
#'   `frequentist_power_at_equivalent_tie` and its interval bounds.
#' @param correlation Correlation between the two estimates, passed to
#'   [mover_difference_ci()].
#'
#' @return A matrix with one row per input row and columns `lower` and `upper`.
#'
#' @keywords internal
power_difference_bounds <- function(results_df, correlation = 0) {
  if (nrow(results_df) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 2,
                  dimnames = list(NULL, c("lower", "upper"))))
  }

  bounds <- mapply(
    mover_difference_ci,
    estimate_1 = results_df$success_proba,
    lower_1 = results_df$conf_int_success_proba_lower,
    upper_1 = results_df$conf_int_success_proba_upper,
    estimate_2 = results_df$frequentist_power_at_equivalent_tie,
    lower_2 = results_df$frequentist_power_at_equivalent_tie_lower,
    upper_2 = results_df$frequentist_power_at_equivalent_tie_upper,
    MoreArgs = list(correlation = correlation)
  )

  return(t(bounds))
}


#' Replace success probabilities by power differences
#'
#' @description Turns the success probability columns into the difference
#'   against the frequentist power at equivalent type I error, with an interval
#'   that accounts for the uncertainty of both estimates rather than shifting
#'   the success probability interval by the comparator's point estimate.
#'
#' @param results_df A data frame of operating characteristics.
#' @param correlation Correlation between the two estimates. The default of `0`
#'   treats them as independent, which is the conservative choice: the two are
#'   estimated from simulations sharing a seed, and any positive correlation
#'   between them would only narrow the interval.
#'
#'   Note separately that the comparator is a deterministic function of a single
#'   type I error estimate shared by every row of a panel, so the comparator
#'   errors are perfectly correlated *across* rows whatever this within-row
#'   correlation is; points within a panel must not be read as independent.
#'
#' @return The data frame with `success_proba` and its interval bounds replaced
#'   by the difference and its interval.
#'
#' @keywords internal
add_power_difference_columns <- function(results_df, correlation = 0) {
  bounds <- power_difference_bounds(results_df, correlation = correlation)

  results_df$success_proba <- results_df$success_proba -
    results_df$frequentist_power_at_equivalent_tie
  results_df$conf_int_success_proba_lower <- bounds[, "lower"]
  results_df$conf_int_success_proba_upper <- bounds[, "upper"]

  return(results_df)
}


#' Flag power gains and power losses against the comparator
#'
#' @description Adds the power difference, its interval, and the gain and loss
#'   flags derived from that single interval. Deciding both from the same
#'   interval keeps gains and losses at the same stringency; they were
#'   previously detected with two different rules -- a z-test built from
#'   interval widths on one side, and non-overlap of two intervals on the other.
#'
#' @param results_df A data frame of operating characteristics.
#' @param correlation Correlation between the two estimates, passed to
#'   [mover_difference_ci()].
#'
#' @return The data frame with `power_difference`, `power_difference_lower`,
#'   `power_difference_upper`, `power_gain` and `power_loss` columns added.
#'
#' @keywords internal
flag_power_differences <- function(results_df, correlation = 0) {
  bounds <- power_difference_bounds(results_df, correlation = correlation)

  results_df$power_difference <- results_df$success_proba -
    results_df$frequentist_power_at_equivalent_tie
  results_df$power_difference_lower <- bounds[, "lower"]
  results_df$power_difference_upper <- bounds[, "upper"]

  results_df$power_gain <- !is.na(results_df$power_difference_lower) &
    results_df$power_difference_lower > 0
  results_df$power_loss <- !is.na(results_df$power_difference_upper) &
    results_df$power_difference_upper < 0

  return(results_df)
}
