#' Whether power can be computed analytically for this target data
#'
#' @description The analytical power formulas apply to a normal summary measure
#'   on a continuous endpoint; the remaining endpoints have to be simulated.
#'   Sharing the predicate keeps the power computation and the propagation of
#'   its uncertainty on the same branch.
#'
#'   Recurrent events are simulated too. Mepolizumab used to be priced in
#'   closed form, but its trials are generated patient by patient from a
#'   negative binomial and the standard error is re-estimated in each one, so
#'   the closed form assumed a test the Bayesian methods were never compared
#'   against.
#'
#' @param target_data Target data object.
#'
#' @return `TRUE` when power has a closed form for this target data.
#'
#' @keywords internal
uses_analytical_power <- function(target_data) {
  return(target_data$endpoint == "normal" ||
    target_data$endpoint == "continuous")
}


#' Simulate the p-values of the frequentist test
#'
#' @description Generates `n_replicates` trials and returns the p-value of the
#'   test in each one. Returning the p-values rather than the decisions lets the
#'   power be read off at any significance level without re-simulating.
#'
#' @param target_data Target data object.
#' @param frequentist_test Type of frequentist test to apply.
#' @param theta_0 Boundary of the null hypothesis space.
#' @param alternative Direction of the alternative hypothesis.
#' @param simulation_config Simulation configuration.
#' @param n_replicates Number of trials to simulate.
#'
#' @return A numeric vector of `n_replicates` p-values.
#'
#' @keywords internal
simulate_test_p_values <- function(target_data,
                                   frequentist_test,
                                   theta_0,
                                   alternative,
                                   simulation_config,
                                   n_replicates) {
  if (frequentist_test != "t-test") {
    stop("Only implemented for a t-test.")
  }

  return(trial_p_values(
    simulated_trials(target_data, simulation_config, n_replicates),
    theta_0 = theta_0,
    alternative = alternative
  ))
}

#' The trials a simulated frequentist power is read from
#'
#' @description Seeds the generator with `simulation_config$seed` and
#'   generates `n_replicates` trials, keeping the estimable ones. Every
#'   simulated power of a design - the separate analysis's at the equivalent
#'   and at the nominal type I error, and the pooled analysis's - reads these
#'   same trials, so they are generated once per design and shared. Only the
#'   three summaries the tests read are kept, which is what makes holding one
#'   set per design affordable.
#'
#' @param target_data Target data object.
#' @param simulation_config Simulation configuration.
#' @param n_replicates Number of trials to simulate.
#'
#' @return A data frame with `treatment_effect_estimate`,
#'   `standard_deviation`, `sample_size_control` and `sample_size_treatment`,
#'   one row per estimable trial.
#' @noRd
simulated_trials <- function(target_data, simulation_config, n_replicates) {
  withr::local_seed(simulation_config$seed)

  samples <- estimable_replicates(target_data$generate(n_replicates))
  return(data.frame(
    treatment_effect_estimate = samples$treatment_effect_estimate,
    standard_deviation = samples$standard_deviation,
    sample_size_control = samples$sample_size_control,
    sample_size_treatment = samples$sample_size_treatment
  ))
}

#' The separate analysis's t-test p-value in each trial
#'
#' @param trials Trials from [simulated_trials()].
#' @param theta_0 Boundary of the null hypothesis space.
#' @param alternative Direction of the alternative hypothesis.
#'
#' @return A numeric vector of p-values, one per trial.
#' @noRd
trial_p_values <- function(trials, theta_0, alternative) {
  return(one_sample_t_test_p_values(
    estimate = trials$treatment_effect_estimate,
    mu = theta_0,
    standard_deviation = trials$standard_deviation,
    n = target_effective_arm_size(trials),
    alternative = alternative
  ))
}

#' P-values of a one-sample t-test from summary statistics, for many trials
#'
#' @description The p-value `BSDA::tsum.test()` returns for a one-sample test,
#'   computed for every replicate at once. Calling `tsum.test()` once per
#'   replicate cost ~50us each in argument checking and the warning it raises
#'   about `var.equal`, which on 10,000 replicates was most of the test. The
#'   statistic is built with the same operations as `tsum.test()`, so the
#'   p-values agree with it to the last bit.
#'
#' @param estimate Sample means, one per replicate.
#' @param mu Value of the mean under the null hypothesis.
#' @param standard_deviation Sample standard deviations.
#' @param n Sample sizes.
#' @param alternative "greater", "less" or "two.sided".
#'
#' @return A numeric vector of p-values, one per replicate.
#' @noRd
one_sample_t_test_p_values <- function(estimate, mu, standard_deviation, n, alternative) {
  statistic <- (estimate - mu) / sqrt(standard_deviation^2 / n)
  degrees_of_freedom <- n - 1

  switch(
    alternative,
    greater = return(1 - stats::pt(statistic, degrees_of_freedom)),
    less = return(stats::pt(statistic, degrees_of_freedom)),
    two.sided = return(2 * stats::pt(-abs(statistic), degrees_of_freedom)),
    stop("alternative must be one of \"greater\", \"less\" or \"two.sided\".")
  )
}
#' Direction of the alternative hypothesis implied by the null space
#'
#' @param null_space Side of the null space, either left or right.
#'
#' @return "greater" or "less".
#' @noRd
alternative_from_null_space <- function(null_space) {
  if (null_space == "left") {
    return("greater")
  } else if (null_space == "right") {
    return("less")
  } else {
    stop("Null space must be either 'left' or 'right'")
  }
}

#' Closed-form power of the frequentist test
#'
#' Vectorised over `alpha`: every pwr function used here accepts a vector of
#' significance levels and returns one power per level, which is what lets
#' compute_power_with_tie_ci() price 1000 sampled type I errors in a single
#' call instead of 1000 (the per-call overhead of pwr - match.arg, argument
#' assembly, simplify2array - costs more than the power calculation itself).
#'
#' @param alpha Significance level(s).
#' @param target_data Target data object.
#' @param frequentist_test Type of frequentist test to apply.
#' @param theta_0 Boundary of the null hypothesis space.
#' @param alternative Direction of the alternative hypothesis.
#'
#' @return A numeric vector of powers, one per element of `alpha`.
#' @noRd
analytical_power <- function(alpha, target_data, frequentist_test, theta_0, alternative) {
  if (target_data$summary_measure_likelihood == "normal") {
    effect_size <- (target_data$treatment_effect - theta_0) / target_data$standard_deviation

    if (frequentist_test == "t-test") {
      return(pwr::pwr.t.test(
        d = effect_size,
        n = target_effective_arm_size(target_data),
        sig.level = alpha,
        type = "one.sample",
        alternative = alternative
      )$power)
    } else if (frequentist_test == "z-test") {
      return(pwr::pwr.norm.test(
        d = effect_size,
        n = target_effective_arm_size(target_data),
        sig.level = alpha,
        alternative = alternative
      )$power)
    } else {
      stop("Unsupported test type.")
    }
  } else if (target_data$summary_measure_likelihood == "binomial") {
    h <- pwr::ES.h(target_data$treatment_rate, target_data$control_rate)

    return(pwr::pwr.2p2n.test(
      h = h,
      n1 = target_data$sample_size_treatment,
      n2 = target_data$sample_size_control,
      sig.level = alpha,
      alternative = alternative
    )$power)
  } else {
    stop("Unsupported likelihood type.")
  }
}

#' Check the target data fields the power computations read
#'
#' Hoisted out of compute_freq_power(): compute_power_with_tie_ci()
#' evaluates power at 1000 sampled alphas per result row and target_data is
#' identical across all of them, so asserting inside that call ran the same
#' three checks 3000 times a row. assert_single_number() costs ~190us a
#' call (it deparses, matches the call and dispatches over a list of
#' assertion functions), which made argument checking ~87% of the analysis.
#'
#' @param target_data Target study data, as returned by load_data().
#'
#' @return `NULL`, invisibly. Called for the error it raises.
#' @noRd
assert_target_data_numbers <- function(target_data) {
  assert_single_number(target_data$treatment_effect)
  assert_single_number(target_data$standard_deviation)
  assert_single_number(target_data$sample_size_control)
  assert_single_number(target_data$sample_size_treatment)

  return(invisible(NULL))
}



#' Compute the frequentist power
#'
#' @description This function computes the power of a test for a given significance level
#' @description This function computes the power of a test for a given significance level
#'
#' @param alpha The significance level.
#' @param target_data The target data containing sample size, treatment effect, standard deviation, and summary measure distribution.
#' @param frequentist_test The type of frequentist test ("t-test" or "z-test").
#' @param theta_0 The null hypothesis value.
#' @param target_data Target data object
#' @param frequentist_test Type of frequentist test to apply, either z-test or t-test
#' @param theta_0 Boundary of the null hypothesis space
#' @param null_space Side of the null space, either left or right.
#' @param simulation_config Simulation configuration.
#' @param case_study Optional case-study name.
#' @param n_replicates Number of Monte Carlo replicates for non-analytical power calculations.
#' @param p_values Optional p-values already simulated for this design by
#'   [simulate_test_p_values()], with the same seed and `n_replicates`. `NULL`
#'   simulates them. Ignored when the power has a closed form.
#' @param null_p_values Optional p-values of the same test on the trials of the
#'   design's null scenario. When given, the test rejects at the threshold its
#'   actual type I error there equals `alpha` at, rather than at `alpha`
#'   itself - see [calibrated_levels()]. Ignored when the power has a closed
#'   form.
#'
#' @return A list containing the power and its confidence interval.
#'
#' @keywords internal
compute_freq_power <- function(alpha,
                               target_data,
                               frequentist_test,
                               theta_0,
                               null_space,
                               simulation_config,
                               case_study = NULL,
                               n_replicates = 1000,
                               p_values = NULL,
                               null_p_values = NULL) {
  alternative <- alternative_from_null_space(null_space)

  if (is.na(alpha)){
    return(list(
      power = NA_real_,
      conf_int_power = rep(NA_real_, 2)
    ))
  }

  power <- NA # Default value in case of an unsupported distribution

  if (target_data$summary_measure_likelihood == "normal") {
    if (uses_analytical_power(target_data)) {
      power <- analytical_power(alpha, target_data, frequentist_test, theta_0, alternative)

      conf_int_power <- c(power, power)
    } else {
      # In this case, we cannot use an analytical computation of power
      if (is.null(p_values)) {
        p_values <- simulate_test_p_values(
          target_data = target_data,
          frequentist_test = frequentist_test,
          theta_0 = theta_0,
          alternative = alternative,
          simulation_config = simulation_config,
          n_replicates = n_replicates
        )
      }

      threshold <- if (is.null(null_p_values)) alpha else calibrated_levels(alpha, null_p_values)
      test_decisions <- p_values < threshold
      power <- mean(test_decisions)
      conf_int_power <- binom.test(sum(test_decisions), length(test_decisions), conf.level = 0.95)$conf.int
    }
  } else if (target_data$summary_measure_likelihood == "binomial") {
    power <- analytical_power(alpha, target_data, frequentist_test, theta_0, alternative)

    conf_int_power <- c(power, power)
  } else {
    stop("Unsupported likelihood type.")
  }

  return(list(power = power, conf_int_power = conf_int_power))
}

#' Compute the frequentist power
#'
#' @description This function computes the power of a test for a pooled analysis at a given significance level
#'
#' @param alpha The significance level.
#' @param target_data The target data containing sample size, treatment effect, standard deviation, and summary measure distribution.
#' @param source_data Source study data
#' @param frequentist_test The type of frequentist test ("t-test" or "z-test").
#' @param theta_0 The null hypothesis value.
#' @param target_data Target data object
#' @param frequentist_test Type of frequentist test to apply, either z-test or t-test
#' @param theta_0 Boundary of the null hypothesis space
#' @param null_space Side of the null space, either left or right.
#' @param simulation_config Simulation configuration.
#' @param case_study Optional case-study name.
#' @param n_replicates Number of Monte Carlo replicates for non-analytical power calculations.
#' @param trials Optional trials already simulated for this design, with the
#'   same seed and `n_replicates`, as the separate analysis's power reads them.
#'   `NULL` simulates them. Ignored when the power has a closed form.
#'
#' @return The power of the test.
#'
#' @keywords internal
compute_freq_power_pooling <- function(alpha,
                                       target_data,
                                       source_data,
                                       frequentist_test,
                                       theta_0,
                                       null_space,
                                       simulation_config,
                                       case_study =  NULL,
                                       n_replicates = 1000,
                                       trials = NULL) {

  if (null_space == "left") {
    alternative <- "greater"
  } else if (null_space == "right") {
    alternative <- "less"
  } else {
    stop("Null space must be either 'left' or 'right'")
  }

  assert_single_number(target_data$treatment_effect)
  assert_single_number(target_data$standard_deviation)
  assert_single_number(target_effective_arm_size(target_data))

  power <- NA # Default value in case of an unsupported distribution

  if (target_data$summary_measure_likelihood == "normal") {
    if (uses_analytical_power(target_data)) {
      target_treatment_effect_standard_error <- target_data$standard_deviation / sqrt(target_effective_arm_size(target_data))

      pooled_treatment_effect <- (
        source_data$treatment_effect_estimate / (
          source_data$standard_error ^ 2 / target_treatment_effect_standard_error ^
            2 + 1
        )
      ) + (
        target_data$treatment_effect / (
          1 + target_treatment_effect_standard_error ^ 2 / source_data$standard_error ^
            2
        )
      )


      pooled_standard_error_2 <- 1 / (1 / source_data$standard_error ^ 2 + 1 / target_treatment_effect_standard_error ^
                                        2)

      # The source estimate is the one observed, fixed across the simulated
      # trials, and only the target estimate varies: the pooled estimate is
      # w_S theta_S + w_T theta_T_hat, with theta_T_hat ~ N(theta_T, SE_T^2).
      # Its spread across trials is therefore w_T SE_T, not the pooled
      # standard error, which also counts the source estimate's sampling
      # variability. The test compares (pooled - theta_0) / SE_pooled with the
      # critical value of the test the simulated branch below applies, so the
      # power is the probability that theta_T_hat clears the matching bound.
      target_weight <- source_data$standard_error ^ 2 /
        (source_data$standard_error ^ 2 + target_treatment_effect_standard_error ^ 2)
      source_part <- pooled_treatment_effect - target_weight * target_data$treatment_effect
      pooled_sample_size <- target_effective_arm_size(target_data) + source_data$equivalent_source_sample_size_per_arm

      critical <- if (frequentist_test == "t-test") {
        stats::qt(1 - alpha, df = pooled_sample_size - 1)
      } else if (frequentist_test == "z-test") {
        stats::qnorm(1 - alpha)
      } else {
        stop("Only implemented for a t-test or z-test.")
      }

      bound <- if (alternative == "greater") {
        theta_0 + critical * sqrt(pooled_standard_error_2)
      } else {
        theta_0 - critical * sqrt(pooled_standard_error_2)
      }
      target_bound <- (bound - source_part) / target_weight
      power <- stats::pnorm(
        (target_bound - target_data$treatment_effect) / target_treatment_effect_standard_error,
        lower.tail = alternative == "less"
      )
      conf_int_power <- c(power, power)
    } else {
      # In this case, we cannot use an analytical computation of power. These
      # are the same trials the separate analysis's power is read from.
      target_data_samples <- if (is.null(trials)) {
        simulated_trials(target_data, simulation_config, n_replicates)
      } else {
        trials
      }

      if (frequentist_test != "t-test") {
        stop("Only implemented for a t-test.")
      }

      # Every replicate at once: the arithmetic is elementwise, so each entry
      # is exactly what the former per-replicate loop computed.
      target_treatment_effect_standard_error <- target_data_samples$standard_deviation / sqrt(target_effective_arm_size(target_data_samples))

      pooled_treatment_effect <- (
        source_data$treatment_effect_estimate / (
          source_data$standard_error ^ 2 / target_treatment_effect_standard_error ^
            2 + 1
        )
      ) + (
        target_data_samples$treatment_effect_estimate / (
          1 + target_treatment_effect_standard_error ^ 2 / source_data$standard_error ^
            2
        )
      )

      pooled_standard_error_2 <- 1 / (1 / source_data$standard_error ^ 2 + 1 / target_treatment_effect_standard_error ^
                                        2)

      pooled_sample_size <- target_effective_arm_size(target_data) + source_data$equivalent_source_sample_size_per_arm

      pooled_variance <- pooled_standard_error_2 * pooled_sample_size

      p_values <- one_sample_t_test_p_values(
        estimate = pooled_treatment_effect,
        mu = theta_0,
        standard_deviation = sqrt(pooled_variance),
        n = pooled_sample_size,
        alternative = alternative
      )
      test_decisions <- as.numeric(p_values < alpha)

      power <- mean(test_decisions)
      conf_int_power <- binom.test(sum(test_decisions), length(test_decisions), conf.level = 0.95)$conf.int
    }
  } else if (target_data$summary_measure_likelihood == "binomial") {
    power <- pooled_binomial_power(
      alpha = alpha,
      target_data = target_data,
      source_data = source_data,
      alternative = alternative
    )
    conf_int_power <- c(power, power)
  } else {
    stop("This likelihood is not supported.")
  }

  assert_single_number(power)
  return(list(power = power, conf_int_power = conf_int_power))
}


#' Draw plausible values of the equivalent type I error
#'
#' @description Draws from the posterior of the type I error implied by the
#'   replicates it was estimated from, rather than from a normal centred on the
#'   estimate whose spread is read off the width of an exact interval.
#'
#' @param alpha A list with the type I error estimate (`mean`), its exact
#'   interval bounds and, when available, its Monte Carlo standard error
#'   (`mcse`) or replicate count (`n_replicates`).
#' @param n_samples Number of draws to return.
#'
#' @return A numeric vector of `n_samples` draws, or `NA` when the replicate
#'   count behind the estimate cannot be recovered.
#'
#' @keywords internal
sample_equivalent_tie <- function(alpha, n_samples) {
  if (alpha$conf_int_upper <= alpha$conf_int_lower) {
    # The type I error carries no uncertainty of its own, so every draw sits at
    # the estimate and the simulated power supplies the only spread.
    return(rep(alpha$mean, n_samples))
  }

  n_replicates <- alpha$n_replicates
  if (is.null(n_replicates) || is.na(n_replicates)) {
    n_replicates <- binomial_replicate_count(
      estimate = alpha$mean,
      mcse = if (is.null(alpha$mcse)) NA_real_ else alpha$mcse,
      conf_int_lower = alpha$conf_int_lower,
      conf_int_upper = alpha$conf_int_upper
    )
  }

  return(sample_binomial_proportion(n_samples, alpha$mean, n_replicates))
}


#' Compute the frequentist power at an estimated type I error
#'
#' @description Propagates the uncertainty of the equivalent type I error, and
#'   the Monte Carlo error of the power itself, into an interval for the power a
#'   separate frequentist analysis would reach at that type I error.
#'
#'   The reported bounds are quantiles of the resulting posterior for the power,
#'   not a frequentist confidence interval; they are stored under the existing
#'   `conf_int_power` name for continuity with the columns downstream.
#'
#' @param alpha A list describing the estimated type I error: its `mean`, its
#'   exact interval bounds, and its `mcse` or `n_replicates` when available.
#' @param target_data Target data object.
#' @param frequentist_test Type of frequentist test to apply, either z-test or t-test.
#' @param theta_0 Boundary of the null hypothesis space.
#' @param null_space Side of the null space, either left or right.
#' @param simulation_config Simulation configuration.
#' @param case_study Optional case-study name.
#' @param n_replicates Number of Monte Carlo replicates for non-analytical power.
#' @param n_samples Number of draws of the type I error.
#' @param p_values Optional p-values of the separate analysis, already
#'   simulated for this design by [simulate_test_p_values()]. They depend on the
#'   design alone, not on the borrowing method, so a caller pricing many rows of
#'   one design simulates them once and passes them here. `NULL` simulates them.
#'   Ignored when the power has a closed form.
#' @param null_p_values Optional p-values of the same test on the trials of the
#'   design's null scenario. When given, each sampled type I error is matched
#'   by the test's actual rejection rate there rather than by its nominal
#'   level - see [calibrated_levels()]. Ignored when the power has a closed
#'   form.
#'
#' @return A list with the power, its interval, and the number of draws used.
#'
#' @keywords internal
compute_power_with_tie_ci <- function(alpha,
                                      target_data,
                                      frequentist_test,
                                      theta_0,
                                      null_space,
                                      simulation_config,
                                      case_study = NULL,
                                      n_replicates = 1000,
                                      n_samples = 1000,
                                      p_values = NULL,
                                      null_p_values = NULL) {
  missing_result <- list(
    power = NA_real_,
    conf_int_power = rep(NA_real_, 2),
    n_effective_samples = 0L
  )

  if (is.na(alpha$mean) || is.na(alpha$conf_int_upper) || is.na(alpha$conf_int_lower)) {
    return(missing_result)
  }

  alpha_samples <- sample_equivalent_tie(alpha, n_samples)
  # A degenerate type I error of zero or one carries no information about the
  # power, so the whole estimate is reported as missing rather than some draws
  # being dropped from an otherwise usable sample.
  if (anyNA(alpha_samples) || any(alpha_samples <= 0) || any(alpha_samples >= 1)) {
    return(missing_result)
  }

  if (uses_analytical_power(target_data)) {
    # Once per result row rather than once per sampled alpha - see
    # assert_target_data_numbers().
    assert_target_data_numbers(target_data)

    # Power is a deterministic function of alpha here, so the type I error is
    # the only source of uncertainty. Every sampled level is priced in one
    # vectorised call - see analytical_power().
    alternative <- alternative_from_null_space(null_space)
    power_samples <- analytical_power(
      alpha_samples, target_data, frequentist_test, theta_0, alternative
    )
  } else {
    # Simulate the trials once and read the rejection count off the same
    # p-values at every sampled alpha, then propagate the Monte Carlo error of
    # each count through its own posterior. The reported interval therefore
    # carries both the type I error uncertainty and the replicate noise, which
    # re-simulating under a fixed seed used to suppress entirely.
    if (is.null(p_values)) {
      alternative <- if (null_space == "left") "greater" else "less"
      p_values <- simulate_test_p_values(
        target_data = target_data,
        frequentist_test = frequentist_test,
        theta_0 = theta_0,
        alternative = alternative,
        simulation_config = simulation_config,
        n_replicates = n_replicates
      )
    }
    # The number of p-values strictly below each sampled alpha, read off the
    # sorted p-values by binary search rather than one full pass per alpha.
    if (is.unsorted(p_values)) {
      p_values <- sort(p_values)
    }
    levels <- if (is.null(null_p_values)) {
      alpha_samples
    } else {
      calibrated_levels(alpha_samples, null_p_values)
    }
    rejections <- as.numeric(
      findInterval(levels, p_values, left.open = TRUE)
    )
    power_samples <- stats::rbeta(
      n_samples,
      shape1 = rejections + 0.5,
      shape2 = length(p_values) - rejections + 0.5
    )
  }

  return(list(
    power = mean(power_samples, na.rm = TRUE),
    conf_int_power = unname(
      stats::quantile(power_samples, probs = c(0.025, 0.975), na.rm = TRUE)
    ),
    n_effective_samples = length(alpha_samples)
  ))
}


#' Key under which one design's simulated trials are cached
#'
#' The trials follow from the design, the seed and the number of replicates,
#' so a key naming all three identifies them.
#'
#' @param design_key Keys from [nominal_tie_design_key()].
#' @param simulation_config The simulation configuration, for its seed.
#' @param n_replicates Number of trials simulated.
#'
#' @return A character vector, one key per design key.
#' @noRd
trial_cache_key <- function(design_key, simulation_config, n_replicates) {
  seed <- simulation_config$seed
  return(paste(design_key,
        if (is.null(seed)) "no seed" else format(seed),
        n_replicates, sep = "\003"))
}


#' The cached trials of each design, where the cache holds them
#'
#' @param design_key Keys from [nominal_tie_design_key()].
#' @param simulation_config The simulation configuration.
#' @param n_replicates Number of trials simulated per design.
#' @param trial_cache An environment of cached trials, or `NULL`.
#'
#' @return A list with one entry per design key: its trials, or `NULL`.
#' @noRd
cached_design_trials <- function(design_key, simulation_config, n_replicates, trial_cache) {
  trials <- vector("list", length(design_key))
  if (is.null(trial_cache)) {
    return(trials)
  }

  keys <- trial_cache_key(design_key, simulation_config, n_replicates)
  for (d in which(vapply(keys, exists, logical(1), envir = trial_cache, inherits = FALSE))) {
    trials[[d]] <- get(keys[d], envir = trial_cache)
  }
  return(trials)
}


#' Simulate the trials of one design
#'
#' @param row A results row describing the design.
#' @param simulation_config The simulation configuration.
#' @param n_replicates Number of trials to simulate.
#'
#' @return The design's trials, as [simulated_trials()] returns them.
#' @noRd
design_trials <- function(row, simulation_config, n_replicates) {
  target_data <- load_data(row, type = "target", reload_data_objects = TRUE)
  return(simulated_trials(target_data, simulation_config, n_replicates))
}


#' The separate analysis's p-values for every design that needs them
#'
#' @description Designs whose power has a closed form get `NULL`. The others
#'   read their trials from `trial_cache` when it holds them, and simulate
#'   them otherwise, on a cluster when there are enough of them to pay for
#'   one; the trials simulated are added to the cache.
#'
#' @param design_rows One results row per design.
#' @param design_keys Their keys from [nominal_tie_design_key()].
#' @param design_target_data Their target data, as `load_data()` returns it.
#' @param frequentist_test The test the p-values are computed for.
#' @param simulation_config The simulation configuration.
#' @param n_replicates Number of trials simulated per design.
#' @param parallelization Whether the caller asked for parallelism.
#' @param trial_cache An environment to read from and add to, or `NULL`.
#' @param cluster A cluster from [new_analysis_cluster()] to run on, or
#'   `NULL` to start one here if one is needed.
#' @param exact Logical, one per design: whether its reference is computed
#'   exactly by enumeration instead, so needs no simulated trials.
#'
#' @return A list with one entry per design: its sorted p-values, or `NULL`.
#' @noRd
equivalent_tie_design_p_values <- function(design_rows,
                                           design_keys,
                                           design_target_data,
                                           frequentist_test,
                                           simulation_config,
                                           n_replicates,
                                           parallelization,
                                           trial_cache = NULL,
                                           cluster = NULL,
                                           exact = rep(FALSE, nrow(design_rows))) {
  p_values <- vector("list", nrow(design_rows))

  # which() drops an undetermined answer along with the closed-form designs.
  needs_simulation <- which(
    !vapply(design_target_data, uses_analytical_power, logical(1)) & !exact
  )
  if (length(needs_simulation) == 0) {
    return(p_values)
  }

  trials <- cached_design_trials(
    design_keys[needs_simulation], simulation_config, n_replicates, trial_cache
  )
  to_simulate <- which(vapply(trials, is.null, logical(1)))

  if (length(to_simulate) > 0) {
    # Only the rows being simulated are sent to the workers.
    rows_to_simulate <- lapply(needs_simulation[to_simulate], function(d) {
      return(design_rows[d, , drop = FALSE])
    })

    if (analysis_uses_cluster(parallelization, length(to_simulate) * n_replicates,
                              min_rows = ANALYSIS_PARALLEL_MIN_TRIALS)) {
      if (is.null(cluster)) {
        cluster <- new_analysis_cluster()
        on.exit(cluster$stop(), add = TRUE)
      }
      cluster$get()
      `%dopar%` <- foreach::`%dopar%`

      # The foreach body is evaluated outside the package namespace, where an
      # internal function is not visible; bound here, it travels to the
      # workers as an exported variable, with the namespace as its environment.
      simulate_design <- design_trials
      trials[to_simulate] <- foreach::foreach(row = rows_to_simulate, .packages = c("dplyr", "yaml")) %dopar% {
        simulate_design(row, simulation_config, n_replicates)
      }
    } else {
      trials[to_simulate] <- lapply(rows_to_simulate, design_trials,
                                    simulation_config = simulation_config,
                                    n_replicates = n_replicates)
    }

    if (!is.null(trial_cache)) {
      keys <- trial_cache_key(
        design_keys[needs_simulation][to_simulate], simulation_config, n_replicates
      )
      for (j in seq_along(to_simulate)) {
        assign(keys[j], trials[[to_simulate[j]]], envir = trial_cache)
      }
    }
  }

  for (j in seq_along(needs_simulation)) {
    row <- design_rows[needs_simulation[j], , drop = FALSE]
    p_values[[needs_simulation[j]]] <- sort(trial_p_values(
      trials[[j]],
      theta_0 = row$theta_0,
      alternative = alternative_from_null_space(row$null_space)
    ))
  }

  return(p_values)
}


#' Compute the frequentist power at equivalent tie
#'
#' @description This function computes the frequentist power at equivalent tie for a given set of results and analysis configuration.
#'
#' @param results The results data frame.
#' @param analysis_config The analysis configuration.
#' @param simulation_config The simulation configuration.
#' @param parallelization Whether to simulate the designs in parallel.
#' @param n_replicates Number of Monte Carlo replicates the simulated power
#'   estimates are built from, for the case studies analytical_power() cannot
#'   be used for.
#' @param trial_cache Optional environment holding the simulated trials of
#'   each design, filled here and read back by
#'   `frequentist_power_at_nominal_tie()`, whose separate and pooled powers
#'   read the same trials. `NULL` keeps them for this call only.
#' @param cluster Optional shared cluster from `new_analysis_cluster()`, as
#'   [simulation_analysis()] passes it. `NULL` starts one for this call if the
#'   work warrants it.
#'
#' @return The final results data frame with power and frequentist test columns added.
#'
#' @keywords internal
frequentist_power_at_equivalent_tie <- function(results, analysis_config, simulation_config, parallelization = FALSE, n_replicates = 1000, trial_cache = NULL, cluster = NULL) {
  if (nrow(results) == 0) {
    stop("The results dataframe is empty.")
  }

  # Remove the tie, mcse_tie, conf_int_tie_lower and conf_int_tie_upper if they exist
  results <- results[, !(
    names(results) %in% c(
      "tie",
      "mcse_tie",
      "conf_int_tie_lower",
      "conf_int_tie_upper",
      "frequentist_power_at_equivalent_tie",
      "frequentist_power_at_equivalent_tie_lower",
      "frequentist_power_at_equivalent_tie_upper",
      "frequentist_test",
      "frequentist_reference_calibration"
    )
  )]

  # The type I error rate is joined back onto the scenarios it was computed
  # for - see scenario_matching_columns().
  matching_columns <- scenario_matching_columns(results)

  for (case_study in unique(results$case_study)){
    if (sum(results[results$case_study == case_study, ]['target_treatment_effect'] == results$theta_0) == 0){
      stop("theta_0 not included among the target study treatment effects considered in the simulation study!")
    }
  }

  # Select the results that correspond to TIE computation.
  results_freq_df_tie <- results[results$target_treatment_effect == results$theta_0, c(
    matching_columns,
    c(
      "success_proba",
      "mcse_success_proba",
      "conf_int_success_proba_lower",
      "conf_int_success_proba_upper"
    )
  )]

  # For these results the probability of success corresponds to TIE, so we rename the columns accordingly.
  results_freq_df_tie <- results_freq_df_tie %>%
    dplyr::rename(
      tie = success_proba,
      mcse_tie = mcse_success_proba,
      conf_int_tie_lower = conf_int_success_proba_lower,
      conf_int_tie_upper = conf_int_success_proba_upper
    )

  # We merge the TIE results onto the corresponding scenarios of the results dataframe.
  results <- dplyr::left_join(results, results_freq_df_tie, by = matching_columns)

  frequentist_test <- analysis_config[["frequentist_test"]]

  has_tie <- !is.na(results$tie)
  if (!all(has_tie)) {
    warning(sprintf(
      "TIE is NA for %d of %d rows, so their power at equivalent TIE is left missing.",
      sum(!has_tie), nrow(results)
    ))
  }

  # The trials the separate analysis is read off depend on the design alone,
  # while the results hold one row per design *and* method-parameter
  # combination (56 rows a design in the paper's environment). The target data
  # and the simulated p-values are therefore built once per design, and only
  # the part that depends on the row's own type I error runs per row. The
  # p-values are exactly the ones a per-row simulation would give, because
  # simulate_test_p_values() reseeds from simulation_config$seed every time.
  design_key <- nominal_tie_design_key(
    results[intersect(NOMINAL_TIE_DESIGN_COLUMNS, names(results))]
  )
  design_keys <- unique(design_key[has_tie])
  design_rows <- results[match(design_keys, design_key), , drop = FALSE]
  row_design <- match(design_key, design_keys)

  design_target_data <- lapply(seq_len(nrow(design_rows)), function(d) {
    return(load_data(design_rows[d, , drop = FALSE],
              type = "target",
              reload_data_objects = TRUE))
  })

  # The reference is calibrated on its actual type I error, read off the
  # design's null scenario, which is one of the designs here: its rows carry a
  # TIE of their own.
  design_exact <- exact_reference_designs(design_rows)
  design_null <- match(null_design_keys(results)[match(design_keys, design_key)],
                       design_keys)

  design_p_values <- equivalent_tie_design_p_values(
    design_rows = design_rows,
    design_keys = design_keys,
    design_target_data = design_target_data,
    frequentist_test = frequentist_test,
    simulation_config = simulation_config,
    n_replicates = n_replicates,
    parallelization = parallelization,
    trial_cache = trial_cache,
    cluster = cluster,
    exact = design_exact
  )

  # The exactly enumerated designs: every outcome with its probability, enough
  # to price a randomised test at any level exactly.
  design_support <- vector("list", nrow(design_rows))
  for (d in which(design_exact)) {
    design_support[[d]] <- design_target_data[[d]]$enumerate_support()
  }

  # The type I error draws and the posterior of each rejection count are the
  # only random part left. Seeding them once makes the columns reproducible,
  # whatever the cache held and however the p-values were computed.
  withr::local_seed(simulation_config$seed)

  power <- rep(NA_real_, nrow(results))
  power_lower <- rep(NA_real_, nrow(results))
  power_upper <- rep(NA_real_, nrow(results))
  test_column <- rep(NA_character_, nrow(results))
  calibration_column <- rep(NA_character_, nrow(results))

  for (i in which(has_tie)) {
    d <- row_design[i]
    d_null <- design_null[d]

    if (design_exact[d] && !is.na(d_null)) {
      # The type I error is exact here, and so is the reference's power at it.
      exact_power <- exact_reference_power(
        alpha = results$tie[i],
        alternative_support = design_support[[d]],
        null_support = design_support[[d_null]],
        theta_0 = results$theta_0[i],
        alternative = alternative_from_null_space(results$null_space[i])
      )
      power_estimation <- list(power = exact_power,
                               conf_int_power = c(exact_power, exact_power))
    } else {
      alpha <- list(
        mean = results$tie[i],
        conf_int_lower = results$conf_int_tie_lower[i],
        conf_int_upper = results$conf_int_tie_upper[i],
        mcse = results$mcse_tie[i]
      )

      power_estimation <- compute_power_with_tie_ci(
        alpha = alpha,
        target_data = design_target_data[[d]],
        frequentist_test = frequentist_test,
        theta_0 = results$theta_0[i],
        null_space = results$null_space[i],
        case_study = results$case_study[i],
        simulation_config = simulation_config,
        n_replicates = n_replicates,
        p_values = design_p_values[[d]],
        null_p_values = if (is.na(d_null)) NULL else design_p_values[[d_null]]
      )
    }

    power[i] <- power_estimation$power
    power_lower[i] <- power_estimation$conf_int_power[1]
    power_upper[i] <- power_estimation$conf_int_power[2]
    test_column[i] <- frequentist_test
    calibration_column[i] <- reference_calibration_label(
      design_target_data[[d]], exact = design_exact[d], has_null = !is.na(d_null)
    )
  }

  results$frequentist_power_at_equivalent_tie <- power
  results$frequentist_power_at_equivalent_tie_lower <- power_lower
  results$frequentist_power_at_equivalent_tie_upper <- power_upper
  results$frequentist_test <- test_column
  results$frequentist_reference_calibration <- calibration_column

  return(results)
}


#' Power of the separate and pooling baselines for one scenario row
#'
#' @description Neither baseline depends on the borrowing method, only on the
#'   design, so both branches of [frequentist_power_at_nominal_tie()] need
#'   exactly this computation. A worker cannot see a copy of it that lives
#'   inside the loop body, so it lives here rather than being written twice.
#'
#' @param row One row of the results data frame.
#' @param nominal_tie The type I error rate the power is evaluated at.
#' @param frequentist_test The test the power is computed for.
#' @param simulation_config The simulation configuration.
#' @param n_replicates Number of Monte Carlo replicates the simulated power
#'   estimates are built from, for the case studies analytical_power() cannot
#'   be used for.
#' @param trials Optional trials already simulated for this design, as
#'   [simulated_trials()] returns them. `NULL` simulates them if the powers
#'   need them.
#' @param null_row The results row of the design's null scenario, which the
#'   separate analysis's test is calibrated on, or `NULL` to leave it at its
#'   nominal level.
#' @param null_trials Optional trials already simulated for the null scenario.
#'
#' @return A named list holding the six power columns for that row, and the
#'   calibration the separate analysis's test received.
#' @noRd
nominal_tie_power_row <- function(row,
                                  nominal_tie,
                                  frequentist_test,
                                  simulation_config,
                                  n_replicates,
                                  trials = NULL,
                                  null_row = NULL,
                                  null_trials = NULL) {
  target_data <- load_data(row, type = "target", reload_data_objects = TRUE)
  source_data <- load_data(row, type = "source", reload_data_objects = TRUE)

  assert_target_data_numbers(target_data)
  exact <- exact_reference_designs(row)
  alternative <- alternative_from_null_space(row$null_space)

  # Where both powers are simulated they read the same trials - each would
  # otherwise reseed and generate them again - so they are generated once
  # here, unless the caller already holds them.
  simulated <- identical(target_data$summary_measure_likelihood, "normal") &&
    isFALSE(uses_analytical_power(target_data)) &&
    identical(frequentist_test, "t-test")
  if (!simulated) {
    trials <- NULL
  } else if (is.null(trials)) {
    trials <- simulated_trials(target_data, simulation_config, n_replicates)
  }
  p_values <- if (is.null(trials)) {
    NULL
  } else {
    trial_p_values(trials, row$theta_0, alternative_from_null_space(row$null_space))
  }

  # The separate analysis's test is calibrated on its actual type I error in
  # the null scenario: exactly where the operating characteristics are
  # enumerated, on simulated null trials where they are simulated.
  null_p_values <- NULL
  if (!is.null(null_row) && !exact && !is.null(p_values)) {
    if (is.null(null_trials)) {
      null_target <- load_data(null_row, type = "target", reload_data_objects = TRUE)
      null_trials <- simulated_trials(null_target, simulation_config, n_replicates)
    }
    null_p_values <- trial_p_values(null_trials, row$theta_0, alternative)
  }

  separate <- if (exact && !is.null(null_row)) {
    null_target <- load_data(null_row, type = "target", reload_data_objects = TRUE)
    exact_power <- exact_reference_power(
      alpha = nominal_tie,
      alternative_support = target_data$enumerate_support(),
      null_support = null_target$enumerate_support(),
      theta_0 = row$theta_0,
      alternative = alternative
    )
    list(power = exact_power, conf_int_power = c(exact_power, exact_power))
  } else {
    compute_freq_power(
      alpha = nominal_tie,
      target_data = target_data,
      frequentist_test = frequentist_test,
      theta_0 = row$theta_0,
      null_space = row$null_space,
      case_study = row$case_study,
      simulation_config = simulation_config,
      n_replicates = n_replicates,
      p_values = p_values,
      null_p_values = null_p_values
    )
  }

  pooling <- compute_freq_power_pooling(
    alpha = nominal_tie,
    target_data = target_data,
    source_data = source_data,
    frequentist_test = frequentist_test,
    theta_0 = row$theta_0,
    null_space = row$null_space,
    case_study = row$case_study,
    simulation_config = simulation_config,
    n_replicates = n_replicates,
    trials = trials
  )

  # The bounds are unnamed so that binding the rows together downstream gives
  # plain numeric columns, whichever of the two intervals came back.
  return(list(
    nominal_frequentist_power_separate = separate$power,
    nominal_frequentist_power_separate_lower = unname(separate$conf_int_power[1]),
    nominal_frequentist_power_separate_upper = unname(separate$conf_int_power[2]),
    nominal_frequentist_power_pooling = pooling$power,
    nominal_frequentist_power_pooling_lower = unname(pooling$conf_int_power[1]),
    nominal_frequentist_power_pooling_upper = unname(pooling$conf_int_power[2]),
    frequentist_reference_calibration = reference_calibration_label(
      target_data, exact = exact, has_null = !is.null(null_row)
    )
  ))
}

#' The columns that identify one baseline computation
#'
#' Neither baseline depends on the borrowing method, so rows that agree on
#' every column here yield the same two numbers however they were analysed.
#'
#' This list must stay a *superset* of what [nominal_tie_power_row()] reads.
#' An extra column only splits a group that could have been shared, which
#' costs a little time; a missing one merges rows that genuinely differ, which
#' is silently wrong. It is therefore every design column of the results
#' frame, not only the six `load_data()` happens to read today.
#' @noRd
NOMINAL_TIE_DESIGN_COLUMNS <- c(
  "case_study", "theta_0", "null_space",
  "target_sample_size_control", "target_sample_size_treatment", "target_treatment_effect",
  "target_standard_deviation", "target_to_source_std_ratio",
  "target_treatment_rate", "target_control_rate",
  "drift", "treatment_drift", "control_drift",
  "source_denominator", "source_denominator_change_factor",
  "source_treatment_effect_estimate", "source_standard_error",
  "source_sample_size_control", "source_sample_size_treatment",
  "source_treatment_rate", "source_control_rate",
  "equivalent_source_sample_size_per_arm",
  "summary_measure_likelihood", "endpoint", "sampling_approximation",
  "dropout_probability", "event_time_distribution", "treatment_delay"
)


#' Group rows of a results frame by the design they describe
#'
#' @param design The design columns of the results frame.
#'
#' @return A character vector, one entry per row, equal exactly for rows
#'   describing the same design.
#' @noRd
nominal_tie_design_key <- function(design) {
  if (length(design) == 0) {
    return(rep("", nrow(design)))
  }

  # A missing entry is given a token no value of its own can take, so that it
  # never collides with a column that literally holds the string "NA".
  tokens <- lapply(design, function(column) {
    column <- as.character(column)
    return(ifelse(is.na(column), "\001missing\001", column))
  })

  return(do.call(paste, c(tokens, list(sep = "\002"))))
}


#' Compute the frequentist power at the nominal type I error rate
#'
#' @description Computes the power of the separate and the pooled analysis at
#'   the nominal type I error rate, which the plots use as the two baselines
#'   every borrowing method is read against.
#'
#'   Both baselines are a property of the design alone, while the results
#'   frame holds one row per design *and* method-parameter combination: the
#'   paper's environment repeats each of its 330 designs 56 times. They are
#'   therefore computed once per design and copied to the rows that share it.
#'   That is exact rather than an approximation, because
#'   `simulate_test_p_values()` reseeds from `simulation_config$seed` on every
#'   call, so the repeats were identical to the last bit anyway.
#'
#' @param results The results data frame.
#' @param analysis_config The analysis configuration.
#' @param simulation_config The simulation configuration.
#' @param parallelization Whether the caller asked for parallelism, as
#'   `analysis_runs_in_parallel()` resolves it. A design costs seconds here -
#'   two power computations, either of which may be a simulation - so a run
#'   with many of them is worth spreading over a cluster.
#' @param n_replicates Number of Monte Carlo replicates the simulated power
#'   estimates are built from, for the case studies analytical_power() cannot
#'   be used for.
#' @param trial_cache Optional environment of simulated trials left by
#'   [frequentist_power_at_equivalent_tie()]. Both baselines read the same
#'   trials, so a design found there is not simulated again.
#' @param cluster Optional shared cluster from `new_analysis_cluster()`.
#'   `NULL` starts one for this call if the work warrants it.
#'
#' @return The results data frame with the six baseline power columns added.
frequentist_power_at_nominal_tie <- function(results, analysis_config, simulation_config, parallelization = FALSE, n_replicates = 1000, trial_cache = NULL, cluster = NULL) {
  if (nrow(results) == 0) {
    stop("The results dataframe is empty.")
  }

  results <- results[, !(
    names(results) %in% c(
     "nominal_frequentist_power_separate",
     "nominal_frequentist_power_separate_lower",
     "nominal_frequentist_power_separate_upper",
     "nominal_frequentist_power_pooling",
     "nominal_frequentist_power_pooling_lower",
     "nominal_frequentist_power_pooling_upper",
     "frequentist_reference_calibration"
    )
  )]

  nominal_tie <- analysis_config[["nominal_tie"]]
  frequentist_test <- analysis_config[["frequentist_test"]]

  # One row per distinct design, and the index that puts each computed value
  # back on every row sharing it.
  design_key <- nominal_tie_design_key(
    results[intersect(NOMINAL_TIE_DESIGN_COLUMNS, names(results))]
  )
  representatives <- !duplicated(design_key)
  design_rows <- results[representatives, , drop = FALSE]
  design_index <- match(design_key, design_key[representatives])

  # Each design's null scenario, among the designs, for calibrating the
  # separate analysis's test on its actual type I error.
  null_index <- match(null_design_keys(results)[representatives],
                      design_key[representatives])

  # Every simulated design's trials, once each: the equivalent-TIE step leaves
  # them in `trial_cache`, and those it did not are simulated here, on a
  # cluster when there are enough. A design's separate power reads both its
  # own trials and its null scenario's - another design here - so holding them
  # all before the per-design work keeps either from being drawn twice.
  if (is.null(trial_cache)) {
    trial_cache <- new.env(parent = emptyenv())
  }
  equivalent_tie_design_p_values(
    design_rows = design_rows,
    design_keys = design_key[representatives],
    design_target_data = lapply(seq_len(nrow(design_rows)), function(d) {
      return(load_data(design_rows[d, , drop = FALSE], type = "target", reload_data_objects = TRUE))
    }),
    frequentist_test = frequentist_test,
    simulation_config = simulation_config,
    n_replicates = n_replicates,
    parallelization = parallelization,
    trial_cache = trial_cache,
    cluster = cluster,
    exact = exact_reference_designs(design_rows)
  )
  cached_trials <- cached_design_trials(
    design_key[representatives], simulation_config, n_replicates, trial_cache
  )

  # Only the distinct designs are handed to the workers, rather than the whole
  # results frame: on the paper's environment that is 330 rows to serialise
  # instead of 18,480, each carrying its parameters as JSON.
  compute_design <- function(i, trials = NULL, null_trials = NULL) {
    return(nominal_tie_power_row(
      row = design_rows[i, ],
      nominal_tie = nominal_tie,
      frequentist_test = frequentist_test,
      simulation_config = simulation_config,
      n_replicates = n_replicates,
      trials = trials,
      null_row = if (is.na(null_index[i])) NULL else design_rows[null_index[i], ],
      null_trials = null_trials
    ))
  }

  # A design whose trials are cached costs milliseconds - two tests on
  # trials already drawn - so those are computed here.
  computed_list <- vector("list", nrow(design_rows))
  is_cached <- !vapply(cached_trials, is.null, logical(1))
  for (i in which(is_cached)) {
    computed_list[[i]] <- compute_design(
      i, cached_trials[[i]],
      if (is.na(null_index[i])) NULL else cached_trials[[null_index[i]]]
    )
  }
  # Dropped before compute_design() goes to the workers: foreach ships the
  # environment it was defined in, which would carry every cached trial to
  # each of them.
  cached_trials <- NULL
  to_compute <- which(!is_cached)

  if (length(to_compute) > 0 &&
      analysis_uses_cluster(parallelization, length(to_compute))) {
    if (is.null(cluster)) {
      cluster <- new_analysis_cluster()
      on.exit(cluster$stop(), add = TRUE)
    }
    cluster$get()
    `%dopar%` <- foreach::`%dopar%`

    computed_list[to_compute] <- foreach::foreach(i = to_compute, .packages = c("dplyr", "yaml", "pwr", "BSDA")) %dopar% {
      compute_design(i)
    }
  } else if (length(to_compute) > 0) {
    # One progress bar for the whole loop. Building it inside the loop, as
    # this did, restarts it on every step, and what gets printed is then
    # thousands of bars that never pass the first few percent.
    pb <- txtProgressBar(min = 0, max = length(to_compute), style = 3)
    on.exit(close(pb), add = TRUE)

    for (k in seq_along(to_compute)) {
      computed_list[[to_compute[k]]] <- compute_design(to_compute[k])

      # Update progress bar
      setTxtProgressBar(pb, k)
    }
  }

  # Assigned by name rather than cbind()ed on, so that a stale column is
  # replaced instead of being shadowed by a second copy of itself.
  computed <- dplyr::bind_rows(computed_list)
  for (column in names(computed)) {
    values <- computed[[column]]
    results[[column]] <- (if (is.character(values)) values else as.numeric(values))[design_index]
  }

  return(results)
}


#' Power of the pooled two-proportion test, the source counts held fixed
#'
#' @description The pooled analysis adds the source study's responders to the
#'   target's and tests the pooled response rates with the arcsine two-sample
#'   test that [pwr::pwr.2p2n.test()] describes. The source counts are the ones
#'   observed, fixed across the simulated trials, so the rejection rate is
#'   summed over the target responder counts alone, weighted by their binomial
#'   probabilities. [pwr::pwr.2p2n.test()] itself would treat the pooled counts
#'   as a fresh sample of the pooled size, and so overstate how much the pooled
#'   rates vary from one trial to the next.
#'
#'   The test is of equal response rates, i.e. a null boundary of zero.
#'
#' @param alpha One-sided significance level.
#' @param target_data Target data, with its true response rates.
#' @param source_data Source data.
#' @param alternative "greater" or "less".
#' @return The power.
#' @keywords internal
pooled_binomial_power <- function(alpha, target_data, source_data, alternative) {
  n_treatment <- as.integer(target_data$sample_size_treatment)
  n_control <- as.integer(target_data$sample_size_control)
  source_treatment <- counts_from_rate(source_data$treatment_rate, source_data$sample_size_treatment)
  source_control <- counts_from_rate(source_data$control_rate, source_data$sample_size_control)
  pooled_treatment_size <- n_treatment + source_data$sample_size_treatment
  pooled_control_size <- n_control + source_data$sample_size_control

  target_treatment_counts <- 0:n_treatment
  target_control_counts <- 0:n_control
  rate_treatment <- (source_treatment + target_treatment_counts) / pooled_treatment_size
  rate_control <- (source_control + target_control_counts) / pooled_control_size

  # Rows are the treatment responders, columns the control responders.
  statistic <- outer(2 * asin(sqrt(rate_treatment)), 2 * asin(sqrt(rate_control)), "-") /
    sqrt(1 / pooled_treatment_size + 1 / pooled_control_size)
  critical <- stats::qnorm(1 - alpha)
  reject <- if (alternative == "greater") statistic > critical else statistic < -critical

  probability <- outer(
    stats::dbinom(target_treatment_counts, n_treatment, target_data$treatment_rate),
    stats::dbinom(target_control_counts, n_control, target_data$control_rate)
  )
  return(sum(probability[reject]))
}


## ---- Calibrating the reference test on its actual type I error -------------
##
## The power of a frequentist test of the target data alone is the reference
## every borrowing method's power is read against, at the method's own type I
## error and at the nominal one. The test is a one-sample t-test, and its
## significance level used to be set to the type I error it should match. That
## level is only the test's actual size when its statistic is t-distributed:
## the Wald statistic of a Cox model, a negative binomial regression or a
## difference in proportions is approximately normal instead, and a t-test on
## it is slightly conservative - at 123 per arm it rejects 2.4% of the time at
## a nominal 2.5%. Matched by its nominal level, the reference then spent less
## type I error than the method it was compared with, and every method looked
## about 1.5 points more powerful than it is.
##
## The reference's critical value is therefore read off the test's actual
## behaviour in the design's null scenario: from the same simulated trials
## where the operating characteristics are simulated, and exactly, with a
## randomised boundary, where they are enumerated. The continuous endpoints
## keep their closed-form t-test, which is exact there.

#' Columns that identify a scenario apart from its treatment effect
#'
#' The type I error is computed in each scenario's null counterpart and joined
#' back onto it on these columns, so a design axis missing here would match
#' one scenario's TIE to several rows.
#' @noRd
SCENARIO_MATCHING_COLUMNS <- c(
  "method",
  "parameters",
  "control_drift",
  "source_denominator",
  "source_denominator_change_factor",
  "case_study",
  "target_to_source_std_ratio",
  "target_sample_size_control", "target_sample_size_treatment",
  "theta_0",
  "null_space",
  "sampling_approximation",
  "summary_measure_likelihood",
  "source_sample_size_treatment",
  "source_sample_size_control",
  "endpoint",
  "source_standard_error",
  "source_treatment_effect_estimate",
  "equivalent_source_sample_size_per_arm"
)

#' The scenario-matching columns a results frame carries
#'
#' Results written before the time-to-event design axes existed do not carry
#' them, and such a run only ever had one design, so they join only when the
#' columns are actually there.
#'
#' @param results A results frame.
#' @return A character vector of column names.
#' @noRd
scenario_matching_columns <- function(results) {
  return(c(SCENARIO_MATCHING_COLUMNS, intersect(time_to_event_design_columns, names(results))))
}

#' The design key of each row's null scenario
#'
#' @description The null scenario of a row is the one that agrees with it on
#'   every scenario-matching column, method aside, and whose treatment effect
#'   is `theta_0`.
#'
#' @param results A results frame.
#' @return A character vector, one per row: the [nominal_tie_design_key()] of
#'   its null scenario, or `NA` when the frame holds none.
#' @noRd
null_design_keys <- function(results) {
  columns <- intersect(setdiff(scenario_matching_columns(results), c("method", "parameters")),
                       names(results))
  scenario <- nominal_tie_design_key(results[columns])
  design <- nominal_tie_design_key(results[intersect(NOMINAL_TIE_DESIGN_COLUMNS, names(results))])
  is_null <- results$target_treatment_effect == results$theta_0
  is_null[is.na(is_null)] <- FALSE
  return(design[is_null][match(scenario, scenario[is_null])])
}

#' Whether a row's operating characteristics were enumerated exactly
#'
#' @param rows Rows of a results frame.
#' @return A logical vector, one per row; `FALSE` where the frame predates the
#'   `exact_ocs` column or leaves it missing.
#' @noRd
exact_reference_designs <- function(rows) {
  if (!"exact_ocs" %in% names(rows)) {
    return(rep(FALSE, nrow(rows)))
  }
  return(as.logical(rows$exact_ocs) %in% TRUE)
}

#' Thresholds at which a test's actual type I error equals given levels
#'
#' @description For each level `alpha`, the p-value threshold `c` such that
#'   the share of the null p-values strictly below `c` is
#'   `floor(alpha * N) / N`, within `1 / N` of `alpha`. A test rejecting when
#'   its p-value is below `c` then has the actual rejection rate `alpha` in the
#'   null scenario, whatever its nominal level there. Null p-values tied at the
#'   threshold all stay unrejected, so the rate never exceeds `alpha`.
#'
#' @param alpha Levels, in `[0, 1]`.
#' @param null_p_values The test's p-values on the null scenario's trials.
#' @return One threshold per level; `Inf` rejects every trial.
#' @keywords internal
calibrated_levels <- function(alpha, null_p_values) {
  null_p_values <- sort(null_p_values[!is.na(null_p_values)])
  n <- length(null_p_values)
  if (n == 0) {
    stop("No null p-value to calibrate the test on.", call. = FALSE)
  }
  # The small offset keeps a level that is a whole number of trials up to
  # rounding - 0.025 * 10000 - on that number rather than one below it.
  rejected <- pmin(pmax(floor(alpha * n + 1e-9), 0), n)
  return(ifelse(rejected >= n, Inf, null_p_values[pmin(rejected + 1, n)]))
}

#' The reference test's p-value and probability for each enumerated outcome
#'
#' @description Outcomes whose summary measure is not estimable are left out
#'   and the rest renormalised, as the enumerated operating characteristics
#'   do. An outcome whose statistic is undefined - both arms at 0 or at n
#'   responders, so no variability - is kept with an infinite p-value: it
#'   never rejects.
#'
#' @param support Output of `enumerate_support()`.
#' @param theta_0 Boundary of the null hypothesis space.
#' @param alternative "greater" or "less".
#' @return A list with the p-values `p` and probabilities `w`.
#' @noRd
reference_outcomes <- function(support, theta_0, alternative) {
  samples <- support$samples
  estimable <- is.finite(samples$treatment_effect_estimate) &
    is.finite(samples$standard_deviation)
  p <- trial_p_values(samples[estimable, , drop = FALSE], theta_0, alternative)
  p[!is.finite(p)] <- Inf
  w <- support$weights[estimable]
  return(list(p = p, w = w / sum(w)))
}

#' Exact power of the reference test at an exact type I error
#'
#' @description The most powerful use of a test statistic at a given size
#'   rejects every outcome below a threshold and, on a discrete sample space,
#'   a fraction of the outcomes at it. The threshold `c` and the fraction
#'   `gamma` are chosen on the null scenario's enumerated outcomes so that the
#'   rejected probability there is exactly `alpha`, and the power is the
#'   probability the same randomised test rejects in the scenario of
#'   interest. Without randomisation a binary endpoint's test could only reach
#'   a few levels, and could not be matched to an arbitrary type I error.
#'
#' @param alpha The type I error to match, in `[0, 1]`.
#' @param alternative_support,null_support Enumerated outcomes of the scenario
#'   of interest and of its null scenario, from `enumerate_support()`.
#' @param theta_0 Boundary of the null hypothesis space.
#' @param alternative "greater" or "less".
#' @return The power, one value per element of `alpha`.
#' @keywords internal
exact_reference_power <- function(alpha, alternative_support, null_support,
                                  theta_0, alternative) {
  null <- reference_outcomes(null_support, theta_0, alternative)
  scenario <- reference_outcomes(alternative_support, theta_0, alternative)

  finite <- is.finite(null$p)
  values <- sort(unique(null$p[finite]))
  mass <- as.vector(rowsum(null$w[finite], match(null$p[finite], values),
                           reorder = TRUE))
  below <- c(0, cumsum(mass))

  return(vapply(alpha, function(level) {
    if (is.na(level)) {
      return(NA_real_)
    }
    j <- findInterval(level, below)
    if (j > length(values)) {
      # Every outcome that can reject does, and the level is not reached.
      return(sum(scenario$w[is.finite(scenario$p)]))
    }
    threshold <- values[j]
    gamma <- (level - below[j]) / mass[j]
    return(sum(scenario$w[scenario$p < threshold]) +
      gamma * sum(scenario$w[scenario$p == threshold]))
  }, numeric(1)))
}

#' How the reference test of a design is calibrated
#'
#' @param target_data The design's target data.
#' @param exact Whether its operating characteristics are enumerated exactly.
#' @param has_null Whether its null scenario was found.
#' @return `"exact, randomised"`, `"simulated null"` or `"none"`.
#' @noRd
reference_calibration_label <- function(target_data, exact, has_null) {
  if (!has_null) {
    return("none")
  }
  if (exact) {
    return("exact, randomised")
  }
  if (isTRUE(uses_analytical_power(target_data))) {
    return("none")
  }
  return("simulated null")
}
