# The binomial power prior on the lattice as it was computed before the block
# and kernel-cache rewrite (R/binomial_lattice_fast.R), copied verbatim from
# R/binomial_quadrature.R and R/binomial_ebpp.R with every function prefixed
# with `reference_`. The tests compare the package's functions against these.
# They use the package's unchanged binomial_lattice_width(),
# binomial_lattice_source_mass() and grid_posterior().

reference_lattice_shift_matrix <- function(values, rows, shifts, fill = 0) {
  N <- length(values)
  padded <- c(rep(fill, N), values, rep(fill, N))
  index <- rep(rows, times = length(shifts)) + rep(shifts, each = length(rows)) + N
  matrix(padded[index], nrow = length(rows))
}


reference_binomial_power_prior_target_terms <- function(n_control_source,
                                              n_successes_control_source,
                                              n_treatment_source,
                                              n_successes_treatment_source,
                                              n_control,
                                              n_successes_control,
                                              n_treatment,
                                              n_successes_treatment,
                                              n_lattice = 1000L) {
  N <- as.integer(n_lattice)
  rates <- (seq_len(N) - 0.5) / N
  arm_log_likelihood <- function(n, successes) {
    log_likelihood <- successes * log(rates) + (n - successes) * log1p(-rates)
    log_likelihood - max(log_likelihood)
  }
  source_control <- arm_log_likelihood(n_control_source, n_successes_control_source)
  source_treatment <- arm_log_likelihood(n_treatment_source, n_successes_treatment_source)
  target_control <- exp(arm_log_likelihood(n_control, n_successes_control))
  target_treatment <- exp(arm_log_likelihood(n_treatment, n_successes_treatment))

  control_points <- which(target_control > 1e-20)
  treatment_points <- which(target_treatment > 1e-20)
  differences <- seq(min(treatment_points) - max(control_points),
                     max(treatment_points) - min(control_points))

  # J x K: the target likelihood, zero where the target treatment rate leaves [0, 1].
  likelihood <- target_control[control_points] *
    reference_lattice_shift_matrix(target_treatment, control_points, differences)

  # N x K: summed over the target control rate against the prior density of
  # the risk difference given the two control rates, 1 / (N - |i - j|).
  target <- binomial_lattice_width(N)[, control_points, drop = FALSE] %*% likelihood

  source_log_likelihood <- source_control +
    reference_lattice_shift_matrix(source_treatment, seq_len(N), differences, fill = -Inf)

  list(
    n_lattice = N,
    differences = differences,
    target = target,
    source_log_likelihood = source_log_likelihood,
    source_control_log_likelihood = source_control,
    source_treatment_log_likelihood = source_treatment
  )
}


reference_binomial_power_prior_log_marginal <- function(terms, power_parameter,
                                              bilinear = reference_binomial_power_prior_bilinear(terms)) {
  mass <- binomial_lattice_source_mass(terms$n_lattice)
  control <- exp(outer(terms$source_control_log_likelihood, power_parameter))
  treatment <- exp(outer(terms$source_treatment_log_likelihood, power_parameter))
  # The numerator is a bilinear form in the two source arms' discounted
  # likelihoods, sum_{i,s} T(i, s - i) a_i^gamma b_s^gamma, so every power
  # parameter costs two matrix-vector products rather than an exponential per
  # lattice point.
  numerator <- log(colSums(control * (bilinear %*% treatment)))
  normalizer <- log(colSums(control * (mass %*% treatment)))
  values <- numerator - normalizer

  # Where the products underflow - only under extreme conflict, at power
  # parameters far from the maximum - the numerator is recomputed on the log
  # scale.
  for (index in which(!is.finite(values))) {
    keep <- terms$target > 0 & is.finite(terms$source_log_likelihood)
    exponent <- log(terms$target[keep]) +
      power_parameter[index] * terms$source_log_likelihood[keep]
    largest <- max(exponent)
    values[index] <- largest + log(sum(exp(exponent - largest))) - normalizer[index]
  }
  values
}


reference_binomial_power_prior_bilinear <- function(terms) {
  N <- terms$n_lattice
  rows <- rep(seq_len(N), times = length(terms$differences))
  columns <- rows + rep(terms$differences, each = N)
  inside <- columns >= 1L & columns <= N
  bilinear <- matrix(0, N, N)
  bilinear[cbind(rows[inside], columns[inside])] <- terms$target[inside]
  bilinear
}


reference_binomial_power_prior_empirical_bayes <- function(terms) {
  bilinear <- reference_binomial_power_prior_bilinear(terms)
  grid <- seq(0, 1, by = 0.05)
  values <- reference_binomial_power_prior_log_marginal(terms, grid, bilinear)
  best <- which.max(values)
  lower <- grid[max(1L, best - 1L)]
  upper <- grid[min(length(grid), best + 1L)]
  refined <- stats::optimize(
    function(gamma) reference_binomial_power_prior_log_marginal(terms, gamma, bilinear),
    interval = c(lower, upper), maximum = TRUE, tol = 1e-6
  )
  if (refined$objective > values[best]) refined$maximum else grid[best]
}


reference_binomial_power_prior_terms_posterior <- function(terms, power_parameter) {
  exponent <- power_parameter * terms$source_log_likelihood
  exponent[!is.finite(exponent)] <- -Inf
  largest <- max(exponent[terms$target > 0])
  density <- colSums(terms$target * exp(exponent - largest))
  N <- terms$n_lattice
  differences <- terms$differences
  support <- c(min(differences) - 1L, differences, max(differences) + 1L)
  grid_posterior(support / N, c(0, density, 0))
}


reference_binomial_power_prior_posterior <- function(power_parameter,
                                           n_control_source,
                                           n_successes_control_source,
                                           n_treatment_source,
                                           n_successes_treatment_source,
                                           n_control,
                                           n_successes_control,
                                           n_treatment,
                                           n_successes_treatment,
                                           n_lattice = 1000L,
                                           control_rate = NULL) {
  if (!is.numeric(power_parameter) || length(power_parameter) != 1 ||
      is.na(power_parameter) || power_parameter < 0 || power_parameter > 1) {
    stop("The power parameter must be a single number in [0, 1].", call. = FALSE)
  }
  gamma <- power_parameter
  N <- as.integer(n_lattice)
  rates <- (seq_len(N) - 0.5) / N

  # Each arm's likelihood on the lattice, scaled to a maximum of one.
  arm_log_likelihood <- function(n, successes) {
    log_likelihood <- successes * log(rates) + (n - successes) * log1p(-rates)
    log_likelihood - max(log_likelihood)
  }
  source_control <- exp(gamma * arm_log_likelihood(n_control_source, n_successes_control_source))
  source_treatment <- exp(gamma * arm_log_likelihood(n_treatment_source, n_successes_treatment_source))
  target_control <- exp(arm_log_likelihood(n_control, n_successes_control))
  target_treatment <- exp(arm_log_likelihood(n_treatment, n_successes_treatment))

  control_points <- if (is.null(control_rate)) {
    which(target_control > 1e-20)
  } else {
    min(N, max(1L, as.integer(ceiling(control_rate * N))))
  }
  if (!is.null(control_rate)) {
    target_control[] <- 0
    target_control[control_points] <- 1
  }
  treatment_points <- which(target_treatment > 1e-20)
  differences <- seq(min(treatment_points) - max(control_points),
                     max(treatment_points) - min(control_points))

  # The source control rates that carry any discounted likelihood: all of them
  # at a power parameter near 0, a narrow band near 1.
  source_rows <- which(source_control > 1e-300)

  # R x K: the discounted source likelihood at every source control rate i and
  # treatment effect k, zero where the source treatment rate leaves [0, 1].
  source_weight <- source_control[source_rows] *
    reference_lattice_shift_matrix(source_treatment, source_rows, differences)

  # J x K: summed over the source control rate against the prior density of
  # theta given the two control rates, 1 / (N - |i - j|).
  width <- binomial_lattice_width(N)[source_rows, control_points, drop = FALSE]
  kernel <- crossprod(width, source_weight)

  # The target likelihood, zero where the target treatment rate leaves [0, 1].
  likelihood <- target_control[control_points] *
    reference_lattice_shift_matrix(target_treatment, control_points, differences)

  density <- colSums(kernel * likelihood)
  support <- c(min(differences) - 1L, differences, max(differences) + 1L)
  grid_posterior(support / N, c(0, density, 0))
}
