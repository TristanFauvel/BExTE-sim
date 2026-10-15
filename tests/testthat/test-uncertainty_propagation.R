test_that("binomial_replicate_count recovers the replicate count from the Monte Carlo standard error", {
  n_replicates <- 1000
  successes <- 25
  estimate <- successes / n_replicates
  mcse <- sqrt(estimate * (1 - estimate) / n_replicates)

  expect_equal(
    binomial_replicate_count(
      estimate = estimate,
      mcse = mcse,
      conf_int_lower = NA_real_,
      conf_int_upper = NA_real_
    ),
    n_replicates
  )
})

test_that("binomial_replicate_count falls back on the exact interval when no success was observed", {
  n_replicates <- 500
  exact_interval <- binom.test(0, n_replicates, conf.level = 0.95)$conf.int

  expect_equal(
    binomial_replicate_count(
      estimate = 0,
      mcse = 0,
      conf_int_lower = exact_interval[1],
      conf_int_upper = exact_interval[2]
    ),
    n_replicates
  )
})

test_that("binomial_replicate_count falls back on the exact interval when every replicate succeeded", {
  n_replicates <- 250
  exact_interval <- binom.test(n_replicates, n_replicates, conf.level = 0.95)$conf.int

  expect_equal(
    binomial_replicate_count(
      estimate = 1,
      mcse = 0,
      conf_int_lower = exact_interval[1],
      conf_int_upper = exact_interval[2]
    ),
    n_replicates
  )
})

test_that("binomial_replicate_count approximates the count from the interval width when no standard error was stored", {
  n_replicates <- 1000
  estimate <- 0.025
  exact_interval <- binom.test(estimate * n_replicates, n_replicates, conf.level = 0.95)$conf.int

  recovered <- binomial_replicate_count(
    estimate = estimate,
    mcse = NA_real_,
    conf_int_lower = exact_interval[1],
    conf_int_upper = exact_interval[2]
  )

  expect_false(is.na(recovered))
  expect_gt(recovered, 0.5 * n_replicates)
  expect_lt(recovered, 2 * n_replicates)
})

test_that("sample_binomial_proportion keeps every draw inside the unit interval without discarding any", {
  set.seed(1)
  draws <- sample_binomial_proportion(
    n_samples = 2000,
    estimate = 0.002,
    n_replicates = 1000
  )

  expect_length(draws, 2000)
  expect_true(all(is.finite(draws)))
  expect_true(all(draws > 0 & draws < 1))
})

test_that("sample_binomial_proportion targets the Jeffreys posterior rather than a symmetric normal", {
  set.seed(2)
  n_replicates <- 1000
  estimate <- 0.025
  draws <- sample_binomial_proportion(
    n_samples = 200000,
    estimate = estimate,
    n_replicates = n_replicates
  )

  successes <- estimate * n_replicates
  expected_mean <- (successes + 0.5) / (n_replicates + 1)

  expect_equal(mean(draws), expected_mean, tolerance = 1e-3)
  # The Jeffreys posterior for a small proportion is right-skewed, unlike the
  # symmetric normal the interval width used to be converted into.
  expect_gt(mean(draws > estimate) , 0.5 - 0.05)
  expect_gt(
    quantile(draws, 0.975) - expected_mean,
    expected_mean - quantile(draws, 0.025)
  )
})

test_that("mover_difference_ci is wider than either of the two intervals it combines", {
  interval <- mover_difference_ci(
    estimate_1 = 0.8, lower_1 = 0.76, upper_1 = 0.84,
    estimate_2 = 0.6, lower_2 = 0.55, upper_2 = 0.65
  )

  expect_equal(unname(interval[2] - interval[1]), sqrt(0.04^2 + 0.05^2) * 2)
  expect_gt(interval[2] - interval[1], 0.84 - 0.76)
  expect_gt(interval[2] - interval[1], 0.65 - 0.55)
})

test_that("mover_difference_ci brackets the difference of the point estimates", {
  interval <- mover_difference_ci(
    estimate_1 = 0.8, lower_1 = 0.76, upper_1 = 0.84,
    estimate_2 = 0.6, lower_2 = 0.55, upper_2 = 0.65
  )

  expect_lt(interval[1], 0.2)
  expect_gt(interval[2], 0.2)
})

test_that("mover_difference_ci preserves the asymmetry of an exact binomial interval", {
  exact_interval <- binom.test(5, 1000, conf.level = 0.95)$conf.int
  estimate <- 0.005

  interval <- mover_difference_ci(
    estimate_1 = estimate,
    lower_1 = exact_interval[1],
    upper_1 = exact_interval[2],
    estimate_2 = 0,
    lower_2 = 0,
    upper_2 = 0
  )

  expect_equal(unname(interval[1]), estimate - (estimate - exact_interval[1]))
  expect_equal(unname(interval[2]), estimate + (exact_interval[2] - estimate))
  # An exact interval is right-skewed at this estimate, and MOVER must not
  # symmetrise it.
  expect_gt(interval[2] - estimate, estimate - interval[1])
})

test_that("mover_difference_ci narrows the interval when the two estimates are positively correlated", {
  independent <- mover_difference_ci(
    estimate_1 = 0.8, lower_1 = 0.76, upper_1 = 0.84,
    estimate_2 = 0.6, lower_2 = 0.55, upper_2 = 0.65
  )
  correlated <- mover_difference_ci(
    estimate_1 = 0.8, lower_1 = 0.76, upper_1 = 0.84,
    estimate_2 = 0.6, lower_2 = 0.55, upper_2 = 0.65,
    correlation = 0.5
  )

  expect_lt(correlated[2] - correlated[1], independent[2] - independent[1])
})

test_that("mover_difference_ci collapses to zero width for perfectly correlated identical estimates", {
  interval <- mover_difference_ci(
    estimate_1 = 0.7, lower_1 = 0.65, upper_1 = 0.75,
    estimate_2 = 0.7, lower_2 = 0.65, upper_2 = 0.75,
    correlation = 1
  )

  expect_equal(unname(interval), c(0, 0))
})

test_that("mover_difference_ci propagates missing inputs instead of inventing a bound", {
  interval <- mover_difference_ci(
    estimate_1 = 0.8, lower_1 = NA_real_, upper_1 = 0.84,
    estimate_2 = 0.6, lower_2 = 0.55, upper_2 = 0.65
  )

  expect_true(all(is.na(interval)))
})
