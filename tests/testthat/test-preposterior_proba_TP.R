
test_that("true-positive probability rejects an invalid null space", {
  expect_error(
    BExTE:::preposterior_proba_TP_MC(
      conditional_proba_success = c(0.1, 0.9),
      design_prior_samples = c(-1, 1),
      theta_0 = 0,
      null_space = "middle"
    ),
    "Null space must be either 'left' or 'right'",
    fixed = TRUE
  )

  expect_error(
    preposterior_proba_TP(
      conditional_proba_success = c(0.1, 0.9),
      treatment_effect_values = c(-1, 1),
      theta_0 = 0,
      null_space = "middle",
      design_prior_pdf = c(0.5, 0.5)
    ),
    "Null space must be either 'left' or 'right'",
    fixed = TRUE
  )
})

test_that("true-positive probability rejects misaligned inputs", {
  expect_error(
    BExTE:::preposterior_proba_TP_MC(
      conditional_proba_success = c(0.1, 0.9),
      design_prior_samples = c(-2, -1, 1, 2),
      theta_0 = 0,
      null_space = "left"
    ),
    "must have the same length",
    fixed = TRUE
  )

  expect_error(
    preposterior_proba_TP(
      conditional_proba_success = c(0.1, 0.5),
      treatment_effect_values = c(-1, 0, 1),
      theta_0 = 0,
      null_space = "left",
      design_prior_pdf = c(0.2, 0.6, 0.2)
    ),
    "must have the same length",
    fixed = TRUE
  )
})

test_that("true-positive probability adds the design-prior mass beyond the grid", {
  treatment_effect_values <- seq(0, 1, length.out = 101)

  result <- preposterior_proba_TP(
    conditional_proba_success = rep(1, length(treatment_effect_values)),
    treatment_effect_values = treatment_effect_values,
    theta_0 = 0,
    null_space = "left",
    design_prior_pdf = stats::dnorm(treatment_effect_values),
    design_prior_cdf = function(x) stats::pnorm(x)
  )

  # The rejection rate is 1 over the whole alternative space, so the integral is
  # the design-prior mass above the first alternative grid point, tail included.
  expect_equal(result, 1 - stats::pnorm(0.01), tolerance = 1e-5)
})

test_that("true-positive probability warns when the grid stops before the rejection rate reaches one", {
  treatment_effect_values <- seq(0, 1, length.out = 101)

  expect_warning(
    preposterior_proba_TP(
      conditional_proba_success = rep(0.4, length(treatment_effect_values)),
      treatment_effect_values = treatment_effect_values,
      theta_0 = 0,
      null_space = "left",
      design_prior_pdf = stats::dnorm(treatment_effect_values),
      design_prior_cdf = function(x) stats::pnorm(x)
    ),
    "does not extend far enough into the alternative space"
  )
})
