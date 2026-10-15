## The heterogeneity-prior labels have to say what the model does: the
## inverse gamma is placed on tau^2 and the half-normal on tau, and the
## smallest inverse-gamma shape is 1/1000. Rounded to two decimal places
## before the label was built, that shape used to print as "alpha = 0".

commensurate_label <- function(...) {
  make_labels_from_parameters(
    data.frame(..., check.names = FALSE),
    "commensurate_power_prior"
  )
}

test_that("the inverse gamma is labelled as a prior on tau squared", {
  expect_identical(
    commensurate_label(
      heterogeneity_prior.family = "inverse_gamma",
      heterogeneity_prior.alpha = 1 / 1000,
      heterogeneity_prior.beta = 1
    ),
    "$\\tau^2 \\sim IG(\\alpha = 0.001, \\beta = 1)$"
  )
  expect_identical(
    commensurate_label(
      heterogeneity_prior.family = "inverse_gamma",
      heterogeneity_prior.alpha = 1 / 7,
      heterogeneity_prior.beta = 1
    ),
    "$\\tau^2 \\sim IG(\\alpha = 0.14, \\beta = 1)$"
  )
})

test_that("the half-normal is labelled as a prior on tau", {
  expect_identical(
    commensurate_label(
      heterogeneity_prior.family = "half_normal",
      heterogeneity_prior.std_dev = 5
    ),
    "$\\tau \\sim HN(5)$"
  )
})

test_that("the family is inferred when only the prior's parameters are given", {
  expect_identical(
    commensurate_label(alpha = 1 / 3, beta = 1),
    "$\\tau^2 \\sim IG(\\alpha = 0.33, \\beta = 1)$"
  )
})

test_that("the Cauchy is labelled as a prior on log tau without reusing gamma", {
  # gamma is the power parameter in the paper, so the Cauchy scale cannot be
  # written as gamma in the same figure.
  expect_identical(
    commensurate_label(
      heterogeneity_prior.family = "cauchy",
      heterogeneity_prior.location = 0,
      heterogeneity_prior.scale = 30
    ),
    "$\\log \\tau \\sim Cauchy(0, 30)$"
  )
})
