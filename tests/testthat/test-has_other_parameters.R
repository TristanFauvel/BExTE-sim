## Every metric-versus-parameters plot and table loops over one of a method's
## parameters at a time and conditions on the rest. Deciding whether there are
## any "rest" left is not the same question as whether the frame has rows: a
## method with a single parameter loses its only column, and a data frame with
## no columns keeps every one of its rows. A nrow()-only check therefore said
## there was still something to condition on, and the loop then asked
## make_labels_from_parameters() to label a row holding nothing, which died
## with "invalid subscript type 'list'" and took the whole figure with it.

test_that("dropping the only parameter leaves nothing to condition on", {
  # separate and pooling carry initial_prior and nothing else.
  one_parameter <- data.frame(initial_prior = rep("noninformative", 3))
  remaining <- unique(one_parameter[, -1])

  # The shape that fooled the old check: rows intact, no columns.
  expect_equal(nrow(remaining), 1)
  expect_equal(ncol(remaining), 0)
  expect_false(has_other_parameters(remaining))
})

test_that("a method with several parameters still has others to condition on", {
  # RMP carries prior_weight, initial_prior and empirical_bayes.
  three <- data.frame(
    prior_weight = c(0, 0.5, 1),
    initial_prior = "noninformative",
    empirical_bayes = TRUE
  )

  expect_true(has_other_parameters(unique(three[, -1])))
})

test_that("a remainder that collapses to a vector counts as nothing", {
  # Two parameters: dropping one leaves a bare vector, which has no nrow().
  # This is the case the original check already handled, and it must keep
  # behaving the same way - it decides whether a figure is split by parameter.
  two <- data.frame(
    power_parameter = c(0, 0.5),
    initial_prior = "noninformative"
  )

  expect_false(has_other_parameters(unique(two[, -1])))
})

test_that("an empty remainder counts as nothing", {
  expect_false(has_other_parameters(data.frame()))
  expect_false(has_other_parameters(NULL))
})
