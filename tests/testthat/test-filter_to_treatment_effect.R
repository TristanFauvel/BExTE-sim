test_that("filter_to_treatment_effect keeps the matching target effects", {
  df <- data.frame(
    target_treatment_effect = c(0, 0.25, 0.5),
    source_treatment_effect_estimate = 0.5
  )
  expect_equal(filter_to_treatment_effect(df, "consistent")$target_treatment_effect, 0.5)
  expect_equal(filter_to_treatment_effect(df, "partially_consistent")$target_treatment_effect, 0.25)
})

test_that("filter_to_treatment_effect rejects other effects", {
  df <- data.frame(target_treatment_effect = 0, source_treatment_effect_estimate = 0.5)
  expect_error(filter_to_treatment_effect(df, "no_effect"), "non-positive")
  expect_error(filter_to_treatment_effect(df, "other"), "must be")
})
