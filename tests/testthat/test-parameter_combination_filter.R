## Selecting the parameter combination a figure plots is done by merging the
## results against a one-row data frame of the wanted parameter values.
##
## The merge has to be handed that data frame itself. Writing
## `parameters_combinations[, ]` drops a SINGLE-column frame to a bare vector,
## and merge() then has no column in common to join on: it either cross-joins,
## keeping every parameter value, or returns nothing when a name happens to
## collide. Figure S3 plots the conditional power prior at gamma = 0.25 and is
## the only manifest entry whose selection has exactly one parameter, which is
## why all five of its gamma values were drawn joined into one zig-zagging
## line while two-parameter figures like S4 were fine.

test_that("a single-parameter selection keeps only the matching rows", {
  results <- data.frame(
    power_parameter = rep(c(0, 0.25, 0.5, 0.75, 1), each = 3),
    drift = rep(1:3, times = 5),
    success_proba = seq_len(15) / 15
  )

  kept <- filter_to_parameter_combinations(
    results, data.frame(power_parameter = 0.25)
  )

  expect_equal(nrow(kept), 3)
  expect_equal(unique(kept$power_parameter), 0.25)
  ## One value per drift point is the whole property: more than one is what
  ## makes the line double back on itself.
  expect_equal(sort(kept$drift), 1:3)
})

test_that("a multi-parameter selection also keeps only the matching rows", {
  results <- expand.grid(
    shape_parameter = c(0.01, 1, 20),
    equivalence_margin = c(0.1, 0.5),
    drift = 1:3
  )

  kept <- filter_to_parameter_combinations(
    results, data.frame(shape_parameter = 20, equivalence_margin = 0.5)
  )

  expect_equal(nrow(kept), 3)
  expect_equal(unique(kept$shape_parameter), 20)
  expect_equal(unique(kept$equivalence_margin), 0.5)
})

test_that("a combination that matches nothing yields no rows rather than everything", {
  results <- data.frame(power_parameter = c(0, 0.5), drift = 1:2)

  kept <- filter_to_parameter_combinations(
    results, data.frame(power_parameter = 0.25)
  )

  expect_equal(nrow(kept), 0)
})
