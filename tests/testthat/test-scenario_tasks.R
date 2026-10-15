# The parallel simulation hands a binary endpoint's scenarios to the workers by
# design, so that one worker's inference cache serves every drift of it. These
# tests pin the grouping: what counts as one design, that every scenario is
# handed out exactly once, and that a method with few designs still keeps every
# worker busy.

task_cases <- function(designs = 2, drifts = 5, weights = 1) {
  grid <- expand.grid(
    drift = seq(-0.2, 0.2, length.out = drifts),
    prior_weight = seq_len(weights),
    target_sample_size_control = 40 + 10 * seq_len(designs)
  )
  grid$target_sample_size_treatment <- 2 * grid$target_sample_size_control
  grid$treatment_drift <- grid$drift
  grid$control_drift <- 0
  grid$target_treatment_effect <- 0.3 + grid$drift
  grid$case_study <- "example"
  grid$method <- "RMP"
  grid
}

test_that("a scenario table not shared across drifts gives one task per scenario", {
  cases <- task_cases()
  expect_identical(scenario_tasks(cases, FALSE, n_workers = 4), as.list(seq_len(nrow(cases))))
})

test_that("the drifts of one design and parameter setting form one task", {
  cases <- task_cases(designs = 3, drifts = 5, weights = 4)
  tasks <- scenario_tasks(cases, TRUE, n_workers = 4)

  expect_length(tasks, 12)
  for (task in tasks) {
    design <- unique(cases[task, c(target_size_columns, "prior_weight")])
    expect_equal(nrow(design), 1)
    expect_length(task, 5)
  }
})

test_that("every scenario is handed out exactly once, in order within its task", {
  cases <- task_cases(designs = 2, drifts = 7, weights = 3)
  cases <- cases[sample(nrow(cases)), ]
  rownames(cases) <- NULL
  tasks <- scenario_tasks(cases, TRUE, n_workers = 11)

  expect_setequal(unlist(tasks), seq_len(nrow(cases)))
  expect_length(unlist(tasks), nrow(cases))
  for (task in tasks) expect_false(is.unsorted(task))
})

test_that("fewer designs than workers are cut into pieces that keep the workers busy", {
  cases <- task_cases(designs = 2, drifts = 33)
  tasks <- scenario_tasks(cases, TRUE, n_workers = 11)

  # Two designs over eleven workers: each is cut in six.
  expect_length(tasks, 12)
  expect_true(all(lengths(tasks) %in% c(5, 6)))

  # A design with fewer drifts than pieces gives one piece per drift.
  few <- task_cases(designs = 1, drifts = 3)
  expect_length(scenario_tasks(few, TRUE, n_workers = 11), 3)
})

test_that("only binary endpoints share analyses across drifts", {
  expect_true(analyses_shared_across_drifts(list(endpoint = "binary")))
  expect_false(analyses_shared_across_drifts(list(endpoint = "continuous")))
  expect_false(analyses_shared_across_drifts(list(endpoint = "time_to_event")))
})
