test_that("get_parallel_worker_count always returns a valid cluster size", {
  expect_identical(
    get_parallel_worker_count(detected_cores = NA_integer_),
    1L
  )
  expect_identical(
    get_parallel_worker_count(detected_cores = 1L),
    1L
  )
  expect_identical(
    get_parallel_worker_count(detected_cores = 8L),
    7L
  )
  expect_identical(
    get_parallel_worker_count(max_workers = 8L, detected_cores = 32L),
    8L
  )
})

test_that("BEXTE_MAX_WORKERS caps the worker count", {
  expect_identical(
    get_parallel_worker_count(detected_cores = 12L, env_max_workers = "6"),
    6L
  )
  # The tighter of the two caps wins.
  expect_identical(
    get_parallel_worker_count(max_workers = 4L, detected_cores = 12L,
                              env_max_workers = "6"),
    4L
  )
  # A cap above the machine's size changes nothing.
  expect_identical(
    get_parallel_worker_count(detected_cores = 12L, env_max_workers = "64"),
    11L
  )
  expect_identical(
    get_parallel_worker_count(detected_cores = 12L, env_max_workers = ""),
    11L
  )
})

test_that("an invalid BEXTE_MAX_WORKERS is an error, not a silent default", {
  expect_error(
    get_parallel_worker_count(detected_cores = 12L, env_max_workers = "six"),
    "BEXTE_MAX_WORKERS"
  )
  expect_error(
    get_parallel_worker_count(detected_cores = 12L, env_max_workers = "0"),
    "BEXTE_MAX_WORKERS"
  )
})

test_that("limit_mcmc_chain_parallelism runs chains sequentially when scenarios are parallel", {
  mcmc_config <- list(num_chains = 4L, parallel_chains = 4L, chain_length = 5000L)

  adjusted <- limit_mcmc_chain_parallelism(mcmc_config, parallelization = TRUE)

  expect_identical(adjusted$parallel_chains, 1L)
  # Every other setting is left alone.
  expect_identical(adjusted$num_chains, 4L)
  expect_identical(adjusted$chain_length, 5000L)
})

test_that("limit_mcmc_chain_parallelism leaves the configuration alone when scenarios are sequential", {
  mcmc_config <- list(num_chains = 4L, parallel_chains = 4L, chain_length = 5000L)

  expect_identical(
    limit_mcmc_chain_parallelism(mcmc_config, parallelization = FALSE),
    mcmc_config
  )
})

test_that("a logical parallelization setting applies to every method", {
  expect_true(method_runs_in_parallel(TRUE, "PDCCPP"))
  expect_true(method_runs_in_parallel(TRUE, "commensurate_power_prior"))
  expect_false(method_runs_in_parallel(FALSE, "PDCCPP"))
  expect_false(method_runs_in_parallel(FALSE, "commensurate_power_prior"))
})


test_that("a list of method names parallelises only those methods", {
  # yaml::read_yaml() returns a bare character vector for a sequence of
  # scalars, and a list once any entry is itself a mapping, so both shapes
  # have to resolve the same way.
  for (setting in list("PDCCPP", c("PDCCPP", "NPP"), list("PDCCPP"))) {
    expect_true(method_runs_in_parallel(setting, "PDCCPP"))
    expect_false(method_runs_in_parallel(setting, "commensurate_power_prior"))
  }
})


test_that("an absent parallelization setting keeps every method sequential", {
  expect_false(method_runs_in_parallel(NULL, "PDCCPP"))
})


test_that("an unusable parallelization setting is rejected rather than guessed", {
  # Silently treating these as FALSE would turn a typo into a run that is
  # thirty times slower than intended, with nothing in the output to say so.
  expect_error(method_runs_in_parallel(NA, "PDCCPP"), "TRUE or FALSE")
  expect_error(method_runs_in_parallel(c(TRUE, FALSE), "PDCCPP"), "TRUE or FALSE")
  expect_error(method_runs_in_parallel(list(), "PDCCPP"), "TRUE or FALSE")
  expect_error(method_runs_in_parallel(42, "PDCCPP"), "TRUE or FALSE")
})


test_that("every shipped configuration resolves to a usable setting", {
  # The per-method list form is exercised above; what matters here is that no
  # shipped configuration carries a value the resolver would reject at the
  # start of a multi-day run.
  configs <- list.files(
    system.file("conf", package = "BExTEsim"),
    pattern = "^scenarios_config[.]yml$", recursive = TRUE, full.names = TRUE
  )
  expect_gt(length(configs), 0)

  for (config in configs) {
    setting <- yaml::read_yaml(config)$parallelization
    resolved <- method_runs_in_parallel(setting, "PDCCPP")
    expect_true(is.logical(resolved) && length(resolved) == 1L && !is.na(resolved))
  }
})


test_that("chain parallelism is capped only for a method that runs in parallel", {
  mcmc_config <- list(num_chains = 4L, parallel_chains = 4L)
  setting <- "PDCCPP"

  capped <- limit_mcmc_chain_parallelism(
    mcmc_config, method_runs_in_parallel(setting, "PDCCPP")
  )
  untouched <- limit_mcmc_chain_parallelism(
    mcmc_config, method_runs_in_parallel(setting, "commensurate_power_prior")
  )

  expect_equal(capped$parallel_chains, 1L)
  expect_equal(untouched$parallel_chains, 4L)
})

## The simulation loops pass a per-scenario callback as `.options.snow$progress`:
## it drives the console progress bar and the progress file BExTE-app's Run
## tab reads. Only doSNOW honours that option - doParallel warns "ignoring
## unrecognized snow option(s): progress" and never calls it, which leaves both
## readouts frozen until a whole case study/method block is finished. Which
## backend gets registered is therefore part of the contract, not an
## implementation detail.
test_that("the registered parallel backend reports each task as it finishes", {
  `%dopar%` <- foreach::`%dopar%`

  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  register_parallel_backend(cl)

  seen <- integer(0)
  opts <- list(progress = function(n) seen <<- c(seen, n))
  invisible(
    foreach::foreach(i = 1:6, .combine = c, .options.snow = opts) %dopar% i
  )

  expect_equal(seen, 1:6)
})
