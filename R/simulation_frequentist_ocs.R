#' Simulate a Scenario
#'
#' @description This function simulates a scenario based on the provided inputs and returns the results of the simulation.
#'
#' @param scenario The scenario to simulate.
#' @param simulation_config Simulation configuration
#' @param scenarios_config Configuration for the simulation.
#' @param freq_filename File name for the frequentist OCs
#' @param config_dir Configuration files directory
#' @param case_studies_config_dir Case studies configurations directory
#' @param frequentist_metrics List of frequentist metrics, as defined by
#'   sourcing `metrics_config.R`.
#' @param inference_metrics List of inference metrics, as defined by sourcing
#'   `metrics_config.R`.
#' @return The results of the simulation.
#' @keywords internal
frequentist_ocs_scenario_simulation <- function(scenario,
                                                simulation_config,
                                                scenarios_config,
                                                freq_filename,
                                                config_dir,
                                                case_studies_config_dir,
                                                frequentist_metrics,
                                                inference_metrics) {
  withr::local_seed(simulation_config$seed)

  target_sample_size_control <- scenario$target_sample_size_control[[1]]

  target_sample_size_treatment <- scenario$target_sample_size_treatment[[1]]
  treatment_drift <- scenario$treatment_drift[[1]]
  control_drift <- scenario$control_drift[[1]]
  case_study <- scenario$case_study[[1]]
  method <- scenario$method[[1]]
  method_parameters <- scenario$parameters[[1]]
  source_denominator <- scenario$source_denominator[[1]]
  target_to_source_std_ratio <- scenario$target_to_source_std_ratio[[1]]
  dropout_probability <- scenario$dropout_probability[[1]]
  event_time_distribution <- scenario$event_time_distribution[[1]]
  # Absent from scenarios written before the delayed-effect axis existed.
  treatment_delay <- if (is.null(scenario$treatment_delay)) 0 else scenario$treatment_delay[[1]]

  case_study_config <- yaml::read_yaml(file.path(case_studies_config_dir, paste0(case_study, ".yml")))
  mcmc_config <- read_config(paste0(config_dir, "/mcmc_config.yml"), mcmc_config_schema)
  # `parallelization` may name the methods that run in parallel rather than
  # apply to all of them, so it has to be resolved for this scenario's method
  # before anything keys off it.
  runs_in_parallel <- method_runs_in_parallel(
    scenarios_config$parallelization,
    method
  )
  mcmc_config <- limit_mcmc_chain_parallelism(mcmc_config, runs_in_parallel)

  source_data <- SourceData$new(case_study_config, source_denominator)

  model <- Model$new()

  model <- model$create(
    case_study_config = case_study_config,
    method = method,
    method_parameters = method_parameters,
    source_data = source_data,
    mcmc_config = mcmc_config
  )

  theta_0 <- case_study_config$theta_0
  critical_value <- simulation_config$critical_value
  n_replicates <- case_study_n_replicates(scenarios_config, case_study)
  exact_enumeration <- case_study_exact_enumeration(scenarios_config, case_study)
  if (exact_enumeration) {
    check_exact_enumeration_supported(case_study_config, case_study)
  }

  summary_measure_likelihood <- case_study_config$summary_measure_likelihood
  endpoint <- case_study_config$endpoint

  sampling_approximation <- case_study_config$sampling_approximation

  # Instantiate a target data object
  target_data <- TargetDataFactory$new()

  target_data <- target_data$create(
    source_data = source_data,
    case_study_config = case_study_config,
    target_sample_size_control = target_sample_size_control,
    target_sample_size_treatment = target_sample_size_treatment,
    control_drift = control_drift,
    treatment_drift = treatment_drift,
    summary_measure_likelihood = summary_measure_likelihood,
    target_to_source_std_ratio = target_to_source_std_ratio,
    dropout_probability = dropout_probability,
    event_time_distribution = event_time_distribution,
    treatment_delay = treatment_delay
  )

  # A method whose hyperparameters follow from the design fixes them here:
  # after the source data and the scenario are known, and before any
  # replicate is generated, so that every replicate of this scenario is
  # analysed under the same prior. Methods taking their hyperparameters from
  # the configuration grid are unaffected - the base class does nothing.
  model$calibrate_for_design(target_data)

  # Get the current RNG state
  computation_state <- list(rng_state = .Random.seed, computation_time = 0)

  start_time <- Sys.time()

  null_space <- case_study_config$null_space

  # Estimate the frequentist OCs for the model in the scenario considered
  sim_outputs <- model$estimate_frequentist_operating_characteristics(
    theta_0 = theta_0,
    target_data = target_data,
    n_replicates = n_replicates,
    critical_value = critical_value,
    confidence_level = simulation_config$confidence_level,
    null_space = null_space,
    n_samples_quantiles_estimation = simulation_config$n_samples_quantiles_estimation,
    case_study = case_study,
    method = method,
    simulation_config = simulation_config,
    exact_enumeration = exact_enumeration
  )

  end_time <- Sys.time()

  computation_state$computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  sim_outputs$posterior_parameters <- format_parameters_to_json(sim_outputs$posterior_parameters)
  scenario$parallelization <- runs_in_parallel

  check_colnames(scenario, expected_colnames_scenario)

  check_colnames(data.frame(sim_outputs), expected_colnames_results)

  source_data_df <- data.frame(source_data$to_dict())

  check_colnames(source_data_df, expected_colnames_source)

  target_data_df <- data.frame(lapply(target_data$to_dict(), function(x) {
    if (is.null(x)) {
      return(NA)
    } else {
      return(x)
    }
  }))

  expected_colnames_target <- c(
    "summary_measure_likelihood",
    "target_sample_size_control",
    "target_sample_size_treatment",
    "target_treatment_effect",
    "target_standard_deviation",
    "target_treatment_rate",
    "target_control_rate"
  )

  check_colnames(target_data_df, expected_colnames_target)

  additional_columns <- c("sampling_approximation", "rng_state", "computation_time")

  results_frequentist_ocs <- cbind(
    scenario,
    sim_outputs,
    source_data_df,
    target_data_df,
    sampling_approximation,
    theta_0,
    null_space,
    rng_state = jsonlite::toJSON(computation_state$rng_state, pretty = TRUE),
    computation_time = computation_state$computation_time
  )
  results_frequentist_ocs$source_denominator <- unlist(results_frequentist_ocs$source_denominator)
  results_frequentist_ocs <- results_frequentist_ocs[!duplicated(names(results_frequentist_ocs), fromLast = TRUE)]

  results_frequentist_ocs$parameters <- format_parameters_to_json(results_frequentist_ocs$parameters, escape = TRUE)


  for (col in names(results_frequentist_ocs)) {
    # Check if the column is a list
    if (is.list(results_frequentist_ocs[[col]])) {
      # Attempt to simplify to either numeric or character, depending on content
      if (is.numeric(unlist(results_frequentist_ocs[[col]][1]))) {
        results_frequentist_ocs[[col]] <- as.numeric(unlist(results_frequentist_ocs[[col]]))
      } else if (is.character(unlist(results_frequentist_ocs[[col]][1]))) {
        results_frequentist_ocs[[col]] <- as.character(unlist(results_frequentist_ocs[[col]]))
      } else {
        results_frequentist_ocs[[col]] <- unlist(results_frequentist_ocs[[col]])
      }
    }
  }

  # Reorder the columns to ease reading
  # Combine all expected column names into a single vector in the desired order
  expected_colnames_order <- c(
    expected_colnames_scenario,
    expected_colnames_target,
    expected_colnames_source,
    additional_columns,
    expected_colnames_results
  )
  expected_colnames_order <- unique(expected_colnames_order)

  results_frequentist_ocs <- results_frequentist_ocs[, expected_colnames_order, drop = FALSE]

  check_confidence_intervals(results_frequentist_ocs, c(names(inference_metrics), names(frequentist_metrics)))


  # When this method's scenarios run sequentially, the csv files containing the results are built in an iterative manner, so that in case of bugs results are not lost (not possible if parallelization is used).
  if (!runs_in_parallel) {
    if (!file.exists(freq_filename)) {
      # Get the column names
      fieldnames_freq <- colnames(results_frequentist_ocs)

      # Create an empty data frame with these column names
      empty_df <- data.frame(matrix(ncol = length(fieldnames_freq), nrow = 0))
      colnames(empty_df) <- fieldnames_freq

      # Create an empty data frame with these column names
      write.table(
        empty_df,
        file = freq_filename,
        sep = ",",
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
      )
    }

    write.table(
      results_frequentist_ocs,
      file = freq_filename,
      sep = ",",
      col.names = FALSE,
      row.names = FALSE,
      append = TRUE,
      quote = TRUE
    )
  }

  return(results_frequentist_ocs)
}

## The sample size factors to simulate a given case study at.
##
## scenarios_config$case_study_sample_size_factors, when present, maps a case
## study to the factors it is wanted at, so a run can avoid the cross product
## of every case study with every factor. A case study the map does not
## mention keeps the run-wide sample_size_factors, which is what every config
## written before the field existed does.
case_study_factors <- function(scenarios_config, case_study) {
  per_case_study <- scenarios_config$case_study_sample_size_factors
  if (is.null(per_case_study) || is.null(per_case_study[[case_study]])) {
    return(scenarios_config$sample_size_factors)
  }
  return(per_case_study[[case_study]])
}

## The ratio of the treatment to the control arm size of a case study's
## target study.
##
## scenarios_config$case_study_target_allocation_ratio, when present, maps a
## case study to its own ratio. Otherwise the run-wide target_allocation_ratio
## applies, and a config that sets neither has equal arms.
target_allocation_ratio <- function(scenarios_config, case_study) {
  per_case_study <- scenarios_config$case_study_target_allocation_ratio
  ratio <- if (!is.null(per_case_study) && !is.null(per_case_study[[case_study]])) {
    as.numeric(unlist(per_case_study[[case_study]]))
  } else if (!is.null(scenarios_config$target_allocation_ratio)) {
    scenarios_config$target_allocation_ratio
  } else {
    1
  }
  if (length(ratio) != 1 || !is.finite(ratio) || ratio <= 0) {
    stop("The target allocation ratio of ", case_study,
         " should be a single positive number.", call. = FALSE)
  }
  return(ratio)
}

## The number of replicates to simulate a given case study with.
##
## scenarios_config$case_study_n_replicates, when present, maps a case study to
## its own replicate count, so one whose analysis is expensive - a binomial
## likelihood fitted by MCMC - can be run at fewer replicates than the rest. A
## case study the map does not mention keeps the run-wide n_replicates.
case_study_n_replicates <- function(scenarios_config, case_study) {
  per_case_study <- scenarios_config$case_study_n_replicates
  if (is.null(per_case_study) || is.null(per_case_study[[case_study]])) {
    return(scenarios_config$n_replicates)
  }
  return(as.numeric(unlist(per_case_study[[case_study]])))
}

#' Run simulations based on the given environment
#'
#' @param env The environment to run the simulations in
#' @param simulation_config Simulation configuration
#' @param scenarios_config Configuration for the simulation.
#' @param analysis_config Analysis configuration
#' @param config_dir Configuration directory
#' @param case_studies_config_dir Case studies configuration directory
#' @param logging_file_path Currently unused.
#' @param frequentist_metrics List of frequentist metrics, as defined by
#'   sourcing `metrics_config.R`.
#' @param inference_metrics List of inference metrics, as defined by sourcing
#'   `metrics_config.R`.
#' @return None
#' @keywords internal
simulation_frequentist_ocs <- function(env,
                                       simulation_config,
                                       scenarios_config,
                                       analysis_config,
                                       config_dir,
                                       case_studies_config_dir,
                                       logging_file_path,
                                       frequentist_metrics,
                                       inference_metrics){
  case_studies <- scenarios_config$case_studies
  methods <- scenarios_config$methods

  if (length(case_studies) == 0 || length(methods) == 0) {
    stop("Select at least one case study and one method to run the simulation.")
  }

  # Create scenarios and create a summary table
  cases <- simulation_scenarios(config_dir = config_dir, scenarios_config = scenarios_config, case_studies_config_dir = case_studies_config_dir)

  scenarios_table_ranges(cases, paste0("./results/", env))

  # The loops below narrow `cases` to one case study and method at a time;
  # the tracker keeps a run-wide total across them, which is what
  # BExTE-app's Run tab reads to draw its progress bar. The total is summed per
  # case study rather than taken from the frame above, because a run that
  # restricts a case study's factors simulates fewer scenarios than the
  # full cross product - counting those would leave the bar short of its
  # own total for the whole run.
  run_total <- sum(vapply(case_studies, function(cs) {
    per_case_study <- scenarios_config
    per_case_study$case_studies <- cs
    per_case_study$sample_size_factors <- case_study_factors(scenarios_config, cs)
    return(nrow(simulation_scenarios(
      config_dir = config_dir, scenarios_config = per_case_study,
      case_studies_config_dir = case_studies_config_dir
    )))
  }, numeric(1)))
  run_progress <- run_progress_tracker(env, total = run_total)

  results_dir <- paste0("./results/", env, "/", "frequentist")

  # Replicates are generated once per scenario and shared by all its methods -
  # see generate_replicates(). The directory is emptied first, so that nothing
  # an earlier run generated, possibly with other code, is ever reused.
  generation_cache_dir <- file.path(
    normalizePath(paste0("./results/", env), mustWork = FALSE), "generated"
  )
  unlink(generation_cache_dir, recursive = TRUE)
  simulation_config$generation_cache_dir <- generation_cache_dir

  if (!(any("separate" %in% methods))) {
    warning("Separate analysis not included in the methods.")
    futile.logger::flog.warn("Separate analysis not included in the methods.")
  }

  for (case_study in case_studies) {
    for (method in methods) {
      scenarios_config$case_studies <- case_study
      scenarios_config$methods <- method
      # Only the factors this case study is actually wanted at. Without
      # this the run takes the cross product of every case study with
      # every factor, and simulates combinations nothing asked for - the
      # paper's figures, for instance, use each case study at two of the
      # three factors, so a third of the grid was being computed and
      # never plotted. Case studies absent from the map are unrestricted.
      scenarios_config$sample_size_factors <- case_study_factors(
        scenarios_config, case_study
      )

      # Scenarios of one case study and method are what share analyses, so the
      # cache is emptied here rather than growing for the whole run. The keys
      # already keep other scenarios apart; this bounds the memory it holds.
      inference_cache_reset()

      cases <- simulation_scenarios(config_dir = config_dir, scenarios_config = scenarios_config, case_studies_config_dir = case_studies_config_dir)

      if (nrow(cases) == 0) {
        stop("No simulation results")
      }

      outputs_config <- yaml::read_yaml(system.file("conf/outputs_config.yml", package = "BExTEsim"))
      freq_filename <- paste0(
        results_dir,
        "/",
        case_study,
        "/",
        method,
        "/",
        outputs_config$frequentist_ocs_results_filename
      )

      if (nrow(cases) == 0) {
        stop("No simulation results after filtering")
      }

      # Progress bar function in R
      pb <- progress::progress_bar$new(format = "Progress : [:bar] :elapsed | eta: :eta",
                                       total = nrow(cases),
                                       width = 60)
      run_progress$starting(case_study, method)
      # foreach calls this with the number of tasks finished so far in this
      # block, the sequential loop below with its own index - the same count
      # either way, which is what the tracker expects.
      progress <- function(n) {
        pb$tick()
        return(run_progress$tick(n))
      }
      opts <- list(progress = progress)

      # Create result folder if does not exist
      dir.create(
        paste0(results_dir, "/", case_study, "/", method),
        showWarnings = FALSE,
        recursive = TRUE
      )

      # Delete existing results file. If they exist and are open, close them and delete them.
      if (file.exists(freq_filename)) {
        if (simulation_config$delete_old_results){
          freq_file <- file(freq_filename, "w")
          if (isOpen(freq_file)) {
            close(freq_file)
          }
          file.remove(freq_filename)
        } else {
          next
        }
      }

      # Temporarily disable the global error handler within the tryCatch block. Otherwise there will be interference between the tryCatch block and the global error handler.
      withr::local_options(error = NULL)
      if (!method_runs_in_parallel(scenarios_config$parallelization, method)) {
        # Loop through cases

        for (i in seq_len(nrow(cases))) {
          scenario <- cases[i, ]

          tryCatch({
            results <- frequentist_ocs_scenario_simulation(
              scenario = scenario,
              simulation_config = simulation_config,
              scenarios_config = scenarios_config,
              freq_filename = freq_filename,
              config_dir = config_dir,
              case_studies_config_dir = case_studies_config_dir,
              frequentist_metrics = frequentist_metrics,
              inference_metrics = inference_metrics
            )
            progress(i)
          }, error = function(e) {
            futile.logger::flog.error(paste("Error at iteration", i, ":", e$message))
            futile.logger::flog.error("Call stack:\n%s", paste(deparse(sys.calls()), collapse = "\n"))

            # Save state
            save_state(iteration = i, scenario = scenario, env = env)
            stop(e)
          })
        }
        # At the end of the loop, close file if open
        freq_file <- file(freq_filename, "a")
        if (isOpen(freq_file)) {
          close(freq_file)
        }
      } else {
        ncores <- get_parallel_worker_count(max_workers = 96L)

        # Set a range of ports for cluster
        withr::local_options(clusterPort = c(11000, 11999))

        # Register parallel backend. It has to be one that calls the
        # `opts` progress callback built above as each scenario comes back -
        # see register_parallel_backend().
        check_parallel_packages()
        `%dopar%` <- foreach::`%dopar%`
        cl <- parallel::makeCluster(ncores)
        register_parallel_backend(cl)

        # Export the library paths to each worker
        paths <- .libPaths()
        # BExTEsim may not be an installed package at all (e.g. BExTE-app's
        # dev-mode background process only ever `pkgload::load_all()`s it),
        # so each worker needs the same source path to fall back to.
        pkg_root <- find.package("BExTEsim")
        parallel::clusterExport(
          cl,
          varlist = c(
            "paths",
            "pkg_root",
            "frequentist_metrics",
            "inference_metrics",
            "simulation_config",
            "freq_filename",
            "config_dir",
            "case_studies_config_dir"
          ),
          envir = environment()
        )

        # Evaluate the expression to load BExTEsim in each worker: library() for
        # an installed package (the common case), falling back to
        # pkgload::load_all() when it is only loaded from source.
        parallel::clusterEvalQ(cl, {
          .libPaths(paths)
          if (!requireNamespace("BExTEsim", quietly = TRUE)) {
            pkgload::load_all(pkg_root, quiet = TRUE)
          } else {
            library(BExTEsim)
          }
        })

        # On a binary endpoint the drifts of one design go to the same worker,
        # whose inference cache then serves all of them - see scenario_tasks().
        # Each task carries its own rows, so the workers need not be sent the
        # whole scenario table.
        case_study_config <- yaml::read_yaml(file.path(case_studies_config_dir, paste0(case_study, ".yml")))
        task_rows <- scenario_tasks(
          cases,
          shares_across_drifts = analyses_shared_across_drifts(case_study_config),
          n_workers = ncores
        )
        tasks <- lapply(task_rows, function(rows) {
          return(list(rows = rows, cases = cases[rows, , drop = FALSE]))
        })

        # doSNOW reports how many tasks have finished, so the scenarios done are
        # counted in proportion; they reach the block's total with the last task.
        task_progress <- function(n) {
          done <- if (n >= length(tasks)) nrow(cases) else floor(n * nrow(cases) / length(tasks))
          pb$update(done / nrow(cases))
          return(run_progress$tick(done))
        }

        # Parallel computation over cases
        results_by_task <- foreach::foreach(
          task = tasks,
          .options.snow = list(progress = task_progress)
        ) %dopar% {
          worker_id <- Sys.getpid()  # Get the process ID of the worker

          lapply(seq_along(task$rows), function(j) {
            i <- task$rows[j]
            scenario <- task$cases[j, ]

            return(tryCatch({
              result <- frequentist_ocs_scenario_simulation(
                scenario = scenario,
                simulation_config = simulation_config,
                scenarios_config = scenarios_config,
                freq_filename = freq_filename,
                config_dir = config_dir,
                case_studies_config_dir = case_studies_config_dir,
                frequentist_metrics = frequentist_metrics,
                inference_metrics = inference_metrics
              )
              futile.logger::flog.info(paste("Iteration", i, " completed by worker", worker_id))
              result  # Return the result (assuming it’s a data frame)
            }, error = function(e) {
              futile.logger::flog.error(paste("Error at iteration", i, "by worker", worker_id, ":", e$message))
              futile.logger::flog.error("Call stack:\n%s", paste(deparse(sys.calls()), collapse = "\n"))

              # Save state even if there's an error
              save_state(iteration = i, scenario = scenario, worker_id = worker_id, env = env)

              stop(e)
            }))
          })
        }
        # Back into scenario order, as the one-scenario-per-task loop combined
        # them, so the results file does not depend on how tasks were formed.
        per_scenario <- unlist(results_by_task, recursive = FALSE)
        results <- do.call(rbind, per_scenario[order(unlist(task_rows))])
        write_csv(results, freq_filename)
        parallel::stopCluster(cl)
      }

      # Restore the global error handler
      withr::local_options(error = global_error_handler)

      if (tolower(case_study) == "aprepitant" ||
          is_commensurate_method(method)) {
        # remove remaining MCMC csv files related to the run
        clear_stan_draws(case_study, method)
      }
    }

    # Every method of this case study has run, so its replicates will not be
    # read again.
    unlink(generation_cache_dir, recursive = TRUE)
  }
  return(invisible(NULL))
}
