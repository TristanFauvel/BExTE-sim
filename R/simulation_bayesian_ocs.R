#' Simulate a Scenario
#'
#' @description This function simulates a scenario based on the provided inputs and returns the results of the simulation.
#'
#' @param scenario The scenario to simulate.
#' @param simulation_config Simulation configuration
#' @param scenarios_config Scenario-generation configuration
#' @param bayes_filename Results were results are stored for Bayesian OCs
#' @param config_dir Configuration files directory
#' @param case_studies_config_dir Case-study configuration directory
#' @return The results of the simulation.
#' @keywords internal
bayesian_ocs_scenario_simulation <- function(scenario,
                                             simulation_config,
                                             scenarios_config,
                                             bayes_filename,
                                             config_dir,
                                             case_studies_config_dir) {
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
  # As in the frequentist driver: `parallelization` may name the methods that
  # run in parallel rather than apply to all of them.
  runs_in_parallel <- method_runs_in_parallel(
    scenarios_config$parallelization,
    method
  )
  mcmc_config <- limit_mcmc_chain_parallelism(mcmc_config, runs_in_parallel)

  source_data <- SourceData$new(case_study_config = case_study_config,
                                source_denominator = source_denominator)

  model <- Model$new()

  model <- model$create(
    case_study_config = case_study_config,
    method = method,
    method_parameters = method_parameters,
    source_data = source_data,
    mcmc_config = mcmc_config
  )

  # The Bayesian operating characteristics draw a target treatment effect
  # from the design prior for each sample, and
  # estimate_bayesian_operating_characteristics() builds a target data object
  # per draw, with a drift of its own. Calibrating from one of those would
  # give a different prior for every draw. The prior belongs to the scenario,
  # so it is calibrated once here, from the scenario's own design, and held
  # fixed across the draws and across the three design prior types.
  model$calibrate_for_design(
    TargetDataFactory$new()$create(
      source_data = source_data,
      case_study_config = case_study_config,
      target_sample_size_control = target_sample_size_control,
      target_sample_size_treatment = target_sample_size_treatment,
      control_drift = control_drift,
      treatment_drift = treatment_drift,
      summary_measure_likelihood = case_study_config$summary_measure_likelihood,
      target_to_source_std_ratio = target_to_source_std_ratio,
      dropout_probability = dropout_probability,
      event_time_distribution = event_time_distribution,
      treatment_delay = treatment_delay
    )
  )

  theta_0 <- case_study_config$theta_0
  critical_value <- simulation_config$critical_value
  n_replicates <- scenarios_config$n_replicates

  summary_measure_likelihood <- case_study_config$summary_measure_likelihood
  endpoint <- case_study_config$endpoint

  json_parameters <- format_parameters_to_json(scenario$parameters)

  # Estimate the Bayesian OCs for the model in the scenario considered
  design_prior_types <- c("ui_design_prior", "analysis_prior", "source_posterior")

  # Get the current RNG state
  computation_state <- list(
    rng_state = jsonlite::toJSON(.Random.seed),
    computation_time = 0
  )

  for (design_prior_type in design_prior_types) {
    results_bayesian_ocs <- estimate_bayesian_ocs(
      scenario = scenario,
      case_study_config = case_study_config,
      source_data = source_data,
      target_sample_size_control = target_sample_size_control,
      target_sample_size_treatment = target_sample_size_treatment,
      model = model,
      method_parameters = method_parameters,
      theta_0 = theta_0,
      n_replicates = n_replicates,
      critical_value = critical_value,
      computation_state = computation_state,
      design_prior_type = design_prior_type,
      json_parameters = json_parameters,
      simulation_config = simulation_config,
      target_to_source_std_ratio = target_to_source_std_ratio,
      mcmc_config = mcmc_config
    )

    results_bayesian_ocs$source_denominator <- unlist(results_bayesian_ocs$source_denominator)

    if (design_prior_type == "ui_design_prior") {
      results_bayesian_ocs_ui_design_prior <- results_bayesian_ocs
    } else if (design_prior_type == "analysis_prior") {
      results_bayesian_ocs_analysis_prior <- results_bayesian_ocs
    } else if (design_prior_type == "source_posterior") {
      results_bayesian_ocs_source_posterior <- results_bayesian_ocs
    }
  }

  combined_results <- rbind(
    results_bayesian_ocs_ui_design_prior,
    results_bayesian_ocs_analysis_prior,
    results_bayesian_ocs_source_posterior
  )

  # When this method's scenarios run sequentially, the csv files containing the results are built in an iterative manner, so that in case of bugs results are not lost (not possible if parallelization is used).
  if (!runs_in_parallel) {
    if (!file.exists(bayes_filename)) {
      # Get the column names
      fieldnames_bayes <- colnames(combined_results)

      # Create an empty data frame with these column names
      empty_df <- data.frame(matrix(ncol = length(fieldnames_bayes), nrow = 0))
      colnames(empty_df) <- fieldnames_bayes

      # Create an empty data frame with these column names
      write.table(
        empty_df,
        file = bayes_filename,
        sep = ",",
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
      )
    }

    write.table(
      combined_results,
      file = bayes_filename,
      sep = ",",
      col.names = FALSE,
      row.names = FALSE,
      append = TRUE,
      quote = TRUE
    )
  }
  return(combined_results)
}

#' Names of the Bayesian operating characteristics estimated by Monte Carlo
#'
#' @description Mirrors the values returned by
#'   `Model$estimate_bayesian_operating_characteristics()`, so that scenarios
#'   whose operating characteristics cannot be estimated can still be reported
#'   with the same columns.
#' @return A character vector of metric names.
#' @noRd
bayesian_ocs_metric_names <- function() {
  return(c(
    "average_tie",
    "average_power",
    "prior_proba_no_benefit",
    "prepost_proba_FP",
    "prepost_proba_TP",
    "upper_bound_proba_FP",
    "prior_proba_success"
  ))
}

#' Simulation Bayesian OCs
#'
#' @description This function estimates the Bayesian operating characteristics for the model in the scenario considered.
#' @param scenario Simulation scenario
#' @param case_study_config The case study configuration.
#' @param source_data The source data.
#' @param target_sample_size_control Sample size of the control arm in the target study.
#' @param target_sample_size_treatment Sample size of the treatment arm in the target study.
#' @param model The model.
#' @param method_parameters The method parameters.
#' @param theta_0 The true parameter value.
#' @param n_replicates The number of replicates.
#' @param critical_value The critical value.
#' @param computation_state The computation state.
#' @param design_prior_type The design prior type.
#' @param json_parameters Method parameters in json format
#' @param simulation_config Simulation configuration
#' @param target_to_source_std_ratio Ratio between the source and target study sampling standard deviation.
#' @param mcmc_config MCMC configuration
#' @return The results of the Bayesian operating characteristics estimation.
estimate_bayesian_ocs <- function(scenario,
                                  case_study_config,
                                  source_data,
                                  target_sample_size_control,
                                  target_sample_size_treatment,
                                  model,
                                  method_parameters,
                                  theta_0,
                                  n_replicates,
                                  critical_value,
                                  computation_state,
                                  design_prior_type,
                                  json_parameters,
                                  simulation_config,
                                  target_to_source_std_ratio,
                                  mcmc_config = NULL) {
  if (design_prior_type == "analysis_prior" && isTRUE(model$empirical_bayes)) {
    # The analysis prior of an empirical Bayes method is only defined once the
    # target data have been observed, so it cannot be used as a design prior.
    # Report the combination as missing, as the deterministic analysis in
    # `compute_bayesian_ocs()` does, rather than aborting the whole run.
    metric_names <- bayesian_ocs_metric_names()
    bayesian_ocs <- stats::setNames(as.list(rep(NA_real_, length(metric_names))), metric_names)
    computation_state$computation_time <- 0
  } else {
    design_prior <- DesignPrior$new()

    design_prior <- design_prior$create(
      design_prior_type = design_prior_type,
      model = model,
      source_data = source_data,
      case_study_config = case_study_config,
      simulation_config = simulation_config,
      mcmc_config = mcmc_config,
      case_study = scenario$case_study[[1]],
      # estimate_bayesian_operating_characteristics() simulates every trial with
      # no control drift, so at the source control rate; under a binomial
      # likelihood the design prior is taken given it, and its draws are
      # effects such a trial can have.
      target_control_rate = source_data$control_rate
    )
    start_time <- Sys.time()

    bayesian_ocs <- model$estimate_bayesian_operating_characteristics(
      design_prior = design_prior,
      theta_0 = theta_0,
      source_data = source_data,
      n_replicates = n_replicates,
      critical_value = critical_value,
      confidence_level = simulation_config$confidence_level,
      null_space = case_study_config$null_space,
      n_samples_design_prior = simulation_config$n_samples_design_prior,
      target_sample_size_control = target_sample_size_control,
      target_sample_size_treatment = target_sample_size_treatment,
      case_study_config = case_study_config,
      target_to_source_std_ratio = target_to_source_std_ratio,
      dropout_probability = dropout_probability,
      event_time_distribution = event_time_distribution,
      treatment_delay = treatment_delay,
      simulation_config = simulation_config,
      case_study = scenario$case_study[[1]],
      method = scenario$method[[1]],
      n_samples_quantiles_estimation = simulation_config$n_samples_quantiles_estimation
    )
    end_time <- Sys.time()
    computation_state$computation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }

  results_bayesian_ocs <- cbind(scenario,
                                design_prior = design_prior_type,
                                bayesian_ocs,
                                source_data$to_dict(),
                                computation_state)

  results_bayesian_ocs$parameters <- json_parameters
  results_bayesian_ocs <- results_bayesian_ocs[!duplicated(names(results_bayesian_ocs), fromLast = TRUE)]
  for (col in names(results_bayesian_ocs)) {
    # Check if the column is a list
    if (is.list(results_bayesian_ocs[[col]])) {
      # Attempt to simplify to either numeric or character, depending on content
      if (col == "parameters") {
        results_bayesian_ocs[[col]] <- as.character(unlist(results_bayesian_ocs[[col]]))
      } else if (is.numeric(unlist(results_bayesian_ocs[[col]][1]))) {
        results_bayesian_ocs[[col]] <- as.numeric(unlist(results_bayesian_ocs[[col]]))
      } else if (is.character(unlist(results_bayesian_ocs[[col]][1]))) {
        results_bayesian_ocs[[col]] <- as.character(unlist(results_bayesian_ocs[[col]]))
      } else {
        results_bayesian_ocs[[col]] <- unlist(results_bayesian_ocs[[col]])
      }
    }
  }

  return(results_bayesian_ocs)
}

#' Run simulations based on the given environment
#'
#' @param env The environment to run the simulations in
#' @param simulation_config Simulation configuration
#' @param scenarios_config Configuration for the simulation.
#' @param analysis_config Analysis configuration
#' @param config_dir Configuration directory
#' @param case_studies_config_dir Case studies configuration directory
#'
#' @return None
simulation_bayesian_ocs <- function(env,
                                    simulation_config,
                                    scenarios_config,
                                    analysis_config,
                                    config_dir,
                                    case_studies_config_dir) {
  results_dir <- paste0("./results/", env, "/", "bayesian")

  case_studies <- scenarios_config$case_studies
  methods <- scenarios_config$methods

  for (case_study in case_studies) {
    for (method in methods) {
      scenarios_config$case_studies <- case_study
      scenarios_config$methods <- method

      cases <- simulation_scenarios(config_dir = config_dir, scenarios_config = scenarios_config, case_studies_config_dir = case_studies_config_dir)

      if (nrow(cases) == 0) {
        stop("No simulation results")
      }

      cases <- dplyr::select(
        cases,
        -c(
          "drift",
          "control_drift",
          "treatment_drift",
          "target_treatment_effect"
        )
      )
      cases <- unique(cases)

      outputs_config <- yaml::read_yaml(system.file("conf/outputs_config.yml", package = "BExTEsim"))

      bayes_filename <- paste0(
        results_dir,
        "/",
        case_study,
        "/",
        method,
        "/",
        outputs_config$bayesian_ocs_mc_results_filename
      )

      if (nrow(cases) == 0) {
        stop("No simulation results after filtering")
      }

      withr::local_seed(simulation_config$seed)

      # Progress bar function in R
      pb <- progress::progress_bar$new(format = "Progress : [:bar] :elapsed | eta: :eta",
                                       total = nrow(cases),
                                       width = 60)
      progress <- function(n) {
        return(pb$tick())
      }
      opts <- list(progress = progress)

      # Create result folder if does not exist
      dir.create(
        paste0(results_dir, "/", case_study, "/", method),
        showWarnings = FALSE,
        recursive = TRUE
      )

      # Delete existing results file. If they exist and are open, close them and delete them.
      if (file.exists(bayes_filename)) {
        bayes_file <<- file(bayes_filename, "w")
        if (isOpen(bayes_file)) {
          close(bayes_file)
        }
        file.remove(bayes_filename)
      }

      if (!method_runs_in_parallel(scenarios_config$parallelization, method)) {
        # Loop through cases
        for (i in seq_len(nrow(cases))) {
          scenario <- cases[i, ]
          results <- bayesian_ocs_scenario_simulation(
            scenario = scenario,
            simulation_config = simulation_config,
            scenarios_config = scenarios_config,
            bayes_filename = bayes_filename,
            config_dir = config_dir,
            case_studies_config_dir = case_studies_config_dir
          )
          progress(i)
        }
        # At the end of the loop, close file if open
        bayes_file <- file(bayes_filename, "a")
        if (isOpen(bayes_file)) {
          close(bayes_file)
        }
      } else {
        ncores <- get_parallel_worker_count(max_workers = 8L)

        # Set a range of ports for cluster
        withr::local_options(clusterPort = c(11000, 11999))

        # Register parallel backend, one that calls the `opts` progress
        # callback built above as each scenario comes back - see
        # register_parallel_backend().
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
            "simulation_config",
            "bayes_filename",
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
        # Parallel computation over cases
        results <- foreach::foreach(
          i = seq_len(nrow(cases)),
          .combine = "rbind",
          .options.snow = opts
        ) %dopar% {
          scenario <- cases[i, ]
          bayesian_ocs_scenario_simulation(scenario,
                                           simulation_config,
                                           scenarios_config,
                                           bayes_filename,
                                           config_dir,
                                           case_studies_config_dir)
        }
        write_csv(results, bayes_filename)
        parallel::stopCluster(cl)
      }
    }
  }
  return(invisible(NULL))
}
