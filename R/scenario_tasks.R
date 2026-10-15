# How the parallel simulation hands its scenarios to the workers.
#
# On a binary endpoint the analysis of a trial depends on its responder counts
# alone, not on the drift that produced them, so the inference cache lets the
# scenarios of one design and method share their analyses across drifts. The
# cache lives in each worker, though, and handing the scenarios out one at a
# time spreads a design's drifts over every worker, each of which then analyses
# the design's outcomes afresh. With exact enumeration those first analyses are
# nearly all of the cost: at Aprepitant's n = 71 a method's first drift took
# 69-275 s and the later ones 3-19 s. Handing out a design's drifts together
# lets one worker's cache serve all of them.


#' Scenario columns that carry the drift
#'
#' The columns that differ between the drift scenarios of one design. Every
#' other column of the scenario table identifies the design and the method.
#'
#' @keywords internal
SCENARIO_DRIFT_COLUMNS <- c("drift", "treatment_drift", "control_drift", "target_treatment_effect")


#' Whether a case study's scenarios share analyses across drifts
#'
#' @param case_study_config The case study configuration.
#' @return `TRUE` for a binary endpoint, whose analyses depend on the responder
#'   counts alone.
#' @keywords internal
analyses_shared_across_drifts <- function(case_study_config) {
  return(identical(case_study_config$endpoint, "binary"))
}


#' Split the scenarios into the tasks handed to the workers
#'
#' @description One task per scenario, unless the scenarios share analyses
#'   across drifts: then a task holds the scenarios that differ only in the
#'   drift, so that one worker analyses the design's outcomes once and serves
#'   every drift from its cache. When there are fewer such groups than workers,
#'   each group is cut into contiguous pieces, as few as keep every worker
#'   busy, since a group on one worker would otherwise leave the others idle.
#'
#' @param cases The scenario table, one row per scenario.
#' @param shares_across_drifts Whether the scenarios share analyses across
#'   drifts - see [analyses_shared_across_drifts()].
#' @param n_workers Number of workers.
#' @return A list of integer vectors, the rows of `cases` in each task. Each
#'   row appears in exactly one task, and in increasing order within it.
#' @keywords internal
scenario_tasks <- function(cases, shares_across_drifts, n_workers) {
  rows <- seq_len(nrow(cases))
  if (!isTRUE(shares_across_drifts) || nrow(cases) == 0) {
    return(as.list(rows))
  }

  design_columns <- setdiff(names(cases), SCENARIO_DRIFT_COLUMNS)
  keys <- vapply(rows, function(i) {
    return(rlang::hash(as.list(cases[i, design_columns, drop = FALSE])))
  }, character(1))
  groups <- unname(split(rows, factor(keys, levels = unique(keys))))

  pieces <- max(1L, as.integer(ceiling(n_workers / length(groups))))
  return(unlist(lapply(groups, function(group) {
    n_pieces <- min(length(group), pieces)
    piece <- ceiling(seq_along(group) * n_pieces / length(group))
    return(unname(split(group, piece)))
  }), recursive = FALSE))
}
