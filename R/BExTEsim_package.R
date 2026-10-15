## BExTEsim imports the whole BExTE namespace: its pipeline calls the R6
## method and data classes, and the BExTE functions exported for extension
## packages (see ?BExTE::effective_arm_size and the other functions with
## keyword "internal").
##
## Plotting, table and report code uses packages from Suggests. Each entry
## point calls check_plot_packages(), check_table_packages() or
## check_parallel_packages() before it uses them.
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import BExTE
#' @importFrom dplyr %>% across all_of any_of arrange case_when distinct
#' @importFrom dplyr everything first group_by left_join mutate mutate_if
#' @importFrom dplyr n_distinct pull rename_with row_number rowwise select
#' @importFrom dplyr summarize ungroup where
#' @importFrom tidyr pivot_wider unite unnest
#' @importFrom rlang := sym
#' @importFrom stats approx binom.test pnorm qnorm rbinom rchisq rnorm setNames
#' @importFrom stats quantile update
#' @importFrom utils read.csv setTxtProgressBar txtProgressBar write.csv
#' @importFrom utils write.table
#' @importFrom grDevices cairo_pdf dev.off
## usethis namespace: end
NULL
