# =============================================================================
# FIXED: R/summarize_bootstrap_results.R - summarize_bootstrap_results() docs
# =============================================================================
# 
# This file contains the COMPLETE corrected roxygen2 documentation block for
# the summarize_bootstrap_results() function. The fix escapes curly braces in
# the sgharm parameter example.
#
# INSTRUCTIONS:
# 1. Open R/summarize_bootstrap_results.R
# 2. Find the summarize_bootstrap_results() function
# 3. Replace the roxygen2 documentation block with the content below
#
# KEY CHANGE:
#   OLD: (e.g., c("{age>=50}", "{nodes>=3}"))
#   NEW: (e.g., c("\{age>=50\}", "\{nodes>=3\}"))
#
# =============================================================================

#' Enhanced Bootstrap Results Summary
#'
#' Creates comprehensive output including formatted table with subgroup footnote,
#' diagnostic plots, bootstrap quality metrics, and detailed timing analysis.
#'
#' @param sgharm The selected subgroup object from forestsearch results. Can be:
#'   \itemize{
#'     \item Character vector of factor definitions (e.g., c("\{age>=50\}", "\{nodes>=3\}"))
#'     \item List with \code{sgharm} element containing factor definitions
#'     \item List with \code{sg.harm_label} element (human-readable labels)
#'   }
#' @param boot_results List. Output from forestsearch_bootstrap_dofuture()
#' @param create_plots Logical. Generate diagnostic plots (default: FALSE)
#' @param est.scale Character. "hr" or "1/hr" for effect scale
#'
#' @return List with components:
#'   \describe{
#'     \item{table}{gt table with treatment effects and subgroup footnote}
#'     \item{diagnostics}{List of bootstrap quality metrics}
#'     \item{diagnostics_table_gt}{gt table of diagnostics}
#'     \item{plots}{List of ggplot2 diagnostic plots (if create_plots=TRUE)}
#'     \item{timing}{List of timing analysis (if timing data available)}
#'     \item{subgroup_summary}{List from summarize_bootstrap_subgroups()}
#'   }
#'
#' @details
#' The \code{table} output includes a footnote displaying the identified subgroup
#' definition, analogous to the \code{tab_estimates} table from \code{\link{sg_tables}}.
#' This is achieved by extracting the subgroup definition from \code{sgharm} and
#' passing it to \code{\link{format_bootstrap_table}}.
#'
#' @seealso
#' \code{\link{format_bootstrap_table}} for table creation
#' \code{\link{sg_tables}} for analogous main analysis tables
#' \code{\link{summarize_bootstrap_subgroups}} for subgroup stability analysis
#'
#' @export
