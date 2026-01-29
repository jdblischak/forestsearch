# =============================================================================
# FIXED: R/bootstrap_summaries_helpers.R - format_bootstrap_table() documentation
# =============================================================================
# 
# This file contains the COMPLETE corrected roxygen2 documentation block for
# the format_bootstrap_table() function. The fix escapes curly braces in the
# sg_definition parameter example.
#
# INSTRUCTIONS:
# 1. Open R/bootstrap_summaries_helpers.R
# 2. Find the format_bootstrap_table() function
# 3. Replace the roxygen2 documentation block with the content below
#
# KEY CHANGE:
#   OLD: (e.g., "{age>=50} & {nodes>=3}")
#   NEW: (e.g., "\{age>=50\} & \{nodes>=3\}")
#
# =============================================================================

#' Format Bootstrap Results Table with gt
#'
#' Creates a publication-ready table from ForestSearch bootstrap results,
#' with bias-corrected confidence intervals, informative formatting, and
#' optional subgroup definition footnote.
#'
#' @param FSsg_tab Data frame or matrix from forestsearch_bootstrap_dofuture()$FSsg_tab
#' @param nb_boots Integer. Number of bootstrap iterations performed
#' @param est.scale Character. "hr" or "1/hr" for effect scale
#' @param boot_success_rate Numeric. Proportion of bootstraps that found subgroups
#' @param sg_definition Character. Subgroup definition string to display as footnote
#'   (e.g., "\{age>=50\} & \{nodes>=3\}"). If NULL, no subgroup footnote is added.
#' @param title Character. Custom title (optional)
#' @param subtitle Character. Custom subtitle (optional)
#'
#' @return A gt table object
#'
#' @importFrom gt gt tab_header tab_spanner tab_footnote tab_source_note md
#'   cols_label tab_style cell_fill cell_text cells_body cells_column_labels
#' @importFrom dplyr all_of
#' @export
