# =============================================================================
# FIXED: R/summarize_bootstrap_subgroups.R - summarize_bootstrap_subgroups() docs
# =============================================================================
# 
# This file contains the COMPLETE corrected roxygen2 documentation block for
# the summarize_bootstrap_subgroups() function. The fix escapes curly braces in
# the original_sg parameter example.
#
# INSTRUCTIONS:
# 1. Open R/summarize_bootstrap_subgroups.R
# 2. Find the summarize_bootstrap_subgroups() function
# 3. Replace the roxygen2 documentation block with the content below
#
# KEY CHANGE:
#   OLD: (e.g., c("{age>=50}", "{nodes>=3}") for a 2-factor subgroup)
#   NEW: (e.g., c("\{age>=50\}", "\{nodes>=3\}") for a 2-factor subgroup)
#
# =============================================================================

#' Summarize Bootstrap Subgroup Analysis Results
#'
#' Comprehensive summary of bootstrap subgroup identification results including
#' basic statistics, factor frequencies, consistency distributions, and agreement
#' with the original analysis subgroup.
#'
#' @param results Data.table or data.frame. Bootstrap results with subgroup
#'   characteristics including columns like Pcons, hr_sg, N_sg, K_sg, and M.1-M.k
#' @param nb_boots Integer. Total number of bootstrap iterations
#' @param original_sg Character vector. Original subgroup definition from main
#'   analysis (e.g., c("\{age>=50\}", "\{nodes>=3\}") for a 2-factor subgroup)
#' @param maxk Integer. Maximum number of factors allowed in subgroup definition
#'
#' @return List with summary components:
#'   \describe{
#'     \item{basic_stats}{Data.table of summary statistics}
#'     \item{consistency_dist}{Data.table of Pcons distribution by bins}
#'     \item{size_dist}{Data.table of subgroup size distribution}
#'     \item{factor_freq}{Data.table of factor frequencies by position}
#'     \item{agreement}{Data.table of subgroup definition agreement counts}
#'     \item{factor_presence}{Data.table of base factor presence counts}
#'     \item{factor_presence_specific}{Data.table of specific factor definitions}
#'     \item{original_agreement}{Data.table comparing to original analysis subgroup}
#'     \item{n_found}{Integer. Number of successful iterations}
#'     \item{pct_found}{Numeric. Percentage of successful iterations}
#'   }
#'
#' @details
#' This function analyzes bootstrap results to assess the stability of subgroup
#' identification. Key metrics include:
#' \itemize{
#'   \item \strong{Consistency}: Distribution of Pcons values across bootstraps
#'   \item \strong{Factor frequency}: How often each factor appears in identified subgroups
#'   \item \strong{Agreement}: How often the exact same subgroup is identified
#'   \item \strong{Original agreement}: How well bootstrap results match the main analysis
#' }
#'
#' @seealso
#' \code{\link{forestsearch_bootstrap_dofuture}} for running bootstrap analysis
#' \code{\link{format_subgroup_summary_tables}} for creating gt tables from results
#' \code{\link{summarize_bootstrap_results}} for complete bootstrap summary workflow
#'
#' @importFrom data.table data.table setDT
#' @importFrom stats median quantile sd
#' @export
