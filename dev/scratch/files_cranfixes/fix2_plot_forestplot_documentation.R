# =============================================================================
# FIXED: R/plot_subgroup_results_forestplot.R - roxygen2 documentation
# =============================================================================
# 
# This file contains the COMPLETE corrected roxygen2 documentation block for
# the plot_subgroup_results_forestplot() function. The fix adds the missing 
# @param xlog.
#
# INSTRUCTIONS:
# 1. Open R/plot_subgroup_results_forestplot.R
# 2. Find the roxygen2 documentation block starting with:
#    #' Plot Subgroup Results Forest Plot
# 3. Replace the entire @param section with the content below
#
# =============================================================================

#' Plot Subgroup Results Forest Plot
#'
#' Creates a publication-ready forest plot displaying identified subgroups with
#' hazard ratios, bias-corrected estimates, and cross-validation metrics.
#' This wrapper integrates ForestSearch results with the forestploter package.
#'
#' @description
#' Generates a comprehensive forest plot showing:
#' - ITT (Intent-to-Treat) population estimate
#' - Reference subgroups (e.g., by biomarker levels)
#' - Post-hoc identified subgroups with bias-corrected estimates
#' - Cross-validation agreement metrics as annotations
#'
#' @param fs_results List. A list containing ForestSearch analysis results with elements:
#'   \itemize{
#'     \item \code{fs.est}: ForestSearch estimation object from \code{\link{forestsearch}}
#'     \item \code{fs_bc}: Bootstrap bias-corrected results from
#'       \code{\link{forestsearch_bootstrap_dofuture}}
#'     \item \code{fs_kfold}: K-fold cross-validation results from
#'       \code{\link{forestsearch_Kfold}} or \code{\link{forestsearch_tenfold}} (optional)
#'     \item \code{fs_OOB}: Out-of-bag cross-validation results (optional, alternative to fs_kfold)
#'   }
#' @param df_analysis Data frame. The analysis dataset with outcome, event, and treatment
#'   variables.
#' @param subgroup_list List. Named list of subgroup definitions to include in the plot.
#'   Each element should be a list with:
#'   \itemize{
#'     \item \code{subset_expr}: Character string for subsetting (e.g., "BM > 1")
#'     \item \code{name}: Display name for the subgroup
#'     \item \code{type}: Either "reference" for pre-specified or "posthoc" for identified
#'   }
#' @param outcome.name Character. Name of the survival time variable.
#' @param event.name Character. Name of the event indicator variable.
#' @param treat.name Character. Name of the treatment variable.
#' @param E.name Character. Label for experimental arm (default: "Experimental").
#' @param C.name Character. Label for control arm (default: "Control").
#' @param est.scale Character. Estimate scale: "hr" or "1/hr" (default: "hr").
#' @param xlog Logical. If TRUE (default), the x-axis is plotted on a
#'   logarithmic scale. This is standard for hazard ratio forest plots
#'   where equal distances represent equal relative effects.
#' @param title_text Character. Plot title (default: NULL).
#' @param arrow_text Character vector of length 2. Arrow labels for forest plot
#'   (default: c("Favors Experimental", "Favors Control")).
#' @param footnote_text Character vector. Footnote text for the plot explaining CV metrics
#'   (default provides CV interpretation guidance; set to NULL to omit).
#' @param xlim Numeric vector of length 2. X-axis limits (default: c(0.25, 1.5)).
#' @param ticks_at Numeric vector. X-axis tick positions
#'   (default: c(0.25, 0.70, 1.0, 1.5)).
#' @param show_cv_metrics Logical. Whether to show cross-validation metrics
#'   (default: TRUE if fs_kfold or fs_OOB available).
#' @param cv_source Character. Source for CV metrics:
#'   "auto" (default, uses both if available, otherwise whichever is present),
#'   "kfold" (use fs_kfold only),
#'   "oob" (use fs_OOB only), or
#'   "both" (explicitly use both fs_kfold and fs_OOB, with K-fold first then OOB).
#' @param posthoc_colors Character vector. Colors for post-hoc subgroup rows
#'   (default: c("powderblue", "beige")).
#' @param reference_colors Character vector. Colors for reference subgroup rows
#'   (default: c("yellow", "powderblue")).
#' @param ci_column_spaces Integer. Number of spaces for the CI plot column width.
#'   More spaces = wider CI column (default: 20).
#' @param conf.level Numeric. Confidence level for intervals (default: 0.95 for 95% CI).
#'   Used to calculate the z-multiplier as qnorm(1 - (1 - conf.level)/2).
#'
#' @return A list containing:
#'   \describe{
#'     \item{plot}{The forestploter grob object (can be rendered with plot())}
#'     \item{data}{The data frame used for the forest plot}
#'     \item{row_types}{Character vector of row types for styling reference}
#'     \item{cv_metrics}{Cross-validation metrics text (if available)}
#'   }
#'
#' @details
#' ## ForestSearch Labeling Convention
#'
#' ForestSearch identifies subgroups based on hazard ratio thresholds:
#' \itemize{
#'   \item \code{sg.harm}: Contains the definition of the "harm" or "questionable"
#'     subgroup (H)
#'   \item \code{treat.recommend == 0}: Patient is IN the harm subgroup (H)
#'   \item \code{treat.recommend == 1}: Patient is in the COMPLEMENT subgroup
#'     (Hc, typically benefit)
#' }
#'
#' For \code{est.scale = "hr"} (searching for harm):
#' \itemize{
#'   \item H (treat.recommend=0): Subgroup defined by sg.harm with elevated HR
#'     (harm/questionable)
#'   \item Hc (treat.recommend=1): Complement of sg.harm (potential benefit)
#' }
#'
#' For \code{est.scale = "1/hr"} (searching for benefit):
#' \itemize{
#'   \item Roles are reversed: H becomes the benefit group
#' }
#'
#' @examples
#' \dontrun{
#' # Load ForestSearch results
#' load("fs_analysis_results.Rdata")  # Contains fs.est, fs_bc, fs_kfold
#'
#' # Define subgroups to display
#' subgroups <- list(
#'   bm_gt1 = list(
#'     subset_expr = "BM > 1",
#'     name = "BM > 1",
#'     type = "reference"
#'   ),
#'   bm_gt1_size_gt19 = list(
#'     subset_expr = "BM > 1 & tmrsize > 19",
#'     name = "BM > 1 & Tumor Size > 19",
#'     type = "posthoc"
#'   )
#' )
#'
#' # Create the forest plot
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(fs.est = fs.est, fs_bc = fs_bc, fs_kfold = fs_kfold),
#'   df_analysis = df_itt,
#'   subgroup_list = subgroups,
#'   outcome.name = "os_time",
#'   event.name = "os_event",
#'   treat.name = "combo",
#'   E.name = "Experimental+CT",
#'   C.name = "CT"
#' )
#'
#' # Display the plot
#' plot(result$plot)
#' }
#'
#' @seealso
#' \code{\link{forestsearch}} for running the subgroup analysis
#' \code{\link{forestsearch_bootstrap_dofuture}} for bootstrap bias correction
#' \code{\link{forestsearch_Kfold}} for cross-validation
#'
#' @importFrom survival coxph Surv
#' @importFrom grid gpar
#' @export
