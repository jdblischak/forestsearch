# =============================================================================
# example_plot_sg_results.R - Usage Examples for Plotting ForestSearch Results
# =============================================================================
#
# This script demonstrates how to visualize subgroup results from ForestSearch
# using the plot_sg_results() function and related visualization tools.
#
# Author: ForestSearch Development Team
# =============================================================================


# =============================================================================
# SECTION 1: BASIC USAGE WITH fs$df.est
# =============================================================================

#' @examples
#' # After running ForestSearch analysis:
#' library(ForestSearch)
#' library(survival)
#'
#' # Assuming fs is a forestsearch object with identified subgroup
#' # fs <- forestsearch(df.analysis = my_data, ...)
#'
#' # Basic visualization - all plots combined
#' result <- plot_sg_results(
#'   fs.est = fs,
#'   outcome.name = "Y",
#'   event.name = "Event",
#'   treat.name = "Treat"
#' )
#'
#' # Print summary
#' print(result)
#'
#' # Access the summary data
#' result$summary
#' result$hr_estimates


# =============================================================================
# SECTION 2: CUSTOMIZED KAPLAN-MEIER PLOTS
# =============================================================================

#' Create Kaplan-Meier Plots Only
#'
#' @examples
#' # Generate only KM survival curves
#' plot_sg_results(
#'   fs.est = fs,
#'   outcome.name = "os_time",
#'   event.name = "os_event",
#'   treat.name = "treatment",
#'   plot_type = "km",
#'   # Customize labels
#'   sg0_name = "High Risk (do not recommend)",
#'   sg1_name = "Standard Risk (recommend treatment)",
#'   treat_labels = c("0" = "Placebo", "1" = "Experimental"),
#'   # Plot options
#'   show_ci = TRUE,
#'   show_logrank = TRUE,
#'   show_hr = TRUE,
#'   title = "Overall Survival by Identified Subgroup"
#' )


# =============================================================================
# SECTION 3: USING sg_consistency_out() RESULTS DIRECTLY
# =============================================================================

#' Working with sg_consistency_out() Output
#'
#' The sg_consistency_out() function returns detailed subgroup information
#' that can be used for custom visualizations.
#'
#' @examples
#' # If you have direct access to consistency results:
#' # sg_out <- sg_consistency_out(df, result_new, sg_focus, index.Z, names.Z, ...)
#'
#' # The output contains:
#' # - result: data.table with ranked subgroup candidates
#' # - sg.harm: character vector of subgroup-defining variables
#' # - sg.harm_label: human-readable labels
#' # - df_flag: data frame with treat.recommend flag
#' # - sg.harm.id: vector indicating subgroup membership
#'
#' # To visualize from sg_consistency_out output:
#' plot_from_sg_consistency <- function(sg_out, df, outcome.name, event.name, treat.name) {

#'   # Merge the flag with original data
#'   df_plot <- merge(df, sg_out$df_flag, by = "id")
#'
#'   # Create a pseudo forestsearch object
#'   fs_pseudo <- list(
#'     df.est = df_plot,
#'     sg.harm = sg_out$sg.harm_label,
#'     grp.consistency = list(out_sg = sg_out)
#'   )
#'
#'   # Plot
#'   plot_sg_results(
#'     fs.est = fs_pseudo,
#'     outcome.name = outcome.name,
#'     event.name = event.name,
#'     treat.name = treat.name
#'   )
#' }


# =============================================================================
# SECTION 4: INTEGRATION WITH weightedsurv PACKAGE
# =============================================================================

#' Plot Weighted Kaplan-Meier Curves (if weightedsurv available)
#'
#' For more sophisticated survival curve visualization with risk tables
#' and confidence intervals, use the weightedsurv package integration.
#'
#' @param fs.est ForestSearch result object
#' @param outcome.name Character. Name of outcome variable
#' @param event.name Character. Name of event variable
#' @param treat.name Character. Name of treatment variable
#' @param by.risk Numeric. Risk table interval
#' @param ... Additional arguments to weightedsurv functions
#'
#' @examples
#' plot_sg_weighted_km <- function(
#'   fs.est,
#'   outcome.name = "Y",
#'   event.name = "Event",
#'   treat.name = "Treat",
#'   by.risk = 12,
#'   ...
#' ) {
#'
#'   if (!requireNamespace("weightedsurv", quietly = TRUE)) {
#'     stop("Package 'weightedsurv' required. Install with: devtools::install_github(...)")
#'   }
#'
#'   df <- fs.est$df.est
#'
#'   # Split into subgroups
#'   df_H <- df[df$treat.recommend == 0, ]
#'   df_Hc <- df[df$treat.recommend == 1, ]
#'
#'   # Create counting process data
#'   dfcount_H <- weightedsurv::df_counting(
#'     df_H,
#'     tte.name = outcome.name,
#'     event.name = event.name,
#'     treat.name = treat.name,
#'     arms = c("treat", "control"),
#'     by.risk = by.risk
#'   )
#'
#'   dfcount_Hc <- weightedsurv::df_counting(
#'     df_Hc,
#'     tte.name = outcome.name,
#'     event.name = event.name,
#'     treat.name = treat.name,
#'     arms = c("treat", "control"),
#'     by.risk = by.risk
#'   )
#'
#'   # Plot side by side
#'   par(mfrow = c(1, 2))
#'
#'   weightedsurv::plot_weighted_km(
#'     dfcount_H,
#'     conf.int = TRUE,
#'     show.logrank = TRUE,
#'     put.legend.lr = "topleft",
#'     ymax = 1.05,
#'     main = paste0("Questionable Subgroup (n=", nrow(df_H), ")")
#'   )
#'
#'   weightedsurv::plot_weighted_km(
#'     dfcount_Hc,
#'     conf.int = TRUE,
#'     show.logrank = TRUE,
#'     put.legend.lr = "topleft",
#'     ymax = 1.05,
#'     main = paste0("Recommend Subgroup (n=", nrow(df_Hc), ")")
#'   )
#'
#'   par(mfrow = c(1, 1))
#' }


# =============================================================================
# SECTION 5: CREATING PUBLICATION-READY FOREST PLOTS
# =============================================================================

#' Create Publication Forest Plot Using forestploter
#'
#' For publication-quality forest plots, use plot_subgroup_results_forestplot()
#' which provides comprehensive formatting options.
#'
#' @examples
#' # Using the existing plot_subgroup_results_forestplot function:
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(
#'     fs.est = fs,
#'     fs_bc = fs_bootstrap,  # Optional: bootstrap results
#'     fs_kfold = fs_cv       # Optional: cross-validation results
#'   ),
#'   df_analysis = df,
#'   subgroup_list = list(
#'     biomarker_high = list(
#'       subset_expr = "biomarker > median(biomarker)",
#'       name = "Biomarker High",
#'       type = "reference"
#'     )
#'   ),
#'   outcome.name = "os_time",
#'   event.name = "os_event",
#'   treat.name = "treatment",
#'   E.name = "Experimental",
#'   C.name = "Control",
#'   show_cv_metrics = TRUE
#' )
#'
#' # Display the plot
#' plot(result)


# =============================================================================
# SECTION 6: EXTRACTING AND FORMATTING SUBGROUP LABELS
# =============================================================================

#' Get Human-Readable Subgroup Definition
#'
#' Formats the subgroup definition for display in plots and tables.
#'
#' @param fs.est ForestSearch result object
#' @return Character string with formatted subgroup definition
#'
#' @examples
format_sg_definition <- function(fs.est) {
  if (is.null(fs.est$sg.harm)) {
    return("No subgroup identified")
  }

  # Get labels if available
  labels <- fs.est$sg.harm

  # If grp.consistency has labeled version, use that
  if (!is.null(fs.est$grp.consistency$out_sg$sg.harm_label)) {
    labels <- fs.est$grp.consistency$out_sg$sg.harm_label
  }

  # Format as readable string
  paste(labels, collapse = " AND ")
}


# =============================================================================
# SECTION 7: COMPLETE WORKFLOW EXAMPLE
# =============================================================================

#' Complete Visualization Workflow
#'
#' This example shows a complete workflow from ForestSearch analysis
#' to visualization of results.
#'
#' @examples
#' \dontrun{
#' library(ForestSearch)
#' library(survival)
#'
#' # Step 1: Prepare data
#' data(gbsg, package = "survival")
#' df <- gbsg
#' df$id <- 1:nrow(df)
#' df$Y <- df$rfstime
#' df$Event <- df$status
#' df$Treat <- df$hormon
#'
#' # Create candidate factors
#' df$age_cat <- ifelse(df$age > 50, "old", "young")
#' df$grade_high <- ifelse(df$grade == 3, 1, 0)
#' df$size_large <- ifelse(df$size > 25, 1, 0)
#'
#' # Step 2: Run ForestSearch
#' fs <- forestsearch(
#'   df.analysis = df,
#'   outcome.name = "Y",
#'   event.name = "Event",
#'   treat.name = "Treat",
#'   id.name = "id",
#'   confounders.name = c("age_cat", "grade_high", "size_large", "pgr", "er"),
#'   n.min = 60,
#'   hr.threshold = 1.25,
#'   details = TRUE
#' )
#'
#' # Step 3: Check if subgroup was found
#' if (!is.null(fs$sg.harm)) {
#'   cat("Subgroup found:", format_sg_definition(fs), "\n")
#'
#'   # Step 4: Visualize results
#'
#'   # 4a: Combined visualization
#'   result <- plot_sg_results(
#'     fs.est = fs,
#'     outcome.name = "Y",
#'     event.name = "Event",
#'     treat.name = "Treat",
#'     plot_type = "combined",
#'     sg0_name = "Questionable Benefit",
#'     sg1_name = "Clear Benefit",
#'     treat_labels = c("0" = "No Hormone", "1" = "Hormone Therapy"),
#'     title = "Treatment Effect by ForestSearch-Identified Subgroup"
#'   )
#'
#'   # 4b: Print summary
#'   print(result)
#'
#'   # 4c: Save HR estimates
#'   write.csv(result$hr_estimates, "subgroup_hr_estimates.csv", row.names = FALSE)
#'
#'   # Step 5: Optional - Run bootstrap for bias correction
#'   fs_bc <- forestsearch_bootstrap_dofuture(
#'     fs.est = fs,
#'     nb_boots = 500,
#'     parallel_args = list(plan = "multisession", workers = 4)
#'   )
#'
#'   # Step 6: Create publication forest plot
#'   forest_result <- plot_subgroup_results_forestplot(
#'     fs_results = list(fs.est = fs, fs_bc = fs_bc),
#'     df_analysis = df,
#'     outcome.name = "Y",
#'     event.name = "Event",
#'     treat.name = "Treat"
#'   )
#'
#' } else {
#'   cat("No subgroup meeting criteria was identified.\n")
#' }
#' }


# =============================================================================
# SECTION 8: UTILITY FUNCTION FOR GGPLOT2 VISUALIZATION
# =============================================================================

#' Create ggplot2 Kaplan-Meier Plots
#'
#' Alternative visualization using ggplot2 and survminer for
#' more customizable output.
#'
#' @param fs.est ForestSearch result object
#' @param outcome.name Character. Outcome variable name
#' @param event.name Character. Event variable name
#' @param treat.name Character. Treatment variable name
#' @param ... Additional arguments to survminer::ggsurvplot
#'
#' @return List of ggsurvplot objects
#'
#' @examples
#' \dontrun{
#' plot_sg_ggsurvplot <- function(
#'   fs.est,
#'   outcome.name = "Y",
#'   event.name = "Event",
#'   treat.name = "Treat",
#'   ...
#' ) {
#'
#'   if (!requireNamespace("survminer", quietly = TRUE)) {
#'     stop("Package 'survminer' required for ggplot2 visualizations")
#'   }
#'
#'   df <- fs.est$df.est
#'   df_H <- df[df$treat.recommend == 0, ]
#'   df_Hc <- df[df$treat.recommend == 1, ]
#'
#'   # Create survival formula
#'   surv_formula <- as.formula(paste0(
#'     "survival::Surv(", outcome.name, ", ", event.name, ") ~ ", treat.name
#'   ))
#'
#'   # Fit models
#'   fit_H <- survival::survfit(surv_formula, data = df_H)
#'   fit_Hc <- survival::survfit(surv_formula, data = df_Hc)
#'
#'   # Create plots
#'   p_H <- survminer::ggsurvplot(
#'     fit_H,
#'     data = df_H,
#'     risk.table = TRUE,
#'     pval = TRUE,
#'     conf.int = TRUE,
#'     title = paste0("Questionable Subgroup (n=", nrow(df_H), ")"),
#'     ...
#'   )
#'
#'   p_Hc <- survminer::ggsurvplot(
#'     fit_Hc,
#'     data = df_Hc,
#'     risk.table = TRUE,
#'     pval = TRUE,
#'     conf.int = TRUE,
#'     title = paste0("Recommend Subgroup (n=", nrow(df_Hc), ")"),
#'     ...
#'   )
#'
#'   list(H = p_H, Hc = p_Hc)
#' }
#' }
