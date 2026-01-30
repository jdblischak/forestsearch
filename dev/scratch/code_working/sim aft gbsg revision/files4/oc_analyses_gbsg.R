# =============================================================================
# Operating Characteristics Analysis Functions for GBSG Simulations (Aligned)
# =============================================================================
#
# Functions for running ForestSearch and GRF simulation analyses using
# GBSG-based data generating mechanisms.
#
# Aligned with sim_aft_gbsg_refactored.R and generate_aft_dgm_flex() output:
#   - Supports new hazard_ratios list structure (AHR, AHR_harm, AHR_no_harm)
#   - Includes use_twostage parameter passthrough
#   - Computes AHR metrics from loghr_po when available
#
# Key functions:
#   - run_simulation_analysis(): Main simulation wrapper
#   - run_forestsearch_analysis(): ForestSearch analysis helper
#   - run_grf_analysis(): GRF analysis helper
#   - extract_fs_estimates(): Extract estimates from ForestSearch results
#   - extract_grf_estimates(): Extract estimates from GRF results
#   - summarize_simulation_results(): Create summary tables
#
# =============================================================================

#' @import survival
#' @import data.table
NULL

# =============================================================================
# Global Variables Declaration
# =============================================================================

utils::globalVariables(c(
  "any.H", "ppv", "npv", "sensitivity", "specificity",
  "size.H", "size.Hc", "hr.H.true", "hr.H.hat", "hr.Hc.true", "hr.Hc.hat",
  "hr.itt", "hr.adj.itt", "p.cens", "taumax",
  "analysis", "sim", "aa", "sg_hat",
  # New aligned variables
  "ahr.H.true", "ahr.Hc.true", "ahr.H.hat", "ahr.Hc.hat",
  "loghr_po", "theta_0", "theta_1"
))


# =============================================================================
# Configuration Defaults
# =============================================================================

#' Default ForestSearch Parameters for GBSG Simulations
#'
#' Returns a list of default parameters for ForestSearch analysis
#' in GBSG-based simulations.
#'
#' @return List of default ForestSearch parameters
#'
#' @details
#' Default parameters are optimized for GBSG simulation scenarios with
#' moderate sample sizes (300-1000) and typical event rates.
#'
#' Variable selection defaults:
#' \itemize{
#'   \item use_lasso = TRUE: LASSO-based variable importance (default for FS)
#'   \item use_grf = FALSE: GRF-based variable importance (enable for FSlg)
#' }
#'
#' The `use_twostage` parameter is set to FALSE by default for backward
#' compatibility. Set to TRUE for faster exploratory analyses.
#'
#' @export
default_fs_params <- function() {
  list(
    outcome.name = "y.sim",
    event.name = "event.sim",
    treat.name = "treat",
    id.name = "id",
    # Variable selection: LASSO by default (FS analysis)
    # For FSlg, override with use_grf = TRUE
    use_lasso = TRUE,
    use_grf = FALSE,
    hr.threshold = 1.25,
    hr.consistency = 1.0,
    pconsistency.threshold = 0.90,
    fs.splits = 400,
    n.min = 60,
    d0.min = 12,
    d1.min = 12,
    maxk = 2,
    max.minutes = 5,
    by.risk = 12,
    vi.grf.min = -0.2,
    # Two-stage parameters (aligned with forestsearch)
    use_twostage = FALSE,
    twostage_args = list()
  )
}


#' Default GRF Parameters for GBSG Simulations
#'
#' Returns a list of default parameters for GRF analysis
#' in GBSG-based simulations.
#'
#' @return List of default GRF parameters
#'
#' @export
default_grf_params <- function() {
  list(
    outcome.name = "y.sim",
    event.name = "event.sim",
    treat.name = "treat",
    id.name = "id",
    n.min = 60,
    dmin.grf = 12,
    frac.tau = 0.60,
    maxdepth = 2
  )
}


# =============================================================================
# Helper Functions
# =============================================================================

#' Extract HR from DGM (Backward Compatible)
#'
#' Extracts hazard ratios from DGM object, supporting both old and new formats.
#'
#' @param dgm DGM object (gbsg_dgm or aft_dgm_flex)
#' @param which Character. Which HR to extract: "hr_H", "hr_Hc", "ahr_H", "ahr_Hc",
#'   "hr_overall", "ahr"
#'
#' @return Numeric hazard ratio value
#'
#' @keywords internal
get_dgm_hr <- function(dgm, which = "hr_H") {
  
  # Try new hazard_ratios list first (aligned format)
  if (!is.null(dgm$hazard_ratios)) {
    hr_list <- dgm$hazard_ratios
    result <- switch(which,
      "hr_H" = hr_list$harm_subgroup,
      "hr_Hc" = hr_list$no_harm_subgroup,
      "ahr_H" = hr_list$AHR_harm,
      "ahr_Hc" = hr_list$AHR_no_harm,
      "hr_overall" = hr_list$overall,
      "ahr" = hr_list$AHR,
      NA_real_
    )
    if (!is.null(result) && !is.na(result)) return(result)
  }
  

  # Fall back to old format (direct properties)
  result <- switch(which,
    "hr_H" = dgm$hr_H_true,
    "hr_Hc" = dgm$hr_Hc_true,
    "ahr_H" = dgm$AHR_H_true,
    "ahr_Hc" = dgm$AHR_Hc_true,
    "hr_overall" = dgm$hr_causal,
    "ahr" = dgm$AHR,
    NA_real_
  )
  
  if (is.null(result)) NA_real_ else result
}


#' Create Subgroup Indicator from Factor Definitions
#'
#' Parses factor definitions (e.g., "v1.1", "grade3.1") and creates
#' a binary indicator for subgroup membership.
#'
#' @param df Data frame containing the variables
#' @param sg_factors Character vector of factor definitions
#'
#' @return Integer vector (1 = in subgroup, 0 = not in subgroup)
#'
#' @keywords internal
create_subgroup_indicator <- function(df, sg_factors) {
  indicator <- rep(TRUE, nrow(df))
  
  for (factor_def in sg_factors) {
    if (is.na(factor_def) || factor_def == "") next
    
    parts <- strsplit(factor_def, "\\.")[[1]]
    if (length(parts) >= 2) {
      var_name <- parts[1]
      level <- parts[2]
      
      if (var_name %in% names(df)) {
        indicator <- indicator & (as.character(df[[var_name]]) == level)
      }
    }
  }
  
  as.integer(indicator)
}


#' Compute AHR from loghr_po
#'
#' Computes Average Hazard Ratio from individual log hazard ratios.
#'
#' @param df Data frame with loghr_po column
#' @param subset_indicator Optional logical/integer vector for subsetting
#'
#' @return Numeric AHR value
#'
#' @keywords internal
compute_ahr <- function(df, subset_indicator = NULL) {
  if (!"loghr_po" %in% names(df)) {
    return(NA_real_)
  }
  
  loghr <- df$loghr_po
  
  if (!is.null(subset_indicator)) {
    loghr <- loghr[subset_indicator == 1]
  }
  
  if (length(loghr) == 0 || all(is.na(loghr))) {
    return(NA_real_)
  }
  
  exp(mean(loghr, na.rm = TRUE))
}


# =============================================================================
# Extract ForestSearch Estimates
# =============================================================================

#' Extract Estimates from ForestSearch Results
#'
#' Extracts operating characteristics (HRs, classification metrics, etc.)
#' from ForestSearch analysis results. Aligned with new DGM output structure.
#'
#' @param df Simulated data frame
#' @param fs_res ForestSearch result table (grp.consistency$out_sg$result, or NULL)
#' @param dgm DGM object containing true HRs (supports both old and new formats)
#' @param cox_formula Cox formula for estimation (optional)
#' @param cox_formula_adj Adjusted Cox formula (optional)
#' @param analysis Analysis label (e.g., "FS", "FSlg")
#' @param fs_full Full forestsearch result object (for df.est access)
#' @param verbose Logical. Print extraction details. Default: FALSE
#'
#' @return data.table with extracted estimates including AHR metrics
#'
#' @importFrom data.table data.table
#' @importFrom survival coxph Surv
#' @keywords internal
extract_fs_estimates <- function(
    df,
    fs_res,
    dgm,
    cox_formula = NULL,
    cox_formula_adj = NULL,
    analysis = "FS",
    fs_full = NULL,
    verbose = FALSE
) {
  
  # Initialize output with both HR and AHR metrics
  out <- data.table::data.table(
    analysis = analysis,
    any.H = 0L,
    size.H = NA_integer_,
    size.Hc = nrow(df),
    # Cox-based HRs
    hr.H.true = get_dgm_hr(dgm, "hr_H"),
    hr.H.hat = NA_real_,
    hr.Hc.true = get_dgm_hr(dgm, "hr_Hc"),
    hr.Hc.hat = NA_real_,
    hr.itt = NA_real_,
    hr.adj.itt = NA_real_,
    # AHR metrics (aligned with generate_aft_dgm_flex)
    ahr.H.true = get_dgm_hr(dgm, "ahr_H"),
    ahr.Hc.true = get_dgm_hr(dgm, "ahr_Hc"),
    ahr.H.hat = NA_real_,
    ahr.Hc.hat = NA_real_,
    # Classification metrics
    ppv = NA_real_,
    npv = NA_real_,
    sensitivity = NA_real_,
    specificity = NA_real_,
    # Data summary
    p.cens = 1 - mean(df$event.sim),
    taumax = max(df$y.sim)
  )
  
  # Compute ITT estimates
  out$hr.itt <- tryCatch({
    exp(survival::coxph(
      survival::Surv(y.sim, event.sim) ~ treat,
      data = df
    )$coefficients)
  }, error = function(e) NA_real_)
  
  # Adjusted ITT
  if (!is.null(cox_formula_adj)) {
    out$hr.adj.itt <- tryCatch({
      exp(survival::coxph(cox_formula_adj, data = df)$coefficients["treat"])
    }, error = function(e) NA_real_)
  }
  
  # If no subgroup found
  if (is.null(fs_res) || nrow(fs_res) == 0) {
    if (verbose) {
      message(sprintf("  [%s] No subgroup identified", analysis))
    }
    out$hr.Hc.hat <- out$hr.itt
    # Compute AHR for complement (whole population when no subgroup)
    out$ahr.Hc.hat <- compute_ahr(df)
    return(out)
  }
  
  # Subgroup found
  out$any.H <- 1L
  
  # Try to use df.est from full result
  df_with_sg <- NULL
  if (!is.null(fs_full) && !is.null(fs_full$df.est)) {
    df_with_sg <- fs_full$df.est
  }
  
  # Create subgroup indicator
  if (!is.null(df_with_sg) && "sg" %in% names(df_with_sg)) {
    df$sg_hat <- df_with_sg$sg
  } else if (!is.null(fs_full) && !is.null(fs_full$sg.harm)) {
    df$sg_hat <- create_subgroup_indicator(df, fs_full$sg.harm)
  } else {
    # Try to extract from result table
    sg_def <- fs_res$subgroup[1]
    if (!is.null(sg_def) && !is.na(sg_def)) {
      sg_factors <- strsplit(sg_def, " & ")[[1]]
      df$sg_hat <- create_subgroup_indicator(df, sg_factors)
    } else {
      df$sg_hat <- 0L
    }
  }
  
  # Compute sizes
  out$size.H <- sum(df$sg_hat == 1)
  out$size.Hc <- sum(df$sg_hat == 0)
  
  if (verbose) {
    message(sprintf("  [%s] Subgroup found: n(H) = %d, n(Hc) = %d",
                    analysis, out$size.H, out$size.Hc))
  }
  
  # Compute HR estimates in identified subgroups
  if (out$size.H > 10 && sum(df$sg_hat == 1 & df$event.sim == 1) >= 5) {
    out$hr.H.hat <- tryCatch({
      df_H <- subset(df, sg_hat == 1)
      if (length(unique(df_H$treat)) == 2) {
        exp(survival::coxph(
          survival::Surv(y.sim, event.sim) ~ treat,
          data = df_H
        )$coefficients)
      } else NA_real_
    }, error = function(e) NA_real_)
  }
  
  if (out$size.Hc > 10 && sum(df$sg_hat == 0 & df$event.sim == 1) >= 5) {
    out$hr.Hc.hat <- tryCatch({
      df_Hc <- subset(df, sg_hat == 0)
      if (length(unique(df_Hc$treat)) == 2) {
        exp(survival::coxph(
          survival::Surv(y.sim, event.sim) ~ treat,
          data = df_Hc
        )$coefficients)
      } else NA_real_
    }, error = function(e) NA_real_)
  }
  
  # Compute AHR estimates from loghr_po (aligned metric)
  out$ahr.H.hat <- compute_ahr(df, df$sg_hat)
  out$ahr.Hc.hat <- compute_ahr(df, 1L - df$sg_hat)
  
  if (verbose) {
    message(sprintf("  [%s] HR estimates: H = %.3f, Hc = %.3f",
                    analysis,
                    ifelse(is.na(out$hr.H.hat), NA, out$hr.H.hat),
                    ifelse(is.na(out$hr.Hc.hat), NA, out$hr.Hc.hat)))
    if (!is.na(out$ahr.H.hat)) {
      message(sprintf("  [%s] AHR estimates: H = %.3f, Hc = %.3f",
                      analysis, out$ahr.H.hat, out$ahr.Hc.hat))
    }
  }
  
  # Compute classification metrics
  if ("flag.harm" %in% names(df)) {
    true_H <- df$flag.harm == 1
    hat_H <- df$sg_hat == 1
    
    tp <- sum(true_H & hat_H)
    fp <- sum(!true_H & hat_H)
    tn <- sum(!true_H & !hat_H)
    fn <- sum(true_H & !hat_H)
    
    out$sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    out$specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    out$ppv <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    out$npv <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
    
    if (verbose) {
      message(sprintf("  [%s] Classification: Sens = %.3f, Spec = %.3f, PPV = %.3f, NPV = %.3f",
                      analysis, out$sensitivity, out$specificity, out$ppv, out$npv))
    }
    
    # Recompute true HR in identified subgroup
    if (out$size.H > 10) {
      out$hr.H.true <- tryCatch({
        df_hat_H <- subset(df, sg_hat == 1)
        if (length(unique(df_hat_H$treat)) == 2 &&
            sum(df_hat_H$event.sim) >= 5) {
          exp(survival::coxph(
            survival::Surv(y.sim, event.sim) ~ treat,
            data = df_hat_H
          )$coefficients)
        } else get_dgm_hr(dgm, "hr_H")
      }, error = function(e) get_dgm_hr(dgm, "hr_H"))
    }
  }
  
  out
}


# =============================================================================
# Run ForestSearch Analysis
# =============================================================================

#' Run ForestSearch Analysis
#'
#' Helper function to run ForestSearch and extract estimates.
#' Aligned with forestsearch() parameters including use_twostage.
#'
#' @param data Data frame with simulated trial data
#' @param confounders_name Character vector of confounder names
#' @param params List of ForestSearch parameters
#' @param dgm DGM object for computing true HRs
#' @param cox_formula Cox formula for estimation
#' @param cox_formula_adj Adjusted Cox formula
#' @param analysis_label Character label for this analysis
#' @param verbose Print details
#'
#' @return data.table with analysis estimates
#'
#' @keywords internal
run_forestsearch_analysis <- function(
    data,
    confounders_name,
    params,
    dgm,
    cox_formula = NULL,
    cox_formula_adj = NULL,
    analysis_label = "FS",
    verbose = FALSE
) {
  
  if (verbose) {
    message(sprintf("\n  [%s] Starting ForestSearch analysis...", analysis_label))
    message(sprintf("  [%s] Data: n = %d, events = %d (%.1f%%)",
                    analysis_label, nrow(data), sum(data$event.sim),
                    100 * mean(data$event.sim)))
    message(sprintf("  [%s] Confounders: %s",
                    analysis_label, paste(confounders_name, collapse = ", ")))
    message(sprintf("  [%s] Parameters: use_lasso = %s, use_grf = %s, hr.threshold = %.2f, use_twostage = %s",
                    analysis_label, params$use_lasso, params$use_grf, 
                    params$hr.threshold, params$use_twostage))
  }
  
  # Build forestsearch call arguments
  fs_args <- list(
    df.analysis = data,
    confounders.name = confounders_name,
    details = verbose,
    plot.sg = FALSE
  )
  
  # Valid forestsearch parameters - INCLUDES use_twostage and twostage_args
  param_names <- c(
    "outcome.name", "event.name", "treat.name", "id.name",
    "use_lasso", "use_grf", "conf_force",
    "hr.threshold", "hr.consistency", "pconsistency.threshold",
    "fs.splits", "n.min", "d0.min", "d1.min",
    "maxk", "max.minutes", "by.risk", "vi.grf.min",
    "frac.tau", "dmin.grf", "grf_depth",
    # Two-stage parameters (aligned)
    "use_twostage", "twostage_args"
  )
  
  for (pn in param_names) {
    if (!is.null(params[[pn]])) {
      fs_args[[pn]] <- params[[pn]]
    }
  }
  
  # Run ForestSearch
  if (verbose) {
    message(sprintf("  [%s] Calling forestsearch()...", analysis_label))
  }
  
  fs_result <- tryCatch({
    do.call(forestsearch, fs_args)
  }, error = function(e) {
    warning(sprintf("%s analysis failed: %s", analysis_label, e$message))
    NULL
  })
  
  # Check result
  if (verbose) {
    if (is.null(fs_result)) {
      message(sprintf("  [%s] forestsearch() returned NULL", analysis_label))
    } else {
      sg_def <- if (!is.null(fs_result$sg.harm)) {
        paste(fs_result$sg.harm, collapse = " & ")
      } else {
        "none"
      }
      message(sprintf("  [%s] forestsearch() completed. Subgroup: %s",
                      analysis_label, sg_def))
    }
  }
  
  # Check for results in grp.consistency$out_sg$result
  has_result <- !is.null(fs_result) &&
                !is.null(fs_result$grp.consistency) &&
                !is.null(fs_result$grp.consistency$out_sg) &&
                !is.null(fs_result$grp.consistency$out_sg$result) &&
                nrow(fs_result$grp.consistency$out_sg$result) > 0
  
  if (has_result) {
    extract_fs_estimates(
      df = data,
      fs_res = fs_result$grp.consistency$out_sg$result,
      fs_full = fs_result,
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis = analysis_label,
      verbose = verbose
    )
  } else {
    extract_fs_estimates(
      df = data,
      fs_res = NULL,
      fs_full = NULL,
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis = analysis_label,
      verbose = verbose
    )
  }
}


# =============================================================================
# Extract GRF Estimates
# =============================================================================

#' Extract Estimates from GRF Results
#'
#' Aligned with new DGM output structure including AHR metrics.
#'
#' @param df Simulated data frame
#' @param grf_est GRF estimation result from grf.subg.harm.survival()
#' @param dgm DGM object
#' @param cox_formula Cox formula
#' @param cox_formula_adj Adjusted Cox formula
#' @param analysis Analysis label
#' @param frac_tau Fraction of tau used
#' @param verbose Print extraction details
#' @param debug Print detailed debugging information about GRF result structure
#'
#' @return data.table with extracted estimates
#'
#' @keywords internal
extract_grf_estimates <- function(
    df,
    grf_est,
    dgm,
    cox_formula = NULL,
    cox_formula_adj = NULL,
    analysis = "GRF",
    frac_tau = 1.0,
    verbose = FALSE,
    debug = FALSE
) {
  
  # Initialize output with both HR and AHR metrics
  out <- data.table::data.table(
    analysis = analysis,
    any.H = 0L,
    size.H = NA_integer_,
    size.Hc = nrow(df),
    # Cox-based HRs
    hr.H.true = get_dgm_hr(dgm, "hr_H"),
    hr.H.hat = NA_real_,
    hr.Hc.true = get_dgm_hr(dgm, "hr_Hc"),
    hr.Hc.hat = NA_real_,
    hr.itt = NA_real_,
    hr.adj.itt = NA_real_,
    # AHR metrics (aligned)
    ahr.H.true = get_dgm_hr(dgm, "ahr_H"),
    ahr.Hc.true = get_dgm_hr(dgm, "ahr_Hc"),
    ahr.H.hat = NA_real_,
    ahr.Hc.hat = NA_real_,
    # Classification metrics
    ppv = NA_real_,
    npv = NA_real_,
    sensitivity = NA_real_,
    specificity = NA_real_,
    # Data summary
    p.cens = 1 - mean(df$event.sim),
    taumax = max(df$y.sim)
  )
  
  # Compute ITT
  out$hr.itt <- tryCatch({
    exp(survival::coxph(
      survival::Surv(y.sim, event.sim) ~ treat,
      data = df
    )$coefficients)
  }, error = function(e) NA_real_)
  
  # Adjusted ITT
  if (!is.null(cox_formula_adj)) {
    out$hr.adj.itt <- tryCatch({
      exp(survival::coxph(cox_formula_adj, data = df)$coefficients["treat"])
    }, error = function(e) NA_real_)
  }
  
  # Debug: print GRF result structure
  if (debug && !is.null(grf_est)) {
    message(sprintf("  [%s] GRF result class: %s", analysis, paste(class(grf_est), collapse = ", ")))
    message(sprintf("  [%s] GRF result names: %s", analysis, paste(names(grf_est), collapse = ", ")))
  }
  
  # If no GRF result or no subgroup
  if (is.null(grf_est)) {
    if (verbose) {
      message(sprintf("  [%s] No GRF result - returning ITT estimates", analysis))
    }
    out$hr.Hc.hat <- out$hr.itt
    out$ahr.Hc.hat <- compute_ahr(df)
    return(out)
  }
  
  # Try to extract subgroup indicator from GRF result
  sg_indicator <- NULL
  
  # Try different possible locations
  if ("sg" %in% names(grf_est)) {
    sg_indicator <- grf_est$sg
  } else if ("df.est" %in% names(grf_est) && "sg" %in% names(grf_est$df.est)) {
    sg_indicator <- grf_est$df.est$sg
  } else if ("predictions" %in% names(grf_est)) {
    # Use CATE predictions to define subgroup (positive = harm)
    sg_indicator <- as.integer(grf_est$predictions > 0)
  }
  
  if (is.null(sg_indicator) || all(sg_indicator == 0)) {
    if (verbose) {
      message(sprintf("  [%s] No subgroup identified from GRF", analysis))
    }
    out$hr.Hc.hat <- out$hr.itt
    out$ahr.Hc.hat <- compute_ahr(df)
    return(out)
  }
  
  # Subgroup found
  out$any.H <- 1L
  df$sg_hat <- sg_indicator
  
  # Compute sizes
  out$size.H <- sum(df$sg_hat == 1)
  out$size.Hc <- sum(df$sg_hat == 0)
  
  if (verbose) {
    message(sprintf("  [%s] GRF subgroup: n(H) = %d, n(Hc) = %d",
                    analysis, out$size.H, out$size.Hc))
  }
  
  # Compute HR estimates
  if (out$size.H > 10 && sum(df$sg_hat == 1 & df$event.sim == 1) >= 5) {
    out$hr.H.hat <- tryCatch({
      df_H <- subset(df, sg_hat == 1)
      if (length(unique(df_H$treat)) == 2) {
        exp(survival::coxph(
          survival::Surv(y.sim, event.sim) ~ treat,
          data = df_H
        )$coefficients)
      } else NA_real_
    }, error = function(e) NA_real_)
  }
  
  if (out$size.Hc > 10 && sum(df$sg_hat == 0 & df$event.sim == 1) >= 5) {
    out$hr.Hc.hat <- tryCatch({
      df_Hc <- subset(df, sg_hat == 0)
      if (length(unique(df_Hc$treat)) == 2) {
        exp(survival::coxph(
          survival::Surv(y.sim, event.sim) ~ treat,
          data = df_Hc
        )$coefficients)
      } else NA_real_
    }, error = function(e) NA_real_)
  }
  
  # Compute AHR estimates
  out$ahr.H.hat <- compute_ahr(df, df$sg_hat)
  out$ahr.Hc.hat <- compute_ahr(df, 1L - df$sg_hat)
  
  # Classification metrics
  if ("flag.harm" %in% names(df)) {
    true_H <- df$flag.harm == 1
    hat_H <- df$sg_hat == 1
    
    tp <- sum(true_H & hat_H)
    fp <- sum(!true_H & hat_H)
    tn <- sum(!true_H & !hat_H)
    fn <- sum(true_H & !hat_H)
    
    out$sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    out$specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    out$ppv <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    out$npv <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
    
    if (verbose) {
      message(sprintf("  [%s] Classification: Sens = %.3f, Spec = %.3f",
                      analysis, out$sensitivity, out$specificity))
    }
  }
  
  out
}


# =============================================================================
# Run GRF Analysis
# =============================================================================

#' Run GRF Analysis
#'
#' Helper function to run standalone GRF analysis.
#'
#' @param data Data frame with simulated trial data
#' @param confounders_name Character vector of confounder names
#' @param params List of GRF parameters
#' @param dgm DGM object for computing true HRs
#' @param cox_formula Cox formula for estimation
#' @param cox_formula_adj Adjusted Cox formula
#' @param analysis_label Character label for this analysis
#' @param verbose Print details
#' @param debug Print detailed debugging information
#'
#' @return data.table with analysis estimates
#'
#' @keywords internal
run_grf_analysis <- function(
    data,
    confounders_name,
    params,
    dgm,
    cox_formula = NULL,
    cox_formula_adj = NULL,
    analysis_label = "GRF",
    verbose = FALSE,
    debug = FALSE
) {
  
  if (verbose) {
    message(sprintf("\n  [%s] Starting GRF analysis...", analysis_label))
  }
  
  # Try to get grf.subg.harm.survival function
  grf_fun <- tryCatch({
    get("grf.subg.harm.survival", mode = "function", envir = parent.frame())
  }, error = function(e) {
    tryCatch({
      get("grf.subg.harm.survival", mode = "function", envir = globalenv())
    }, error = function(e2) NULL)
  })
  
  if (is.null(grf_fun)) {
    if (verbose) {
      message(sprintf("  [%s] grf.subg.harm.survival not found. Skipping.", analysis_label))
    }
    warning("grf.subg.harm.survival function not found. Skipping GRF analysis.")
    return(extract_grf_estimates(
      df = data,
      grf_est = NULL,
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis = analysis_label,
      verbose = verbose,
      debug = debug
    ))
  }
  
  # Run GRF analysis
  grf_result <- tryCatch({
    grf_fun(
      data = data,
      confounders.name = confounders_name,
      outcome.name = params$outcome.name,
      event.name = params$event.name,
      id.name = params$id.name,
      treat.name = params$treat.name,
      n.min = params$n.min,
      dmin.grf = params$dmin.grf,
      frac.tau = params$frac.tau,
      maxdepth = params$maxdepth,
      details = verbose
    )
  }, error = function(e) {
    if (verbose) {
      message(sprintf("  [%s] GRF analysis failed: %s", analysis_label, e$message))
    }
    warning(sprintf("%s analysis failed: %s", analysis_label, e$message))
    NULL
  })
  
  extract_grf_estimates(
    df = data,
    grf_est = grf_result,
    dgm = dgm,
    cox_formula = cox_formula,
    cox_formula_adj = cox_formula_adj,
    analysis = analysis_label,
    frac_tau = params$frac.tau,
    verbose = verbose,
    debug = debug
  )
}


# =============================================================================
# Main Simulation Analysis Function
# =============================================================================

#' Run Single Simulation Analysis
#'
#' Executes ForestSearch and/or GRF analysis on a single simulated dataset.
#' This is the core function called within a simulation loop.
#'
#' Aligned with create_gbsg_dgm() and generate_aft_dgm_flex() output structures.
#'
#' @param sim_id Integer. Simulation index for seed offset and tracking
#' @param dgm A DGM object from \code{\link{create_gbsg_dgm}} or similar
#' @param n_sample Integer. Sample size for simulation
#' @param max_follow Numeric. Maximum follow-up time. Default: Inf
#' @param muC_adj Numeric. Censoring adjustment. Default: 0
#' @param confounders_base Character vector. Base confounder names
#' @param n_add_noise Integer. Number of noise variables to add. Default: 0
#' @param run_fs Logical. Run ForestSearch with LASSO variable selection.
#'   Default: TRUE. Analysis label: "FS"
#' @param run_fs_grf Logical. Run ForestSearch with LASSO + GRF variable selection.
#'   Default: TRUE. Analysis label: "FSlg"
#' @param run_grf Logical. Run standalone GRF analysis (grf.subg.harm.survival).
#'   Default: TRUE. Analysis label: "GRF"
#' @param fs_params List. ForestSearch parameters (overrides defaults)
#' @param grf_params List. GRF parameters (overrides defaults)
#' @param cox_formula Formula. Cox model formula for estimation
#' @param cox_formula_adj Formula. Adjusted Cox model formula
#' @param n_sims_total Integer. Total simulations (for progress display)
#' @param seed_base Integer. Base random seed. Default: 8316951
#' @param verbose Logical. Print progress. Default: FALSE
#' @param verbose_n Integer. Only print verbose output for first N simulations.
#'   Default: NULL (print for all simulations when verbose = TRUE)
#' @param debug Logical. Print detailed debugging information. Default: FALSE
#'
#' @return A data.table with analysis results for all requested methods,
#'   including both HR and AHR metrics. Contains columns:
#'   \describe{
#'     \item{sim}{Simulation ID}
#'     \item{sizeH_true, propH_true}{True harm subgroup size/proportion in sample}
#'     \item{analysis}{Analysis method: "FS", "FSlg", or "GRF"}
#'     \item{any.H}{1 if subgroup identified, 0 otherwise}
#'     \item{size.H, size.Hc}{Size of identified H and complement}
#'     \item{hr.H.true, hr.H.hat}{True and estimated HR in identified H}
#'     \item{hr.Hc.true, hr.Hc.hat}{True and estimated HR in identified Hc}
#'     \item{ahr.H.true, ahr.H.hat}{True and estimated AHR in identified H}
#'     \item{sensitivity, specificity, ppv, npv}{Classification metrics}
#'   }
#'
#' @details
#' ## Analysis Methods
#'
#' The function can run up to three analysis types:
#' \itemize{
#'   \item \strong{FS}: ForestSearch with LASSO variable selection only
#'     (use_lasso = TRUE, use_grf = FALSE)
#'   \item \strong{FSlg}: ForestSearch with LASSO + GRF variable selection
#'     (use_lasso = TRUE, use_grf = TRUE)
#'   \item \strong{GRF}: Standalone GRF-based subgroup identification using
#'     grf.subg.harm.survival()
#' }
#'
#' ## Output Column Naming Convention
#'
#' The output distinguishes between:
#' \itemize{
#'   \item \strong{True subgroup from DGM}: sizeH_true, propH_true (known from data generation)
#'   \item \strong{Identified subgroup}: size.H, hr.H.hat (estimated by analysis)
#' }
#'
#' @examples
#' \dontrun{
#' # Create DGM (aligned version)
#' dgm <- create_gbsg_dgm(model = "alt", k_inter = 2, verbose = TRUE)
#'
#' # Run single simulation with LASSO only
#' result <- run_simulation_analysis(
#'   sim_id = 1,
#'   dgm = dgm,
#'   n_sample = 500,
#'   confounders_base = c("v1", "v2", "v3", "v4", "v5", "v6", "v7"),
#'   run_fs = TRUE,      # LASSO only
#'   run_fs_grf = FALSE, # Skip LASSO+GRF
#'   run_grf = FALSE,    # Skip standalone GRF
#'   verbose = TRUE
#' )
#'
#' # Run all three analysis types
#' result_all <- run_simulation_analysis(
#'   sim_id = 1,
#'   dgm = dgm,
#'   n_sample = 500,
#'   run_fs = TRUE,
#'   run_fs_grf = TRUE,
#'   run_grf = TRUE,
#'   verbose = TRUE
#' )
#' # result_all has 3 rows: one for FS, one for FSlg, one for GRF
#'
#' # With use_twostage = TRUE for faster analysis
#' result_fast <- run_simulation_analysis(
#'   sim_id = 1,
#'   dgm = dgm,
#'   n_sample = 500,
#'   fs_params = list(use_twostage = TRUE),
#'   verbose = TRUE
#' )
#' }
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom stats rnorm
#' @export
run_simulation_analysis <- function(
    sim_id,
    dgm,
    n_sample,
    max_follow = Inf,
    muC_adj = 0,
    confounders_base = c("v1", "v2", "v3", "v4", "v5", "v6", "v7"),
    n_add_noise = 0L,
    run_fs = TRUE,
    run_fs_grf = TRUE,
    run_grf = TRUE,
    fs_params = list(),
    grf_params = list(),
    cox_formula = NULL,
    cox_formula_adj = NULL,
    n_sims_total = NULL,
    seed_base = 8316951L,
    verbose = FALSE,
    verbose_n = NULL,
    debug = FALSE
) {
  
  # -------------------------------------------------------------------------
  # Determine effective verbosity based on verbose_n
  # -------------------------------------------------------------------------
  show_verbose <- verbose
  if (verbose && !is.null(verbose_n)) {
    show_verbose <- sim_id <= verbose_n
  }
  
  # -------------------------------------------------------------------------
  # Verbose: Header
  # -------------------------------------------------------------------------
  if (show_verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message(sprintf("Simulation %d", sim_id))
    if (!is.null(n_sims_total)) {
      message(sprintf("  Progress: %d / %d (%.1f%%)",
                      sim_id, n_sims_total, 100 * sim_id / n_sims_total))
    }
    message(paste(rep("=", 60), collapse = ""))
  }
  
  # -------------------------------------------------------------------------
  # Simulate Data
  # -------------------------------------------------------------------------
  if (show_verbose) {
    message("\n[1] Simulating data...")
    message(sprintf("    n_sample = %d, max_follow = %s",
                    n_sample, ifelse(is.infinite(max_follow), "Inf", max_follow)))
  }
  
  sim_data <- simulate_from_gbsg_dgm(
    dgm = dgm,
    n = n_sample,
    sim_id = sim_id,
    max_follow = max_follow,
    muC_adj = muC_adj
  )
  
  if (show_verbose) {
    message(sprintf("    Simulated: n = %d, events = %d (%.1f%%)",
                    nrow(sim_data), sum(sim_data$event.sim),
                    100 * mean(sim_data$event.sim)))
  }
  
  # -------------------------------------------------------------------------
  # Add Noise Variables if Requested
  # -------------------------------------------------------------------------
  confounders_name <- confounders_base
  
  if (n_add_noise > 0) {
    set.seed(seed_base + 1000L * sim_id)
    noise_names <- paste0("noise", seq_len(n_add_noise))
    
    for (nm in noise_names) {
      sim_data[[nm]] <- stats::rnorm(nrow(sim_data))
    }
    confounders_name <- c(confounders_base, noise_names)
    
    if (show_verbose) {
      message(sprintf("    Added %d noise variables", n_add_noise))
    }
  }
  
  # -------------------------------------------------------------------------
  # Compute True Subgroup Properties
  # -------------------------------------------------------------------------
  size_H_true <- sum(sim_data$flag.harm)
  prop_H_true <- mean(sim_data$flag.harm)
  size_Hc_true <- sum(!sim_data$flag.harm)
  prop_Hc_true <- mean(!sim_data$flag.harm)
  
  # Extract HRs using helper (backward compatible)
  hr_H_dgm <- get_dgm_hr(dgm, "hr_H")
  hr_Hc_dgm <- get_dgm_hr(dgm, "hr_Hc")
  ahr_H_dgm <- get_dgm_hr(dgm, "ahr_H")
  ahr_Hc_dgm <- get_dgm_hr(dgm, "ahr_Hc")
  
  if (show_verbose) {
    message("\n[2] True subgroup properties:")
    message(sprintf("    True H:  n = %d (%.1f%%)", size_H_true, 100 * prop_H_true))
    message(sprintf("    True Hc: n = %d (%.1f%%)", size_Hc_true, 100 * prop_Hc_true))
    message(sprintf("    DGM HR_H = %.3f, HR_Hc = %.3f", hr_H_dgm, hr_Hc_dgm))
    if (!is.na(ahr_H_dgm)) {
      message(sprintf("    DGM AHR_H = %.3f, AHR_Hc = %.3f", ahr_H_dgm, ahr_Hc_dgm))
    }
  }
  
  df_pop <- data.table::data.table(
    sim = sim_id,
    sizeH_true = size_H_true,
    propH_true = prop_H_true,
    sizeHc_true = size_Hc_true,
    propHc_true = prop_Hc_true
  )
  
  # -------------------------------------------------------------------------
  # Merge Parameters with Defaults
  # -------------------------------------------------------------------------
  fs_defaults <- default_fs_params()
  fs_merged <- modifyList(fs_defaults, fs_params)
  
  grf_defaults <- default_grf_params()
  grf_merged <- modifyList(grf_defaults, grf_params)
  
  # -------------------------------------------------------------------------
  # Run Analyses
  # -------------------------------------------------------------------------
  results_list <- list()
  
  if (show_verbose) {
    analyses_to_run <- c(
      if (run_fs) "FS (LASSO)" else NULL,
      if (run_fs_grf) "FSlg (LASSO+GRF)" else NULL,
      if (run_grf) "GRF (standalone)" else NULL
    )
    message(sprintf("\n[3] Running analyses: %s", paste(analyses_to_run, collapse = ", ")))
    message(sprintf("    use_twostage = %s", fs_merged$use_twostage))
  }
  
  # ForestSearch with LASSO variable selection only
  if (run_fs) {
    # Explicitly set LASSO only: use_lasso = TRUE, use_grf = FALSE
    fs_params_lasso <- modifyList(fs_merged, list(use_lasso = TRUE, use_grf = FALSE))
    
    fs_result <- run_forestsearch_analysis(
      data = sim_data,
      confounders_name = confounders_name,
      params = fs_params_lasso,
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis_label = "FS",
      verbose = show_verbose
    )
    results_list[["FS"]] <- cbind(df_pop, fs_result)
  }
  
  # ForestSearch with LASSO + GRF variable selection
  if (run_fs_grf) {
    # Explicitly set both: use_lasso = TRUE, use_grf = TRUE
    fs_params_lasso_grf <- modifyList(fs_merged, list(use_lasso = TRUE, use_grf = TRUE))
    
    fs_grf_result <- run_forestsearch_analysis(
      data = sim_data,
      confounders_name = confounders_name,
      params = fs_params_lasso_grf,
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis_label = "FSlg",
      verbose = show_verbose
    )
    results_list[["FSlg"]] <- cbind(df_pop, fs_grf_result)
  }
  
  # Standalone GRF analysis (uses grf.subg.harm.survival)
  if (run_grf) {
    grf_result <- run_grf_analysis(
      data = sim_data,
      confounders_name = confounders_name,
      params = grf_merged,
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis_label = "GRF",
      verbose = show_verbose,
      debug = debug
    )
    results_list[["GRF"]] <- cbind(df_pop, grf_result)
  }
  
  # -------------------------------------------------------------------------
  # Combine Results
  # -------------------------------------------------------------------------
  if (length(results_list) == 0) {
    warning("No analyses were run. Check run_fs, run_fs_grf, run_grf settings.")
    return(NULL)
  }
  
  result <- data.table::rbindlist(results_list, fill = TRUE)
  
  if (show_verbose) {
    message("\n[4] Simulation complete.")
    message(sprintf("    Results: %d rows x %d columns", nrow(result), ncol(result)))
    message(paste(rep("=", 60), collapse = ""), "\n")
  }
  
  result
}


# =============================================================================
# Summary Functions
# =============================================================================

#' Summarize Simulation Results
#'
#' Creates a summary table of operating characteristics across all simulations.
#' Includes both HR and AHR metrics.
#'
#' @param results data.table with simulation results from run_simulation_analysis
#' @param analyses Character vector. Analysis methods to include. Default: all
#' @param digits Integer. Decimal places for proportions. Default: 2
#' @param digits_hr Integer. Decimal places for hazard ratios. Default: 3
#'
#' @return Data frame with summary statistics
#'
#' @export
summarize_simulation_results <- function(
    results,
    analyses = NULL,
    digits = 2,
    digits_hr = 3
) {
  
  if (is.null(analyses)) {
    analyses <- unique(results$analysis)
  }
  
  summaries <- lapply(analyses, function(a) {
    res <- subset(results, analysis == a)
    summarize_single_analysis(res, digits = digits, digits_hr = digits_hr)
  })
  
  summary_df <- do.call(cbind, summaries)
  colnames(summary_df) <- analyses
  
  summary_df
}


#' Summarize Single Analysis Results
#'
#' @param result data.table with results for a single analysis method
#' @param digits Integer. Decimal places for proportions
#' @param digits_hr Integer. Decimal places for hazard ratios
#'
#' @return Data frame with summary statistics
#'
#' @keywords internal
summarize_single_analysis <- function(result, digits = 2, digits_hr = 3) {
  
  # Classification metrics
  class_cols <- c("any.H", "sensitivity", "specificity", "ppv", "npv")
  class_cols <- intersect(class_cols, names(result))
  class_means <- sapply(result[, class_cols, with = FALSE], mean, na.rm = TRUE)
  class_means <- round(class_means, digits)
  
  # Subgroup sizes (only when found)
  res_found <- subset(result, any.H == 1)
  
  if (nrow(res_found) > 0) {
    avg_H <- round(mean(res_found$size.H, na.rm = TRUE), 0)
    min_H <- round(min(res_found$size.H, na.rm = TRUE), 0)
    max_H <- round(max(res_found$size.H, na.rm = TRUE), 0)
  } else {
    avg_H <- min_H <- max_H <- NA
  }
  
  avg_Hc <- round(mean(result$size.Hc, na.rm = TRUE), 0)
  min_Hc <- round(min(result$size.Hc, na.rm = TRUE), 0)
  max_Hc <- round(max(result$size.Hc, na.rm = TRUE), 0)
  
  # HR estimates (only when found)
  if (nrow(res_found) > 0) {
    hr_H_true <- round(mean(res_found$hr.H.true, na.rm = TRUE), digits_hr)
    hr_H_hat <- round(mean(res_found$hr.H.hat, na.rm = TRUE), digits_hr)
    hr_Hc_true <- round(mean(res_found$hr.Hc.true, na.rm = TRUE), digits_hr)
    hr_Hc_hat <- round(mean(res_found$hr.Hc.hat, na.rm = TRUE), digits_hr)
    
    # AHR estimates (aligned)
    ahr_H_hat <- if ("ahr.H.hat" %in% names(res_found)) {
      round(mean(res_found$ahr.H.hat, na.rm = TRUE), digits_hr)
    } else NA
    ahr_Hc_hat <- if ("ahr.Hc.hat" %in% names(res_found)) {
      round(mean(res_found$ahr.Hc.hat, na.rm = TRUE), digits_hr)
    } else NA
  } else {
    hr_H_true <- hr_H_hat <- hr_Hc_true <- hr_Hc_hat <- NA
    ahr_H_hat <- ahr_Hc_hat <- NA
  }
  
  # Overall HR estimates (all simulations)
  hr_H_true_all <- round(mean(result$hr.H.true, na.rm = TRUE), digits_hr)
  hr_Hc_true_all <- round(mean(result$hr.Hc.true, na.rm = TRUE), digits_hr)
  hr_itt_all <- round(mean(result$hr.itt, na.rm = TRUE), digits_hr)
  hr_adj_itt_all <- if ("hr.adj.itt" %in% names(result)) {
    round(mean(result$hr.adj.itt, na.rm = TRUE), digits_hr)
  } else NA
  
  # Combine
  values <- c(
    class_means,
    avg_H, min_H, max_H, avg_Hc, min_Hc, max_Hc,
    hr_H_true, hr_H_hat, hr_Hc_true, hr_Hc_hat,
    ahr_H_hat, ahr_Hc_hat,
    hr_H_true_all, hr_Hc_true_all, hr_itt_all, hr_adj_itt_all
  )
  
  row_names <- c(
    names(class_means),
    "Avg(#H)", "minH", "maxH", "Avg(#Hc)", "minHc", "maxHc",
    "hat(H*)", "hat(hat[H])", "hat(Hc*)", "hat(hat[Hc])",
    "hat(AHR_H)", "hat(AHR_Hc)",
    "hat(H*)all", "hat(Hc*)all", "hat(ITT)all", "hat(ITTadj)all"
  )
  
  out <- data.frame(value = values, stringsAsFactors = FALSE)
  rownames(out) <- row_names
  
  out
}


# =============================================================================
# Format Operating Characteristics Table
# =============================================================================

#' Format Operating Characteristics Results as GT Table
#'
#' Creates a formatted gt table from simulation operating characteristics results.
#'
#' @param results data.table or data.frame. Simulation results from
#'   \code{\link{run_simulation_analysis}} or combined results from multiple simulations.
#' @param analyses Character vector. Analysis methods to include.
#'   Default: NULL (all analyses in results)
#' @param metrics Character vector. Metrics to display. Options include:
#'   "detection", "classification", "hr_estimates", "subgroup_size", "all".
#'   Default: "all"
#' @param digits Integer. Decimal places for proportions. Default: 3
#' @param digits_hr Integer. Decimal places for hazard ratios. Default: 3
#' @param title Character. Table title. Default: "Operating Characteristics Summary"
#' @param subtitle Character. Table subtitle. Default: NULL
#' @param use_gt Logical. Return gt table if TRUE, data.frame if FALSE. Default: TRUE
#'
#' @return A gt table object (if use_gt = TRUE and gt package available) or data.frame
#'
#' @details
#' The function summarizes simulation results across multiple metrics:
#' \itemize{
#'   \item \strong{Detection}: Proportion of simulations finding a subgroup (any.H)
#'   \item \strong{Classification}: Sensitivity, specificity, PPV, NPV
#'   \item \strong{HR Estimates}: Mean hazard ratios in H and Hc subgroups
#'   \item \strong{Subgroup Size}: Average, min, max sizes
#' }
#'
#' @importFrom data.table is.data.table as.data.table
#' @export
format_oc_results <- function(
    results,
    analyses = NULL,
    metrics = "all",
    digits = 3,
    digits_hr = 3,
    title = "Operating Characteristics Summary",
    subtitle = NULL,
    use_gt = TRUE
) {
  
  # Convert to data.table if needed
  if (!data.table::is.data.table(results)) {
    results <- data.table::as.data.table(results)
  }
  
  # Get analyses if not specified
  if (is.null(analyses)) {
    analyses <- unique(results$analysis)
  }
  
  # Filter to requested analyses
  results <- results[results$analysis %in% analyses, ]
  
  # Compute summary statistics for each analysis
  summary_list <- lapply(analyses, function(a) {
    res <- results[results$analysis == a, ]
    n_sims <- nrow(res)
    
    # Detection rate
    detection_rate <- mean(res$any.H, na.rm = TRUE)
    
    # Classification metrics (averaged across all sims)
    sensitivity <- mean(res$sensitivity, na.rm = TRUE)
    specificity <- mean(res$specificity, na.rm = TRUE)
    ppv <- mean(res$ppv, na.rm = TRUE)
    npv <- mean(res$npv, na.rm = TRUE)
    
    # HR estimates (only when subgroup found)
    res_found <- res[res$any.H == 1, ]
    if (nrow(res_found) > 0) {
      hr_H_hat <- mean(res_found$hr.H.hat, na.rm = TRUE)
      hr_Hc_hat <- mean(res_found$hr.Hc.hat, na.rm = TRUE)
      hr_H_true <- mean(res_found$hr.H.true, na.rm = TRUE)
      hr_Hc_true <- mean(res_found$hr.Hc.true, na.rm = TRUE)
      size_H_mean <- mean(res_found$size.H, na.rm = TRUE)
      size_H_min <- min(res_found$size.H, na.rm = TRUE)
      size_H_max <- max(res_found$size.H, na.rm = TRUE)
    } else {
      hr_H_hat <- hr_Hc_hat <- hr_H_true <- hr_Hc_true <- NA
      size_H_mean <- size_H_min <- size_H_max <- NA
    }
    
    # ITT estimate (all simulations)
    hr_itt <- mean(res$hr.itt, na.rm = TRUE)
    
    data.frame(
      Analysis = a,
      N_sims = n_sims,
      Detection = detection_rate,
      Sensitivity = sensitivity,
      Specificity = specificity,
      PPV = ppv,
      NPV = npv,
      HR_H_hat = hr_H_hat,
      HR_Hc_hat = hr_Hc_hat,
      HR_H_true = hr_H_true,
      HR_Hc_true = hr_Hc_true,
      HR_ITT = hr_itt,
      Size_H_mean = size_H_mean,
      Size_H_min = size_H_min,
      Size_H_max = size_H_max,
      stringsAsFactors = FALSE
    )
  })
  
  summary_df <- do.call(rbind, summary_list)
  
  # Filter columns based on metrics
  if (!"all" %in% metrics) {
    cols_to_keep <- c("Analysis", "N_sims")
    if ("detection" %in% metrics) {
      cols_to_keep <- c(cols_to_keep, "Detection")
    }
    if ("classification" %in% metrics) {
      cols_to_keep <- c(cols_to_keep, "Sensitivity", "Specificity", "PPV", "NPV")
    }
    if ("hr_estimates" %in% metrics) {
      cols_to_keep <- c(cols_to_keep, "HR_H_hat", "HR_Hc_hat", "HR_H_true", "HR_Hc_true", "HR_ITT")
    }
    if ("subgroup_size" %in% metrics) {
      cols_to_keep <- c(cols_to_keep, "Size_H_mean", "Size_H_min", "Size_H_max")
    }
    summary_df <- summary_df[, cols_to_keep, drop = FALSE]
  }
  
  # Format as gt table if requested and available
  if (use_gt && requireNamespace("gt", quietly = TRUE)) {
    gt_table <- gt::gt(summary_df)
    
    # Add title
    gt_table <- gt::tab_header(
      gt_table,
      title = title,
      subtitle = subtitle
    )
    
    # Format numeric columns
    numeric_cols <- setdiff(names(summary_df), c("Analysis", "N_sims"))
    
    # Proportion columns (0-1 scale)
    prop_cols <- intersect(c("Detection", "Sensitivity", "Specificity", "PPV", "NPV"),
                           names(summary_df))
    if (length(prop_cols) > 0) {
      gt_table <- gt::fmt_number(
        gt_table,
        columns = gt::all_of(prop_cols),
        decimals = digits
      )
    }
    
    # HR columns
    hr_cols <- intersect(c("HR_H_hat", "HR_Hc_hat", "HR_H_true", "HR_Hc_true", "HR_ITT"),
                         names(summary_df))
    if (length(hr_cols) > 0) {
      gt_table <- gt::fmt_number(
        gt_table,
        columns = gt::all_of(hr_cols),
        decimals = digits_hr
      )
    }
    
    # Size columns
    size_cols <- intersect(c("Size_H_mean", "Size_H_min", "Size_H_max"),
                           names(summary_df))
    if (length(size_cols) > 0) {
      gt_table <- gt::fmt_number(
        gt_table,
        columns = gt::all_of(size_cols),
        decimals = 0
      )
    }
    
    # Rename columns for display
    gt_table <- gt::cols_label(
      gt_table,
      Analysis = "Method",
      N_sims = "N Sims"
    )
    
    # Add column spanners
    if ("all" %in% metrics || "classification" %in% metrics) {
      class_cols <- intersect(c("Sensitivity", "Specificity", "PPV", "NPV"), names(summary_df))
      if (length(class_cols) > 0) {
        gt_table <- gt::tab_spanner(
          gt_table,
          label = "Classification",
          columns = gt::all_of(class_cols)
        )
      }
    }
    
    if ("all" %in% metrics || "hr_estimates" %in% metrics) {
      hr_cols <- intersect(c("HR_H_hat", "HR_Hc_hat", "HR_H_true", "HR_Hc_true", "HR_ITT"),
                           names(summary_df))
      if (length(hr_cols) > 0) {
        gt_table <- gt::tab_spanner(
          gt_table,
          label = "Hazard Ratios",
          columns = gt::all_of(hr_cols)
        )
      }
    }
    
    # Style
    gt_table <- gt::tab_style(
      gt_table,
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels()
    )
    
    # Handle missing values
    if (utils::packageVersion("gt") >= "0.6.0") {
      gt_table <- gt::sub_missing(
        gt_table,
        columns = gt::everything(),
        missing_text = "-"
      )
    }
    
    return(gt_table)
  }
  
  # Return data.frame if gt not available or not requested
  summary_df
}
