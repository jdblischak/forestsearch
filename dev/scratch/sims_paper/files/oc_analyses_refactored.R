# =============================================================================
# Operating Characteristics Analysis Functions (Refactored)
# =============================================================================
#
# Refactored simulation analysis functions for evaluating ForestSearch and
# GRF subgroup identification methods. Aligned with forestsearch package.
#
# Key functions:
#   - run_simulation_analysis(): Main simulation wrapper
#   - run_forestsearch_analysis(): ForestSearch analysis helper
#   - run_grf_analysis(): GRF analysis helper
#   - extract_fs_estimates(): Extract estimates from ForestSearch results
#   - extract_grf_estimates(): Extract estimates from GRF results
#   - summarize_simulation_results(): Create summary tables
#   - compute_theoretical_power(): Analytical power approximation
#
# =============================================================================

# =============================================================================
# Global Variables Declaration
# =============================================================================

utils::globalVariables(c(
  # Simulation results columns
  "any.H", "ppv", "npv", "sensitivity", "specificity",
  "size.H", "size.Hc", "hr.H.true", "hr.H.hat", "hr.Hc.true", "hr.Hc.hat",
  "hr.itt", "hr.adj.itt", "p.cens", "taumax",
  
  # Analysis labels
  "analysis",
  
  # Loop variables
  "sim", "aa"
))


# =============================================================================
# Configuration Defaults
# =============================================================================

#' Default ForestSearch Parameters
#' @keywords internal
default_fs_params <- function() {
  list(
    outcome.name = "y.sim",
    event.name = "event.sim",
    treat.name = "treat",
    id.name = "id",
    use_lasso = TRUE,
    use_grf = FALSE,
    hr.threshold = 1.25,
    hr.consistency = 1.0,
    pconsistency.threshold = 0.90,
    stop.threshold = 0.90,
    fs.splits = 300,
    n.min = 60,
    d0.min = 20,
    d1.min = 20,
    maxk = 2,
    max.minutes = 5,
    pstop_futile = 0.10,
    by.risk = 12,
    vi.grf.min = NULL
  )
}


#' Default GRF Parameters
#' @keywords internal
default_grf_params <- function() {
  list(
    outcome.name = "y.sim",
    event.name = "event.sim",
    treat.name = "treat",
    id.name = "id",
    n.min = 60,
    dmin.grf = 20,
    frac.tau = 1.0,
    maxdepth = 2
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
#' @param sim_id Integer. Simulation index for seed offset and tracking
#' @param dgm A DGM object from \code{\link{create_gbsg_dgm}} or similar
#' @param n_sample Integer. Sample size for simulation
#' @param max_follow Numeric. Maximum follow-up time. Default: Inf
#' @param muC_adj Numeric. Censoring adjustment. Default: 0
#' @param confounders_base Character vector. Base confounder names
#' @param n_add_noise Integer. Number of noise variables to add. Default: 0
#' @param run_fs Logical. Run ForestSearch analysis. Default: TRUE
#' @param run_fs_grf Logical. Run ForestSearch with GRF. Default: TRUE
#' @param run_grf Logical. Run standalone GRF analysis. Default: TRUE
#' @param fs_params List. ForestSearch parameters (overrides defaults)
#' @param grf_params List. GRF parameters (overrides defaults)
#' @param cox_formula Formula. Cox model formula for estimation
#' @param cox_formula_adj Formula. Adjusted Cox model formula
#' @param n_sims_total Integer. Total simulations (for progress display)
#' @param seed_base Integer. Base random seed. Default: 8316951
#' @param verbose Logical. Print progress. Default: FALSE
#'
#' @return A data.table with analysis results for all requested methods
#'
#' @details
#' ## Analysis Methods
#'
#' The function can run up to three analysis types:
#' \itemize{
#'   \item \strong{FS}: ForestSearch with LASSO only
#'   \item \strong{FSlg}: ForestSearch with LASSO + GRF variable selection
#'   \item \strong{GRF}: Standalone GRF-based subgroup identification
#' }
#'
#' ## Noise Variables
#'
#' If \code{n_add_noise > 0}, random normal noise variables are added to test
#' robustness against irrelevant covariates.
#'
#' @examples
#' \dontrun{
#' # Create DGM
#' dgm <- create_gbsg_dgm(model = "alt", k_inter = 2)
#'
#' # Run single simulation
#' result <- run_simulation_analysis(
#'   sim_id = 1,
#'   dgm = dgm,
#'   n_sample = 500,
#'   confounders_base = c("z1", "z2", "z3", "z4", "z5", "size", "grade3"),
#'   run_fs = TRUE,
#'   run_grf = TRUE,
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
    confounders_base = c("z1", "z2", "z3", "z4", "z5", "size", "grade3"),
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
    verbose = FALSE
) {
  
  # -------------------------------------------------------------------------
  # Simulate Data
  # -------------------------------------------------------------------------
  sim_data <- simulate_from_gbsg_dgm(
    dgm = dgm,
    n = n_sample,
    sim_id = sim_id,
    max_follow = max_follow,
    muC_adj = muC_adj
  )
  
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
  }
  
  # -------------------------------------------------------------------------
  # Compute True Subgroup Properties
  # -------------------------------------------------------------------------
  size_H_true <- sum(sim_data$flag.harm)
  prop_H_true <- mean(sim_data$flag.harm)
  size_Hc_true <- sum(!sim_data$flag.harm)
  prop_Hc_true <- mean(!sim_data$flag.harm)
  
  df_pop <- data.table::data.table(
    sim = sim_id,
    sizeH_true = size_H_true,
    propH_true = prop_H_true,
    sizeHc_true = size_Hc_true,
    propHc_true = prop_Hc_true
  )
  
  # -------------------------------------------------------------------------
  # Determine Verbosity for This Simulation
  # -------------------------------------------------------------------------
  show_details <- FALSE
  if (verbose && !is.null(n_sims_total)) {
    if (sim_id <= 1 || sim_id >= (n_sims_total - 1)) {
      show_details <- TRUE
    }
  }
  
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
  
  # ForestSearch (LASSO only)
  if (run_fs) {
    fs_result <- run_forestsearch_analysis(
      data = sim_data,
      confounders_name = confounders_name,
      params = modifyList(fs_merged, list(use_lasso = TRUE, use_grf = FALSE)),
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis_label = "FS",
      verbose = show_details
    )
    results_list[["FS"]] <- cbind(df_pop, fs_result)
  }
  
  # ForestSearch (LASSO + GRF)
  if (run_fs_grf) {
    fs_grf_result <- run_forestsearch_analysis(
      data = sim_data,
      confounders_name = confounders_name,
      params = modifyList(fs_merged, list(use_lasso = TRUE, use_grf = TRUE)),
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis_label = "FSlg",
      verbose = show_details
    )
    results_list[["FSlg"]] <- cbind(df_pop, fs_grf_result)
  }
  
  # Standalone GRF
  if (run_grf) {
    grf_result <- run_grf_analysis(
      data = sim_data,
      confounders_name = confounders_name,
      params = grf_merged,
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis_label = "GRF",
      verbose = show_details
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
  
  data.table::rbindlist(results_list, fill = TRUE)
}


# =============================================================================
# ForestSearch Analysis Helper
# =============================================================================

#' Run ForestSearch Analysis
#'
#' Helper function to run ForestSearch and extract estimates.
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
  
  # Build forestsearch call arguments
  fs_args <- list(
    df.analysis = data,
    Allconfounders.name = confounders_name,
    details = verbose,
    plot.sg = verbose
  )
  
  # Add parameters
  param_names <- c(
    "outcome.name", "event.name", "treat.name", "id.name",
    "use_lasso", "use_grf", "conf_force",
    "hr.threshold", "hr.consistency", "pconsistency.threshold",
    "stop.threshold", "fs.splits", "n.min", "d0.min", "d1.min",
    "maxk", "max.minutes", "pstop_futile", "by.risk", "vi.grf.min",
    "frac.tau", "dmin.grf", "grf_depth"
  )
  
  for (pn in param_names) {
    if (!is.null(params[[pn]])) {
      fs_args[[pn]] <- params[[pn]]
    }
  }
  
  # Run ForestSearch
  fs_result <- tryCatch({
    do.call(forestsearch, fs_args)
  }, error = function(e) {
    if (verbose) {
      warning(sprintf("%s analysis failed: %s", analysis_label, e$message))
    }
    NULL
  })
  
  # Extract estimates
  if (!is.null(fs_result) && !is.null(fs_result$grp.consistency$out_hr)) {
    extract_fs_estimates(
      df = data,
      fs_res = fs_result$grp.consistency$out_hr,
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis = analysis_label
    )
  } else {
    extract_fs_estimates(
      df = data,
      fs_res = NULL,
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis = analysis_label
    )
  }
}


# =============================================================================
# GRF Analysis Helper
# =============================================================================

#' Run GRF Analysis
#'
#' Helper function to run GRF-based subgroup identification and extract estimates.
#'
#' @param data Data frame with simulated trial data
#' @param confounders_name Character vector of confounder names
#' @param params List of GRF parameters
#' @param dgm DGM object for computing true HRs
#' @param cox_formula Cox formula for estimation
#' @param cox_formula_adj Adjusted Cox formula
#' @param analysis_label Character label for this analysis
#' @param verbose Print details
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
    verbose = FALSE
) {
  
  # Check if grf_subg_harm_survival exists
  if (!exists("grf_subg_harm_survival", mode = "function")) {
    warning("grf_subg_harm_survival function not found. Skipping GRF analysis.")
    return(extract_grf_estimates(
      df = data,
      grf_est = NULL,
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis = analysis_label
    ))
  }
  
  # Run GRF analysis
  grf_result <- tryCatch({
    grf_subg_harm_survival(
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
      warning(sprintf("%s analysis failed: %s", analysis_label, e$message))
    }
    NULL
  })
  
  # Extract estimates
  extract_grf_estimates(
    df = data,
    grf_est = grf_result,
    dgm = dgm,
    cox_formula = cox_formula,
    cox_formula_adj = cox_formula_adj,
    analysis = analysis_label,
    frac_tau = params$frac.tau
  )
}


# =============================================================================
# Estimate Extraction Functions
# =============================================================================

#' Extract Estimates from ForestSearch Results
#'
#' @param df Simulated data frame
#' @param fs_res ForestSearch results object (or NULL if failed)
#' @param dgm DGM object
#' @param cox_formula Cox formula
#' @param cox_formula_adj Adjusted Cox formula
#' @param analysis Analysis label
#'
#' @return data.table with extracted estimates
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
    analysis = "FS"
) {
  
  # Initialize output
  out <- data.table::data.table(
    analysis = analysis,
    any.H = 0L,
    size.H = NA_integer_,
    size.Hc = nrow(df),
    hr.H.true = NA_real_,
    hr.H.hat = NA_real_,
    hr.Hc.true = dgm$hr_Hc_true,
    hr.Hc.hat = NA_real_,
    hr.itt = NA_real_,
    hr.adj.itt = NA_real_,
    ppv = NA_real_,
    npv = NA_real_,
    sensitivity = NA_real_,
    specificity = NA_real_,
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
    out$hr.Hc.hat <- out$hr.itt
    return(out)
  }
  
  # Subgroup found - extract details
  out$any.H <- 1L
  
  # Get subgroup definition
  sg_def <- fs_res[1, ]
  
  # Create subgroup indicator in data
  # This would need the actual subgroup definition from fs_res
  # For now, placeholder logic
  sg_factors <- character(0)
  for (i in seq_len(min(7, ncol(fs_res)))) {
    col_name <- paste0("M.", i)
    if (col_name %in% names(fs_res) && !is.na(fs_res[[col_name]][1])) {
      sg_factors <- c(sg_factors, fs_res[[col_name]][1])
    }
  }
  
  # Apply subgroup definition
  if (length(sg_factors) > 0) {
    # Convert factor definitions to indicator
    df$sg_hat <- create_subgroup_indicator(df, sg_factors)
    
    out$size.H <- sum(df$sg_hat, na.rm = TRUE)
    out$size.Hc <- sum(!df$sg_hat, na.rm = TRUE)
    
    # Compute HRs
    if (out$size.H > 10) {
      out$hr.H.hat <- tryCatch({
        exp(survival::coxph(
          survival::Surv(y.sim, event.sim) ~ treat,
          data = subset(df, sg_hat == 1)
        )$coefficients)
      }, error = function(e) NA_real_)
    }
    
    if (out$size.Hc > 10) {
      out$hr.Hc.hat <- tryCatch({
        exp(survival::coxph(
          survival::Surv(y.sim, event.sim) ~ treat,
          data = subset(df, sg_hat == 0)
        )$coefficients)
      }, error = function(e) NA_real_)
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
      
      # True HR in estimated subgroup
      out$hr.H.true <- tryCatch({
        exp(survival::coxph(
          survival::Surv(y.sim, event.sim) ~ treat,
          data = subset(df, sg_hat == 1)
        )$coefficients)
      }, error = function(e) NA_real_)
    }
  }
  
  out
}


#' Extract Estimates from GRF Results
#'
#' @param df Simulated data frame
#' @param grf_est GRF results object (or NULL if failed)
#' @param dgm DGM object
#' @param cox_formula Cox formula
#' @param cox_formula_adj Adjusted Cox formula
#' @param analysis Analysis label
#' @param frac_tau Fraction of tau used
#'
#' @return data.table with extracted estimates
#'
#' @importFrom data.table data.table
#' @keywords internal

extract_grf_estimates <- function(
    df,
    grf_est,
    dgm,
    cox_formula = NULL,
    cox_formula_adj = NULL,
    analysis = "GRF",
    frac_tau = 1.0
) {
  
  # Use same structure as FS estimates
  out <- extract_fs_estimates(
    df = df,
    fs_res = NULL,
    dgm = dgm,
    cox_formula = cox_formula,
    cox_formula_adj = cox_formula_adj,
    analysis = analysis
  )
  
  # If GRF result available, extract additional info
  if (!is.null(grf_est) && !is.null(grf_est$harm_indicator)) {
    out$any.H <- 1L
    out$size.H <- sum(grf_est$harm_indicator, na.rm = TRUE)
    out$size.Hc <- sum(!grf_est$harm_indicator, na.rm = TRUE)
    
    # Add GRF-specific estimates if available
    if (!is.null(grf_est$hr_harm)) {
      out$hr.H.hat <- grf_est$hr_harm
    }
    if (!is.null(grf_est$hr_no_harm)) {
      out$hr.Hc.hat <- grf_est$hr_no_harm
    }
  }
  
  out
}


#' Create Subgroup Indicator from Factor Definitions
#'
#' @param df Data frame
#' @param sg_factors Character vector of factor definitions (e.g., "v1.1", "v3.1")
#'
#' @return Logical vector indicating subgroup membership
#'
#' @keywords internal

create_subgroup_indicator <- function(df, sg_factors) {
  
  # Initialize all TRUE
  indicator <- rep(TRUE, nrow(df))
  
  for (factor_def in sg_factors) {
    if (is.na(factor_def) || factor_def == "") next
    
    # Parse factor definition (e.g., "v1.1" means v1 == 1)
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


# =============================================================================
# Summary Functions
# =============================================================================

#' Summarize Simulation Results
#'
#' Creates summary statistics from simulation results across multiple methods.
#'
#' @param results data.table with combined simulation results
#' @param analyses Character vector of analysis labels to summarize.
#'   If NULL, uses all unique values in results$analysis
#' @param digits Integer. Decimal places for rates. Default: 2
#' @param digits_hr Integer. Decimal places for hazard ratios. Default: 3
#'
#' @return data.frame with summary statistics
#'
#' @examples
#' \dontrun{
#' summary_table <- summarize_simulation_results(results)
#' print(summary_table)
#' }
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
  
  # Process each analysis
  summaries <- lapply(analyses, function(a) {
    res <- subset(results, analysis == a)
    summarize_single_analysis(res, digits = digits, digits_hr = digits_hr)
  })
  
  # Combine into single data frame
  summary_df <- do.call(cbind, summaries)
  colnames(summary_df) <- analyses
  
  summary_df
}


#' Summarize Single Analysis Results
#'
#' @param result data.table with results for a single analysis type
#' @param digits Decimal places for rates
#' @param digits_hr Decimal places for HRs
#'
#' @return data.frame with summary statistics
#'
#' @keywords internal

summarize_single_analysis <- function(result, digits = 2, digits_hr = 3) {
  
  # Classification metrics
  class_cols <- c("any.H", "ppv", "npv", "sensitivity", "specificity")
  class_means <- sapply(result[, class_cols, with = FALSE], mean, na.rm = TRUE)
  class_means <- round(class_means, digits)
  
  # Size statistics (among found subgroups)
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
  
  # HR statistics (among found)
  if (nrow(res_found) > 0) {
    hr_H_true <- round(mean(res_found$hr.H.true, na.rm = TRUE), digits_hr)
    hr_H_hat <- round(mean(res_found$hr.H.hat, na.rm = TRUE), digits_hr)
    hr_Hc_true <- round(mean(res_found$hr.Hc.true, na.rm = TRUE), digits_hr)
    hr_Hc_hat <- round(mean(res_found$hr.Hc.hat, na.rm = TRUE), digits_hr)
  } else {
    hr_H_true <- hr_H_hat <- hr_Hc_true <- hr_Hc_hat <- NA
  }
  
  # Overall HR statistics
  hr_H_true_all <- round(mean(result$hr.H.true, na.rm = TRUE), digits_hr)
  hr_Hc_true_all <- round(mean(result$hr.Hc.true, na.rm = TRUE), digits_hr)
  hr_itt_all <- round(mean(result$hr.itt, na.rm = TRUE), digits_hr)
  hr_adj_itt_all <- round(mean(result$hr.adj.itt, na.rm = TRUE), digits_hr)
  
  # Combine
  out <- data.frame(
    value = c(
      class_means,
      avg_H, min_H, max_H, avg_Hc, min_Hc, max_Hc,
      hr_H_true, hr_H_hat, hr_Hc_true, hr_Hc_hat,
      hr_H_true_all, hr_Hc_true_all, hr_itt_all, hr_adj_itt_all
    ),
    stringsAsFactors = FALSE
  )
  
  rownames(out) <- c(
    "any.H", "sensH", "sensHc", "ppH", "ppHc",
    "Avg(#H)", "minH", "maxH", "Avg(#Hc)", "minHc", "maxHc",
    "hat(H*)", "hat(hat[H])", "hat(Hc*)", "hat(hat[Hc])",
    "hat(H*)all", "hat(Hc*)all", "hat(ITT)all", "hat(ITTadj)all"
  )
  
  out
}


# =============================================================================
# Theoretical Power Approximation
# =============================================================================

#' Compute Theoretical Power for Subgroup Detection
#'
#' Calculates the approximate probability of detecting a subgroup with
#' hazard ratio below threshold based on asymptotic normal approximation.
#'
#' @param n_sg Integer. Subgroup sample size
#' @param prop_cens Numeric. Proportion censored (0-1). Default: 0.3
#' @param theta Numeric. True hazard ratio in subgroup
#' @param hr_threshold Numeric. Detection threshold for HR
#' @param hr_consistency Numeric. Consistency threshold for HR. Default: same as hr_threshold
#' @param alpha Numeric. One-sided significance level. Default: 0.025
#' @param verbose Logical. Print result. Default: FALSE
#'
#' @return Numeric probability of detection
#'
#' @details
#' Uses the approximation based on the asymptotic variance of the log hazard
#' ratio estimator: Var(log(HR)) ≈ 4/d, where d is the number of events.
#'
#' @examples
#' \dontrun
#' # Power to detect HR = 1.5 with n = 100 in subgroup
#' compute_theoretical_power(n_sg = 100, theta = 1.5, hr_threshold = 1.25)
#' }
#'
#' @importFrom stats pnorm qnorm
#' @export

compute_theoretical_power <- function(
    n_sg,
    prop_cens = 0.3,
    theta,
    hr_threshold,
    hr_consistency = hr_threshold,
    alpha = 0.025,
    verbose = FALSE
) {
  
  # Number of events
  d_sg <- n_sg * (1 - prop_cens)
  
  # Power calculation
  power <- 1 - stats::pnorm(
    sqrt(d_sg / 4) * (log(hr_threshold) - log(theta) - 
                       (2 / sqrt(d_sg)) * stats::qnorm(1 - alpha))
  )
  
  if (verbose) {
    message(sprintf("Power approximation: %.3f", power))
    message(sprintf("  n_sg = %d, events = %.0f, theta = %.2f, threshold = %.2f",
                    n_sg, d_sg, theta, hr_threshold))
  }
  
  power
}


#' Joint Density for Two-Split Detection
#'
#' Computes the joint density for detecting a subgroup in both splits
#' of a cross-validation procedure.
#'
#' @param x Numeric vector of length 2 (log HRs in two splits)
#' @param theta Numeric. True hazard ratio
#' @param prop_cens Numeric. Proportion censored
#' @param n_sg Integer. Subgroup sample size
#' @param k1 Numeric. Log of mean HR threshold
#' @param k2 Numeric. Log of consistency HR threshold
#'
#' @return Numeric density value
#'
#' @keywords internal

joint_density_two_splits <- function(x, theta, prop_cens = 0.3, n_sg, k1, k2) {
  
  d_sg <- n_sg * (1 - prop_cens)
  mu_split <- log(theta)
  sig2_split <- 8 / d_sg
  sig_split <- sqrt(sig2_split)
  
  dens_1 <- stats::dnorm(x[1], mean = mu_split, sd = sig_split)
  dens_2 <- stats::dnorm(x[2], mean = mu_split, sd = sig_split)
  
  # Both splits must pass thresholds
  indicator <- ((x[1] + x[2]) >= 2 * k1) * (x[1] >= k2) * (x[2] >= k2)
  
  indicator * dens_1 * dens_2
}


# =============================================================================
# Output Functions
# =============================================================================

#' Save Simulation Results with Metadata
#'
#' @param results data.table with simulation results
#' @param dgm DGM object used for simulations
#' @param summary_table Summary statistics table
#' @param runtime_hours Numeric. Total runtime in hours
#' @param output_file Character. Path to output .Rdata file
#' @param power_approx List with theoretical power approximations (optional
#'
#' @export

save_simulation_results <- function(
    results,
    dgm,
    summary_table,
    runtime_hours,
    output_file,
    power_approx = NULL
) {
  
  save(
    results,
    dgm,
    summary_table,
    runtime_hours,
    power_approx,
    file = output_file
  )
  
  message(sprintf("Results saved to: %s", output_file))
}


# =============================================================================
# Kable/gt Table Generation
# =============================================================================

#' Create Formatted Summary Table
#'
#' @param summary_df Summary data frame from summarize_simulation_results
#' @param n_sims Integer. Number of simulations
#' @param model_type Character. "alt" or "null"
#' @param power_approx Numeric. Theoretical power approximation (optional)
#' @param row_subset Character vector. Rows to include in table
#'
#' @return gt table or kable object
#'
#' @export

create_summary_kable <- function(
    summary_df,
    n_sims,
    model_type = "alt",
    power_approx = NULL,
    row_subset = c("any.H", "sensH", "sensHc", "ppH", "ppHc",
                   "Avg(#H)", "minH", "maxH", "Avg(#Hc)", "minHc", "maxHc")
) {
  
  # Subset rows
  summary_df <- summary_df[row_subset, , drop = FALSE]
  
  # Format for display
  summary_df <- format(summary_df, digits = 3, drop0trailing = TRUE)
  
  # Set sensitivity to NA for null model
  if (model_type == "null") {
    summary_df["sensH", ] <- NA
  }
  
  # Create table
  if (requireNamespace("gt", quietly = TRUE)) {
    
    tbl <- gt::gt(summary_df, rownames_to_stub = TRUE) |>
      gt::tab_header(
        title = "Simulation Results Summary",
        subtitle = sprintf("Based on %d simulations", n_sims)
      ) |>
      gt::tab_row_group(
        label = "Finding H",
        rows = 1:5
      ) |>
      gt::tab_row_group(
        label = "Size of H and H-complement",
        rows = 6:11
      )
    
    if (!is.null(power_approx)) {
      tbl <- gt::tab_footnote(
        tbl,
        footnote = sprintf("Theoretical power approximation: %.3f", power_approx)
      )
    }
    
    tbl <- gt::tab_footnote(
      tbl,
      footnote = sprintf("Number of simulations: %d", n_sims)
    )
    
    return(tbl)
    
  } else if (requireNamespace("kableExtra", quietly = TRUE)) {
    
    kableExtra::kbl(
      summary_df,
      format = "latex",
      booktabs = TRUE,
      caption = sprintf("Simulation Results (n = %d)", n_sims)
    ) |>
      kableExtra::kable_styling("striped", full_width = FALSE) |>
      kableExtra::group_rows("Finding H", 1, 5) |>
      kableExtra::group_rows("Size of H and H-complement", 6, 11)
    
  } else {
    
    summary_df
  }
}
