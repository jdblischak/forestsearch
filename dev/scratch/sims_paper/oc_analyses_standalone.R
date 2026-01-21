# =============================================================================
# STANDALONE Operating Characteristics Analysis Functions
# =============================================================================
#
# This file contains ALL functions needed for simulation analysis.
# It does NOT depend on any package version of these functions.
# Only requires: forestsearch package for forestsearch() function
#
# Usage:
#   1. Load forestsearch package: devtools::load_all("path/to/forestsearch")
#   2. Source THIS file: source("oc_analyses_standalone.R")
#   3. Run simulations
#
# =============================================================================

library(survival)
library(data.table)

cat("Loading STANDALONE Operating Characteristics Analysis Functions...\n")

# =============================================================================
# Global Variables Declaration
# =============================================================================

utils::globalVariables(c(
  "any.H", "ppv", "npv", "sensitivity", "specificity",
  "size.H", "size.Hc", "hr.H.true", "hr.H.hat", "hr.Hc.true", "hr.Hc.hat",
  "hr.itt", "hr.adj.itt", "p.cens", "taumax",
  "analysis", "sim", "aa", "sg_hat"
))


# =============================================================================
# Configuration Defaults - NO stop.threshold or pstop_futile!
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
    fs.splits = 300,
    n.min = 60,
    d0.min = 20,
    d1.min = 20,
    maxk = 2,
    max.minutes = 5,
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
# Helper: Create Subgroup Indicator
# =============================================================================

#' Create Subgroup Indicator from Factor Definitions
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


# =============================================================================
# Extract ForestSearch Estimates
# =============================================================================

#' Extract Estimates from ForestSearch Results
#' @keywords internal
extract_fs_estimates <- function(
    df,
    fs_res,
    dgm,
    cox_formula = NULL,
    cox_formula_adj = NULL,
    analysis = "FS",
    fs_full = NULL
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

  # Subgroup found
  out$any.H <- 1L

  # Try to use df.est from full result
  df_with_sg <- NULL
  if (!is.null(fs_full) && !is.null(fs_full$df.est)) {
    df_with_sg <- fs_full$df.est
  }

  # Create subgroup indicator
  if (!is.null(df_with_sg) && "treat.recommend" %in% names(df_with_sg)) {
    df$sg_hat <- as.integer(df_with_sg$treat.recommend == 0)
  } else {
    sg_factors <- character(0)
    for (i in seq_len(min(7, ncol(fs_res)))) {
      col_name <- paste0("M.", i)
      if (col_name %in% names(fs_res) && !is.na(fs_res[[col_name]][1])) {
        sg_factors <- c(sg_factors, fs_res[[col_name]][1])
      }
    }

    if (length(sg_factors) > 0) {
      df$sg_hat <- create_subgroup_indicator(df, sg_factors)
    } else {
      out$hr.Hc.hat <- out$hr.itt
      return(out)
    }
  }

  out$size.H <- sum(df$sg_hat, na.rm = TRUE)
  out$size.Hc <- sum(!df$sg_hat, na.rm = TRUE)

  # Compute HRs in identified subgroups
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

    if (out$size.H > 10) {
      out$hr.H.true <- tryCatch({
        exp(survival::coxph(
          survival::Surv(y.sim, event.sim) ~ treat,
          data = subset(df, sg_hat == 1)
        )$coefficients)
      }, error = function(e) dgm$hr_H_true)
    }
  }

  out
}


# =============================================================================
# Run ForestSearch Analysis - FIXED VERSION
# =============================================================================

#' Run ForestSearch Analysis
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
    confounders.name = confounders_name,
    details = verbose,
    plot.sg = verbose
  )

  # FIXED: Only valid forestsearch parameters - NO stop.threshold or pstop_futile!
  param_names <- c(
    "outcome.name", "event.name", "treat.name", "id.name",
    "use_lasso", "use_grf", "conf_force",
    "hr.threshold", "hr.consistency", "pconsistency.threshold",
    "fs.splits", "n.min", "d0.min", "d1.min",
    "maxk", "max.minutes", "by.risk", "vi.grf.min",
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
    warning(sprintf("%s analysis failed: %s", analysis_label, e$message))
    NULL
  })

  # Check for results
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
      analysis = analysis_label
    )
  } else {
    extract_fs_estimates(
      df = data,
      fs_res = NULL,
      fs_full = NULL,
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis = analysis_label
    )
  }
}


# =============================================================================
# Run GRF Analysis
# =============================================================================

#' Extract GRF Estimates
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
  out <- extract_fs_estimates(
    df = df,
    fs_res = NULL,
    fs_full = NULL,
    dgm = dgm,
    cox_formula = cox_formula,
    cox_formula_adj = cox_formula_adj,
    analysis = analysis
  )

  if (!is.null(grf_est) && !is.null(grf_est$harm_indicator)) {
    out$any.H <- 1L
    out$size.H <- sum(grf_est$harm_indicator, na.rm = TRUE)
    out$size.Hc <- sum(!grf_est$harm_indicator, na.rm = TRUE)

    if (!is.null(grf_est$hr_harm)) {
      out$hr.H.hat <- grf_est$hr_harm
    }
    if (!is.null(grf_est$hr_no_harm)) {
      out$hr.Hc.hat <- grf_est$hr_no_harm
    }
  }

  out
}


#' Run GRF Analysis
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

  grf_fun <- tryCatch({
    get("grf_subg_harm_survival", mode = "function", envir = parent.frame())
  }, error = function(e) {
    tryCatch({
      get("grf_subg_harm_survival", mode = "function", envir = globalenv())
    }, error = function(e2) NULL)
  })

  if (is.null(grf_fun)) {
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
      warning(sprintf("%s analysis failed: %s", analysis_label, e$message))
    }
    NULL
  })

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
# Main Simulation Analysis Function - FIXED VERSION
# =============================================================================

#' Run Single Simulation Analysis
#'
#' @param sim_id Integer. Simulation index
#' @param dgm A DGM object from create_gbsg_dgm()
#' @param n_sample Integer. Sample size
#' @param max_follow Numeric. Maximum follow-up time. Default: Inf
#' @param muC_adj Numeric. Censoring adjustment. Default: 0
#' @param confounders_base Character vector. Base confounder names
#' @param n_add_noise Integer. Number of noise variables to add. Default: 0
#' @param run_fs Logical. Run ForestSearch analysis. Default: TRUE
#' @param run_fs_grf Logical. Run ForestSearch with GRF. Default: TRUE
#' @param run_grf Logical. Run standalone GRF analysis. Default: TRUE
#' @param fs_params List. ForestSearch parameters (overrides defaults)
#' @param grf_params List. GRF parameters (overrides defaults)
#' @param cox_formula Formula. Cox model formula
#' @param cox_formula_adj Formula. Adjusted Cox formula
#' @param n_sims_total Integer. Total simulations (for progress)
#' @param seed_base Integer. Base random seed. Default: 8316951
#' @param verbose Logical. Print progress. Default: FALSE
#'
#' @return A data.table with analysis results
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

  # Simulate Data
  sim_data <- simulate_from_gbsg_dgm(
    dgm = dgm,
    n = n_sample,
    sim_id = sim_id,
    max_follow = max_follow,
    muC_adj = muC_adj
  )

  # Add Noise Variables if Requested
  confounders_name <- confounders_base

  if (n_add_noise > 0) {
    set.seed(seed_base + 1000L * sim_id)
    noise_names <- paste0("noise", seq_len(n_add_noise))

    for (nm in noise_names) {
      sim_data[[nm]] <- stats::rnorm(nrow(sim_data))
    }
    confounders_name <- c(confounders_base, noise_names)
  }

  # Compute True Subgroup Properties
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

  # Determine Verbosity
  show_details <- FALSE
  if (verbose && !is.null(n_sims_total)) {
    if (sim_id <= 1 || sim_id >= (n_sims_total - 1)) {
      show_details <- TRUE
    }
  }

  # Merge Parameters with Defaults - FIXED: no stop.threshold or pstop_futile!
  fs_defaults <- default_fs_params()
  fs_merged <- modifyList(fs_defaults, fs_params)

  grf_defaults <- default_grf_params()
  grf_merged <- modifyList(grf_defaults, grf_params)

  # Run Analyses
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

  # Combine Results
  if (length(results_list) == 0) {
    warning("No analyses were run. Check run_fs, run_fs_grf, run_grf settings.")
    return(NULL)
  }

  data.table::rbindlist(results_list, fill = TRUE)
}


# =============================================================================
# Summary Functions
# =============================================================================

#' Summarize Simulation Results
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
#' @keywords internal
summarize_single_analysis <- function(result, digits = 2, digits_hr = 3) {

  class_cols <- c("any.H", "ppv", "npv", "sensitivity", "specificity")
  class_means <- sapply(result[, class_cols, with = FALSE], mean, na.rm = TRUE)
  class_means <- round(class_means, digits)

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

  if (nrow(res_found) > 0) {
    hr_H_true <- round(mean(res_found$hr.H.true, na.rm = TRUE), digits_hr)
    hr_H_hat <- round(mean(res_found$hr.H.hat, na.rm = TRUE), digits_hr)
    hr_Hc_true <- round(mean(res_found$hr.Hc.true, na.rm = TRUE), digits_hr)
    hr_Hc_hat <- round(mean(res_found$hr.Hc.hat, na.rm = TRUE), digits_hr)
  } else {
    hr_H_true <- hr_H_hat <- hr_Hc_true <- hr_Hc_hat <- NA
  }

  hr_H_true_all <- round(mean(result$hr.H.true, na.rm = TRUE), digits_hr)
  hr_Hc_true_all <- round(mean(result$hr.Hc.true, na.rm = TRUE), digits_hr)
  hr_itt_all <- round(mean(result$hr.itt, na.rm = TRUE), digits_hr)
  hr_adj_itt_all <- round(mean(result$hr.adj.itt, na.rm = TRUE), digits_hr)

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
# Confirmation
# =============================================================================

cat("STANDALONE functions loaded successfully!\n")
cat("\nAvailable functions:\n")
cat("  - run_simulation_analysis()\n")
cat("  - run_forestsearch_analysis()\n")
cat("  - run_grf_analysis()\n")
cat("  - extract_fs_estimates()\n")
cat("  - extract_grf_estimates()\n")
cat("  - summarize_simulation_results()\n")
cat("  - default_fs_params()\n")
cat("  - default_grf_params()\n")

cat("\nVerifying NO bad parameters:\n")
defs <- default_fs_params()
if ("stop.threshold" %in% names(defs) || "pstop_futile" %in% names(defs)) {
  cat("  ERROR: Bad parameters found!\n")
} else {
  cat("  OK: No stop.threshold or pstop_futile\n")
}

cat("\nUsage:\n")
cat("  dgm <- create_gbsg_dgm(model = 'alt', k_inter = 2)\n")
cat("  result <- run_simulation_analysis(\n")
cat("    sim_id = 1,\n")
cat("    dgm = dgm,\n")
cat("    n_sample = 500,\n")
cat("    confounders_base = c('z1', 'z2', 'z3', 'z4', 'z5', 'size', 'grade3'),\n")
cat("    run_fs = TRUE,\n")
cat("    run_grf = FALSE,\n")
cat("    verbose = TRUE\n")
cat("  )\n")


# 3. Run your simulation
dgm <- create_gbsg_dgm(model = "alt", k_inter = 2)
result <- run_simulation_analysis(
  sim_id = 1,
  dgm = dgm,
  n_sample = 500,
  confounders_base = c("z1", "z2", "z3", "z4", "z5", "size", "grade3"),
  run_fs = TRUE,
  run_fs_grf = FALSE,
  run_grf = FALSE,
  verbose = TRUE
)
