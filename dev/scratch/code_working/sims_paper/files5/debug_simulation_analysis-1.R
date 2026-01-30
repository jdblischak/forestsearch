# =============================================================================
# Debug Script: run_simulation_analysis()
# =============================================================================
#
# Purpose: Step-by-step debugging of the simulation analysis workflow
# without parallel processing (dofuture).
#
# This script:
#   1. Creates synthetic trial data mimicking simulate_from_gbsg_dgm() output
#   2. Creates a mock DGM object with required elements
#   3. Tests each component function individually
#   4. Runs forestsearch() directly to verify it works
#   5. Tests the full run_simulation_analysis() pipeline
#
# =============================================================================

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

library(survival)
library(data.table)

# Load forestsearch package (adjust path as needed)
# devtools::load_all("path/to/forestsearch")

cat("=== Debug Script for run_simulation_analysis() ===\n\n")

# -----------------------------------------------------------------------------
# SECTION 0: Define Internal Helper Functions
# -----------------------------------------------------------------------------
# These functions are from oc_analyses_refactored.R and are needed for testing

cat("SECTION 0: Defining helper functions...\n")

#' Create Subgroup Indicator from Factor Definitions
#' @keywords internal
create_subgroup_indicator <- function(df, sg_factors) {
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

  # Try to use df.est from full result (has treat.recommend column)
  df_with_sg <- NULL
  if (!is.null(fs_full) && !is.null(fs_full$df.est)) {
    df_with_sg <- fs_full$df.est
  }

  # Create subgroup indicator
  if (!is.null(df_with_sg) && "treat.recommend" %in% names(df_with_sg)) {
    # Use treat.recommend: 0 = harm subgroup, 1 = complement
    df$sg_hat <- as.integer(df_with_sg$treat.recommend == 0)
  } else {
    # Fallback: try to parse M.1, M.2 factor definitions
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
      # Cannot determine subgroup
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

  # Compute classification metrics if true subgroup is known
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


#' Run ForestSearch Analysis Helper
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
  # NOTE: forestsearch() uses confounders.name, not Allconfounders.name
  fs_args <- list(
    df.analysis = data,
    confounders.name = confounders_name,
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
    warning(sprintf("%s analysis failed: %s", analysis_label, e$message))
    NULL
  })

  # Extract estimates
  # Check for results in grp.consistency$out_sg$result (correct path)
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

cat("  Helper functions defined:\n")
cat("    - create_subgroup_indicator()\n")
cat("    - extract_fs_estimates()\n")
cat("    - run_forestsearch_analysis()\n\n")

# -----------------------------------------------------------------------------
# SECTION 1: Create Synthetic Trial Data
# -----------------------------------------------------------------------------

cat("SECTION 1: Creating synthetic trial data...\n")

set.seed(12345)
n <- 500

# Create covariates that match expected column names
synthetic_data <- data.frame(
 id = 1:n,

 # Continuous covariates (standardized)
 z1 = rnorm(n),           # e.g., age
 z2 = rnorm(n),           # e.g., tumor size
 z3 = rnorm(n),           # e.g., nodes
 z4 = rnorm(n),           # e.g., progesterone receptor
 z5 = rnorm(n),           # e.g., estrogen receptor

 # Binary/categorical covariates
 size = rbinom(n, 1, 0.4),
 grade3 = rbinom(n, 1, 0.3),

 # Treatment assignment (1:1 randomization)
 treat = rbinom(n, 1, 0.5)
)

# Define true harm subgroup: z1 > 0 AND grade3 == 1
# This gives ~15% in harm subgroup
synthetic_data$flag.harm <- as.integer(synthetic_data$z1 > 0 & synthetic_data$grade3 == 1)

cat(sprintf("  Sample size: %d\n", n))
cat(sprintf("  True harm subgroup size: %d (%.1f%%)\n",
           sum(synthetic_data$flag.harm),
           100 * mean(synthetic_data$flag.harm)))

# Generate survival times with treatment effect heterogeneity
# Baseline hazard ~ Weibull(shape=1.5, scale=50)
# HR = 0.7 in complement (treatment benefit)
# HR = 1.5 in harm subgroup (treatment harm)

baseline_time <- rweibull(n, shape = 1.5, scale = 50)

# Apply treatment effect
hr_complement <- 0.7
hr_harm <- 1.5

treatment_effect <- ifelse(
 synthetic_data$flag.harm == 1,
 ifelse(synthetic_data$treat == 1, hr_harm, 1),
 ifelse(synthetic_data$treat == 1, hr_complement, 1)
)

# Transform times by treatment effect (AFT-style)
event_time <- baseline_time / treatment_effect^0.5

# Administrative censoring at 60 months
cens_time <- runif(n, 30, 80)
admin_cens <- 60

synthetic_data$y.sim <- pmin(event_time, cens_time, admin_cens)
synthetic_data$event.sim <- as.integer(event_time <= pmin(cens_time, admin_cens))
synthetic_data$t.sim <- event_time  # True (latent) event time

cat(sprintf("  Event rate: %.1f%%\n", 100 * mean(synthetic_data$event.sim)))
cat(sprintf("  Max follow-up: %.1f\n", max(synthetic_data$y.sim)))

# Verify treatment effects in subgroups
cat("\n  Verifying treatment effects:\n")

# ITT
fit_itt <- coxph(Surv(y.sim, event.sim) ~ treat, data = synthetic_data)
cat(sprintf("    ITT HR: %.3f\n", exp(coef(fit_itt))))

# Harm subgroup
fit_harm <- coxph(Surv(y.sim, event.sim) ~ treat,
                 data = subset(synthetic_data, flag.harm == 1))
cat(sprintf("    Harm subgroup HR: %.3f (target: %.1f)\n",
           exp(coef(fit_harm)), hr_harm))

# Complement
fit_comp <- coxph(Surv(y.sim, event.sim) ~ treat,
                 data = subset(synthetic_data, flag.harm == 0))
cat(sprintf("    Complement HR: %.3f (target: %.1f)\n",
           exp(coef(fit_comp)), hr_complement))

cat("\n  Synthetic data created successfully.\n")
cat("  Column names:", paste(names(synthetic_data), collapse = ", "), "\n\n")


# -----------------------------------------------------------------------------
# SECTION 2: Create Mock DGM Object
# -----------------------------------------------------------------------------

cat("SECTION 2: Creating mock DGM object...\n")

# The DGM object needs specific elements for extract_fs_estimates()
mock_dgm <- list(
 model = "alt",
 hr_H_true = hr_harm,
 hr_Hc_true = hr_complement,
 hr_itt_true = NA,  # Will be computed

 # Additional metadata
 subgroup_definition = "z1 > 0 & grade3 == 1",
 prop_harm = mean(synthetic_data$flag.harm),

 # Placeholders for compatibility
 model_params = list(mu = 3.5, sigma = 0.8),
 cens_params = list(type = "weibull", mu = 4.0, sigma = 0.5)
)

class(mock_dgm) <- c("gbsg_dgm", "list")

cat("  DGM object created with:\n")
cat(sprintf("    hr_H_true: %.2f\n", mock_dgm$hr_H_true))
cat(sprintf("    hr_Hc_true: %.2f\n", mock_dgm$hr_Hc_true))
cat(sprintf("    prop_harm: %.3f\n", mock_dgm$prop_harm))
cat("\n")


# -----------------------------------------------------------------------------
# SECTION 3: Test forestsearch() Directly
# -----------------------------------------------------------------------------

cat("SECTION 3: Testing forestsearch() directly...\n")

# Define confounders
confounders_base <- c("z1", "z2", "z3", "z4", "z5", "size", "grade3")

cat("  Confounders:", paste(confounders_base, collapse = ", "), "\n")

# Check that confounders exist and are appropriate types
cat("  Checking confounder types:\n")
for (conf in confounders_base) {
 if (conf %in% names(synthetic_data)) {
   val <- synthetic_data[[conf]]
   cat(sprintf("    %s: %s (range: %.2f to %.2f)\n",
               conf, class(val), min(val), max(val)))
 } else {
   cat(sprintf("    %s: MISSING!\n", conf))
 }
}

# Test forestsearch with verbose output
cat("\n  Calling forestsearch()...\n")
cat("  " , rep("-", 50), "\n", sep = "")

fs_result <- tryCatch({
 forestsearch(
   df.analysis = synthetic_data,
   confounders.name = confounders_base,  # CORRECT parameter name
   outcome.name = "y.sim",
   event.name = "event.sim",
   treat.name = "treat",
   id.name = "id",
   use_lasso = TRUE,
   use_grf = FALSE,
   hr.threshold = 1.25,
   hr.consistency = 1.0,
   pconsistency.threshold = 0.90,
   fs.splits = 100,  # Reduced for debugging
   n.min = 30,
   d0.min = 10,
   d1.min = 10,
   maxk = 2,
   max.minutes = 2,
   details = TRUE,
   plot.sg = FALSE
 )
}, error = function(e) {
 cat("  ERROR in forestsearch():\n")
 cat("  ", e$message, "\n")
 return(NULL)
})

cat("  ", rep("-", 50), "\n", sep = "")

# Examine result structure
if (!is.null(fs_result)) {
 cat("\n  ForestSearch result structure:\n")
 cat("    Class:", paste(class(fs_result), collapse = ", "), "\n")
 cat("    Top-level elements:", paste(names(fs_result), collapse = ", "), "\n")

 # Check for subgroup
 if (!is.null(fs_result$sg.harm)) {
   cat("\n  Subgroup found!\n")
   cat("    sg.harm:", paste(fs_result$sg.harm, collapse = " & "), "\n")
 } else {
   cat("\n  No subgroup found (sg.harm is NULL)\n")
 }

 # Check grp.consistency structure
 if (!is.null(fs_result$grp.consistency)) {
   gc <- fs_result$grp.consistency
   cat("\n  grp.consistency structure:\n")
   cat("    Elements:", paste(names(gc), collapse = ", "), "\n")

   # Check out_sg
   if (!is.null(gc$out_sg)) {
     cat("    out_sg elements:", paste(names(gc$out_sg), collapse = ", "), "\n")

     if (!is.null(gc$out_sg$result)) {
       cat("    out_sg$result dimensions:",
           nrow(gc$out_sg$result), "x", ncol(gc$out_sg$result), "\n")
       cat("    out_sg$result columns:",
           paste(names(gc$out_sg$result), collapse = ", "), "\n")

       # Show first row
       if (nrow(gc$out_sg$result) > 0) {
         cat("\n    First result row:\n")
         print(gc$out_sg$result[1, ])
       }
     }
   } else {
     cat("    out_sg is NULL\n")
   }
 } else {
   cat("\n  grp.consistency is NULL\n")
 }

 # Check df.est
 if (!is.null(fs_result$df.est)) {
   cat("\n  df.est structure:\n")
   cat("    Dimensions:", nrow(fs_result$df.est), "x", ncol(fs_result$df.est), "\n")
   cat("    Columns:", paste(names(fs_result$df.est), collapse = ", "), "\n")

   if ("treat.recommend" %in% names(fs_result$df.est)) {
     cat("    treat.recommend distribution:\n")
     print(table(fs_result$df.est$treat.recommend, useNA = "ifany"))
   }
 } else {
   cat("\n  df.est is NULL\n")
 }

} else {
 cat("\n  forestsearch() returned NULL - check error above\n")
}

cat("\n")


# -----------------------------------------------------------------------------
# SECTION 4: Test extract_fs_estimates() Function
# -----------------------------------------------------------------------------

cat("SECTION 4: Testing extract_fs_estimates()...\n")

# Test with NULL result (no subgroup found)
cat("  Testing with NULL fs_res (no subgroup):\n")
est_null <- extract_fs_estimates(
  df = synthetic_data,
  fs_res = NULL,
  fs_full = NULL,
  dgm = mock_dgm,
  cox_formula = NULL,
  cox_formula_adj = NULL,
  analysis = "TEST_NULL"
)
print(est_null)

# Test with actual result if available
if (!is.null(fs_result) &&
    !is.null(fs_result$grp.consistency) &&
    !is.null(fs_result$grp.consistency$out_sg) &&
    !is.null(fs_result$grp.consistency$out_sg$result)) {

  cat("\n  Testing with actual forestsearch result:\n")
  est_real <- extract_fs_estimates(
    df = synthetic_data,
    fs_res = fs_result$grp.consistency$out_sg$result,
    fs_full = fs_result,
    dgm = mock_dgm,
    cox_formula = NULL,
    cox_formula_adj = NULL,
    analysis = "TEST_REAL"
  )
  print(est_real)

  # Check classification metrics
  cat("\n  Classification metrics interpretation:\n")
  cat(sprintf("    any.H: %d (subgroup detected)\n", est_real$any.H))
  cat(sprintf("    size.H: %d (identified harm subgroup size)\n", est_real$size.H))
  cat(sprintf("    size.Hc: %d (identified complement size)\n", est_real$size.Hc))
  cat(sprintf("    sensitivity: %.3f (true harm correctly identified)\n",
              est_real$sensitivity))
  cat(sprintf("    specificity: %.3f (true complement correctly identified)\n",
              est_real$specificity))
  cat(sprintf("    ppv: %.3f (positive predictive value)\n", est_real$ppv))
  cat(sprintf("    npv: %.3f (negative predictive value)\n", est_real$npv))
} else {
  cat("\n  No forestsearch result available to test extract_fs_estimates() with real data\n")
}

cat("\n")


# -----------------------------------------------------------------------------
# SECTION 5: Test run_forestsearch_analysis() Helper
# -----------------------------------------------------------------------------

cat("SECTION 5: Testing run_forestsearch_analysis() helper...\n")

# Define parameters
fs_params <- list(
  outcome.name = "y.sim",
  event.name = "event.sim",
  treat.name = "treat",
  id.name = "id",
  use_lasso = TRUE,
  use_grf = FALSE,
  hr.threshold = 1.25,
  hr.consistency = 1.0,
  pconsistency.threshold = 0.90,
  fs.splits = 100,
  n.min = 30,
  d0.min = 10,
  d1.min = 10,
  maxk = 2,
  max.minutes = 2
)

cat("  Calling run_forestsearch_analysis()...\n")

fs_helper_result <- tryCatch({
  run_forestsearch_analysis(
    data = synthetic_data,
    confounders_name = confounders_base,
    params = fs_params,
    dgm = mock_dgm,
    cox_formula = NULL,
    cox_formula_adj = NULL,
    analysis_label = "FS_TEST",
    verbose = TRUE
  )
}, error = function(e) {
  cat("  ERROR:\n")
  cat("  ", e$message, "\n")
  traceback()
  return(NULL)
})

if (!is.null(fs_helper_result)) {
  cat("\n  run_forestsearch_analysis() result:\n")
  print(fs_helper_result)
}

cat("\n")


# -----------------------------------------------------------------------------
# SECTION 6: Test Full run_simulation_analysis() Pipeline
# -----------------------------------------------------------------------------

cat("SECTION 6: Testing full run_simulation_analysis() pipeline...\n")

cat("  Note: run_simulation_analysis() calls simulate_from_gbsg_dgm() internally.\n")
cat("  For this test, we need a proper DGM object.\n\n")

# Check if the full function exists (from oc_analyses_refactored.R)
if (exists("run_simulation_analysis")) {
  cat("  run_simulation_analysis() is available.\n")
  cat("  To test the full pipeline:\n")
  cat("    1. Create a real DGM with create_gbsg_dgm()\n")
  cat("    2. Call run_simulation_analysis() with that DGM\n")
} else {
  cat("  run_simulation_analysis() not found.\n")
  cat("  This function is in oc_analyses_refactored.R and requires simulate_from_gbsg_dgm().\n")
  cat("  The helper functions tested above are the core components.\n")
}

cat("\n")


# -----------------------------------------------------------------------------
# SECTION 7: Manual Pipeline Test (Step by Step)
# -----------------------------------------------------------------------------

cat("SECTION 7: Manual pipeline walkthrough...\n\n")

cat("Step 1: Data simulation (using synthetic_data created above)\n")
cat("  ✓ synthetic_data has", nrow(synthetic_data), "rows\n")
cat("  ✓ Required columns present:",
   all(c("y.sim", "event.sim", "treat", "id", "flag.harm") %in% names(synthetic_data)), "\n")

cat("\nStep 2: Compute true subgroup properties\n")
size_H_true <- sum(synthetic_data$flag.harm)
prop_H_true <- mean(synthetic_data$flag.harm)
cat(sprintf("  True H size: %d (%.1f%%)\n", size_H_true, 100 * prop_H_true))

cat("\nStep 3: Create population data.table\n")
df_pop <- data.table(
 sim = 1,
 sizeH_true = size_H_true,
 propH_true = prop_H_true,
 sizeHc_true = sum(!synthetic_data$flag.harm),
 propHc_true = mean(!synthetic_data$flag.harm)
)
print(df_pop)

cat("\nStep 4: Build forestsearch arguments\n")
fs_args <- list(
 df.analysis = synthetic_data,
 confounders.name = confounders_base,  # CRITICAL: correct parameter name
 outcome.name = "y.sim",
 event.name = "event.sim",
 treat.name = "treat",
 id.name = "id",
 use_lasso = TRUE,
 use_grf = FALSE,
 hr.threshold = 1.25,
 hr.consistency = 1.0,
 pconsistency.threshold = 0.90,
 fs.splits = 100,
 n.min = 30,
 maxk = 2,
 details = FALSE,
 plot.sg = FALSE
)
cat("  Arguments prepared:\n")
cat("    df.analysis: data.frame with", nrow(fs_args$df.analysis), "rows\n")
cat("    confounders.name:", paste(fs_args$confounders.name, collapse = ", "), "\n")
cat("    outcome.name:", fs_args$outcome.name, "\n")
cat("    event.name:", fs_args$event.name, "\n")

cat("\nStep 5: Call forestsearch via do.call()\n")
fs_result_manual <- tryCatch({
 do.call(forestsearch, fs_args)
}, error = function(e) {
 cat("  ERROR:", e$message, "\n")
 NULL
})

if (!is.null(fs_result_manual)) {
 cat("  ✓ forestsearch() completed successfully\n")

 # Check result path
 has_result <- !is.null(fs_result_manual$grp.consistency) &&
               !is.null(fs_result_manual$grp.consistency$out_sg) &&
               !is.null(fs_result_manual$grp.consistency$out_sg$result) &&
               nrow(fs_result_manual$grp.consistency$out_sg$result) > 0

 cat("  Result path exists:", has_result, "\n")

 if (has_result) {
   cat("  Subgroup found:", !is.null(fs_result_manual$sg.harm), "\n")
   if (!is.null(fs_result_manual$sg.harm)) {
     cat("  sg.harm:", paste(fs_result_manual$sg.harm, collapse = " & "), "\n")
   }
 }
} else {
 cat("  ✗ forestsearch() failed\n")
}

cat("\n")


# -----------------------------------------------------------------------------
# SECTION 8: Summary and Recommendations
# -----------------------------------------------------------------------------

cat("=== SUMMARY ===\n\n")

cat("Key findings:\n")
cat("1. Parameter name: forestsearch() uses 'confounders.name', NOT 'Allconfounders.name'\n")
cat("2. Result path: fs_result$grp.consistency$out_sg$result (not $out_hr)\n")
cat("3. Subgroup indicator: fs_result$df.est$treat.recommend (0=harm, 1=complement)\n")
cat("4. DGM needs: hr_H_true, hr_Hc_true elements for extract_fs_estimates()\n")

cat("\nDebugging checklist:\n")
cat("[ ] Verify forestsearch is loaded: exists('forestsearch')\n")
cat("[ ] Check confounders exist in data: all(confounders %in% names(data))\n")
cat("[ ] Verify outcome/event columns: 'y.sim', 'event.sim' in names(data)\n")
cat("[ ] Check for sufficient events: sum(data$event.sim) > 50\n")
cat("[ ] Verify DGM structure: has hr_H_true, hr_Hc_true\n")

cat("\nTo run this debug script:\n")
cat("  1. Load forestsearch package: devtools::load_all('path/to/forestsearch')\n")
cat("  2. Source this script: source('debug_simulation_analysis.R')\n")
cat("  3. Check output for errors at each section\n")

cat("\n=== END DEBUG SCRIPT ===\n")


# Create DGM
# dgm <- create_gbsg_dgm(model = "alt", k_inter = 2)
# # Run single simulation
# result <- run_simulation_analysis(
#   sim_id = 1,
#   dgm = dgm,
#   n_sample = 500,
#   confounders_base = c("z1", "z2", "z3", "z4", "z5", "size", "grade3"),
#   run_fs = TRUE,
#   run_grf = TRUE,
#   verbose = TRUE
# )

