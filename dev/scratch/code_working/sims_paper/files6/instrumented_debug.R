# =============================================================================
# INSTRUMENTED Debug Version - Traces Every Step
# =============================================================================
#
# This version adds print statements at EVERY step to identify exactly
# where run_simulation_analysis() stops executing.
#
# =============================================================================

library(survival)
library(data.table)

cat("=== INSTRUMENTED DEBUG: run_simulation_analysis() ===\n\n")

# =============================================================================
# INSTRUMENTED: run_forestsearch_analysis
# =============================================================================

run_forestsearch_analysis_DEBUG <- function(
    data,
    confounders_name,
    params,
    dgm,
    cox_formula = NULL,
    cox_formula_adj = NULL,
    analysis_label = "FS",
    verbose = FALSE
) {

  cat("    [DEBUG] Entered run_forestsearch_analysis_DEBUG()\n")
  cat("    [DEBUG] analysis_label:", analysis_label, "\n")
  cat("    [DEBUG] data dimensions:", nrow(data), "x", ncol(data), "\n")
  cat("    [DEBUG] confounders_name:", paste(confounders_name, collapse = ", "), "\n")
  cat("    [DEBUG] verbose:", verbose, "\n")

  # Check data columns
  cat("    [DEBUG] Checking required columns in data...\n")
  required <- c("y.sim", "event.sim", "treat", "id")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    cat("    [DEBUG] ERROR: Missing columns:", paste(missing, collapse = ", "), "\n")
    return(NULL)
  }
  cat("    [DEBUG] All required columns present\n")

  # Check confounders
  cat("    [DEBUG] Checking confounders in data...\n")
  missing_conf <- setdiff(confounders_name, names(data))
  if (length(missing_conf) > 0) {
    cat("    [DEBUG] ERROR: Missing confounders:", paste(missing_conf, collapse = ", "), "\n")
    return(NULL)
  }
  cat("    [DEBUG] All confounders present\n")

  # Build arguments
  cat("    [DEBUG] Building fs_args list...\n")
  fs_args <- list(
    df.analysis = data,
    confounders.name = confounders_name,
    details = verbose,
    plot.sg = verbose
  )
  cat("    [DEBUG] Base fs_args created\n")

  # Add parameters
  param_names <- c(
    "outcome.name", "event.name", "treat.name", "id.name",
    "use_lasso", "use_grf", "conf_force",
    "hr.threshold", "hr.consistency", "pconsistency.threshold",
    "fs.splits", "n.min", "d0.min", "d1.min",
    "maxk", "max.minutes", "by.risk", "vi.grf.min",
    "frac.tau", "dmin.grf", "grf_depth"
  )

  params_added <- character(0)
  for (pn in param_names) {
    if (!is.null(params[[pn]])) {
      fs_args[[pn]] <- params[[pn]]
      params_added <- c(params_added, pn)
    }
  }
  cat("    [DEBUG] Parameters added:", paste(params_added, collapse = ", "), "\n")

  # Show key parameters
  cat("    [DEBUG] Key fs_args values:\n")
  cat("      outcome.name:", fs_args$outcome.name, "\n")
  cat("      event.name:", fs_args$event.name, "\n")
  cat("      treat.name:", fs_args$treat.name, "\n")
  cat("      id.name:", fs_args$id.name, "\n")
  cat("      use_lasso:", fs_args$use_lasso, "\n")
  cat("      use_grf:", fs_args$use_grf, "\n")
  cat("      n.min:", fs_args$n.min, "\n")
  cat("      fs.splits:", fs_args$fs.splits, "\n")

  # Check if forestsearch exists
  cat("    [DEBUG] Checking if forestsearch() exists...\n")
  if (!exists("forestsearch", mode = "function")) {
    cat("    [DEBUG] ERROR: forestsearch() function NOT FOUND!\n")
    cat("    [DEBUG] Available functions with 'forest' in name:\n")
    all_fns <- ls(envir = .GlobalEnv)
    forest_fns <- grep("forest", all_fns, value = TRUE, ignore.case = TRUE)
    cat("      ", paste(forest_fns, collapse = ", "), "\n")
    return(NULL)
  }
  cat("    [DEBUG] forestsearch() exists\n")

  # Run ForestSearch
  cat("    [DEBUG] Calling do.call(forestsearch, fs_args)...\n")
  cat("    [DEBUG] ========== FORESTSEARCH START ==========\n")

  fs_result <- tryCatch({
    result <- do.call(forestsearch, fs_args)
    cat("    [DEBUG] ========== FORESTSEARCH END ==========\n")
    cat("    [DEBUG] forestsearch() returned successfully\n")
    result
  }, error = function(e) {
    cat("    [DEBUG] ========== FORESTSEARCH ERROR ==========\n")
    cat("    [DEBUG] ERROR in forestsearch():\n")
    cat("    [DEBUG]   Message:", e$message, "\n")
    cat("    [DEBUG]   Call:", deparse(e$call), "\n")
    NULL
  }, warning = function(w) {
    cat("    [DEBUG] WARNING in forestsearch():", w$message, "\n")
    invokeRestart("muffleWarning")
  })

  # Check result
  cat("    [DEBUG] Checking fs_result...\n")
  if (is.null(fs_result)) {
    cat("    [DEBUG] fs_result is NULL\n")
  } else {
    cat("    [DEBUG] fs_result class:", paste(class(fs_result), collapse = ", "), "\n")
    cat("    [DEBUG] fs_result names:", paste(names(fs_result), collapse = ", "), "\n")

    if (!is.null(fs_result$sg.harm)) {
      cat("    [DEBUG] sg.harm:", paste(fs_result$sg.harm, collapse = " & "), "\n")
    } else {
      cat("    [DEBUG] sg.harm is NULL (no subgroup found)\n")
    }
  }

  # Check result path
  cat("    [DEBUG] Checking result extraction path...\n")
  has_result <- !is.null(fs_result) &&
                !is.null(fs_result$grp.consistency) &&
                !is.null(fs_result$grp.consistency$out_sg) &&
                !is.null(fs_result$grp.consistency$out_sg$result) &&
                nrow(fs_result$grp.consistency$out_sg$result) > 0

  cat("    [DEBUG] has_result:", has_result, "\n")

  if (!is.null(fs_result)) {
    cat("    [DEBUG] grp.consistency exists:", !is.null(fs_result$grp.consistency), "\n")
    if (!is.null(fs_result$grp.consistency)) {
      cat("    [DEBUG] out_sg exists:", !is.null(fs_result$grp.consistency$out_sg), "\n")
      if (!is.null(fs_result$grp.consistency$out_sg)) {
        cat("    [DEBUG] out_sg names:", paste(names(fs_result$grp.consistency$out_sg), collapse = ", "), "\n")
        cat("    [DEBUG] result exists:", !is.null(fs_result$grp.consistency$out_sg$result), "\n")
        if (!is.null(fs_result$grp.consistency$out_sg$result)) {
          cat("    [DEBUG] result rows:", nrow(fs_result$grp.consistency$out_sg$result), "\n")
        }
      }
    }
  }

  # Return the raw result for inspection
  cat("    [DEBUG] Returning fs_result for inspection\n")
  return(fs_result)
}


# =============================================================================
# INSTRUMENTED: run_simulation_analysis
# =============================================================================

run_simulation_analysis_DEBUG <- function(
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

  cat("\n[DEBUG] ============================================\n")
  cat("[DEBUG] run_simulation_analysis_DEBUG() STARTED\n")
  cat("[DEBUG] ============================================\n")
  cat("[DEBUG] sim_id:", sim_id, "\n")
  cat("[DEBUG] n_sample:", n_sample, "\n")
  cat("[DEBUG] run_fs:", run_fs, "\n")
  cat("[DEBUG] run_fs_grf:", run_fs_grf, "\n")
  cat("[DEBUG] run_grf:", run_grf, "\n")
  cat("[DEBUG] verbose:", verbose, "\n")

  # Check DGM
  cat("\n[DEBUG] Checking DGM object...\n")
  cat("[DEBUG] dgm class:", paste(class(dgm), collapse = ", "), "\n")
  cat("[DEBUG] dgm names:", paste(names(dgm), collapse = ", "), "\n")
  cat("[DEBUG] dgm$model:", dgm$model, "\n")
  cat("[DEBUG] dgm$hr_H_true:", dgm$hr_H_true, "\n")
  cat("[DEBUG] dgm$hr_Hc_true:", dgm$hr_Hc_true, "\n")

  # -------------------------------------------------------------------------
  # Simulate Data
  # -------------------------------------------------------------------------
  cat("\n[DEBUG] Calling simulate_from_gbsg_dgm()...\n")

  sim_data <- tryCatch({
    result <- simulate_from_gbsg_dgm(
      dgm = dgm,
      n = n_sample,
      sim_id = sim_id,
      max_follow = max_follow,
      muC_adj = muC_adj
    )
    cat("[DEBUG] simulate_from_gbsg_dgm() SUCCESS\n")
    result
  }, error = function(e) {
    cat("[DEBUG] ERROR in simulate_from_gbsg_dgm():\n")
    cat("[DEBUG]   Message:", e$message, "\n")
    NULL
  })

  if (is.null(sim_data)) {
    cat("[DEBUG] sim_data is NULL - ABORTING\n")
    return(NULL)
  }

  cat("[DEBUG] sim_data dimensions:", nrow(sim_data), "x", ncol(sim_data), "\n")
  cat("[DEBUG] sim_data columns:", paste(names(sim_data), collapse = ", "), "\n")
  cat("[DEBUG] Events:", sum(sim_data$event.sim), "/", nrow(sim_data), "\n")
  cat("[DEBUG] Harm subgroup:", sum(sim_data$flag.harm), "/", nrow(sim_data), "\n")

  # -------------------------------------------------------------------------
  # Add Noise Variables if Requested
  # -------------------------------------------------------------------------
  confounders_name <- confounders_base

  if (n_add_noise > 0) {
    cat("[DEBUG] Adding", n_add_noise, "noise variables\n")
    set.seed(seed_base + 1000L * sim_id)
    noise_names <- paste0("noise", seq_len(n_add_noise))

    for (nm in noise_names) {
      sim_data[[nm]] <- stats::rnorm(nrow(sim_data))
    }
    confounders_name <- c(confounders_base, noise_names)
  }

  cat("[DEBUG] confounders_name:", paste(confounders_name, collapse = ", "), "\n")

  # -------------------------------------------------------------------------
  # Compute True Subgroup Properties
  # -------------------------------------------------------------------------
  cat("\n[DEBUG] Computing true subgroup properties...\n")
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
  cat("[DEBUG] df_pop created\n")
  print(df_pop)

  # -------------------------------------------------------------------------
  # Merge Parameters with Defaults
  # -------------------------------------------------------------------------
  cat("\n[DEBUG] Merging parameters with defaults...\n")

  fs_defaults <- list(
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
    n.min = 60,
    d0.min = 20,
    d1.min = 20,
    maxk = 2,
    max.minutes = 2,  # Reduced for debugging
    by.risk = 12,
    vi.grf.min = NULL
  )

  fs_merged <- modifyList(fs_defaults, fs_params)
  cat("[DEBUG] fs_merged created\n")

  # -------------------------------------------------------------------------
  # Run Analyses
  # -------------------------------------------------------------------------
  results_list <- list()

  # ForestSearch (LASSO only)
  if (run_fs) {
    cat("\n[DEBUG] ============================================\n")
    cat("[DEBUG] RUNNING ForestSearch (LASSO only)\n")
    cat("[DEBUG] ============================================\n")

    fs_params_run <- modifyList(fs_merged, list(use_lasso = TRUE, use_grf = FALSE))

    fs_result <- run_forestsearch_analysis_DEBUG(
      data = sim_data,
      confounders_name = confounders_name,
      params = fs_params_run,
      dgm = dgm,
      cox_formula = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis_label = "FS",
      verbose = TRUE  # Force verbose for debugging
    )

    cat("\n[DEBUG] fs_result returned\n")
    if (!is.null(fs_result)) {
      cat("[DEBUG] fs_result class:", paste(class(fs_result), collapse = ", "), "\n")
    } else {
      cat("[DEBUG] fs_result is NULL\n")
    }
  }

  cat("\n[DEBUG] ============================================\n")
  cat("[DEBUG] run_simulation_analysis_DEBUG() COMPLETED\n")
  cat("[DEBUG] ============================================\n")

  return(fs_result)
}


# =============================================================================
# RUN THE DEBUG TEST
# =============================================================================

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("RUNNING DEBUG TEST\n")
cat("=" , rep("=", 60), "\n", sep = "")

cat("\nStep 1: Create DGM\n")
dgm <- tryCatch({
  create_gbsg_dgm(model = "alt", k_inter = 2)
}, error = function(e) {
  cat("ERROR creating DGM:", e$message, "\n")
  NULL
})

if (is.null(dgm)) {
  stop("Cannot continue without DGM")
}

cat("DGM created successfully\n")

cat("\nStep 2: Run instrumented simulation analysis\n")
result <- run_simulation_analysis_DEBUG(
  sim_id = 1,
  dgm = dgm,
  n_sample = 500,
  confounders_base = c("z1", "z2", "z3", "z4", "z5", "size", "grade3"),
  run_fs = TRUE,
  run_fs_grf = FALSE,  # Disable to simplify
  run_grf = FALSE,     # Disable to simplify
  verbose = TRUE
)

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("DEBUG TEST COMPLETED\n")
cat("=" , rep("=", 60), "\n", sep = "")

if (!is.null(result)) {
  cat("\nResult returned:\n")
  cat("Class:", paste(class(result), collapse = ", "), "\n")
  if ("sg.harm" %in% names(result)) {
    cat("sg.harm:", ifelse(is.null(result$sg.harm), "NULL",
                          paste(result$sg.harm, collapse = " & ")), "\n")
  }
} else {
  cat("\nResult is NULL - check debug output above for where it failed\n")
}

library(forestsearch)

# Now it should work
dgm <- create_gbsg_dgm(model = "alt", k_inter = 2)

result <- run_simulation_analysis_DEBUG(
  sim_id = 1,
  dgm = dgm,
  n_sample = 500,
  confounders_base = c("z1", "z2", "z3", "z4", "z5", "size", "grade3"),
  run_fs = TRUE,
  run_grf = FALSE,
  verbose = TRUE
)

