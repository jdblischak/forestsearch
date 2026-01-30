# =============================================================================
# Diagnostic Script: run_simulation_analysis() with Real DGM
# =============================================================================
#
# This script diagnoses why run_simulation_analysis() fails to execute
# forestsearch() when using a real DGM from create_gbsg_dgm()
#
# =============================================================================

library(survival)
library(data.table)

library(forestsearch)

# Load forestsearch package
# devtools::load_all("path/to/forestsearch")

cat("=== Diagnostic: run_simulation_analysis() with Real DGM ===\n\n")

# =============================================================================
# STEP 1: Create DGM and Examine Structure
# =============================================================================

cat("STEP 1: Creating DGM...\n")

dgm <- create_gbsg_dgm(model = "alt", k_inter = 2)

cat("  DGM class:", paste(class(dgm), collapse = ", "), "\n")
cat("  DGM elements:", paste(names(dgm), collapse = ", "), "\n")

# Check for required elements
cat("\n  Checking required DGM elements:\n")
cat("    hr_H_true:", ifelse(is.null(dgm$hr_H_true), "MISSING", dgm$hr_H_true), "\n")
cat("    hr_Hc_true:", ifelse(is.null(dgm$hr_Hc_true), "MISSING", dgm$hr_Hc_true), "\n")
cat("    model:", ifelse(is.null(dgm$model), "MISSING", dgm$model), "\n")

# =============================================================================
# STEP 2: Test simulate_from_gbsg_dgm() Directly
# =============================================================================

cat("\nSTEP 2: Testing simulate_from_gbsg_dgm()...\n")

sim_data <- tryCatch({
  simulate_from_gbsg_dgm(
    dgm = dgm,
    n = 500,
    sim_id = 1,
    max_follow = Inf,
    muC_adj = 0
  )
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  NULL
})

if (!is.null(sim_data)) {
  cat("  Simulated data dimensions:", nrow(sim_data), "x", ncol(sim_data), "\n")
  cat("  Column names:\n    ", paste(names(sim_data), collapse = ", "), "\n")

  # Check for required columns
  required_cols <- c("id", "y.sim", "event.sim", "treat", "flag.harm",
                     "z1", "z2", "z3", "z4", "z5", "size", "grade3")
  missing_cols <- setdiff(required_cols, names(sim_data))

  if (length(missing_cols) > 0) {
    cat("\n  WARNING: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
  } else {
    cat("\n  All required columns present.\n")
  }

  # Show data summary
  cat("\n  Data summary:\n")
  cat("    Events:", sum(sim_data$event.sim), "/", nrow(sim_data),
      sprintf("(%.1f%%)\n", 100 * mean(sim_data$event.sim)))
  cat("    Treated:", sum(sim_data$treat), "/", nrow(sim_data), "\n")
  cat("    Harm subgroup:", sum(sim_data$flag.harm), "/", nrow(sim_data),
      sprintf("(%.1f%%)\n", 100 * mean(sim_data$flag.harm)))
  cat("    Max follow-up:", max(sim_data$y.sim), "\n")

  # Check confounder types
  cat("\n  Confounder types:\n")
  confounders_base <- c("z1", "z2", "z3", "z4", "z5", "size", "grade3")
  for (conf in confounders_base) {
    if (conf %in% names(sim_data)) {
      val <- sim_data[[conf]]
      cat(sprintf("    %s: %s, range [%.2f, %.2f]\n",
                  conf, class(val)[1], min(val, na.rm = TRUE), max(val, na.rm = TRUE)))
    }
  }
} else {
  cat("  simulate_from_gbsg_dgm() failed - cannot continue\n")
  stop("Simulation failed")
}

# =============================================================================
# STEP 3: Test forestsearch() Directly with Simulated Data
# =============================================================================

cat("\nSTEP 3: Testing forestsearch() directly with simulated data...\n")

confounders_base <- c("z1", "z2", "z3", "z4", "z5", "size", "grade3")

# Verify confounders exist
missing_conf <- setdiff(confounders_base, names(sim_data))
if (length(missing_conf) > 0) {
  cat("  ERROR: Missing confounders:", paste(missing_conf, collapse = ", "), "\n")
  cat("  Available columns:", paste(names(sim_data), collapse = ", "), "\n")
  stop("Missing confounders")
}

cat("  Calling forestsearch() with:\n")
cat("    confounders.name:", paste(confounders_base, collapse = ", "), "\n")
cat("    outcome.name: y.sim\n")
cat("    event.name: event.sim\n")
cat("    treat.name: treat\n")
cat("    id.name: id\n")

fs_direct <- tryCatch({
  forestsearch(
    df.analysis = sim_data,
    confounders.name = confounders_base,  # CORRECT parameter
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
    n.min = 60,
    d0.min = 20,
    d1.min = 20,
    maxk = 2,
    max.minutes = 2,
    details = TRUE,
    plot.sg = FALSE
  )
}, error = function(e) {
  cat("  ERROR in forestsearch():", e$message, "\n")
  NULL
})

if (!is.null(fs_direct)) {
  cat("\n  forestsearch() completed successfully!\n")
  cat("  sg.harm:", ifelse(is.null(fs_direct$sg.harm), "NULL (no subgroup)",
                          paste(fs_direct$sg.harm, collapse = " & ")), "\n")

  # Check result structure
  if (!is.null(fs_direct$grp.consistency)) {
    cat("  grp.consistency elements:", paste(names(fs_direct$grp.consistency), collapse = ", "), "\n")

    if (!is.null(fs_direct$grp.consistency$out_sg)) {
      cat("  out_sg elements:", paste(names(fs_direct$grp.consistency$out_sg), collapse = ", "), "\n")

      if (!is.null(fs_direct$grp.consistency$out_sg$result)) {
        cat("  out_sg$result rows:", nrow(fs_direct$grp.consistency$out_sg$result), "\n")
      }
    }
  }
}

# =============================================================================
# STEP 4: Check run_simulation_analysis() Function Source
# =============================================================================

cat("\nSTEP 4: Checking run_simulation_analysis() function...\n")

if (exists("run_simulation_analysis")) {
  # Get the function body and check for the bug
  fn_body <- deparse(body(run_simulation_analysis))
  fn_text <- paste(fn_body, collapse = "\n")

  # Check for the bug: Allconfounders.name vs confounders.name
  if (grepl("Allconfounders.name", fn_text)) {
    cat("  BUG FOUND: Function uses 'Allconfounders.name' (WRONG)\n")
    cat("  This is the OLD buggy version!\n")
  } else if (grepl("confounders.name", fn_text)) {
    cat("  Function uses 'confounders.name' (CORRECT)\n")
  } else {
    cat("  Could not determine parameter usage\n")
  }

  # Check environment
  cat("  Function environment:", environmentName(environment(run_simulation_analysis)), "\n")

} else {
  cat("  run_simulation_analysis() not found!\n")
}

# =============================================================================
# STEP 5: Check run_forestsearch_analysis() Helper
# =============================================================================

cat("\nSTEP 5: Checking run_forestsearch_analysis() helper...\n")

if (exists("run_forestsearch_analysis")) {
  fn_body <- deparse(body(run_forestsearch_analysis))
  fn_text <- paste(fn_body, collapse = "\n")

  # Check for bugs
  if (grepl("Allconfounders.name", fn_text)) {
    cat("  BUG FOUND: Uses 'Allconfounders.name' (WRONG)\n")
  } else if (grepl("confounders.name\\s*=", fn_text)) {
    cat("  Uses 'confounders.name' (CORRECT)\n")
  }

  # Check result extraction path
  if (grepl("out_hr", fn_text)) {
    cat("  BUG FOUND: Uses 'out_hr' path (WRONG)\n")
  }
  if (grepl("out_sg\\$result", fn_text) || grepl("out_sg\\[\\[\"result\"\\]\\]", fn_text)) {
    cat("  Uses 'out_sg$result' path (CORRECT)\n")
  }

} else {
  cat("  run_forestsearch_analysis() not found - may be internal to package\n")
}

# =============================================================================
# STEP 6: Manual Test of run_simulation_analysis Components
# =============================================================================

cat("\nSTEP 6: Manual walkthrough of run_simulation_analysis()...\n")

# Simulate what run_simulation_analysis does step by step

cat("\n  6a. Simulating data (already done in Step 2)\n")

cat("\n  6b. Computing true subgroup properties...\n")
size_H_true <- sum(sim_data$flag.harm)
prop_H_true <- mean(sim_data$flag.harm)
cat(sprintf("    True H size: %d (%.1f%%)\n", size_H_true, 100 * prop_H_true))

cat("\n  6c. Creating df_pop...\n")
df_pop <- data.table(
  sim = 1,
  sizeH_true = size_H_true,
  propH_true = prop_H_true,
  sizeHc_true = sum(!sim_data$flag.harm),
  propHc_true = mean(!sim_data$flag.harm)
)
print(df_pop)

cat("\n  6d. Merging parameters with defaults...\n")
if (exists("default_fs_params")) {
  fs_defaults <- default_fs_params()
  cat("    default_fs_params() found\n")
  cat("    Default outcome.name:", fs_defaults$outcome.name, "\n")
  cat("    Default event.name:", fs_defaults$event.name, "\n")
} else {
  cat("    default_fs_params() NOT FOUND\n")
}

cat("\n  6e. Building forestsearch arguments...\n")
cat("    This is where the bug occurs if using 'Allconfounders.name'\n")

# Correct way:
fs_args_correct <- list(
  df.analysis = sim_data,
  confounders.name = confounders_base,  # CORRECT
  outcome.name = "y.sim",
  event.name = "event.sim",
  treat.name = "treat",
  id.name = "id",
  details = TRUE,
  plot.sg = FALSE
)

# Show what the arguments should be
cat("    Correct arguments:\n")
cat("      df.analysis: [", nrow(fs_args_correct$df.analysis), " rows]\n")
cat("      confounders.name:", paste(fs_args_correct$confounders.name, collapse = ", "), "\n")

# =============================================================================
# STEP 7: Test with Fixed run_forestsearch_analysis()
# =============================================================================

cat("\nSTEP 7: Testing with FIXED run_forestsearch_analysis()...\n")

# Define the fixed version locally
run_forestsearch_analysis_FIXED <- function(
    data,
    confounders_name,
    params,
    dgm,
    cox_formula = NULL,
    cox_formula_adj = NULL,
    analysis_label = "FS",
    verbose = FALSE
) {

  cat("    [FIXED] Building forestsearch arguments...\n")

  # CORRECT: use confounders.name
  fs_args <- list(
    df.analysis = data,
    confounders.name = confounders_name,  # FIXED!
    details = verbose,
    plot.sg = verbose
  )

  # Add parameters
  param_names <- c(
    "outcome.name", "event.name", "treat.name", "id.name",
    "use_lasso", "use_grf", "conf_force",
    "hr.threshold", "hr.consistency", "pconsistency.threshold",
    "stop.threshold", "fs.splits", "n.min", "d0.min", "d1.min",
    "maxk", "max.minutes", "pstop_futile", "by.risk", "vi.grf.min"
  )

  for (pn in param_names) {
    if (!is.null(params[[pn]])) {
      fs_args[[pn]] <- params[[pn]]
    }
  }

  cat("    [FIXED] Calling forestsearch()...\n")

  # Run ForestSearch
  fs_result <- tryCatch({
    do.call(forestsearch, fs_args)
  }, error = function(e) {
    cat("    [FIXED] ERROR:", e$message, "\n")
    NULL
  })

  # Check result
  if (!is.null(fs_result)) {
    cat("    [FIXED] forestsearch() completed\n")
    cat("    [FIXED] sg.harm:", ifelse(is.null(fs_result$sg.harm), "NULL",
                                       paste(fs_result$sg.harm, collapse = " & ")), "\n")

    # CORRECT: Check out_sg$result path
    has_result <- !is.null(fs_result$grp.consistency) &&
                  !is.null(fs_result$grp.consistency$out_sg) &&
                  !is.null(fs_result$grp.consistency$out_sg$result) &&
                  nrow(fs_result$grp.consistency$out_sg$result) > 0

    cat("    [FIXED] Has valid result:", has_result, "\n")

    return(list(
      fs_result = fs_result,
      has_result = has_result
    ))
  }

  return(NULL)
}

# Test the fixed version
fs_params_test <- list(
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
  n.min = 60,
  d0.min = 20,
  d1.min = 20,
  maxk = 2,
  max.minutes = 2
)

fixed_result <- run_forestsearch_analysis_FIXED(
  data = sim_data,
  confounders_name = confounders_base,
  params = fs_params_test,
  dgm = dgm,
  analysis_label = "FS_FIXED",
  verbose = TRUE
)

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("SUMMARY\n")
cat("=" , rep("=", 60), "\n", sep = "")

cat("\nDiagnostic results:\n")
cat("  1. DGM creation: OK\n")
cat("  2. simulate_from_gbsg_dgm(): ", ifelse(!is.null(sim_data), "OK", "FAILED"), "\n")
cat("  3. forestsearch() direct call: ", ifelse(!is.null(fs_direct), "OK", "FAILED"), "\n")
cat("  7. Fixed helper function: ", ifelse(!is.null(fixed_result), "OK", "FAILED"), "\n")

cat("\nLikely issues:\n")
cat("  - The PACKAGE version of run_simulation_analysis() has bugs\n")
cat("  - Need to source the FIXED oc_analyses_refactored.R file\n")
cat("  - Or update the package with the corrected code\n")

cat("\nTO FIX:\n")
cat("  Option 1: Source the fixed file before running:\n")
cat("    source('oc_analyses_refactored.R')\n")
cat("    result <- run_simulation_analysis(...)\n")
cat("\n")
cat("  Option 2: Copy fixed functions to package R/ directory and reinstall\n")

cat("\n=== END DIAGNOSTIC ===\n")


dgm <- create_gbsg_dgm(model = "alt", k_inter = 2)

result <- run_simulation_analysis(
  sim_id = 1,
  dgm = dgm,
  n_sample = 500,
  confounders_base = c("z1", "z2", "z3", "z4", "z5", "size", "grade3"),
  run_fs = TRUE,
  run_grf = FALSE,
  verbose = TRUE
)



# 1. Load forestsearch package
devtools::load_all("path/to/forestsearch")

# 2. SOURCE THIS FIXED FILE to override package functions
source("oc_analyses_refactored.R")

# 3. Verify the fix is loaded
cat("Checking run_forestsearch_analysis:\n")
fn_body <- deparse(body(run_forestsearch_analysis))
if (any(grepl("confounders.name =", fn_body))) {
  cat("  ✓ FIXED version loaded\n")
} else {
  cat("  ✗ OLD version still active\n")
}

# 4. Now run
dgm <- create_gbsg_dgm(model = "alt", k_inter = 2)
result <- run_simulation_analysis(
  sim_id = 1,
  dgm = dgm,
  n_sample = 500,
  confounders_base = c("z1", "z2", "z3", "z4", "z5", "size", "grade3"),
  run_fs = TRUE,
  run_grf = FALSE,
  verbose = TRUE
)
