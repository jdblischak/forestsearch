# =============================================================================
# apply_cran_fixes.R - ForestSearch CRAN Check Fixes
# =============================================================================
#
# Run this script from your ForestSearch package root directory to fix
# the R CMD check WARNING and NOTE issues.
#
# Issues resolved:
#   WARNING: Undocumented arguments 'stop_threshold' and 'xlog'
#   NOTE: Lost braces in Rd files (3 occurrences)
#
# Usage:
#   setwd("path/to/ForestSearch")
#   source("apply_cran_fixes.R")
#
# =============================================================================

cat(paste(rep("=", 77), collapse = ""), "\n")
cat("Applying CRAN check fixes to ForestSearch\n")
cat(paste(rep("=", 77), collapse = ""), "\n\n")

# -----------------------------------------------------------------------------
# Helper function for safe string replacement
# -----------------------------------------------------------------------------
safe_replace <- function(file_path, old_text, new_text, description) {
  if (!file.exists(file_path)) {
    cat(sprintf("  [SKIP] File not found: %s\n", file_path))
    return(FALSE)
  }
 
  content <- readLines(file_path, warn = FALSE)
  content_str <- paste(content, collapse = "\n")
 
  if (!grepl(old_text, content_str, fixed = TRUE)) {
    cat(sprintf("  [SKIP] Pattern not found in %s\n", basename(file_path)))
    return(FALSE)
  }
 
  new_content <- gsub(old_text, new_text, content_str, fixed = TRUE)
  writeLines(new_content, file_path)
  cat(sprintf("  [OK] %s: %s\n", basename(file_path), description))
  return(TRUE)
}

# =============================================================================
# FIX 1: Add @param stop_threshold to forestsearch()
# =============================================================================
cat("FIX 1: Adding @param stop_threshold documentation...\n")

# Try primary pattern (pconsistency.threshold -> showten_subgroups)
fix1_done <- safe_replace(
  "R/forest_search_revised.R",
  "#' @param pconsistency.threshold Numeric. Minimum consistency proportion. Default 0.90.\n#' @param showten_subgroups",
  "#' @param pconsistency.threshold Numeric. Minimum consistency proportion. Default 0.90.
#' @param stop_threshold Numeric. Early stopping threshold for consistency
#'   evaluation. When a candidate subgroup's estimated consistency probability
#'   exceeds this threshold, evaluation stops early. Default 0.95.
#' @param showten_subgroups",
  "Added stop_threshold parameter"
)

# =============================================================================
# FIX 2: Add @param xlog to plot_subgroup_results_forestplot()
# =============================================================================
cat("\nFIX 2: Adding @param xlog documentation...\n")

fix2_done <- safe_replace(
  "R/plot_subgroup_results_forestplot.R",
  "#' @param est.scale Character. Estimate scale: \"hr\" or \"1/hr\" (default: \"hr\").\n#' @param title_text",
  "#' @param est.scale Character. Estimate scale: \"hr\" or \"1/hr\" (default: \"hr\").
#' @param xlog Logical. If TRUE (default), the x-axis is plotted on a
#'   logarithmic scale. This is standard for hazard ratio forest plots
#'   where equal distances represent equal relative effects.
#' @param title_text",
  "Added xlog parameter"
)

# =============================================================================
# FIX 3: Escape braces in format_bootstrap_table() - sg_definition
# =============================================================================
cat("\nFIX 3: Escaping braces in format_bootstrap_table documentation...\n")

safe_replace(
  "R/bootstrap_summaries_helpers.R",
  "(e.g., \"{age>=50} & {nodes>=3}\")",
  "(e.g., \"\\{age>=50\\} & \\{nodes>=3\\}\")",
  "Escaped braces in sg_definition"
)

# =============================================================================
# FIX 4: Escape braces in summarize_bootstrap_results() - sgharm
# =============================================================================
cat("\nFIX 4: Escaping braces in summarize_bootstrap_results documentation...\n")

safe_replace(
  "R/summarize_bootstrap_results.R",
  "(e.g., c(\"{age>=50}\", \"{nodes>=3}\"))",
  "(e.g., c(\"\\{age>=50\\}\", \"\\{nodes>=3\\}\"))",
  "Escaped braces in sgharm"
)

# =============================================================================
# FIX 5: Escape braces in summarize_bootstrap_subgroups() - original_sg
# =============================================================================
cat("\nFIX 5: Escaping braces in summarize_bootstrap_subgroups documentation...\n")

safe_replace(
  "R/summarize_bootstrap_subgroups.R",
  "(e.g., c(\"{age>=50}\", \"{nodes>=3}\") for a 2-factor subgroup)",
  "(e.g., c(\"\\{age>=50\\}\", \"\\{nodes>=3\\}\") for a 2-factor subgroup)",
  "Escaped braces in original_sg"
)

# =============================================================================
# Final instructions
# =============================================================================
cat("\n")
cat(paste(rep("=", 77), collapse = ""), "\n")
cat("FIXES COMPLETE\n")
cat(paste(rep("=", 77), collapse = ""), "\n")
cat("\n")
cat("Next steps:\n")
cat("  1. Review changes in R/ files\n")
cat("  2. Run: devtools::document()\n")
cat("  3. Run: devtools::check()\n")
cat("\n")
cat("Expected result: 0 errors, 0 warnings, 0 notes\n")
cat("\n")
cat("If any fixes showed [SKIP], see CRAN_FIXES_MANUAL.md for manual instructions.\n")
cat(paste(rep("=", 77), collapse = ""), "\n")
