#!/usr/bin/env Rscript
# fix_all_nonascii.R
# Script to fix non-ASCII characters in all R package files

# List of files with non-ASCII characters (from R CMD check)
files_to_fix <- c(
  "R/bootstrap_summaries_helpers.R",
  "R/cox_ahr_cde_wrapper.R",
  "R/cox_spline_fit.R",
  "R/format_subgroup_summary_tables.R",
  "R/generate_aft_dgm_helpers.R",
  "R/get_FSdata_refactored.r",
  "R/summarize_bootstrap_results.R",
  "R/summarize_bootstrap_subgroups.R",
  "R/summary_utility_functions.R"
)

fix_nonascii_in_file <- function(file_path) {
  
  if (!file.exists(file_path)) {
    cat("File not found:", file_path, "\n")
    return(FALSE)
  }
  
  cat("Processing:", file_path, "\n")
  
  # Read file
  lines <- readLines(file_path, encoding = "UTF-8", warn = FALSE)
  original_lines <- lines
  
  # Replace dagger symbols with Unicode escapes (for footnotes)
  lines <- gsub("\u2020", "\\\\u2020", lines)  # † (dagger)
  lines <- gsub("\u2021", "\\\\u2021", lines)  # ‡ (double dagger)
  
  # Replace em dash with double hyphen
  lines <- gsub("\u2014", "--", lines)         # — (em dash)
  lines <- gsub("\u2013", "--", lines)         # – (en dash)
  
  # Replace checkmark with Unicode escape
  lines <- gsub("\u2713", "\\\\u2713", lines)  # ✓ (checkmark)
  
  # Replace warning symbol with Unicode escape
  lines <- gsub("\u26A0", "\\\\u26A0", lines)  # ⚠ (warning)
  
  # Replace bullet with hyphen
  lines <- gsub("\u2022", "-", lines)          # • (bullet)
  
  # Replace arrows with ASCII equivalents
  lines <- gsub("\u2192", "->", lines)         # → (right arrow)
  lines <- gsub("\u2190", "<-", lines)         # ← (left arrow)
  
  # Replace comparison symbols with ASCII equivalents
  lines <- gsub("\u2265", ">=", lines)         # ≥ (greater than or equal)
  lines <- gsub("\u2264", "<=", lines)         # ≤ (less than or equal)
  
  # Replace multiplication and division symbols
  lines <- gsub("\u00D7", "*", lines)          # × (multiplication)
  lines <- gsub("\u00F7", "/", lines)          # ÷ (division)
  
  # Replace minus sign with hyphen
  lines <- gsub("\u2212", "-", lines)          # − (minus sign)
  
  # Replace ellipsis with three periods
  lines <- gsub("\u2026", "...", lines)        # … (ellipsis)
  
  # Replace smart quotes with straight quotes
  lines <- gsub("[\u201C\u201D]", '"', lines)  # " " (curly double quotes)
  lines <- gsub("[\u2018\u2019]", "'", lines)  # ' ' (curly single quotes)
  
  # Check if any changes were made
  if (identical(lines, original_lines)) {
    cat("  No changes needed\n")
    return(TRUE)
  }
  
  # Write back to file
  writeLines(lines, file_path, useBytes = TRUE)
  
  # Verify no non-ASCII characters remain
  remaining <- grep("[^\x01-\x7F]", lines)
  if (length(remaining) > 0) {
    cat("  WARNING: Still has non-ASCII characters on lines:", 
        paste(remaining, collapse = ", "), "\n")
    cat("  Run tools::showNonASCIIfile('", file_path, "') to investigate\n", sep = "")
    return(FALSE)
  }
  
  cat("  Fixed successfully\n")
  return(TRUE)
}

# Main execution
cat("\n=== Fixing Non-ASCII Characters in R Package Files ===\n\n")

# Make sure we're in the package root
if (!file.exists("DESCRIPTION")) {
  cat("ERROR: Not in package root directory. Please run this from your package root.\n")
  cat("Current directory:", getwd(), "\n")
  quit(status = 1)
}

# Process all files
results <- sapply(files_to_fix, fix_nonascii_in_file)

# Summary
cat("\n=== Summary ===\n")
cat("Files processed:", length(results), "\n")
cat("Files fixed:", sum(results), "\n")
cat("Files with issues:", sum(!results), "\n")

if (all(results)) {
  cat("\nAll files fixed successfully!\n")
  cat("\nNext steps:\n")
  cat("1. Run: devtools::document()\n")
  cat("2. Run: devtools::check()\n")
  cat("3. Verify no non-ASCII warnings remain\n")
} else {
  cat("\nSome files still have issues. Review warnings above.\n")
}

cat("\n")
