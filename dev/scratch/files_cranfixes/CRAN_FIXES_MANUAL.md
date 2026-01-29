# ForestSearch CRAN Check Fixes - Complete Guide

## Summary

This document provides **complete fixes** for the following R CMD check issues:

| Type | File | Issue |
|------|------|-------|
| **WARNING** | `forestsearch.Rd` | Undocumented argument `stop_threshold` |
| **WARNING** | `plot_subgroup_results_forestplot.Rd` | Undocumented argument `xlog` |
| **NOTE** | `format_bootstrap_table.Rd` | Lost braces in Rd (line 27) |
| **NOTE** | `summarize_bootstrap_results.Rd` | Lost braces in Rd (line 17) |
| **NOTE** | `summarize_bootstrap_subgroups.Rd` | Lost braces in Rd (line 16) |

---

## Fix 1: Add `@param stop_threshold` to `forestsearch()`

**File:** `R/forest_search_revised.R`

### BEFORE (around line with pconsistency.threshold):
```r
#' @param pconsistency.threshold Numeric. Minimum consistency proportion. Default 0.90.
#' @param showten_subgroups Logical. Show top 10 subgroups. Default FALSE.
```

### AFTER:
```r
#' @param pconsistency.threshold Numeric. Minimum consistency proportion. Default 0.90.
#' @param stop_threshold Numeric. Early stopping threshold for consistency
#'   evaluation. When a candidate subgroup's estimated consistency probability
#'   exceeds this threshold, evaluation stops early. Default 0.95.
#' @param showten_subgroups Logical. Show top 10 subgroups. Default FALSE.
```

---

## Fix 2: Add `@param xlog` to `plot_subgroup_results_forestplot()`

**File:** `R/plot_subgroup_results_forestplot.R`

### BEFORE:
```r
#' @param est.scale Character. Estimate scale: "hr" or "1/hr" (default: "hr").
#' @param title_text Character. Plot title (default: NULL).
```

### AFTER:
```r
#' @param est.scale Character. Estimate scale: "hr" or "1/hr" (default: "hr").
#' @param xlog Logical. If TRUE (default), the x-axis is plotted on a
#'   logarithmic scale. This is standard for hazard ratio forest plots
#'   where equal distances represent equal relative effects.
#' @param title_text Character. Plot title (default: NULL).
```

---

## Fix 3: Escape curly braces in `format_bootstrap_table()`

**File:** `R/bootstrap_summaries_helpers.R`

### BEFORE:
```r
#' @param sg_definition Character. Subgroup definition string to display as footnote
#'   (e.g., "{age>=50} & {nodes>=3}"). If NULL, no subgroup footnote is added.
```

### AFTER:
```r
#' @param sg_definition Character. Subgroup definition string to display as footnote
#'   (e.g., "\{age>=50\} & \{nodes>=3\}"). If NULL, no subgroup footnote is added.
```

---

## Fix 4: Escape curly braces in `summarize_bootstrap_results()`

**File:** `R/summarize_bootstrap_results.R`

### BEFORE:
```r
#' @param sgharm The selected subgroup object from forestsearch results. Can be:
#'   \itemize{
#'     \item Character vector of factor definitions (e.g., c("{age>=50}", "{nodes>=3}"))
```

### AFTER:
```r
#' @param sgharm The selected subgroup object from forestsearch results. Can be:
#'   \itemize{
#'     \item Character vector of factor definitions (e.g., c("\{age>=50\}", "\{nodes>=3\}"))
```

---

## Fix 5: Escape curly braces in `summarize_bootstrap_subgroups()`

**File:** `R/summarize_bootstrap_subgroups.R`

### BEFORE:
```r
#' @param original_sg Character vector. Original subgroup definition from main
#'   analysis (e.g., c("{age>=50}", "{nodes>=3}") for a 2-factor subgroup)
```

### AFTER:
```r
#' @param original_sg Character vector. Original subgroup definition from main
#'   analysis (e.g., c("\{age>=50\}", "\{nodes>=3\}") for a 2-factor subgroup)
```

---

## Quick Reference: Find & Replace Commands

For users who prefer find-and-replace:

| File | Find | Replace |
|------|------|---------|
| `bootstrap_summaries_helpers.R` | `"{age>=50} & {nodes>=3}"` | `"\{age>=50\} & \{nodes>=3\}"` |
| `summarize_bootstrap_results.R` | `c("{age>=50}", "{nodes>=3}")` | `c("\{age>=50\}", "\{nodes>=3\}")` |
| `summarize_bootstrap_subgroups.R` | `c("{age>=50}", "{nodes>=3}") for a 2-factor` | `c("\{age>=50\}", "\{nodes>=3\}") for a 2-factor` |

---

## Verification Steps

After making all changes:

```r
# Step 1: Regenerate documentation
devtools::document()

# Step 2: Run full check
devtools::check()
```

**Expected result:** 
```
── R CMD check results ─────────────────────────────────
0 errors ✔ | 0 warnings ✔ | 0 notes ✔
```

---

## Technical Notes

### Why escape braces?

In Rd (R documentation) format, curly braces `{}` are used for macro arguments (like `\code{...}`). When you want **literal** curly braces to appear in the output, you must escape them as `\{` and `\}`.

In R source files (roxygen2 comments), you write `\{` which roxygen2 converts to the properly escaped form in the generated `.Rd` file.

### Parameter documentation order

The `@param` entries should match the order of parameters in the function signature:

- `stop_threshold` appears after `pconsistency.threshold` in `forestsearch()`
- `xlog` appears after `est.scale` in `plot_subgroup_results_forestplot()`

### Files included in this fix package

| File | Description |
|------|-------------|
| `CRAN_FIXES_MANUAL.md` | This guide with complete before/after examples |
| `apply_cran_fixes.R` | Automated R script to apply fixes |
| `fix1_forestsearch_documentation.R` | Complete fixed roxygen2 block for forestsearch() |
| `fix2_plot_forestplot_documentation.R` | Complete fixed roxygen2 block for plot function |
| `fix3_bootstrap_helpers_documentation.R` | Fixed roxygen2 for format_bootstrap_table() |
| `fix4_summarize_bootstrap_results_documentation.R` | Fixed roxygen2 for summarize_bootstrap_results() |
| `fix5_summarize_bootstrap_subgroups_documentation.R` | Fixed roxygen2 for summarize_bootstrap_subgroups() |
