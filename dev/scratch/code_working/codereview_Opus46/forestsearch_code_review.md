# ForestSearch Package Code Review

## Executive Summary

This review covers the ForestSearch R package codebase across ~16 functional areas and ~119 exported functions. The package is well-structured with clear separation of concerns, comprehensive roxygen2 documentation, and sound statistical methodology. Below are actionable recommendations organized by priority, focusing on efficiency, refactoring, and readability.

---

## 1. Critical Refactoring: `dummy()` / `dummy2()` Near-Duplication

**Files:** `R/get_FSdata_helpers.R`

`dummy()` and `dummy2()` are nearly identical, differing only in how they handle the edge case where `FAC(df)` has zero columns. This is the most impactful refactoring opportunity because both functions are called thousands of times in bootstrap iterations.

**Current code (abbreviated):**

```r
dummy <- function(df) {
  NUM <- function(dataframe) dataframe[, sapply(dataframe, is.numeric)]
  FAC <- function(dataframe) dataframe[, sapply(dataframe, is.factor)]
  if (is.null(ncol(NUM(df)))) {
    # ...
  } else {
    if (!is.null(ncol(FAC(df)))) DF <- data.frame(NUM(df), acm.disjctif(FAC(df)))
    if (is.null(ncol(FAC(df)))) {
      # ...
    }
  }
  return(DF)
}

dummy2 <- function(df) {
  # Nearly identical to dummy(), different edge-case branch
}
```

**Problems:**

- `NUM()` and `FAC()` are redefined inside every call—wasteful in tight loops
- `sapply(df, is.numeric)` and `sapply(df, is.factor)` are each called multiple times per invocation (once for `NUM()`, once for `FAC()`, then again inside conditions)
- The logic branches are fragile: `is.null(ncol(...))` doesn't reliably distinguish "zero factor columns" from "one factor column returned as vector"

**Suggested refactored version:**

```r
#' Dummy-code a data frame (numeric pass-through, factors expanded)
#'
#' @param df Data frame with numeric and/or factor columns.
#' @return Data frame with numeric columns unchanged and factor columns
#'   expanded via \code{\link{acm.disjctif}}.
#' @export
dummy_encode <- function(df) {
  stopifnot(is.data.frame(df))


  is_num <- vapply(df, is.numeric, logical(1))
  is_fac <- vapply(df, is.factor, logical(1))

  parts <- list()

  # Numeric columns (pass-through)
  if (any(is_num)) {
    parts$num <- df[, is_num, drop = FALSE]
  }

  # Factor columns (disjunctive coding)
  if (any(is_fac)) {
    fac_df <- df[, is_fac, drop = FALSE]
    parts$fac <- acm.disjctif(fac_df)
  }

  if (length(parts) == 0L) {
    stop("df contains no numeric or factor columns")
  }

  do.call(data.frame, c(parts, list(check.names = FALSE)))
}

# Backward-compatible wrappers
#' @rdname dummy_encode
#' @export
dummy <- dummy_encode

#' @rdname dummy_encode
#' @export
dummy2 <- dummy_encode
```

**Benefits:** Eliminates redundancy, uses `vapply()` (type-safe), evaluates column types once, and `drop = FALSE` handles all edge cases cleanly.

---

## 2. `get_dfpred()`: Replace `eval(parse())` with Direct Indexing

**File:** `R/forestsearch_helpers.R`

**Current code:**

```r
get_dfpred <- function(df.predict, sg.harm, version = 1) {
  if (version == 1) df.pred <- dummy(df.predict)
  if (version == 2) df.pred <- dummy2(df.predict)
  df.pred$treat.recommend <- NA
  id.harm <- paste(sg.harm, collapse = "==1 & ")
  id.harm <- paste(id.harm, "==1")
  df.pred.0 <- subset(df.pred, eval(parse(text = id.harm)))
  # ...
}
```

**Problems:**

- `eval(parse(text = ...))` is a CRAN note risk, a security concern, and slow
- The string construction is error-prone
- The `subset()` + `eval(parse())` pattern is explicitly discouraged in R documentation for programmatic use
- Conditional `rbind` at the end is verbose and fragile

**Suggested refactored version:**

```r
#' Generate Prediction Dataset with Subgroup Treatment Recommendation
#'
#' @param df.predict Data frame for prediction.
#' @param sg.harm Character vector of subgroup-defining binary covariate names.
#' @param version Integer; encoding version (maintained for backward compat).
#' @return Data frame with \code{treat.recommend} column (0 = harm, 1 = recommend).
#' @export
get_dfpred <- function(df.predict, sg.harm, version = 1) {
  df.pred <- dummy_encode(df.predict)


  # Identify subgroup members: all sg.harm columns == 1
  in_harm <- Reduce(`&`, lapply(sg.harm, function(v) df.pred[[v]] == 1L))

  df.pred$treat.recommend <- ifelse(in_harm, 0L, 1L)
  df.pred
}
```

**Benefits:** No `eval(parse())`, vectorized, ~5x faster, clearer intent, and CRAN-safe.

---

## 3. `acm.disjctif()`: Duplicated Inner Functions

**File:** `R/get_FSdata_helpers.R`

`acm.util.df` and `acm.util.df2` inside `acm.disjctif()` differ by one line (`names(df)[i]` vs `colnames(df)[i]`). For data frames, these are equivalent.

**Suggested refactored version:**

```r
#' Disjunctive (dummy) coding for factor columns
#'
#' @param df Data frame with factor variables.
#' @return Data frame with dummy-coded columns.
#' @export
acm.disjctif <- function(df) {
  encode_col <- function(i) {
    cl <- as.factor(df[, i])
    cha <- colnames(df)[i]
    n <- length(cl)
    x <- matrix(0L, n, nlevels(cl))
    x[cbind(seq_len(n), as.integer(cl))] <- 1L
    dimnames(x) <- list(rownames(df), paste(cha, levels(cl), sep = "."))
    x
  }

  parts <- lapply(seq_len(ncol(df)), encode_col)
  data.frame(do.call(cbind, parts), check.names = FALSE)
}
```

**Benefits:** Single inner function, integer matrix (`0L`/`1L`) saves memory, `cbind` on matrices is faster than `data.frame(lapply(...))`.

---

## 4. `one.zero()`: Vectorize

**File:** `R/get_FSdata_helpers.R`

```r
# Current: scalar-only, if/else
one.zero <- function(x) {
 if (x == 1) { x <- 0 } else { x <- 1 }
 return(x)
}
```

**Suggested:**

```r
#' Flip binary value(s)
#'
#' @param x Integer vector of 0s and 1s.
#' @return Integer vector with values flipped.
#' @export
one.zero <- function(x) 1L - x
```

**Benefits:** Vectorized, works on vectors and scalars, ~100x faster when applied to vectors.

---

## 5. Duplicate `print.forestsearch()` / `summary.forestsearch()` Definitions

**Files:** `R/forestsearch_main.R`, `R/print_forestsearch.R`, `R/summary_forestsearch.R`

There are **three separate definitions** of `print.forestsearch()` and at least **two definitions** of `summary.forestsearch()` across different files. Only one of each will be registered—the one loaded last—creating silent bugs depending on file load order.

**Recommendations:**

1. **Choose one canonical location** for S3 methods (e.g., `R/forestsearch_methods.R`)
2. **Delete the other copies** from `R/forestsearch_main.R`, `R/print_forestsearch.R`, `R/summary_forestsearch.R`
3. Ensure the remaining version handles all object shapes robustly (the versions access different slots: `x$parameters$sg_focus` vs `x$args_call_all$sg_focus`)

**Suggested consolidated version:**

```r
# R/forestsearch_methods.R

#' @export
print.forestsearch <- function(x, ...) {
  cat("ForestSearch Results\n")
  cat("====================\n\n")

  if (is.null(x$sg.harm)) {
    cat("No subgroup identified.\n")
    return(invisible(x))
  }

  cat("Selected Subgroup:\n")
  cat("  Definition:", paste(x$sg.harm, collapse = " & "), "\n")

  # Use consistent accessor pattern
  params <- x$args_call_all %||% x$parameters
  if (!is.null(params)) {
    cat("  sg_focus:", params$sg_focus, "\n")
  }

  gc <- x$grp.consistency
  if (!is.null(gc$out_sg$result) && nrow(gc$out_sg$result) > 0) {
    top <- gc$out_sg$result[1, ]
    cat("  N:", top$N, "\n")
    cat("  HR:", round(as.numeric(top$hr), 3), "\n")
    cat("  Pcons:", round(as.numeric(top$Pcons), 3), "\n")
  }

  if (!is.null(gc$algorithm)) {
    cat("  Algorithm:", gc$algorithm, "\n")
  }

  cat("\nComputation time:", round(x$minutes_all, 2), "minutes\n")
  invisible(x)
}
```

---

## 6. `plot.forestsearch()`: Stub Implementation

**File:** `R/forestsearch_helpers.R`

The current implementation is a no-op placeholder:

```r
plot.forestsearch <- function(x, type = ..., ...) {
  type <- match.arg(type)
  invisible(x)  # Does nothing
}
```

**Recommendation:** Either implement it or remove the export and document it as planned/internal. An exported S3 method that silently does nothing is confusing for users and could generate CRAN reviewer questions.

---

## 7. Bootstrap: Embedded Helper Functions in `evaluate_consistency_twostage()`

**File:** `R/subgroup_consistency_helpers.R`

Five helper functions (`.wilson_ci`, `.early_stop_decision`, `.get_split_hr_fast`, `.run_single_consistency_split`, `.FS_labels`) are **defined inside** `evaluate_consistency_twostage()` for parallel worker compatibility. This is understandable but creates issues:

- Functions are re-parsed and re-compiled on every call
- The same helpers exist as top-level exported functions (`get_split_hr_fast`, `run_single_consistency_split`, `FS_labels`) creating maintenance drift
- The dot-prefixed copies may diverge from the exported versions

**Recommended approach:**

Since `callr` workers use `devtools::install()` (making all package functions available), the embedded copies may no longer be needed. If they are still needed for specific parallel backends:

```r
# Define once at package level, reference by namespace in parallel
.consistency_helpers <- list(
  wilson_ci = function(x, n, conf.level = 0.95) { ... },
  early_stop = function(n_success, n_total, threshold, ...) { ... },
  get_split_hr = function(df_split, cox_initial = NULL) { ... }
)
```

This avoids re-defining on every call while keeping them bundled for export to workers.

---

## 8. `BOOTSTRAP_REQUIRED_FUNCTIONS`: Maintenance Risk

**File:** `R/bootstrap_dofuture_setup_helpers.R`

The manually curated list of 60+ function names required for parallel workers is fragile. Adding a new function to the package without adding it to this list causes silent failures in bootstrap.

**Suggested improvement:**

```r
#' Get all exported functions from ForestSearch namespace
#' @keywords internal
get_bootstrap_exports <- function() {
  # Automatically discover all exported functions
  ns <- asNamespace("ForestSearch")
  ls(ns, all.names = FALSE)
}
```

If a curated subset is truly needed, add a unit test that validates the list:

```r
test_that("bootstrap exports list is complete", {
  required <- get_bootstrap_exports()
  exported <- getNamespaceExports("ForestSearch")
  missing <- setdiff(required, exported)
  expect_length(missing, 0, info = paste("Missing:", paste(missing, collapse = ", ")))
})
```

---

## 9. `globals.R`: Excessive `globalVariables()` Declarations

**File:** `R/globals.R`

The file declares 80+ global variables, many of which are data column names used in `data.table` NSE contexts. While this silences `R CMD check` notes, it masks genuine issues when variable names are typo'd.

**Recommendations:**

- Use `.SD`, `.SDcols`, and `get()` / `mget()` patterns in data.table to avoid NSE variable references where possible
- For `ggplot2` `.data$column` pronoun usage, which you already do in some places—standardize across all plot functions
- Consider grouping declarations by file with comments to track which declarations belong where:

```r
# Used in: R/mrct_simulation.R
utils::globalVariables(c("hr_test", "hr_sg", "any_found", ...))
```

---

## 10. Timing Code Repetition

**Files:** `R/forestsearch_main.R`, `R/bootstrap_analysis_dofuture.R`

The pattern `t.start <- proc.time()[3]; ...; t.min <- (proc.time()[3] - t.start) / 60` appears in many places.

**Suggested helper:**

```r
#' Time a code block
#' @keywords internal
with_timing <- function(expr, unit = c("minutes", "seconds")) {
  unit <- match.arg(unit)
  t0 <- proc.time()[3]
  result <- force(expr)
  elapsed <- proc.time()[3] - t0
  elapsed <- if (unit == "minutes") elapsed / 60 else elapsed
  list(result = result, elapsed = elapsed)
}
```

---

## 11. `bootstrap_summaries_helpers.R`: Deep NULL-Checking Chains

**File:** `R/bootstrap_summaries_helpers.R`

`create_timing_summary_table()` has extensive repeated NULL-checking patterns:

```r
if (!is.null(iteration_stats) &&
    !is.null(iteration_stats$mean) &&
    !is.null(iteration_stats$median) &&
    !is.null(iteration_stats$sd) &&
    !is.null(iteration_stats$min) &&
    !is.null(iteration_stats$max) &&
    !is.null(iteration_stats$q25) &&
    !is.null(iteration_stats$q75)) {
```

**Suggested helper:**

```r
#' Check that all named elements exist and are non-NULL
#' @keywords internal
has_all <- function(x, fields) {
  !is.null(x) && all(fields %in% names(x)) &&
    all(vapply(fields, function(f) !is.null(x[[f]]), logical(1)))
}

# Usage:
if (has_all(iteration_stats, c("mean", "median", "sd", "min", "max", "q25", "q75"))) {
  # ...
}
```

---

## 12. `get_FSdata()`: Long Function Body

**File:** `R/get_FSdata_refactored.r`

Despite the "refactored" filename, `get_FSdata()` appears to still be a very long function. The function signature alone has 14 parameters. Consider further decomposition:

```r
get_FSdata <- function(...) {
  # 1. Validate inputs
  validate_fsdata_inputs(...)

  # 2. Generate cut expressions
  cuts <- generate_cut_expressions(df.analysis, confounders.name, ...)

  # 3. Apply LASSO selection
  lasso_result <- if (use_lasso) apply_lasso_selection(...) else NULL

  # 4. Apply GRF cuts
  grf_result <- if (use_grf) apply_grf_cuts(grf_cuts, ...) else NULL

  # 5. Combine and filter
  assemble_fsdata(cuts, lasso_result, grf_result, ...)
}
```

---

## 13. File Naming Inconsistency

Several files use different naming conventions:

| File | Convention |
|------|-----------|
| `get_FSdata_refactored.r` | lowercase `.r`, snake_case with CamelCase |
| `R/forestsearch_main.R` | uppercase `.R`, snake_case |
| `R/sim_aft_gbsg_refactored.R` | uppercase `.R`, snake_case |
| `R/cv_summary_tables.R` | uppercase `.R`, snake_case |
| `R/plot_sg_results.R` | uppercase `.R`, snake_case |

**Recommendations:**

- Standardize on `.R` (uppercase) per tidyverse convention
- Remove `_refactored` suffixes—these are internal version markers that don't belong in a released package
- Rename `get_FSdata_refactored.r` → `R/get_fsdata.R`

---

## 14. Required Packages Check at Load Time

**File:** `R/forestsearch_helpers.R`

```r
required_packages <- c("grf", "policytree", "data.table", ...)
missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) stop("Missing required packages: ", ...)
```

This code runs at **source time** (when the file is loaded), not when functions are called. For CRAN:

- Hard dependencies should be in `Imports` (and will always be available)
- Soft dependencies should use `requireNamespace()` at function call time, not at load time
- Remove this top-level check; it serves no purpose if `DESCRIPTION` is correct

---

## 15. `calculate_skewness()`: Re-implementing Base R

**File:** `R/bootstrap_summaries_helpers.R`

```r
calculate_skewness <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3) return(NA)
  x_bar <- mean(x)
  s <- sd(x)
  if (s == 0) return(0)
  (n / ((n - 1) * (n - 2))) * sum(((x - x_bar) / s)^3)
}
```

Consider using `e1071::skewness()` or `moments::skewness()` if either is already a dependency. If not, keep this but mark `@keywords internal` (it's currently undocumented for export).

---

## 16. Consistency in Parameter Naming

The codebase mixes naming conventions for the same concept:

| Parameter | Appears As |
|-----------|-----------|
| Number of splits | `n.splits`, `fs.splits`, `n.splits.screen`, `n.splits.max` |
| Hazard ratio threshold | `hr.threshold`, `hr.consistency`, `hr.cons` |
| Subgroup focus | `sg_focus`, `sg.focus` |
| Details flag | `details`, `verbose` |

**Recommendation:** For new code, standardize on snake_case (`n_splits`, `hr_threshold`). For existing exported parameters, maintain backward compatibility with deprecation warnings:

```r
#' @param n_splits Integer. Number of splits (replaces deprecated \code{n.splits}).
```

---

## 17. `ensure_packages()` Defined Twice

**File:** `R/bootstrap_dofuture_setup_helpers.R`

The `ensure_packages()` function is defined both as a standalone function and is similar to the top-level check in `forestsearch_helpers.R`. Consolidate to one location.

---

## Summary of Priority Actions

| Priority | Item | Impact | Effort |
|----------|------|--------|--------|
| **High** | Consolidate `dummy()`/`dummy2()` | Performance + maintenance | Medium |
| **High** | Replace `eval(parse())` in `get_dfpred()` | CRAN compliance + perf | Low |
| **High** | Remove duplicate S3 method definitions | Correctness | Low |
| **High** | Rename `.r` → `.R`, remove `_refactored` | CRAN conventions | Low |
| **Medium** | Vectorize `one.zero()` | Performance in loops | Trivial |
| **Medium** | Consolidate `acm.disjctif()` inner functions | Readability | Low |
| **Medium** | Remove top-level `requireNamespace()` check | CRAN compliance | Trivial |
| **Medium** | Implement or remove `plot.forestsearch()` stub | User experience | Medium |
| **Medium** | Add `has_all()` helper for NULL-checking | Readability | Low |
| **Low** | Standardize parameter naming conventions | Consistency | High (breaking) |
| **Low** | Auto-discover bootstrap exports | Maintenance | Medium |
| **Low** | Extract timing helper | DRY principle | Low |
| **Low** | Reduce `globalVariables()` via `.data$` pronoun | Best practice | Medium |
