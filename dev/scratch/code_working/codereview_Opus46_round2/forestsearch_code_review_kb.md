# ForestSearch Package — Code Review (Knowledge Base)

**Date:** 2026-02-08
**Scope:** Full codebase from project Knowledge Base
**Focus:** CRAN compliance, tidyverse style, correctness, performance

---

## Executive Summary

The ForestSearch package is well-structured with strong modular design, comprehensive roxygen2 documentation, and sound statistical methodology. Many issues from the original review document (`forestsearch_code_review.md`) have already been addressed, including:

- ✅ `dummy()` / `dummy2()` consolidated into `dummy_encode()` with aliases
- ✅ `acm.disjctif()` refactored with single `encode_col()` inner function
- ✅ `one.zero()` vectorized to `1L - x`
- ✅ `get_dfpred()` refactored with `evaluate_comparison()` — no `eval(parse())`
- ✅ `plot.forestsearch()` fully implemented, dispatches to `plot_sg_results()`
- ✅ `get_bootstrap_exports()` auto-discovers namespace exports
- ✅ `evaluate_cuts_once()` caches cut evaluations
- ✅ Rich S3 ecosystem: `print`/`plot` methods for `fs_kfold`, `fs_tenfold`, `fs_forestplot`, `fs_forest_theme`, `fs_sg_plot`

This review identifies **15 remaining findings** organized by priority.

---

## Critical (3 items — fix before CRAN submission)

### C1. Duplicate S3 Methods: `print.forestsearch()` / `summary.forestsearch()`

**File:** `R/forestsearch_main.R` (lines visible in KB)

Both `print.forestsearch()` and `summary.forestsearch()` are defined in `R/forestsearch_main.R`. The old review document also references definitions in `R/print_forestsearch.R` and `R/summary_forestsearch.R`. If those files still exist, only one definition per method gets registered (last loaded alphabetically wins), creating silent bugs.

**Verified in KB:** The `forestsearch_main.R` versions access `x$grp.consistency$sg_focus` and `x$args_call_all`. If other files access `x$parameters$sg_focus`, users could get different output depending on load order.

**Recommendation:** Audit for duplicate files. Consolidate into a single `R/forestsearch_methods.R`:

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

  # Robust accessor: try both possible slot locations
  gc <- x$grp.consistency
  if (!is.null(gc)) {
    cat("  sg_focus:", gc$sg_focus %||% x$args_call_all$sg_focus, "\n")

    if (!is.null(gc$out_sg$result) && nrow(gc$out_sg$result) > 0) {
      top <- gc$out_sg$result[1, ]
      cat("  N:", top$N, "\n")
      cat("  HR:", round(as.numeric(top$hr), 3), "\n")
      cat("  Pcons:", round(as.numeric(top$Pcons), 3), "\n")
    }

    if (!is.null(gc$algorithm)) {
      cat("  Algorithm:", gc$algorithm, "\n")
    }

    if (!is.null(gc$n_candidates_evaluated)) {
      cat("  Candidates evaluated:", gc$n_candidates_evaluated, "\n")
      cat("  Candidates passed:", gc$n_passed, "\n")
    }
  }

  cat("\nComputation time:", round(x$minutes_all, 2), "minutes\n")
  invisible(x)
}

#' @export
summary.forestsearch <- function(object, ...) {
  # ... single canonical version here ...
}
```

Then delete the method definitions from `R/forestsearch_main.R` (and `R/print_forestsearch.R`, `R/summary_forestsearch.R` if they exist).

---

### C2. Remaining `eval(parse())` in Three Locations

While `get_dfpred()` has been fixed, `eval(parse())` remains in at least **three places**:

| Location | Usage | Risk |
|----------|-------|------|
| `R/subgroup_consistency_helpers.R` — `sg_consistency_out()` | `subset(df, eval(parse(text = id.harm)))` for plotting | CRAN note |
| `R/plot_subgroup_results_forestplot.R` — reference & posthoc subgroups | `subset(df_analysis, eval(parse(text = sg$subset_expr)))` | CRAN note |
| `R/plot_km_band_forestsearch.R` — `create_reference_subgroup_columns()` | `eval(parse(text = subset_expr), envir = df)` | CRAN note |

**Recommendation for `sg_consistency_out()` — replace with vectorized indexing:**

```r
# Instead of:
id.harm <- paste(paste(sg.harm, collapse = "==1 & "), "==1")
df.sub <- subset(df, eval(parse(text = id.harm)))

# Use:
in_harm <- Reduce(`&`, lapply(sg.harm, function(v) df[[v]] == 1L))
df.sub  <- df[in_harm, , drop = FALSE]
df.subC <- df[!in_harm, , drop = FALSE]
```

**Recommendation for forest plot and KM band functions:** These accept user-provided `subset_expr` strings, so they're harder to eliminate. Create a safe evaluation helper:

```r
#' Safely subset data using an expression string
#'
#' @param df Data frame.
#' @param expr Character. Subset expression (e.g., "BM > 1 & size > 19").
#' @return Subset of df, or NULL on failure.
#' @keywords internal
safe_subset <- function(df, expr) {
  tryCatch({
    # Parse in a controlled environment containing only df columns
    env <- list2env(df, parent = baseenv())
    idx <- eval(parse(text = expr), envir = env)
    df[idx, , drop = FALSE]
  }, error = function(e) {
    warning("Failed to evaluate subset expression: ", expr,
            " - ", e$message, call. = FALSE)
    NULL
  })
}
```

This isolates evaluation from the global environment and is clearer about intent. Even if CRAN still flags `eval(parse())`, having it wrapped in a documented, controlled helper with `tryCatch` is defensible for user-supplied expressions.

---

### C3. `_pkgdown.yml` UTF-8 Encoding (Mojibake)

**Line 86:**
```
infinitesimal jackknife variance estimation (LeÃ³n et al., 2024).
```

"LeÃ³n" is mojibake for "León" — the file was likely saved with a non-UTF-8 encoding or double-encoded.

**Fix:** Replace with proper UTF-8:
```yaml
infinitesimal jackknife variance estimation (León et al., 2024).
```

Ensure the file is saved as UTF-8 without BOM. In RStudio: *File → Save with Encoding → UTF-8*.

---

## High (4 items — should fix before CRAN)

### H1. `BOOTSTRAP_REQUIRED_FUNCTIONS` List Is Redundant

**File:** `R/bootstrap_dofuture_setup_helpers.R`

The manually curated `BOOTSTRAP_REQUIRED_FUNCTIONS` list (60+ entries organized into 6 categories) is never actually used for the `foreach` globals specification. Instead, `bootstrap_results()` calls `get_bootstrap_exports()` which auto-discovers all namespace exports:

```r
# In bootstrap_results():
.options.future = list(
  globals = structure(TRUE, add = c(
    get_bootstrap_exports(),       # <-- auto-discovers everything
    "fs.est", "df_boot_analysis", ...
  ))
)
```

The curated list exists as documentation but risks becoming stale and misleading.

**Recommendation:** Either:
1. **Remove `BOOTSTRAP_REQUIRED_FUNCTIONS`** entirely (since `get_bootstrap_exports()` handles it), or
2. **Add a unit test** that validates the curated list stays in sync:

```r
test_that("BOOTSTRAP_REQUIRED_FUNCTIONS covers all exports", {
  curated <- unlist(BOOTSTRAP_REQUIRED_FUNCTIONS)
  auto <- get_bootstrap_exports()
  missing <- setdiff(curated, auto)
  expect_length(missing, 0)
})
```

---

### H2. Commented-Out `on.exit()` Block in Bootstrap Main

**File:** `R/bootstrap_dofuture_main.r`

~20 lines of commented-out `on.exit()` parallel cleanup code remain:

```r
# on.exit({
#   # Ensure workers are properly shut down
#   if (exists(".Last.future.plan")) {
#     future::plan(.Last.future.plan)
#   } else {
#     future::plan("sequential")  # Safe fallback
#   }
#   ...
# })
```

**Risk:** CRAN reviewers flag commented-out code. More importantly, this suggests parallel workers may not be cleaned up properly.

**Recommendation:** Either:
- **Uncomment and activate** a proper cleanup strategy, or
- **Delete entirely** and document in the function's `@details` why cleanup is handled elsewhere

---

### H3. File Extension: `.r` vs `.R`

**File:** `R/bootstrap_dofuture_main.r`

This file uses lowercase `.r`. While R itself is case-insensitive on most platforms, CRAN checks on Linux are case-sensitive, and tidyverse convention mandates `.R`.

**Fix:** `git mv R/bootstrap_dofuture_main.r R/bootstrap_dofuture_main.R`

Also check for `R/get_FSdata_refactored.r` (referenced in the old review — if it still exists, rename it too and consider dropping the `_refactored` suffix).

---

### H4. Duplicate `@importFrom` Declarations

**Files:** `R/globals.R` and `R/forestsearch-package.R`

Both files declare overlapping `@importFrom` tags. For example:

- `R/globals.R`: `@importFrom stats AIC BIC model.frame vcov sd median quantile ...`
- `R/forestsearch-package.R`: `@importFrom stats as.formula coef na.exclude na.omit predict setNames`
- `R/globals.R`: `@importFrom graphics hist barplot mtext axis title rug text rect plot.window`
- `R/forestsearch-package.R`: `@importFrom graphics par plot plot.new plot.window lines points abline ...`

Plus `R/globals.R` has `@import data.table` and `@import survival` (blanket imports) while individual functions also have targeted `@importFrom data.table` declarations.

**Risks:**
- `@import data.table` + `@import survival` import *all* symbols, potentially masking base R functions
- Scattered `@importFrom` across multiple files makes namespace hard to audit

**Recommendation:**
1. Replace `@import data.table` with targeted `@importFrom data.table` for the specific symbols used (`:=`, `.SD`, `.N`, `.I`, `.GRP`, `data.table`, `rbindlist`, `setnames`, `setcolorder`, `is.data.table`, `copy`, etc.)
2. Consolidate all `@importFrom` declarations into `R/forestsearch-package.R`
3. Remove the duplicate set from `R/globals.R`

---

## Medium (5 items)

### M1. `evaluate_cuts_once()` Still Uses `eval(parse())`

**File:** `R/get_FSdata_helpers.R`

The caching function evaluates cuts via:
```r
result <- eval(parse(text = thiscut), envir = df)
```

This is called at startup to cache results (not in a hot loop), so performance isn't the concern — CRAN compliance is. The cut expressions are internally generated strings like `"q1 == 1"`, not user input.

**Recommendation:** Since these are generated from known confounders, consider pre-building the expressions as calls:

```r
evaluate_cuts_once <- function(confs, df, details = FALSE) {
  n_confs <- length(confs)
  evaluations <- vector("list", n_confs)
  is_valid <- logical(n_confs)
  has_error <- logical(n_confs)

  for (i in seq_along(confs)) {
    tryCatch({
      # Use evaluate_comparison() instead of eval(parse())
      result <- evaluate_comparison(confs[i], df)
      evaluations[[i]] <- as.logical(result)
      is_valid[i] <- length(unique(result)) > 1
    }, error = function(e) {
      has_error[i] <<- TRUE
      is_valid[i] <<- FALSE
    })
  }

  list(evaluations = evaluations, is_valid = is_valid, has_error = has_error)
}
```

This reuses the already-CRAN-safe `evaluate_comparison()` you built for `get_dfpred()`.

---

### M2. `ensure_packages()` Defined in Multiple Locations

**File:** `R/bootstrap_dofuture_setup_helpers.R`

The `ensure_packages()` utility is defined here and likely exists in other files too (the old review flagged this). Consolidate to one location — either in this file (since bootstrap is its primary consumer) or in a shared `R/utils.R`.

---

### M3. `globals.R` — 100+ Declarations with Potential for Masking

The file declares 100+ global variables spanning data.table symbols, column names, parameters, and simulation variables. While necessary to silence `R CMD check`, this creates risk:
- Typos in column names (`hr` vs `Hr`) become invisible
- Variables intended for one function's scope leak into the global declaration

**Recommendation:** Add `# Source:` comments to each section indicating which file(s) use each set:

```r
# ============================================================================
# Used in: R/bootstrap_analysis_dofuture.R
# ============================================================================
  "boot", "boot_id", "iteration", "converged",
```

This makes auditing much easier and helps identify variables that can be removed when code changes.

---

### M4. Deep NULL-Checking Chains

**File:** `R/bootstrap_summaries_helpers.R` (referenced in old review)

The pattern of checking 7+ nested NULL fields is fragile:

```r
if (!is.null(x) && !is.null(x$a) && !is.null(x$a$b) && ...) { }
```

**Recommendation:** Add a small helper:

```r
#' Check that all named elements exist and are non-NULL
#' @keywords internal
has_all <- function(x, fields) {
  !is.null(x) && all(vapply(fields, function(f) !is.null(x[[f]]), logical(1)))
}

# Usage:
if (has_all(iteration_stats, c("mean", "median", "sd", "min", "max"))) { ... }
```

---

### M5. Timing Code Duplication

**Files:** `R/forestsearch_main.R`, `R/bootstrap_analysis_dofuture.R`, `R/bootstrap_dofuture_main.r`

The pattern:
```r
t.start <- proc.time()[3]
# ... work ...
elapsed <- (proc.time()[3] - t.start) / 60
```

appears at least 5 times across files.

**Recommendation:**

```r
#' Time a code block
#' @param expr Expression to time.
#' @param unit Character. "minutes" or "seconds".
#' @return List with `result` and `elapsed`.
#' @keywords internal
with_timing <- function(expr, unit = c("minutes", "seconds")) {
  unit <- match.arg(unit)
  t0 <- proc.time()[3]
  result <- force(expr)
  elapsed <- proc.time()[3] - t0
  if (unit == "minutes") elapsed <- elapsed / 60
  list(result = result, elapsed = elapsed)
}
```

---

## Low (3 items)

### L1. Parameter Naming Inconsistency

The codebase mixes dot-case and snake_case for the same concepts:

| Concept | Variants Found |
|---------|---------------|
| Number of splits | `n.splits`, `fs.splits`, `n.splits.screen` |
| HR threshold | `hr.threshold`, `hr.consistency` |
| Subgroup focus | `sg_focus`, `sg.focus` |
| Verbosity | `details`, `verbose` |

**Recommendation:** For new/internal code, standardize on snake_case. For existing exported parameters, maintain backward compatibility. Document the mapping in a package-level comment.

---

### L2. `stringr` Dependency in Cross-Validation

**File:** `R/forestsearch_cross-validation.R`

Uses `stringr::str_length()` and `stringr::str_sub()` in `find_covariate_any_match()`:

```r
lc <- stringr::str_length(confs2)
ctoget <- stringr::str_sub(sg_target, 2, max(lc))
```

These have direct base R equivalents: `nchar()` and `substr()`. Removing the `stringr` dependency simplifies installation.

---

### L3. `_pkgdown.yml` Placeholder URL

**Line 13:**
```yaml
url: https://github.com/yourusername/ForestSearch
```

Update to the actual GitHub URL before deploying the site.

---

## Architectural Strengths (Verified in KB)

These are confirmed as well-designed patterns worth maintaining:

1. **`evaluate_comparison()` + `get_dfpred()`**: Clean operator-dispatch pattern replacing `eval(parse())`. Supports `<=`, `>=`, `!=`, `==`, `<`, `>` with bare column name fallback and negation via `!{...}`.

2. **`dummy_encode()` consolidation**: Single implementation with `vapply()` type-safety, `drop = FALSE` for edge cases, and backward-compatible `dummy`/`dummy2` aliases via `@rdname`.

3. **GRF pipeline**: `create_grf_config()` → `grf.subg.harm.survival()` → `grf.estimates.out()` with clean config/result separation and `create_success_result()`/`create_null_result()` factory functions.

4. **`filter_call_args()`**: Elegant utility reducing `do.call()` boilerplate throughout the package by matching argument lists to function signatures.

5. **Bootstrap infrastructure**: `get_bootstrap_exports()` auto-discovers exports, `resolve_bootstrap_parallel_args()` provides clean fallback chain, `BOOTSTRAP_REQUIRED_PACKAGES` documents dependencies explicitly.

6. **S3 method ecosystem**: Comprehensive `print`/`plot` methods for 6+ object types (`forestsearch`, `fs_kfold`, `fs_tenfold`, `fs_forestplot`, `fs_forest_theme`, `fs_sg_plot`) with helpful usage instructions in print output.

7. **Publication-ready output**: `plot_subgroup_results_forestplot()` with `create_forest_theme()` customization, `sens_text()` for CV metrics, and `create_subgroup_summary_df()` for tabular summaries.

---

## Summary Table

| Priority | # | Key Items | Effort |
|----------|---|-----------|--------|
| **Critical** | 3 | Duplicate S3 methods, remaining `eval(parse())` ×3, `_pkgdown.yml` encoding | Low–Med |
| **High** | 4 | Dead `BOOTSTRAP_REQUIRED_FUNCTIONS`, commented-out code, `.r` extension, duplicate imports | Low |
| **Medium** | 5 | `evaluate_cuts_once()` eval/parse, duplicate `ensure_packages`, globals audit, NULL-check helper, timing helper | Low–Med |
| **Low** | 3 | Parameter naming, `stringr` dependency, placeholder URL | Low |

**Total: 15 findings.** Critical + High items estimated at ~2–3 hours focused work.
