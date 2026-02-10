# ForestSearch Package — Fresh Independent Code Review

**Date:** February 9, 2026
**Scope:** Complete codebase as attached to project knowledge base
**Reviewer approach:** All findings derived solely from current source files

---

## Executive Summary

ForestSearch is a well-architected R package for exploratory subgroup identification in survival-endpoint clinical trials. The codebase demonstrates strong domain expertise and a thoughtful modular design spanning ~16 functional areas. Key strengths include comprehensive roxygen2 documentation, clean separation of bootstrap/CV/consistency evaluation logic, a robust `safe_eval_expr()` pattern that replaces many `eval(parse())` calls, and well-designed S3 methods in a canonical `forestsearch_methods.R` file.

This review identifies **28 findings** organized by priority. The most impactful issues center on: duplicate function definitions across files (creating load-order bugs), residual `eval(parse())` usage, over-exported internal helpers, inconsistent file naming, and a brittle manually-curated function list for parallel workers.

### Severity Counts

| Priority | Count | CRAN Risk |
|----------|-------|-----------|
| Critical | 4 | Would likely trigger NOTE or rejection |
| High | 7 | Best-practice violations affecting reliability |
| Medium | 10 | Code quality, consistency, maintainability |
| Low | 7 | Style, DRY, future-proofing |

---

## Critical Priority (CRAN Blockers / Correctness)

### C1. Duplicate `forestsearch_KfoldOut()` Definition

**Files:** `R/forestsearch_cross-validation.R` AND `R/crossvalidation_helpers.R`

Both files export `forestsearch_KfoldOut()` with the same signature. Only one will be registered depending on file load order (alphabetical in R). The two implementations diverge:

- **`crossvalidation_helpers.R`** version has a hard `requireNamespace("weightedsurv")` dependency
- **`forestsearch_cross-validation.R`** version does NOT require `weightedsurv`
- The older version uses raw `subset(df_CV, cvindex==ks)` patterns; the newer version uses cleaner data.table operations

**Impact:** Silent correctness bug — which version runs depends on file collation order. The `weightedsurv` dependency may cause installation failures if that package isn't available.

**Fix:**

```r
# 1. Delete the definition from R/crossvalidation_helpers.R
# 2. Keep only the R/forestsearch_cross-validation.R version
# 3. If weightedsurv functionality is needed, gate it:
if (requireNamespace("weightedsurv", quietly = TRUE)) {
  # weightedsurv-specific logic
}
```

### C2. Duplicate `CV_sgs()` Definition

**Files:** `R/forestsearch_cross-validation.R` (contains `find_covariate_any_match()` — a cleaner refactored version) AND `R/crossvalidation_helpers.R` (contains `CV_sgs()` with inline covariate matching logic)

The `CV_sgs()` function in `crossvalidation_helpers.R` contains extensive inline brace-parsing and `charmatch()` logic that duplicates the cleaner `find_covariate_any_match()` in `forestsearch_cross-validation.R`. The inline version has fragile variable scoping (`rm("bb1","bb2","cfs")`) and repeated copy-paste blocks for `sg1a` vs `sg2a`.

**Fix:** Consolidate to one file. The `forestsearch_cross-validation.R` version with `find_covariate_any_match()` as a helper is cleaner — adopt it and delete the `crossvalidation_helpers.R` copy.

### C3. Residual `eval(parse())` in `evaluate_cuts_once()` and `get_FSdata()`

**File:** `R/get_FSdata_helpers.R`, `R/get_FSdata_refactored.r`

While `get_dfpred()` and `plot_subgroup_results_forestplot.R` have been successfully migrated to `safe_eval_expr()` / `evaluate_comparison()`, the cut-caching function still uses bare `eval(parse())`:

```r
# In evaluate_cuts_once():
result <- eval(parse(text = thiscut), envir = df)

# In get_FSdata_refactored.r (fallback path):
result <- eval(parse(text = thiscut), envir = df.FS)
```

These cut expressions are internally generated (e.g., `"age <= 50"`) so the security risk is low, but CRAN reviewers flag `eval(parse())` patterns. The `safe_eval_expr()` helper already exists and provides sandboxed evaluation via `list2env(as.list(df), parent = baseenv())`.

**Fix:** Replace both with `safe_eval_expr()`:

```r
# In evaluate_cuts_once():
result <- safe_eval_expr(df, thiscut)
if (is.null(result)) {
  has_error[i] <- TRUE
  is_valid[i] <- FALSE
  next
}
```

### C4. File Extension Inconsistency: `.r` vs `.R`

**Files:** `R/get_FSdata_refactored.r`, `R/bootstrap_dofuture_main.r`, `R/subgroup_consistency_main.r`, `R/bootstrap_analysis_dofuture.R`

Three files use lowercase `.r` while the rest use uppercase `.R`. Per tidyverse style guide and CRAN conventions, `.R` is standard. Some build tools on case-sensitive filesystems may handle these differently.

Additionally, `_refactored` suffixes in filenames are internal version markers that should not appear in a released package.

**Fix:**
- `get_FSdata_refactored.r` → `get_fsdata.R`
- `bootstrap_dofuture_main.r` → `bootstrap_dofuture_main.R`
- `subgroup_consistency_main.r` → `subgroup_consistency_main.R`

---

## High Priority (Reliability / Best Practice)

### H1. `stop_threshold` Documented Twice in `forestsearch()` Signature

**File:** `R/forestsearch_main.R`

The `@param stop_threshold` roxygen tag appears **twice** in the documentation block:

```r
#' @param stop_threshold Numeric. Early stopping threshold for consistency
#'   evaluation. When a candidate subgroup's estimated consistency probability
#'   exceeds this threshold, evaluation stops early. Default 0.95.
...
#' @param stop_threshold Numeric. Early stopping threshold for consistency
#'   evaluation. ...
```

This will cause a roxygen2 warning during `devtools::document()` and may result in garbled help pages.

**Fix:** Remove the duplicate `@param stop_threshold` block.

### H2. Over-Exported Helper Functions

Multiple small utility functions are `@export`ed that should be `@keywords internal` or `@noRd`:

| Function | File | Why internal? |
|----------|------|---------------|
| `one.zero()` | `get_FSdata_helpers.R` | Single-line binary flip |
| `ztrail()` | `get_FSdata_helpers.R` | Binary trailing-zeros counter |
| `acm.disjctif()` | `get_FSdata_helpers.R` | Low-level encoding detail |
| `is_flag_continuous()` | `get_FSdata_helpers.R` | Cut-expression classifier |
| `is_flag_drop()` | `get_FSdata_helpers.R` | Cut-expression validator |
| `get_cut_name()` | `get_FSdata_helpers.R` | Cut-expression parser |
| `evaluate_cuts_once()` | `get_FSdata_helpers.R` | Internal caching |
| `process_conf_force_expr()` | `get_FSdata_helpers.R` | Expression pre-processor |
| `filter_by_lassokeep()` | `get_FSdata_helpers.R` | Simple vector filter |
| `run_single_consistency_split()` | `subgroup_consistency_helpers.R` | Per-split worker |

Exporting these inflates the public API, creates backward-compatibility obligations, and clutters the namespace. Users should never need to call `ztrail()` or `one.zero()` directly.

**Fix:** Change `@export` → `@keywords internal` for all pure-internal helpers. For those called by parallel workers, they remain available in the package namespace via `:::` or through proper `.export` lists in `foreach`.

### H3. `BOOTSTRAP_REQUIRED_FUNCTIONS` — Brittle Manual Curation

**File:** `R/bootstrap_dofuture_setup_helpers.R`

A manually maintained list of 60+ function names is used to ensure parallel workers have access to needed functions. This includes functions that may not exist (`clean_data`, `prepare_data`, `run_bootstrap`, `run_grf`, `evaluate_subgroups`, `summarize_results`, `policy_tree`) and misses functions that DO exist.

**Impact:** Adding any new function without updating this list causes silent bootstrap failures. Conversely, listing nonexistent functions causes warnings or errors during worker initialization.

**Fix — Option A (preferred):** If using `future` with `plan(callr)` or `plan(multisession)`, the entire installed package is available. Remove the manual list and rely on package namespace:

```r
# Workers automatically have access to all exported functions
# No manual list needed when package is installed
```

**Fix — Option B:** Auto-discover exports:

```r
get_fs_exports <- function() {
  getNamespaceExports("ForestSearch")
}
```

### H4. `crossvalidation_helpers.R` — Legacy Code with Style Debt

**File:** `R/crossvalidation_helpers.R`

This file contains the older implementation of `CV_sgs()` and `forestsearch_KfoldOut()` with significant style issues:

- Inconsistent spacing (`if(dda ==1)` vs `if (dda == 1)`)
- `rm()` calls inside functions to clean up local variables (unnecessary; they're GC'd on exit)
- Repeated copy-paste blocks for processing `sg1a` and `sg2a`
- No input validation
- Uses `subset()` in programmatic context (warned against in `?subset` documentation)

Since cleaner versions of these functions exist in `forestsearch_cross-validation.R`, this file is a consolidation target.

**Fix:** Delete `R/crossvalidation_helpers.R` entirely after verifying `forestsearch_cross-validation.R` covers all functionality. If any logic is unique to the old file, migrate it first.

### H5. `plot.forestsearch()` — Missing from `forestsearch_methods.R`

**File:** `R/forestsearch_methods.R` header comments reference it, and `_pkgdown.yml` lists `plot.forestsearch` in the Core Algorithm reference section, but the actual method is not defined in the canonical methods file. If it exists elsewhere as a stub (the previous review noted a no-op version), it should either be:

1. Implemented meaningfully in `forestsearch_methods.R`, or
2. Removed from `_pkgdown.yml` and any `@export` tags

An exported S3 method that silently does nothing will confuse users and draw CRAN reviewer questions.

### H6. `lasso_selection()` Uses `set.seed()` at Function Level

**File:** `R/get_FSdata_helpers.R`

```r
lasso_selection <- function(df, confounders.name, outcome.name, event.name,
                            seedit = 8316951) {
  set.seed(seedit)
  ...
}
```

Calling `set.seed()` inside a package function is discouraged because it silently resets the user's RNG state. CRAN reviewers specifically look for this. The `glmnet::cv.glmnet()` call is the only stochastic element.

**Fix:** Use `withr::with_seed()` or document the side effect:

```r
lasso_selection <- function(..., seedit = 8316951) {
  withr::local_seed(seedit)
  # ... rest of function
}
```

If `withr` is not a dependency, save and restore the seed:

```r
old_seed <- .Random.seed
on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
set.seed(seedit)
```

### H7. `globals.R` — Overly Broad `globalVariables()` Declarations

**File:** `R/globals.R`

The file declares 100+ global variables including extremely generic names (`y`, `id`, `age`, `size`, `desc`, `z`, `E`, `N`, `K`, `Y`). This suppresses R CMD check NOTEs but masks genuine typos — if someone mistypes `hr_individua` instead of `hr_individual`, the check won't catch it.

**Fix (incremental):**
1. Group declarations by source file with comments (partially done)
2. Replace data.table NSE with `.SD`/`.SDcols` patterns where feasible
3. Use `ggplot2::aes(.data$column)` instead of bare `aes(column)`
4. Consider removing declarations for variables that are only used in examples/tests

---

## Medium Priority (Code Quality / Consistency)

### M1. Parameter Naming Inconsistency

The codebase mixes dot-separated and snake_case naming:

| Concept | Variants Found |
|---------|---------------|
| Number of splits | `n.splits`, `fs.splits`, `n.splits.screen` |
| HR threshold | `hr.threshold`, `hr.consistency`, `hr.cons` |
| Subgroup focus | `sg_focus` (✓ consistent) |
| Outcome variable | `outcome.name` (dot) |
| Treat variable | `treat.name` (dot) |
| Minimum size | `n.min` (dot) |
| Two-stage | `use_twostage` (snake) |
| Show details | `details`, `verbose` |

The existing dot-separated parameters (`outcome.name`, `hr.threshold`, etc.) are part of the exported API and should be preserved for backward compatibility. New parameters correctly use snake_case (`use_twostage`, `sg_focus`).

**Recommendation:** Document the convention: "Legacy parameters use dot-separation for backward compatibility; new parameters use snake_case." Avoid introducing new dot-separated parameters.

### M2. `crossvalidation_helpers.R` Covariate Matching — Fragile String Parsing

**File:** `R/crossvalidation_helpers.R` (and partially `R/forestsearch_cross-validation.R`)

The `CV_sgs()` function uses `charmatch()` for prefix matching of brace-wrapped confounder names. This is fragile for names that are prefixes of each other (e.g., `z1` vs `z11`). The code has explicit workarounds:

```r
# Find exact match when multiple confounders match (e.g., "z1" and "z11")
if (length(index_name) > 1) {
  confs2 <- confs[index_name]
  lc <- stringr::str_length(confs2)
  ...
}
```

**Fix:** Replace `charmatch()` with exact matching using word boundaries:

```r
# Use exact matching with delimiters
match_confounder <- function(sg_label, confs) {
  for (i in seq_along(confs)) {
    pattern <- paste0("\\b", confs[i], "\\b")
    if (grepl(pattern, sg_label)) return(i)
  }
  NA_integer_
}
```

### M3. Redundant `dummy()` / `dummy2()` Aliases

**File:** `R/get_FSdata_helpers.R`

Both `dummy` and `dummy2` are now aliases for `dummy_encode`:

```r
#' @rdname dummy_encode
#' @export
dummy <- dummy_encode

#' @rdname dummy_encode
#' @export
dummy2 <- dummy_encode
```

This is a good consolidation, but both aliases are still exported AND listed separately in `_pkgdown.yml`. For CRAN, having three exports that are identical functions is unusual.

**Fix:** Keep `dummy_encode` as the canonical export. Deprecate `dummy` and `dummy2`:

```r
#' @rdname dummy_encode
#' @export
dummy <- function(df) {
  .Deprecated("dummy_encode")
  dummy_encode(df)
}
```

Or, if backward compatibility is essential, at minimum update `_pkgdown.yml` to list only `dummy_encode`.

### M4. Inconsistent Error Handling Patterns

Three different error-handling styles appear across the codebase:

```r
# Pattern 1: try() with inherits check
run_bootstrap <- try(do.call(forestsearch, args_FS_boot), TRUE)
if (inherits(run_bootstrap, "try-error")) { ... }

# Pattern 2: tryCatch() with explicit handler
result <- tryCatch(
  do.call(subgroup.search, search_args),
  error = function(e) { warning(...); return(NULL) }
)

# Pattern 3: tryCatch() + suppressWarnings
resCV <- suppressWarnings({foreach::foreach(...) %dofuture% { ... }})
```

**Recommendation:** Standardize on `tryCatch()` with explicit handlers. The `try()` pattern is less informative and requires a separate `inherits()` check. Reserve `suppressWarnings()` for known benign warnings only, with a comment explaining what's being suppressed.

### M5. `safe_eval_expr()` vs `evaluate_comparison()` — Two Approaches to Expression Evaluation

**Files:** `R/forestsearch_helpers.R`

The codebase now has two expression evaluation mechanisms:

1. **`safe_eval_expr()`** — sandboxed `eval(parse())` with `baseenv()` parent. Used for arbitrary expressions like `"BM > 1 & tmrsize > 19"`.
2. **`evaluate_comparison()`** — operator dispatch without `eval(parse())`. Used for single comparisons like `"er <= 0"`.

Both are good solutions. However, `get_dfpred()` uses `evaluate_comparison()` per-factor and then `Reduce(&, ...)` to combine, while `create_reference_subgroup_columns()` uses `safe_eval_expr()` for compound expressions.

**Recommendation:** Document the intended division of labor:
- `evaluate_comparison()`: Single variable comparisons (no `&` or `|`)
- `safe_eval_expr()`: Compound expressions
- `evaluate_cuts_once()`: Should be migrated to use `safe_eval_expr()` (see C3)

### M6. `forestsearch_main.R` — Very Long Function Body

The `forestsearch()` function spans what appears to be 10+ sections with inline processing. While the section headers and comments are excellent, the function is doing too much:

1. Capture arguments (Section 1)
2. Validate inputs (Section 2)
3. Data preparation via `get_FSdata` (Section 3-4)
4. GRF variable importance screening (Section 5)
5. Subgroup search via `subgroup.search` (Section 6)
6. Check for valid subgroups (Section 8)
7. Consistency evaluation (Section 9)
8. Result assembly and prediction (Section 10+)

Sections 5 (GRF screening) and 9 (consistency evaluation) are particularly amenable to extraction as helper functions.

**Fix (incremental):** Extract the GRF screening logic into `screen_variables_grf()`:

```r
screen_variables_grf <- function(df, FSconfounders.name, Y, Treat, Event,
                                  is.RCT, vi.grf.min, max_n_confounders, ...) {
  X <- as.matrix(df[, FSconfounders.name])
  X <- apply(X, 2, as.numeric)
  tau.rmst <- min(c(max(Y[Treat == 1 & Event == 1]),
                     max(Y[Treat == 0 & Event == 1])))
  cs.forest <- grf::causal_survival_forest(...)
  # ... screening logic ...
  conf.screen
}
```

### M7. `_pkgdown.yml` — UTF-8 Encoding Issue in Description

**File:** `_pkgdown.yml`, line 86

```yaml
    desc: >
      Bootstrap methods for bias-corrected hazard ratio estimation using
      infinitesimal jackknife variance estimation (LeÃ³n et al., 2024).
```

The `ó` in "León" has been corrupted to `Ã³` — a classic UTF-8/Latin-1 mismatch. This will render incorrectly on the pkgdown site.

**Fix:** Ensure the file is saved as UTF-8 and replace with the correct character:

```yaml
      infinitesimal jackknife variance estimation (León et al., 2024).
```

### M8. `run_grf_analysis()` — Fragile Function Lookup

**File:** `R/oc_analyses_gbsg_refactored.R`

```r
grf_fun <- tryCatch({
  get("grf.subg.harm.survival", mode = "function", envir = parent.frame())
}, error = function(e) {
  tryCatch({
    get("grf.subg.harm.survival", mode = "function", envir = globalenv())
  }, error = function(e2) NULL)
})
```

This searches parent frames and the global environment for the function — a fragile pattern that breaks in non-interactive contexts and parallel workers. Since `grf.subg.harm.survival` is exported from this package, it should be called directly.

**Fix:**

```r
# Direct call — the function is in our own package namespace
grf_result <- tryCatch(
  do.call(grf.subg.harm.survival, grf_args),
  error = function(e) {
    warning(sprintf("%s analysis failed: %s", analysis_label, e$message))
    NULL
  }
)
```

### M9. `create_summary_table.R` — No Early Return for Missing `gt`

**File:** `R/create_summary_table.R`

The function has a fallback at the very end:

```r
} else {
  message("Note: gt package not available. Returning data frame instead.")
```

But all the gt-specific styling code (30+ lines) runs before hitting this branch. If `gt` is unavailable, none of that code is reachable. A better pattern:

```r
create_summary_table <- function(..., use_gt = TRUE) {
  # Build the data first
  results <- build_summary_data(...)

  # Format output
  if (use_gt && requireNamespace("gt", quietly = TRUE)) {
    return(format_as_gt(results, ...))
  }

  message("gt package not available. Returning data frame.")
  results
}
```

### M10. Timing Pattern Repetition

The pattern `t_start <- proc.time()[3]; ...; t_min <- (proc.time()[3] - t_start) / 60` appears in at least 5 locations across `forestsearch_main.R`, `bootstrap_dofuture_main.r`, `forestsearch_cross-validation.R`, and `oc_analyses_gbsg_refactored.R`.

**Fix:** Create a simple timing utility:

```r
#' @keywords internal
fs_timer <- function() {
  t0 <- proc.time()[3]
  function() (proc.time()[3] - t0) / 60
}

# Usage:
elapsed <- fs_timer()
# ... computation ...
cat("Time:", round(elapsed(), 2), "minutes\n")
```

---

## Low Priority (Style / DRY / Future-Proofing)

### L1. `_pkgdown.yml` URL Placeholder

**File:** `_pkgdown.yml`, lines 13, 63

```yaml
url: https://github.com/yourusername/ForestSearch
href: https://github.com/yourusername/ForestSearch
```

These placeholders should be updated to the actual repository URL before deployment.

### L2. `print_cv_params()` and Bootstrap Parameter Printing — Repeated Patterns

**Files:** `R/forestsearch_cross-validation.R`, `R/bootstrap_dofuture_main.r`

The parameter-printing blocks in CV and bootstrap are nearly identical:

```r
# CV version:
cat("  - sg_focus:", cv_args$sg_focus, "\n")
cat("  - maxk:", cv_args$maxk, "\n")
...

# Bootstrap version:
cat("  - sg_focus:", args_forestsearch_call$sg_focus, "\n")
cat("  - maxk:", args_forestsearch_call$maxk, "\n")
...
```

**Fix:** Extract a shared `print_fs_params(args, context = "bootstrap")` utility that both can call.

### L3. `calculate_skewness()` — Custom Implementation

If `e1071` or `moments` is already in Suggests/Imports, use their `skewness()`. If not, the custom implementation is fine but should be `@keywords internal` rather than exported (if it currently is).

### L4. `BOOTSTRAP_REQUIRED_PACKAGES` — Top-Level Constant

**File:** `R/bootstrap_dofuture_setup_helpers.R`

```r
BOOTSTRAP_REQUIRED_PACKAGES <- c("data.table", "foreach", ...)
```

Top-level code in R package files executes at load time. While this particular constant is harmless, CRAN prefers that package files contain only function definitions and `globalVariables()` calls. Wrap it:

```r
bootstrap_required_packages <- function() {
  c("data.table", "foreach", "doFuture", "doRNG", "survival")
}
```

### L5. Comment Artifacts — `# DEAD` Code Markers

**File:** `R/plot_subgroup_results_forestplot.R`

Several commented-out blocks are marked with `# DEAD`:

```r
# DEAD
# df_sg <- tryCatch(
#   subset(df_analysis, eval(parse(text = sg$subset_expr))),
#   ...
# )

df_sg <- safe_subset(df_analysis, sg$subset_expr)
```

These are helpful during development but should be removed before CRAN submission. The git history preserves the old code.

### L6. `simulate_from_gbsg_dgm()` Function Name vs `simulate_from_dgm()` in pkgdown

**File:** `R/sim_aft_gbsg_refactored.R` exports `simulate_from_gbsg_dgm()`, but `_pkgdown.yml` lists `simulate_from_dgm`. If both exist, they should be cross-referenced. If only one exists, the pkgdown config needs updating.

### L7. Verbose Messaging — Inconsistent Use of `message()` vs `cat()`

The codebase mixes `cat()` and `message()` for user-facing output:

- `forestsearch_main.R` uses `cat()` with `\n`
- `bootstrap_dofuture_setup_helpers.R` uses `message()`
- `oc_analyses_gbsg_refactored.R` uses `message()` with `sprintf()`

Per R best practices, `message()` is preferred because:
- It respects `suppressMessages()`
- It goes to `stderr` (not captured by `capture.output()`)
- It's the standard for informational output

**Fix (incremental):** Standardize new code on `message()`. Migrate `cat()` calls over time.

---

## Architecture Assessment

### What's Working Well

1. **`forestsearch_methods.R`** — Clean canonical S3 methods with `.fs_get()` and `.fs_sg_labels()` helpers for robust slot access across object versions.

2. **`evaluate_comparison()`** — Excellent CRAN-safe replacement for `eval(parse())` with operator dispatch. The longest-first operator matching prevents partial-match bugs.

3. **`safe_eval_expr()` / `safe_subset()`** — Good sandboxed evaluation pattern using `list2env(as.list(df), parent = baseenv())`.

4. **Parallel architecture** — Clear separation of concerns: outer parallelization in bootstrap/CV, forced sequential in inner `forestsearch()` calls. Prevents nested parallelism issues.

5. **Two-stage consistency algorithm** — Well-designed adaptive algorithm with proper early stopping, Wilson CI for proportion estimation, and configurable batch sizes.

6. **`filter_call_args()`** — Elegant utility for managing the complex argument-passing patterns inherent in the bootstrap/CV workflows.

7. **`dummy_encode()` consolidation** — Clean refactor with `vapply()` type-checking, `drop = FALSE` edge-case handling, and backward-compatible aliases.

8. **Comprehensive `@param` / `@return` / `@examples` documentation** — The roxygen2 blocks are thorough with contextual `\dontrun{}` examples.

### Structural Recommendations

1. **File consolidation priority:** Delete `R/crossvalidation_helpers.R` after verifying coverage in `R/forestsearch_cross-validation.R`.

2. **Export audit:** Run `devtools::check()` and count the number of documented exports. Consider whether users need 119 exported functions or whether 30-40 well-designed entry points would suffice.

3. **Testing infrastructure:** No test files were visible in the KB. Adding `testthat` tests for at minimum `get_dfpred()`, `evaluate_comparison()`, `dummy_encode()`, and `forestsearch_KfoldOut()` would catch regressions during the consolidation work.

---

## Summary Action Plan

| Phase | Actions | Files Affected |
|-------|---------|----------------|
| **Phase 1: CRAN blockers** | Delete `crossvalidation_helpers.R`; fix `.r` → `.R`; replace `eval(parse())` in `evaluate_cuts_once()`; remove duplicate `@param` | 4 files |
| **Phase 2: Cleanup** | Mark internal helpers with `@keywords internal`; remove `# DEAD` comments; fix UTF-8 in `_pkgdown.yml`; update placeholder URLs | 6 files |
| **Phase 3: Consolidation** | Replace `BOOTSTRAP_REQUIRED_FUNCTIONS` manual list; standardize on `message()`; extract timing utility; deprecate `dummy()`/`dummy2()` | 4 files |
| **Phase 4: Documentation** | Verify `plot.forestsearch()` exists or remove from pkgdown; verify `simulate_from_dgm` vs `simulate_from_gbsg_dgm`; parameter naming convention note | 2 files |
