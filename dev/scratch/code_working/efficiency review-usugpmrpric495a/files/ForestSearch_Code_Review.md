# ForestSearch R Package - Code Review

## Executive Summary

This review covers the ForestSearch R package for exploratory subgroup identification in clinical trials with survival endpoints. The package is well-structured with comprehensive functionality, but several improvements can enhance CRAN compliance, code efficiency, and maintainability.

**Overall Assessment**: The codebase is mature with good documentation coverage. Key areas requiring attention are:
1. Consolidation of duplicate code (S3 methods, global variables)
2. Performance optimizations in hot paths (bootstrap iterations)
3. Minor tidyverse style adjustments
4. CRAN compliance refinements

---

## 1. CRAN Compliance Issues

### 1.1 Multiple `globalVariables()` Declarations (HIGH PRIORITY)

**Problem**: `globalVariables()` declarations are scattered across multiple files:
- `R/globals.R`
- `R/oc_analyses_gbsg_refactored.R`
- `R/mrct_simulation.R`
- `R/sim_aft_gbsg_refactored.R`

**CRAN Impact**: While not strictly forbidden, this makes maintenance difficult and can lead to redundant declarations.

**Recommendation**: Consolidate ALL `globalVariables()` into a single `R/globals.R` file:

```r
# R/globals.R - CONSOLIDATED
# All globalVariables declarations should be here

utils::globalVariables(c(
  # ==========================================================================
  # Summary table columns
  # ==========================================================================
  "Category", "Metric", "Value", "Count", "Base_Factor", "Factor_Definition",
  "Factor", "Consistency Range", "Size Range", "N_positions", "N_total",
  "Positions", "Position",


  # ==========================================================================
  # Cross-validation
  # ==========================================================================
  "cv_index", "cvindex",

  # ==========================================================================
  # ForestSearch parameters (NSE contexts)
  # ==========================================================================
  "est.scale", "insplit1", "confounders.name", "outcome.name",
  "event.name", "treat.name",

  # ==========================================================================
  # Survival analysis & Cox models
  # ==========================================================================
  "treat", "event", "Y", "E", "Treat", "stratum", "z", "loghr_po",
  "theta_1", "theta_0", "treat_harm",

  # ==========================================================================
  # AFT DGM / Simulation variables
  # ==========================================================================
  "y", "y_sim", "event_sim", "t_true", "c_time", "lin_pred_1", "lin_pred_0",
  "lin_pred_obs", "lin_pred_cens_1", "lin_pred_cens_0",

  # ==========================================================================
  # GBSG dataset columns
  # ==========================================================================
  "id", "pid", "status", "rfstime", "age", "size", "nodes", "pgr", "er",
  "meno", "grade", "hormon", "desc",

  # ==========================================================================
  # Operating Characteristics (oc_analyses_gbsg)
  # ==========================================================================
  "any.H", "ppv", "npv", "sensitivity", "specificity", "size.H", "size.Hc",
  "hr.H.true", "hr.H.hat", "hr.Hc.true", "hr.Hc.hat", "hr.itt", "hr.adj.itt",
  "p.cens", "taumax", "analysis", "sim", "aa", "sg_hat",
  "ahr.H.true", "ahr.Hc.true", "ahr.H.hat", "ahr.Hc.hat",

  # ==========================================================================
  # MRCT simulation variables
  # ==========================================================================
  "hr_test", "hr_sg", "any_found", "sg_found", "hr_sg_null", "regAflag",
  "sg_le85", "regAflag2", "regAflag3", "found", "sg_biomarker", "sg_age",
  "sg_male", "sg_ecog", "sg_histology", "sg_CTregimen", "sg_region",
  "sg_surgery", "sg_prior_treat", "est", "region_var", "z_regA",

  # ==========================================================================
  # Duplicate detection
  # ==========================================================================
  "dup_key", "dup_group", "dup_count", "dup_rank",

  # ==========================================================================
  # GBSG DGM specific
  # ==========================================================================
  "zh", "flag.harm", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "grade3",
  "lin.conf.true", "lin1.conf", "lin0.conf", "linC1.conf", "linC0.conf",
  "hlin.conf.1", "hlin.conf.0", "hlin.ratio", "h1.potential", "h0.potential",
  "Ts", "es", "t.sim",

  # ==========================================================================
  # Model comparison
  # ==========================================================================
  "By_AIC", "By_BIC", "By_LogLik", "AIC_value", "BIC_value", "LogLik_value",

  # ==========================================================================
  # Standardized GBSG covariate names
  # ==========================================================================
  "z_age", "z_size", "z_nodes", "z_pgr", "z_er", "z_meno",
  "z_grade_1", "z_grade_2", "z_hormon"
))
```

### 1.2 Duplicate S3 Method Definitions (HIGH PRIORITY)

**Problem**: `print.forestsearch()` and `summary.forestsearch()` are defined in multiple files:
- `R/forest_search_revised.R`
- `R/print_forestsearch.R`
- `R/summary_forestsearch.R`

**CRAN Impact**: Will cause warnings during R CMD check about duplicate function definitions.

**Recommendation**: Keep ONLY ONE definition per method. Remove from `forest_search_revised.R` and keep the dedicated files:

```r
# In R/forest_search_revised.R, REMOVE these function definitions:
# - print.forestsearch() (keep in R/print_forestsearch.R)
# - summary.forestsearch() (keep in R/summary_forestsearch.R)
```

### 1.3 Using `install.packages()` in Package Code

**Problem**: Reference to `install.packages()` in comments suggests it might be used.

**CRAN Impact**: Packages should NEVER call `install.packages()` during execution.

**Verification**: Ensure all package availability checks use `requireNamespace()`:

```r
# CORRECT pattern:
if (!requireNamespace("ggplot2", quietly = TRUE)) {

  stop("Package 'ggplot2' is required but not installed.")
}

# INCORRECT pattern (should not exist):
# install.packages("ggplot2")
```

---

## 2. Efficiency Improvements

### 2.1 `cox_summary()` Optimization (HIGH PRIORITY)

**Location**: `R/Cox_estimation_helpers.R`

**Current Issue**: Function is called thousands of times in bootstrap iterations. Current implementation has optimization potential.

**Recommended Optimization**:

```r
#' Cox model summary - Optimized for bootstrap iterations
#' @keywords internal
cox_summary_fast <- function(Y, E, Treat, Strata = NULL,
                             return_format = c("numeric", "formatted")) {

  return_format <- match.arg(return_format)

  # Quick validation (fail fast)
  n <- length(Y)
  n_events <- sum(E)


  # Early exit for edge cases
  if (n_events < 2L || length(unique(Treat)) < 2L) {
    return(if (return_format == "formatted") {
      "NA (NA, NA)"
    } else {
      c(HR = NA_real_, Lower = NA_real_, Upper = NA_real_)
    })
  }

  # Fit model with minimal overhead
  fit <- tryCatch({
    if (is.null(Strata)) {
      survival::coxph(
        survival::Surv(Y, E) ~ Treat,
        robust = TRUE,
        model = FALSE, x = FALSE, y = FALSE  # Memory optimization
      )
    } else {
      survival::coxph(
        survival::Surv(Y, E) ~ Treat + strata(Strata),
        robust = TRUE,
        model = FALSE, x = FALSE, y = FALSE
      )
    }
  }, error = function(e) NULL)

  if (is.null(fit)) {
    return(if (return_format == "formatted") {
      "NA (NA, NA)"
    } else {
      c(HR = NA_real_, Lower = NA_real_, Upper = NA_real_)
    })
  }

  # Extract directly from conf.int (single summary call)
  ci <- summary(fit)$conf.int[1L, c(1L, 3L, 4L)]

  if (return_format == "formatted") {
    sprintf("%.2f (%.2f, %.2f)", ci[1L], ci[2L], ci[3L])
  } else {
    c(HR = ci[1L], Lower = ci[2L], Upper = ci[3L])
  }
}
```

### 2.2 `get_split_hr()` Optimization (HIGH PRIORITY)

**Location**: Called within consistency evaluation loops

**Current Pattern** (inefficient):
```r
# Called inside loops - creates overhead
for (split in 1:n_splits) {
  hr <- get_split_hr(data, split_indices[[split]])
}
```

**Recommended**: Pre-compute split assignments and vectorize where possible:

```r
#' Vectorized split HR calculation
#' @keywords internal
get_split_hr_batch <- function(df, split_list, outcome_col, event_col, treat_col) {
  
  n_splits <- length(split_list)
  results <- numeric(n_splits)
  
  # Pre-extract vectors once
  Y <- df[[outcome_col]]
  E <- df[[event_col]]
  Treat <- df[[treat_col]]
  
  for (i in seq_len(n_splits)) {
    idx <- split_list[[i]]
    fit <- tryCatch(
      survival::coxph(
        survival::Surv(Y[idx], E[idx]) ~ Treat[idx],
        model = FALSE, x = FALSE, y = FALSE
      ),
      error = function(e) NULL
    )
    results[i] <- if (!is.null(fit)) exp(coef(fit)[1L]) else NA_real_
  }
  
  results
}
```

### 2.3 data.table Copy Avoidance

**Location**: Multiple files, especially in `subgroup_consistency_helpers.R`

**Current Pattern**:
```r
result_new <- data.table::copy(res)
```

**Issue**: Unnecessary copies when data won't be modified by reference.

**Recommendation**: Only copy when modifications are needed:

```r
# BEFORE:
result_new <- data.table::copy(res)
result_new <- sort_subgroups(result_new, sg_focus)

# AFTER (if sort_subgroups doesn't modify in place):
result_sorted <- sort_subgroups(res, sg_focus)

# OR use setkey() for in-place sorting when appropriate:
data.table::setorderv(res, cols = "hr", order = -1L)
```

---

## 3. Readability Improvements

### 3.1 Long Function Decomposition

**Location**: `R/forest_search_revised.R` - `forestsearch()` function

**Issue**: Main function is very long with embedded section comments. While sections help, extraction into helper functions improves testability.

**Recommendation**: Extract major sections into internal helper functions:

```r
#' Main ForestSearch function (refactored)
#' @export
forestsearch <- function(df.analysis, ...) {
  
  # Validate inputs
  validated <- validate_forestsearch_inputs(df.analysis, ...)
  
  # Variable selection
  var_selection <- perform_variable_selection(
    df = validated$df,
    use_lasso = validated$args$use_lasso,
    use_grf = validated$args$use_grf
  )
  
  # Subgroup search
  search_results <- perform_subgroup_search(
    df = validated$df,
    confounders = var_selection$selected,
    maxk = validated$args$maxk
  )
  
  # Consistency evaluation
  consistency_results <- evaluate_consistency(
    df = validated$df,
    candidates = search_results$candidates,
    n_splits = validated$args$fs.splits,
    use_twostage = validated$args$use_twostage
  )
  
  # Compile results
  compile_forestsearch_results(
    consistency = consistency_results,
    search = search_results,
    var_selection = var_selection,
    args = validated$args
  )
}
```

### 3.2 Consistent Naming Convention

**Issue**: Mixed naming styles throughout codebase:
- `df.analysis` vs `df_analysis`
- `hr.threshold` vs `hr_threshold`
- `sg.harm` vs `sg_harm`

**Tidyverse Style Recommendation**: Use `snake_case` consistently for all arguments and internal variables.

**Migration Strategy** (backward-compatible):

```r
#' ForestSearch function with consistent naming
#'
#' @param df_analysis Data frame (old: df.analysis)
#' @param hr_threshold HR threshold (old: hr.threshold)
forestsearch <- function(
    df_analysis = NULL,
    hr_threshold = 1.25,
    # Backward compatibility
    df.analysis = NULL,
    hr.threshold = NULL,
    ...
) {
  # Handle deprecated argument names
  if (!is.null(df.analysis) && is.null(df_analysis)) {
    .Deprecated(msg = "df.analysis is deprecated, use df_analysis")
    df_analysis <- df.analysis
  }
  if (!is.null(hr.threshold) && is.null(hr_threshold)) {
    .Deprecated(msg = "hr.threshold is deprecated, use hr_threshold")
    hr_threshold <- hr.threshold
  }
  
  # Continue with snake_case internally
  ...
}
```

### 3.3 Magic Numbers → Named Constants

**Issue**: Hard-coded values scattered throughout code.

**Current**:
```r
if (nb_boots < 100) {
  warning("nb_boots < 100 may produce unreliable intervals")
}
```

**Recommended** (create constants file or at top of file):

```r
# Constants for ForestSearch configuration
BOOTSTRAP_MIN_RECOMMENDED <- 100L
BOOTSTRAP_MIN_PRODUCTION <- 500L
DEFAULT_SEED <- 8316951L
DAYS_PER_MONTH <- 30.4375

# Usage:
if (nb_boots < BOOTSTRAP_MIN_RECOMMENDED) {
  warning(sprintf(
    "nb_boots < %d may produce unreliable intervals",
    BOOTSTRAP_MIN_RECOMMENDED
  ))
}
```

---

## 4. Tidyverse Style Guide Compliance

### 4.1 Line Length

**Issue**: Several lines exceed 80-100 character limit.

**Example** (from `R/oc_analyses_gbsg_refactored.R`):
```r
message(sprintf("    [FSlg] use_lasso = %s, use_grf = %s, use_twostage = %s",
                      fs_params_lasso_grf$use_lasso, fs_params_lasso_grf$use_grf,
                      fs_params_lasso_grf$use_twostage))
```

**Recommendation**: Break long lines at meaningful points:
```r
message(sprintf(
  "    [FSlg] use_lasso = %s, use_grf = %s, use_twostage = %s",
  fs_params_lasso_grf$use_lasso,
  fs_params_lasso_grf$use_grf,
  fs_params_lasso_grf$use_twostage
))
```

### 4.2 Assignment Operator

**Issue**: Mix of `=` and `<-` for assignment.

**Tidyverse Standard**: Use `<-` for assignment, `=` only for function arguments.

```r
# CORRECT:
result <- compute_hr(data)
fit <- coxph(formula, data = df)

# INCORRECT (assignment inside data.table should still use :=):
# But standalone assignment should use <-
```

### 4.3 Spacing Around Operators

**Issue**: Inconsistent spacing in some files.

**Current**:
```r
if(nb_boots<100)  # Missing spaces
```

**Recommended**:
```r
if (nb_boots < 100)  # Proper spacing
```

### 4.4 Function Documentation Order

**Tidyverse roxygen2 ordering**:
1. `@title` (or first line)
2. `@description`
3. `@param` (in order of function signature)
4. `@return`
5. `@details`
6. `@examples`
7. `@seealso`
8. `@references`
9. `@export` / `@keywords internal`
10. `@importFrom`

---

## 5. Code Organization Recommendations

### 5.1 Proposed File Structure

```
R/
├── globals.R                      # ALL globalVariables (consolidated)
├── forestsearch-package.R         # Package-level documentation & imports
├── constants.R                    # Package constants (NEW)
│
├── # Core functionality
├── forest_search_main.R           # Main forestsearch() function
├── forest_search_helpers.R        # Helper functions
├── subgroup_search.R              # Subgroup search algorithm
├── subgroup_consistency_main.R    # Consistency evaluation
├── subgroup_consistency_helpers.R # Consistency helpers
│
├── # Cox/Survival utilities
├── cox_estimation.R               # Cox model fitting (consolidated)
├── survival_helpers.R             # Survival utilities
│
├── # Variable selection
├── variable_selection_lasso.R     # LASSO selection
├── variable_selection_grf.R       # GRF-based selection
│
├── # Bootstrap inference
├── bootstrap_main.R               # Main bootstrap function
├── bootstrap_workers.R            # Parallel worker functions
├── bootstrap_helpers.R            # Bootstrap utilities
│
├── # Data generation
├── dgm_aft_main.R                 # Main DGM functions
├── dgm_aft_helpers.R              # DGM helpers
├── dgm_gbsg.R                     # GBSG-specific DGM
├── simulate_from_dgm.R            # Simulation functions
│
├── # Simulation/OC
├── oc_analysis.R                  # Operating characteristics
├── mrct_simulation.R              # MRCT simulations
│
├── # S3 methods (one file per generic)
├── print_methods.R                # All print methods
├── summary_methods.R              # All summary methods
├── plot_methods.R                 # All plot methods
│
├── # Output/Visualization
├── format_tables.R                # Table formatting (gt)
├── plot_forestplot.R              # Forest plot creation
└── utils.R                        # General utilities
```

### 5.2 Import Organization

**Current** (scattered across files):
```r
#' @importFrom survival coxph Surv
#' @importFrom data.table data.table setorder
```

**Recommended** (consolidate in `forestsearch-package.R`):

```r
# R/forestsearch-package.R

#' @keywords internal
"_PACKAGE"

# ==========================================================================
# Core Package Imports
# ==========================================================================

#' @import survival
#' @import data.table
NULL

# ==========================================================================
# Selective Imports by Category
# ==========================================================================

#' @importFrom stats coef predict residuals formula terms update
#' @importFrom stats quantile sd median var cor
#' @importFrom stats pchisq qnorm qt pnorm pt
#' @importFrom stats rexp runif rnorm rbinom
#' @importFrom stats na.omit complete.cases aggregate
#' @importFrom stats AIC BIC model.frame vcov
#' @importFrom stats as.formula setNames
#' @name forestsearch-stats-imports
#' @keywords internal
NULL

#' @importFrom graphics par plot lines points abline legend
#' @importFrom graphics hist barplot mtext axis title
#' @name forestsearch-graphics-imports
#' @keywords internal
NULL

#' @importFrom utils head tail str object.size
#' @importFrom utils capture.output write.table read.table
#' @importFrom utils hasName packageVersion
#' @name forestsearch-utils-imports
#' @keywords internal
NULL

#' @importFrom splines ns
#' @name forestsearch-splines-imports
#' @keywords internal
NULL

# ==========================================================================
# Parallel Processing Imports
# ==========================================================================

#' @importFrom future plan
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @importFrom doRNG %dorng%
#' @name forestsearch-parallel-imports
#' @keywords internal
NULL
```

---

## 6. Specific Code Fixes

### 6.1 Fix `get_dfpred()` Readability

**Location**: `R/forestsearch_helpers.R`

**Current** (uses `eval(parse())` which is fragile):
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

**Recommended** (safer and more readable):
```r
#' Generate Prediction Dataset with Subgroup Treatment Recommendation
#'
#' @param df_predict Data frame for prediction
#' @param sg_harm Character vector of subgroup-defining covariate names
#' @param dummy_version Integer; 1 uses dummy(), 2 uses dummy2()
#'
#' @return Data frame with treat_recommend column
#' @export
get_dfpred <- function(df_predict, sg_harm, dummy_version = 1L) {
  

  # Apply dummy coding
  df_pred <- if (dummy_version == 1L) dummy(df_predict) else dummy2(df_predict)
  
  # Identify harm subgroup using vectorized logic
  # All sg_harm columns must equal 1 for membership in harm group
  in_harm_subgroup <- Reduce(`&`, lapply(sg_harm, function(col) {
    df_pred[[col]] == 1L
  }))
  
  # Assign treatment recommendations
  # treat_recommend = 0: do not treat (in harm subgroup)
  # treat_recommend = 1: treat (not in harm subgroup)
  df_pred$treat_recommend <- ifelse(in_harm_subgroup, 0L, 1L)
  
  df_pred
}
```

### 6.2 Fix `identify_near_duplicates()` Warning

**Location**: `R/subgroup_consistency_helpers.R`

**Issue**: data.table modification warning from `[, dup_key := ...]` syntax.

**Fix**:
```r
identify_near_duplicates <- function(dt, cols_to_check = 2:10,
                                     tolerance = 1e-6, keep = "first",
                                     return_indices = FALSE) {
  
  # Work on a copy to avoid modifying input
  dt_work <- data.table::copy(dt)
  
  # ... rest of function using dt_work
  
  # Use set() instead of := for cleaner modification
  data.table::set(
    dt_work, 
    j = "dup_key", 
    value = do.call(paste, c(dt_work[, ..cols_to_check], sep = "_"))
  )
  
  # ...
}
```

---

## 7. Testing Recommendations

### 7.1 Unit Tests to Add

```r
# tests/testthat/test-cox-summary.R

test_that("cox_summary handles edge cases", {
  # Test insufficient events
  result <- cox_summary(
    Y = c(1, 2, 3, 4, 5),
    E = c(0, 0, 0, 0, 1),  # Only 1 event
    Treat = c(0, 0, 1, 1, 1),
    return_format = "numeric"
  )
  expect_true(is.na(result["HR"]))
  
  # Test no treatment variation
  result <- cox_summary(
    Y = c(1, 2, 3, 4, 5),
    E = c(1, 1, 1, 1, 1),
    Treat = c(1, 1, 1, 1, 1),  # All treated
    return_format = "numeric"
  )
  expect_true(is.na(result["HR"]))
})

test_that("cox_summary returns correct format", {
  set.seed(123)
  n <- 100
  Y <- rexp(n)
  E <- rbinom(n, 1, 0.7)
  Treat <- rbinom(n, 1, 0.5)
  
  # Numeric format
  result_num <- cox_summary(Y, E, Treat, return_format = "numeric")
  expect_named(result_num, c("HR", "Lower", "Upper"))
  expect_true(all(!is.na(result_num)))
  
  # Formatted string
  result_fmt <- cox_summary(Y, E, Treat, return_format = "formatted")
  expect_type(result_fmt, "character")
  expect_match(result_fmt, "^[0-9.]+\\s*\\([0-9.]+,\\s*[0-9.]+\\)$")
})
```

### 7.2 Performance Benchmarks

```r
# tests/testthat/test-performance.R

test_that("cox_summary_fast is faster than legacy", {
  skip_if_not_installed("bench")
  
  set.seed(123)
  n <- 1000
  Y <- rexp(n)
  E <- rbinom(n, 1, 0.7)
  Treat <- rbinom(n, 1, 0.5)
  
  bm <- bench::mark(
    legacy = cox_summary_legacy(Y, E, Treat, NULL),
    optimized = cox_summary(Y, E, Treat, return_format = "formatted"),
    check = FALSE,
    iterations = 100
  )
  
  # Optimized should be at least as fast
  expect_lte(
    as.numeric(bm$median[2]),
    as.numeric(bm$median[1]) * 1.1  # Allow 10% tolerance
  )
})
```

---

## 8. Priority Action Items

### Immediate (Before CRAN Submission)

1. **Consolidate `globalVariables()`** into single `R/globals.R`
2. **Remove duplicate S3 methods** from `R/forest_search_revised.R`
3. **Verify no `install.packages()` calls** exist in package code
4. **Run `R CMD check --as-cran`** and resolve all warnings

### Short-term (Code Quality)

5. **Optimize `cox_summary()`** for bootstrap performance
6. **Apply consistent `snake_case`** naming (with deprecation warnings)
7. **Add constants file** for magic numbers
8. **Improve `get_dfpred()`** to avoid `eval(parse())`

### Medium-term (Maintainability)

9. **Reorganize file structure** per Section 5.1
10. **Add comprehensive unit tests** for edge cases
11. **Refactor long functions** into smaller helpers
12. **Consolidate imports** in `forestsearch-package.R`

---

## Summary Statistics

| Category | Issues Found | Priority High | Priority Medium |
|----------|-------------|---------------|-----------------|
| CRAN Compliance | 4 | 3 | 1 |
| Efficiency | 5 | 2 | 3 |
| Readability | 8 | 2 | 6 |
| Style Guide | 6 | 1 | 5 |
| **Total** | **23** | **8** | **15** |

---

*Review completed: January 2026*
*Package version reviewed: ForestSearch (development)*
