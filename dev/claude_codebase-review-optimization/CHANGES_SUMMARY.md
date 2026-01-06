# ForestSearch CRAN Compliance Fixes - Summary of Changes

## Overview

This document summarizes all critical and CRAN compliance fixes applied to the ForestSearch package files. These files are drop-in replacements for the original versions.

---

## Files Provided

| File | Description |
|------|-------------|
| `generate_bootstrap_synthetic_general.R` | Bootstrap synthetic data generation |
| `generate_aft_dgm_improved.R` | AFT model DGM (basic) |
| `generate_aft_dgm_flexible.R` | AFT model DGM (flexible cutpoints) |
| `bootstrap_summaries_helpers_fixed.R` | Bootstrap result summaries |
| `process_factors_corrected.R` | Factor variable processing |
| `NAMESPACE` | Package namespace template |

---

## Critical Bug Fixes

### 1. `enumerate()` Syntax Error (CRITICAL)

**Files affected:** `generate_aft_dgm_improved.R`, `generate_aft_dgm_flexible.R`

**Original (BROKEN):**
```r
for (i, var in enumerate(subgroup_vars)) {
```

**Fixed:**
```r
for (var in subgroup_vars) {
```

**Reason:** R does not have an `enumerate()` function. This was Python syntax that would cause immediate runtime errors.

---

### 2. `library()` Calls Inside Functions (CRAN REJECTION)

**Files affected:** `generate_bootstrap_synthetic_general.R`

**Original (REJECTED BY CRAN):**
```r
generate_gbsg_bootstrap_general <- function(...) {
  library(survival)  # NOT ALLOWED
  data(cancer)
  ...
}
```

**Fixed:**
```r
generate_gbsg_bootstrap_general <- function(...) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }
  gbsg <- survival::gbsg  # Access via namespace
  ...
}
```

**Reason:** CRAN policy prohibits `library()` or `require()` calls inside package functions. Use `requireNamespace()` and namespace-qualified access instead.

---

### 3. Corrupted Unicode Characters (CRAN WARNING)

**Files affected:** `bootstrap_summaries_helpers_fixed.R`

**Original (CORRUPTED):**
```r
labels = c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.9", "0.9-0.95", "â‰¥0.95"),
size_labels <- c("<50", "50-99", "100-149", "150-199", "200-299", "â‰¥300")
```

**Fixed:**
```r
labels = c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.9", "0.9-0.95", ">=0.95"),
size_labels <- c("<50", "50-99", "100-149", "150-199", "200-299", ">=300")
```

**Reason:** The `≥` symbol was corrupted during file transfer. CRAN requires either ASCII or properly escaped Unicode (`\u2265`).

---

## CRAN Compliance Additions

### 4. Added `@export` Tags

All public functions now have explicit `@export` tags in their roxygen2 documentation:

```r
#' @export
generate_bootstrap_synthetic <- function(...) { }
```

**Functions exported:**
- `generate_bootstrap_synthetic`
- `detect_variable_types`
- `generate_gbsg_bootstrap_general`
- `compare_datasets_general`
- `generate_aft_dgm`
- `generate_aft_dgm_flex`
- `simulate_from_dgm`
- `summarize_simulation`
- `summarize_bootstrap_subgroups`
- `process_factor_variables`
- `process_covariates_corrected`
- And others...

---

### 5. Added `@importFrom` Declarations

Each file now includes proper import declarations:

```r
#' @importFrom stats rnorm runif sd median quantile
#' @importFrom survival survreg coxph Surv
#' @importFrom data.table data.table as.data.table .N
```

---

### 6. Enhanced Documentation

All functions now have:
- Complete `@param` descriptions
- `@return` value documentation
- Working `@examples` (wrapped in `\donttest{}` where appropriate)
- `@seealso` cross-references
- `@details` sections for complex functions

---

## Efficiency Improvements

### 7. Vectorized Categorical Perturbation

**Original:** Loop-based perturbation (slow for large datasets)

**Improved:** Pre-computed replacement maps with vectorized operations

```r
# Pre-compute replacement options
replacement_map <- lapply(unique_vals, function(v) setdiff(unique_vals, v))
names(replacement_map) <- as.character(unique_vals)

# Vectorized replacement
new_vals <- vapply(current_vals, function(cv) {
  sample(replacement_map[[cv]], 1)
}, FUN.VALUE = character(1))
```

---

### 8. Cached Integer Detection

**Original:** Checked `all(x == round(x))` inside loop for every variable

**Improved:** Pre-computed before loop

```r
is_integer_var <- vapply(continuous_vars, function(var) {
  all(data[[var]] == round(data[[var]]), na.rm = TRUE)
}, FUN.VALUE = logical(1))
```

---

### 9. Fixed `original_agreement` Calculation

**Original (`bootstrap_summaries_helpers_fixed.R`):**
```r
original_agreement = NULL,  # Simplified for now  ← ALWAYS NULL!
```

**Fixed:** Actually calculates agreement with original subgroup:
```r
if (!is.null(original_sg) && n_found > 0) {
  # Find subgroup column and calculate matches
  matches <- sum(as.character(sg_found[[subgroup_col]]) == orig_sg_char, na.rm = TRUE)
  
  original_agreement <- data.table::data.table(
    Metric = c("Total bootstrap iterations", "Exact match with original", ...),
    Value = c(as.character(nb_boots), sprintf("%d (%.1f%%)", matches, ...), ...)
  )
}
```

---

## Installation Instructions

### Option 1: Replace Individual Files

1. Copy each `.R` file to your package's `R/` directory
2. Run `devtools::document()` to regenerate NAMESPACE and documentation
3. Run `R CMD check` to verify

### Option 2: Use the NAMESPACE Template

1. Copy all `.R` files to `R/`
2. Replace your existing `NAMESPACE` with the provided template (or run `devtools::document()`)
3. Verify with `R CMD check`

---

## Verification Checklist

After installing the fixed files:

- [ ] `R CMD check` produces no ERRORs
- [ ] `R CMD check` produces no WARNINGs about:
  - [ ] Non-standard things in the check directory
  - [ ] Missing exports
  - [ ] Undocumented arguments
  - [ ] Non-ASCII characters
- [ ] All examples run without error
- [ ] `generate_aft_dgm()` works (no `enumerate` error)
- [ ] `generate_gbsg_bootstrap_general()` works (no `library` error)
- [ ] Bootstrap summaries show Table 5 (original_agreement not NULL)

---

## Testing the Fixes

```r
# Test 1: enumerate() fix
library(survival)
data(cancer)

dgm <- generate_aft_dgm(
  data = gbsg,
  continuous_vars = c("age", "size"),
  factor_vars = c("grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("age"),
  subgroup_cuts = list(age = 50),
  model = "alt",
  verbose = TRUE
)
# Should complete without error

# Test 2: library() fix
synth <- generate_gbsg_bootstrap_general(n = 100)
# Should work without loading survival into global environment

# Test 3: Unicode fix (check labels in output)
boot_results <- data.frame(
  Pcons = runif(50, 0.5, 1),
  M.1 = rep("age<=50", 50),
  N_sg = sample(50:200, 50, replace = TRUE)
)

summary <- summarize_bootstrap_subgroups(boot_results, nb_boots = 50)
print(summary$consistency_dist)
# Should show ">=0.95" not corrupted characters
```

---

## Summary Table

| Issue | Severity | Files | Status |
|-------|----------|-------|--------|
| `enumerate()` syntax | CRITICAL | 2 | ✅ Fixed |
| `library()` in functions | CRAN REJECT | 1 | ✅ Fixed |
| Unicode corruption | CRAN WARNING | 1 | ✅ Fixed |
| Missing `@export` | CRAN WARNING | 5 | ✅ Added |
| Missing `@importFrom` | CRAN NOTE | 5 | ✅ Added |
| `original_agreement = NULL` | BUG | 1 | ✅ Fixed |
| Loop inefficiencies | PERFORMANCE | 2 | ✅ Improved |

---

## Contact

These fixes were prepared based on code review of the ForestSearch package. If you encounter issues after applying these fixes, the most likely causes are:

1. **Merge conflicts** with other changes you've made
2. **Missing dependencies** (ensure survival and data.table are in Imports)
3. **roxygen2 regeneration** needed (`devtools::document()`)

Run `R CMD check --as-cran` to verify full CRAN compliance.
