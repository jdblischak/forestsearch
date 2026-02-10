# Consolidate globalVariables() — Removal Guide

## What changed in `R/globals.R`

Three new sections were added before the closing `))`:

1. **GBSG DGM internal variables** — 20 vars from `sim_aft_gbsg_refactored.R`
2. **Operating characteristics** — 20 vars from `oc_analyses_gbsg_refactored.R`  
3. **MRCT simulation variables** — 24 vars from `mrct_simulation.R`
4. **`"ksim"`** added to the existing Cross-validation section

Duplicates already present in `globals.R` (e.g., `"y"`, `"event"`, `"treat"`,
`"y_sim"`, `"event_sim"`, `"theta_0"`, `"theta_1"`, `"loghr_po"`, `"rfstime"`,
`"status"`, `"age"`, etc.) were NOT re-added — they're only in their original
section.

---

## Files to edit: Remove these exact blocks

### 1. `R/forestsearch_cross-validation.R`

**DELETE these 2 lines** (near top of file, before the first function):

```r
# Declare global variables to avoid R CMD check NOTEs for foreach loop variables
utils::globalVariables(c("ksim"))
```

### 2. `R/mrct_simulation.R`

**DELETE this entire block** (lines ~19–37, between the file header and the
first `#'` roxygen block):

```r
# =============================================================================
# Global Variables Declaration (for R CMD check)
# =============================================================================
# These variables are used in data.table NSE contexts

utils::globalVariables(c(
  # foreach loop variable
 "sim",
 
  # MRCT simulation variables
  "hr_test", "hr_sg", "any_found", "sg_found", "hr_sg_null",
  "regAflag", "sg_le85", "regAflag2", "regAflag3", "found",
  "sg_biomarker", "sg_age", "sg_male", "sg_ecog", "sg_histology",
  "sg_CTregimen", "sg_region", "sg_surgery", "sg_prior_treat",
  "analysis", "est",
  

  # Data splitting variables
  "region_var", "z_regA"
))
```

### 3. `R/sim_aft_gbsg_refactored.R`

**DELETE this entire block** (after the file header comment, before `# Constants`):

```r
# =============================================================================
# Global Variables Declaration
# =============================================================================

utils::globalVariables(c(
  # GBSG dataset variables
  "rfstime", "status", "hormon", "er", "age", "pgr", "meno", "nodes",
  "grade", "size",
  
  # Derived variables
  "y", "event", "treat", "z1", "z2", "z3", "z4", "z5", "zh", "flag.harm",
  "v1", "v2", "v3", "v4", "v5", "v6", "v7", "grade3",
  "lin.conf.true", "lin1.conf", "lin0.conf", "linC1.conf", "linC0.conf",
  "hlin.conf.1", "hlin.conf.0", "hlin.ratio", "h1.potential", "h0.potential",
  "Ts", "es", "t.sim", "y.sim", "event.sim",
  
  # New aligned variables (matching generate_aft_dgm_flex)
  "theta_0", "theta_1", "loghr_po", "lin_pred_0", "lin_pred_1"
))
```

### 4. `R/oc_analyses_gbsg_refactored.R`

**DELETE this entire block** (after the `#' @import` directives and before
`# Configuration Defaults`):

```r
# =============================================================================
# Global Variables Declaration
# =============================================================================

utils::globalVariables(c(
  "any.H", "ppv", "npv", "sens", "spec",
  "size.H", "size.Hc", "hr.H.true", "hr.H.hat", "hr.Hc.true", "hr.Hc.hat",
  "hr.itt", "hr.adj.itt", "p.cens", "taumax",
  "analysis", "sim", "aa", "sg_hat",
  # New aligned variables
  "ahr.H.true", "ahr.Hc.true", "ahr.H.hat", "ahr.Hc.hat",
  "loghr_po", "theta_0", "theta_1"
))
```

**KEEP** the `@import` directives that precede this block:
```r
#' @import survival
#' @import data.table
NULL
```

---

## Verification

After making these changes, run:

```r
# Confirm no scattered declarations remain
grep -rn "globalVariables" R/ | grep -v "globals.R"
# Should return ZERO results

# Full check
devtools::check()
# Should still show 0 errors, 0 warnings, 0 notes
```
