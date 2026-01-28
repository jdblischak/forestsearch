# GBSG DGM Architecture: Bug Fix Documentation

## Overview

This document describes the bug fix implemented in `sim_aft_gbsg_refactored.R` to correct the `k_inter` parameter handling in `create_gbsg_dgm()`.

---

## The Bug

### Problem Statement

When using `k_inter` to modulate the hazard ratio in the harm subgroup, the "Target HR(H)" computed by `create_gbsg_dgm()` did not correctly reflect the interaction effect.

### Root Cause

The interaction term `zh = treat * z1 * z3` was not properly updated in the potential outcome design matrices. Specifically:

1. **Original code created `zh` once** based on observed treatment in the source data
2. **Potential outcome matrices (`z_true_1`, `z_true_0`) were created** by copying and modifying `treat` column
3. **BUT the `zh` column was NOT recalculated** to reflect the counterfactual treatment assignment

```r
# ORIGINAL PROBLEMATIC CODE:
z_true_1 <- z_true
z_true_1[, 1] <- 1  # Set treat = 1
z_true_1[, loc_inter] <- z_true[, "z1"] * z_true[, "z3"]  # Used original z_true values

z_true_0 <- z_true
z_true_0[, 1] <- 0  # Set treat = 0
z_true_0[, loc_inter] <- 0  # This was correct (zh=0 when treat=0)
```

The issue was that `z_true_1[, loc_inter]` used values from `z_true` which contained the **original observed treatment**, not the counterfactual treatment=1 scenario.

---

## The Fix

### Key Changes

#### 1. Explicit Potential Outcome Design Matrices

```r
# FIXED CODE:
# Get z1 and z3 as numeric vectors
z1_vec <- as.numeric(dfa$z1)
z3_vec <- as.numeric(dfa$z3)

if (model == "alt") {
  # Potential outcome under treatment (treat = 1)
  z_true_1 <- z_true
  z_true_1[, col_treat] <- 1           # Set treat = 1
  z_true_1[, col_zh] <- z1_vec * z3_vec  # zh = 1 * z1 * z3 = z1 * z3
  
  # Potential outcome under control (treat = 0)
  z_true_0 <- z_true
  z_true_0[, col_treat] <- 0           # Set treat = 0
  z_true_0[, col_zh] <- 0              # zh = 0 * z1 * z3 = 0 (ALWAYS)
}
```

#### 2. Update `zh` After Super-Population Treatment Assignment

```r
# Treatment group
df_treat$treat <- 1L
df_treat$lin.conf.true <- df_treat$lin1.conf
df_treat$zh <- df_treat$z1 * df_treat$z3  # CRITICAL: Recalculate zh

# Control group
df_control$treat <- 0L
df_control$lin.conf.true <- df_control$lin0.conf
df_control$zh <- 0L  # CRITICAL: zh = 0 for control
```

#### 3. Same Fix in `simulate_from_gbsg_dgm()`

```r
# When drawing treatment
df_treat$zh <- df_treat$z1 * df_treat$z3
df_control$zh <- 0L
```

---

## Function Relationship Diagram

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         CALIBRATION WORKFLOW                                 │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────┐
│    calibrate_k_inter()      │
│                             │
│  Parameters:                │
│  • target_hr_harm           │
│  • model = "alt"            │
│  • k_treat = 1              │
│  • k_inter_range            │
│  • tol                      │
└──────────────┬──────────────┘
               │
               │  Uses uniroot() to solve:
               │  HR_H(k_inter) - target_hr_harm = 0
               │
               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                          create_gbsg_dgm()                                   │
│                                                                             │
│  Key Parameters:                                                            │
│  • k_treat: Multiplier for gamma["treat"]                                   │
│  • k_inter: Multiplier for gamma["zh"] (interaction)                        │
│  • k_z3: Multiplier for gamma["z3"]                                         │
└─────────────────────────────────────────────────────────────────────────────┘
               │
               │
               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       INTERNAL WORKFLOW (FIXED)                              │
└─────────────────────────────────────────────────────────────────────────────┘

    Step 1: Prepare GBSG Data
    ─────────────────────────
    • Load survival::gbsg dataset
    • Create z1 (low ER), z2 (age), z3 (premenopausal), z4 (PgR), z5 (nodes)
    • Create zh = treat * z1 * z3 (observed interaction)
    • Create flag.harm = (z1 == 1 & z3 == 1)

    Step 2: Create Design Matrices (FIXED)
    ──────────────────────────────────────
    ┌─────────────────────────────────────────────────────────────────────┐
    │  z_true_1 (Potential outcome under TREATMENT):                      │
    │  • treat column = 1                                                 │
    │  • zh column = z1 * z3  ◀── FIXED: Use z1, z3 directly              │
    │                                                                     │
    │  z_true_0 (Potential outcome under CONTROL):                        │
    │  • treat column = 0                                                 │
    │  • zh column = 0        ◀── FIXED: Always 0 when treat=0            │
    └─────────────────────────────────────────────────────────────────────┘

    Step 3: Fit AFT Model & Apply Effect Modifiers
    ───────────────────────────────────────────────
    • Fit survreg(Surv(y, event) ~ z_true)
    • Extract gamma coefficients
    • Apply modifiers:
        gamma["treat"] <- k_treat * gamma["treat"]
        gamma["zh"]    <- k_inter * gamma["zh"]
        gamma["z3"]    <- k_z3 * gamma["z3"]

    Step 4: Compute Linear Predictors
    ──────────────────────────────────
    • lin1_conf = z_true_1 %*% b_true  (under treatment)
    • lin0_conf = z_true_0 %*% b_true  (under control)

    Step 5: Generate Super-Population (FIXED)
    ─────────────────────────────────────────
    ┌─────────────────────────────────────────────────────────────────────┐
    │  df_treat:                                                          │
    │  • treat = 1                                                        │
    │  • lin.conf.true = lin1.conf                                        │
    │  • zh = z1 * z3           ◀── FIXED: Recalculate                    │
    │                                                                     │
    │  df_control:                                                        │
    │  • treat = 0                                                        │
    │  • lin.conf.true = lin0.conf                                        │
    │  • zh = 0                 ◀── FIXED: Always 0                       │
    └─────────────────────────────────────────────────────────────────────┘

    Step 6: Compute Empirical HRs
    ─────────────────────────────
    • Generate survival times: log_Ts = mu + sigma*epsilon + lin.conf.true
    • Fit Cox model in harm subgroup (flag.harm == 1) → hr_H_true
    • Fit Cox model in complement (flag.harm == 0) → hr_Hc_true


┌─────────────────────────────────────────────────────────────────────────────┐
│                         HELPER FUNCTIONS                                     │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────┐     ┌─────────────────────────────┐
│     cut_numeric()           │     │      cut_size()             │
│ Discretize to quartiles     │     │ Discretize tumor size       │
└─────────────────────────────┘     └─────────────────────────────┘

┌─────────────────────────────┐     ┌─────────────────────────────┐
│   get_dgm_with_output()     │     │  simulate_from_gbsg_dgm()   │
│ Wrapper with calibration    │     │ Generate trial data         │
│ and file path generation    │     │ (also fixed zh handling)    │
└─────────────────────────────┘     └─────────────────────────────┘

┌─────────────────────────────┐     ┌─────────────────────────────┐
│   print.gbsg_dgm()          │     │ validate_k_inter_effect()   │
│ S3 print method             │     │ NEW: Test k_inter behavior  │
└─────────────────────────────┘     └─────────────────────────────┘
```

---

## Mathematical Background

### AFT Model Specification

```
log(T) = μ + σ * ε + X * γ
```

Where:
- `T` = survival time
- `μ` = intercept
- `σ` = scale parameter
- `ε` = standard extreme value error
- `X` = design matrix
- `γ` = regression coefficients

### Hazard Ratio Derivation

The hazard scale coefficients are: `b0 = -γ / σ`

For a subject in the **harm subgroup** (z1=1, z3=1) under **treatment** (treat=1):
- zh = treat × z1 × z3 = 1 × 1 × 1 = 1
- Log-hazard contribution: `b0["treat"] + b0["zh"]`
- HR(H) = exp(b0["treat"] + b0["zh"])

For a subject in the **complement** (z1=0 OR z3=0) under **treatment**:
- zh = treat × z1 × z3 = 0 (since z1*z3 = 0)
- Log-hazard contribution: `b0["treat"]`
- HR(Hc) = exp(b0["treat"])

### Effect of k_inter

```
γ["zh"] ← k_inter × γ_orig["zh"]
b0["zh"] = -γ["zh"] / σ = -k_inter × γ_orig["zh"] / σ

HR(H) = exp(b0["treat"] + k_inter × b0_orig["zh"])
HR(Hc) = exp(b0["treat"])

Ratio = HR(H) / HR(Hc) = exp(k_inter × b0_orig["zh"])
```

When **k_inter = 0**:
- Ratio = exp(0) = 1
- HR(H) = HR(Hc) (no heterogeneity)

When **k_inter > 0** (assuming b0_orig["zh"] > 0):
- Ratio > 1
- HR(H) > HR(Hc)

---

## Verification Tests

### Test 1: k_inter = 0 Should Give No Heterogeneity

```r
dgm_k0 <- create_gbsg_dgm(model = "alt", k_inter = 0, verbose = TRUE)
# Expected: hr_H_true ≈ hr_Hc_true (ratio ≈ 1.0)
```

### Test 2: Increasing k_inter Should Increase Heterogeneity

```r
dgm_k1 <- create_gbsg_dgm(model = "alt", k_inter = 1, verbose = TRUE)
dgm_k2 <- create_gbsg_dgm(model = "alt", k_inter = 2, verbose = TRUE)
dgm_k3 <- create_gbsg_dgm(model = "alt", k_inter = 3, verbose = TRUE)
# Expected: Increasing separation between hr_H_true and hr_Hc_true
```

### Test 3: Calibration Round-Trip

```r
target_hr <- 1.5
k_found <- calibrate_k_inter(target_hr_harm = target_hr, verbose = TRUE)
dgm_verify <- create_gbsg_dgm(model = "alt", k_inter = k_found, verbose = TRUE)
abs(dgm_verify$hr_H_true - target_hr) < 0.01  # Should be TRUE
```

### Test 4: Use Validation Function

```r
results <- validate_k_inter_effect(k_inter_values = c(-2, -1, 0, 1, 2, 3))
# Check that k_inter = 0 row has ratio ≈ 1.0
```

---

## New Features Added

### 1. `validate_k_inter_effect()` Function

New exported function to verify k_inter behavior:

```r
validate_k_inter_effect <- function(
    k_inter_values = c(-2, -1, 0, 1, 2, 3),
    verbose = TRUE,
    ...
)
```

### 2. Enhanced `print.gbsg_dgm()` Method

Now displays:
- Effect modifiers (k_treat, k_inter, k_z3)
- HR ratio (HR(H) / HR(Hc))
- ER threshold value

### 3. Additional Output in DGM Object

```r
dgm$effect_modifiers <- list(
  k_treat = k_treat,
  k_inter = k_inter,
  k_z3 = k_z3
)

dgm$model_params$gamma_orig  # Original coefficients before modification
dgm$subgroup_info$er_threshold  # Actual ER cutoff value
```

---

## Summary of Changes

| Location | Change | Purpose |
|----------|--------|---------|
| `create_gbsg_dgm()` lines ~250-265 | Explicit zh calculation in z_true_1, z_true_0 | Fix potential outcome matrices |
| `create_gbsg_dgm()` lines ~295-305 | Recalculate zh after treatment assignment | Fix super-population zh values |
| `simulate_from_gbsg_dgm()` lines ~470-485 | Recalculate zh when drawing treatment | Consistent fix in simulation |
| New function | `validate_k_inter_effect()` | Verification testing |
| `print.gbsg_dgm()` | Enhanced output | Better diagnostics |

---

## Files Provided

1. **sim_aft_gbsg_refactored.R** - Complete fixed R file with:
   - `create_gbsg_dgm()` - Fixed
   - `simulate_from_gbsg_dgm()` - Fixed
   - `calibrate_k_inter()` - Unchanged (works correctly with fixed DGM)
   - `get_dgm_with_output()` - Minor improvements
   - `validate_k_inter_effect()` - New validation function
   - `print.gbsg_dgm()` - Enhanced

2. **gbsg_dgm_bugfix_documentation.md** - This documentation file
