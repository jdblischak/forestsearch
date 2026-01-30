# GBSG DGM Architecture: create_gbsg_dgm() and calibrate_k_inter()

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
               │  Uses uniroot() to find k_inter
               │  where: dgm$hr_H_true - target_hr_harm = 0
               │
               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                     OBJECTIVE FUNCTION (internal)                            │
│                                                                             │
│   objective <- function(k_val) {                                            │
│       dgm <- create_gbsg_dgm(                                               │
│           model = model,                                                    │
│           k_treat = k_treat,                                                │
│           k_inter = k_val,    ◀─── Varies this to find root                 │
│           ...                                                               │
│       )                                                                     │
│       dgm$hr_H_true - target_hr_harm   ◀─── Returns difference              │
│   }                                                                         │
└─────────────────────────────────────────────────────────────────────────────┘
               │
               │  Calls repeatedly with different k_inter values
               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                          create_gbsg_dgm()                                   │
│                                                                             │
│  Parameters:                                                                │
│  • model = c("alt", "null")                                                 │
│  • k_treat = 1          ◀─── Treatment effect multiplier                    │
│  • k_inter = 1          ◀─── Interaction effect multiplier (KEY PARAM)      │
│  • k_z3 = 1             ◀─── z3 (meno) effect multiplier                    │
│  • z1_quantile = 0.25   ◀─── ER quantile threshold                          │
│  • n_super = 5000       ◀─── Super-population size                          │
│  • cens_type            ◀─── Censoring distribution                         │
│  • seed                 ◀─── Random seed                                    │
│  • verbose                                                                  │
└─────────────────────────────────────────────────────────────────────────────┘
               │
               │
               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                    create_gbsg_dgm() INTERNAL WORKFLOW                       │
└─────────────────────────────────────────────────────────────────────────────┘

    STEP 1: Load & Prepare GBSG Data
    ─────────────────────────────────
    ┌─────────────────────────────────────────────────────────────────────────┐
    │  dfa <- survival::gbsg                                                  │
    │                                                                         │
    │  Create derived variables:                                              │
    │  • y = rfstime / 30.4375 (months)                                       │
    │  • event = status                                                       │
    │  • treat = hormon                                                       │
    │                                                                         │
    │  Create subgroup-defining variables:                                    │
    │  • z1 = ifelse(er <= quantile(er, z1_quantile), 1, 0)  # Low ER         │
    │  • z2 = cut_numeric(age)                               # Age quartiles  │
    │  • z3 = ifelse(meno == 0, 1, 0)                        # Premenopausal  │
    │  • z4 = cut_numeric(pgr)                               # PgR quartiles  │
    │  • z5 = cut_numeric(nodes)                             # Nodes quartiles│
    │                                                                         │
    │  Create interaction term:                                               │
    │  • zh = treat * z1 * z3    ◀─── Three-way interaction                   │
    │                                                                         │
    │  Create harm subgroup indicator:                                        │
    │  • flag.harm = ifelse(z1 == 1 & z3 == 1, 1, 0)                          │
    └─────────────────────────────────────────────────────────────────────────┘
               │
               ▼
    STEP 2: Fit AFT Model
    ──────────────────────
    ┌─────────────────────────────────────────────────────────────────────────┐
    │  Formula: Surv(y, event) ~ treat + z1 + z2 + z3 + z4 + z5 + zh          │
    │                                                                         │
    │  fit_aft <- survreg(..., dist = "weibull")                              │
    │                                                                         │
    │  Extract parameters:                                                    │
    │  • sigma = fit_aft$scale                                                │
    │  • mu = coef(fit_aft)[1]  (intercept)                                   │
    │  • gamma = coef(fit_aft)[-1]  (regression coefficients)                 │
    └─────────────────────────────────────────────────────────────────────────┘
               │
               ▼
    STEP 3: Apply Effect Modifiers  ◀═══════════════════════════════════════╗
    ─────────────────────────────────                                        ║
    ┌─────────────────────────────────────────────────────────────────────────┐
    │                                                                         │
    │  gamma["z3"] <- k_z3 * gamma["z3"]                                      │
    │                                                                         │
    │  gamma["treat"] <- k_treat * gamma["treat"]                             │
    │                                                                         │
    │  if (model == "alt") {                                                  │
    │      gamma["zh"] <- k_inter * gamma["zh"]   ◀─── KEY: Modifies zh coef  ║
    │  }                                                                      │
    │                                                                         │
    │  Convert to hazard scale:                                               │
    │  • b_true = gamma                                                       │
    │  • b0 = -gamma / sigma                                                  │
    └─────────────────────────────────────────────────────────────────────────┘
               │                                                              ║
               ▼                                                              ║
    STEP 4: Generate Randomized Super-Population                              ║
    ────────────────────────────────────────────                              ║
    ┌─────────────────────────────────────────────────────────────────────────┐
    │  Sample n_super observations with replacement                           │
    │                                                                         │
    │  Split into:                                                            │
    │  • df_treat (n_super/2): treat = 1, lin.conf.true = lin1.conf           │
    │  • df_control (n_super/2): treat = 0, lin.conf.true = lin0.conf         │
    │                                                                         │
    │  Combine: df_big = rbind(df_treat, df_control)                          │
    │                                                                         │
    │  ⚠️  POTENTIAL ISSUE: Does zh get recalculated for the super-pop?       ║
    │      zh = treat * z1 * z3                                               │
    │      The zh variable may not be updated after treatment assignment!     ║
    └─────────────────────────────────────────────────────────────────────────┘
               │                                                              ║
               ▼                                                              ║
    STEP 5: Compute Empirical Hazard Ratios                                   ║
    ───────────────────────────────────────                                   ║
    ┌─────────────────────────────────────────────────────────────────────────┐
    │  Generate survival times:                                               │
    │  • epsilon = log(rexp(n_super))                                         │
    │  • log_Ts = mu + sigma * epsilon + df_big$lin.conf.true                 │
    │  • Ts = exp(log_Ts)                                                     │
    │  • es = 1 (all events)                                                  │
    │                                                                         │
    │  Compute HR for harm subgroup:                                          │
    │  hr_H_true <- exp(coxph(                                                │
    │      Surv(Ts, es) ~ treat,                                              │
    │      data = subset(df_big, flag.harm == 1)     ◀─── Uses flag.harm      ║
    │  )$coefficients)                                                        │
    │                                                                         │
    │  ⚠️  BUG LIKELY HERE: The survival times Ts are generated using         ║
    │      lin.conf.true which includes the modified gamma["zh"] coefficient, ║
    │      BUT the flag.harm indicator was set BEFORE resampling and doesn't  ║
    │      account for the actual treatment assignment in the super-pop!      ║
    └─────────────────────────────────────────────────────────────────────────┘
               │
               ▼
    STEP 6: Return Results
    ──────────────────────
    ┌─────────────────────────────────────────────────────────────────────────┐
    │  result <- list(                                                        │
    │      df_super_rand = df_big,                                            │
    │      hr_H_true = hr_H_true,      ◀─── This is returned to calibration   │
    │      hr_Hc_true = hr_Hc_true,                                           │
    │      hr_causal = hr_causal,                                             │
    │      model_params = list(mu, sigma, gamma, b_true, b_hr),               │
    │      cens_params = list(type, mu, sigma),                               │
    │      subgroup_info = list(fs_harm_true, grf_harm_true, ...),            │
    │      analysis_vars = c("v1", "v2", "v3", "v4", "v5", "v6", "v7"),        │
    │      model_type = model,                                                │
    │      n_super = n_super,                                                 │
    │      seed = seed                                                        │
    │  )                                                                      │
    └─────────────────────────────────────────────────────────────────────────┘


┌─────────────────────────────────────────────────────────────────────────────┐
│                           HELPER FUNCTIONS                                   │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────┐     ┌─────────────────────────────┐
│     cut_numeric()           │     │      cut_size()             │
│                             │     │                             │
│  Discretizes continuous     │     │  Discretizes tumor size     │
│  variable into quantile-    │     │  using fixed breakpoints    │
│  based categories           │     │  Default: c(20, 50)         │
│                             │     │                             │
│  Default probs:             │     │  Returns: 1, 2, or 3        │
│  c(0.25, 0.5, 0.75)         │     │                             │
│  Returns: 1, 2, 3, or 4     │     │                             │
└─────────────────────────────┘     └─────────────────────────────┘


┌─────────────────────────────────────────────────────────────────────────────┐
│                         DOWNSTREAM FUNCTIONS                                 │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────┐     ┌─────────────────────────────┐
│   get_dgm_with_output()     │     │  simulate_from_gbsg_dgm()   │
│                             │     │                             │
│  Wrapper that:              │     │  Generates simulated trial  │
│  1. Calls calibrate_k_inter │     │  data from the DGM          │
│     if target_hr provided   │     │                             │
│  2. Creates DGM             │     │  Parameters:                │
│  3. Generates output path   │     │  • dgm (gbsg_dgm object)    │
│                             │     │  • n (sample size)          │
│                             │     │  • rand_ratio               │
│                             │     │  • max_follow               │
│                             │     │  • muC_adj (censoring)      │
└─────────────────────────────┘     └─────────────────────────────┘


┌─────────────────────────────────────────────────────────────────────────────┐
│                        RELATED FLEXIBLE DGM FUNCTIONS                        │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────┐     ┌─────────────────────────────┐
│  find_k_inter_for_target_hr │     │ sensitivity_analysis_k_inter│
│                             │     │                             │
│  Similar to calibrate_k_inter│    │  Evaluates HRs across       │
│  but for generate_aft_dgm_  │     │  range of k_inter values    │
│  flex() DGMs                │     │                             │
│                             │     │  Creates 4-panel plot:      │
│  Uses same uniroot approach │     │  1. HR_harm vs k_inter      │
│                             │     │  2. All HRs vs k_inter      │
│                             │     │  3. HR ratio                │
│                             │     │  4. Key values table        │
└─────────────────────────────┘     └─────────────────────────────┘
```

---

## Potential Bug Analysis: HR(H) Calculation Issue

### The Problem

When `k_inter` is passed to `create_gbsg_dgm()`, the "Target HR(H)" computed may not correctly reflect the effect of the interaction term.

### Root Cause Analysis

Looking at the code flow, there are **two potential issues**:

#### Issue 1: zh Variable Not Updated After Resampling

```r
# ORIGINAL DATA PREPARATION (before resampling)
dfa$zh <- dfa$treat * dfa$z1 * dfa$z3    # Set based on ORIGINAL treatment

# SUPER-POPULATION GENERATION (treatment reassigned)
df_treat$treat <- 1L                      # Treatment CHANGED
df_control$treat <- 0L                    # Treatment CHANGED

# BUT zh is NOT recalculated!
# zh still reflects the ORIGINAL treatment assignment, not the new one
```

**The interaction term `zh = treat * z1 * z3` is computed on the original data, but when the super-population is created with new treatment assignments, `zh` is not updated.**

#### Issue 2: Linear Predictor Construction

The linear predictor `lin.conf.true` needs to properly incorporate the interaction effect:

```r
# The gamma["zh"] coefficient is modified:
gamma["zh"] <- k_inter * gamma["zh"]

# But the linear predictor is computed using the ORIGINAL zh values:
lin1.conf = z_true_1 %*% b_true   # Uses design matrix with treat=1
lin0.conf = z_true_0 %*% b_true   # Uses design matrix with treat=0
```

**If the design matrices `z_true_1` and `z_true_0` don't properly update the `zh` column to reflect the treatment assignment, the k_inter effect won't propagate correctly.**

---

## Recommended Fix

### Option A: Update zh in Design Matrices

When creating `z_true_1` and `z_true_0`, explicitly set the interaction term:

```r
# For treatment group (treat = 1)
z_true_1 <- z_true
z_true_1[, "treat"] <- 1
z_true_1[, "zh"] <- 1 * z_true_1[, "z1"] * z_true_1[, "z3"]  # UPDATE zh

# For control group (treat = 0)  
z_true_0 <- z_true
z_true_0[, "treat"] <- 0
z_true_0[, "zh"] <- 0 * z_true_0[, "z1"] * z_true_0[, "z3"]  # UPDATE zh (= 0)
```

### Option B: Recalculate zh After Treatment Assignment

In the super-population generation step:

```r
# After assigning treatment
df_treat$treat <- 1L
df_treat$zh <- df_treat$treat * df_treat$z1 * df_treat$z3  # Recalculate

df_control$treat <- 0L
df_control$zh <- df_control$treat * df_control$z1 * df_control$z3  # = 0
```

### Option C: Verify Design Matrix Construction

Check that the model formula and design matrix extraction properly handles the interaction:

```r
# Ensure covs_true includes zh
covs_true <- c("treat", "z1", "z2", "z3", "z4", "z5", "zh")

# Verify z_true matrix has correct zh values
z_true <- as.matrix(dfa[, covs_true])

# When creating potential outcome matrices, ensure:
# - z_true_1 has zh = 1 * z1 * z3 (for subjects where z1=1 & z3=1)
# - z_true_0 has zh = 0 (always, since treat=0)
```

---

## Verification Test

To verify the fix works correctly:

```r
# Test 1: Boundary conditions
dgm_k0 <- create_gbsg_dgm(model = "alt", k_inter = 0, verbose = TRUE)
# HR(H) should equal HR(Hc) when k_inter = 0 (no interaction effect)

dgm_k1 <- create_gbsg_dgm(model = "alt", k_inter = 1, verbose = TRUE)
# HR(H) should show some heterogeneity

dgm_k2 <- create_gbsg_dgm(model = "alt", k_inter = 2, verbose = TRUE)
# HR(H) should show MORE heterogeneity (larger departure from HR(Hc))

# Test 2: Calibration round-trip
target_hr <- 1.5
k_found <- calibrate_k_inter(target_hr_harm = target_hr, verbose = TRUE)
dgm_verify <- create_gbsg_dgm(model = "alt", k_inter = k_found, verbose = TRUE)
print(paste("Target:", target_hr, "Achieved:", dgm_verify$hr_H_true))
# Should be very close to target
```

---

## Summary Table

| Function | Location | Purpose | Calls/Called By |
|----------|----------|---------|-----------------|
| `calibrate_k_inter()` | `sim_aft_gbsg_refactored.R` | Find k_inter for target HR | Calls `create_gbsg_dgm()` |
| `create_gbsg_dgm()` | `sim_aft_gbsg_refactored.R` | Create GBSG-based DGM | Called by `calibrate_k_inter()`, `get_dgm_with_output()` |
| `cut_numeric()` | `sim_aft_gbsg_refactored.R` | Discretize continuous vars | Called by `create_gbsg_dgm()` |
| `cut_size()` | `sim_aft_gbsg_refactored.R` | Discretize tumor size | Called by `create_gbsg_dgm()` |
| `get_dgm_with_output()` | `sim_aft_gbsg_refactored.R` | Wrapper with file output | Calls `calibrate_k_inter()`, `create_gbsg_dgm()` |
| `simulate_from_gbsg_dgm()` | `sim_aft_gbsg_refactored.R` | Generate trial data | Called after DGM creation |
| `print.gbsg_dgm()` | `sim_aft_gbsg_refactored.R` | S3 print method | N/A |

---

## AFT Model Coefficients and HR Relationship

The AFT model specifies:
```
log(T) = μ + σ * ε + X * γ
```

Where:
- `γ["treat"]`: Main treatment effect (modified by `k_treat`)
- `γ["zh"]`: Interaction effect treat × z1 × z3 (modified by `k_inter`)

The hazard ratio for the harm subgroup (z1=1, z3=1) should be:
```
HR_H = exp(-γ["treat"]/σ - γ["zh"]/σ)
     = exp(-k_treat * γ_orig["treat"]/σ - k_inter * γ_orig["zh"]/σ)
```

And for the complement (z1=0 OR z3=0):
```
HR_Hc = exp(-γ["treat"]/σ)
      = exp(-k_treat * γ_orig["treat"]/σ)
```

The difference between HR_H and HR_Hc should be driven entirely by `k_inter * γ["zh"]`.
