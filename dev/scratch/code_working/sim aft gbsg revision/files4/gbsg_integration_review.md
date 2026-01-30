# Integration Review: GBSG DGM Codebase with forestsearch Package Infrastructure

## Executive Summary

This document reviews the alignment and integration potential between:
1. **`sim_aft_gbsg_refactored.R`** (`create_gbsg_dgm()`, `simulate_from_gbsg_dgm()`) 
2. **`generate_aft_dgm_main.R`** (`generate_aft_dgm_flex()`) + **`simulate_from_dgm.R`** (`simulate_from_dgm()`)

**Recommendation**: The GBSG-specific functions can be reimplemented as **thin wrappers** around the generic `generate_aft_dgm_flex()` and `simulate_from_dgm()` functions, with GBSG-specific defaults and column name mappings.

---

## 1. DGM Creation: `create_gbsg_dgm()` vs `generate_aft_dgm_flex()`

### 1.1 Structural Comparison

| Aspect | `create_gbsg_dgm()` | `generate_aft_dgm_flex()` |
|--------|---------------------|---------------------------|
| **Input Data** | Hardcoded GBSG dataset via `survival::gbsg` | Any user-provided data.frame |
| **Continuous Vars** | Hardcoded: `age`, `er`, `pgr`, `nodes`, `size` | User-specified via `continuous_vars` |
| **Factor Vars** | Hardcoded: `meno`, `grade` | User-specified via `factor_vars` |
| **Subgroup Definition** | Hardcoded: `z1 == 1 & z3 == 1` (low ER & premenopausal) | Flexible via `subgroup_vars` + `subgroup_cuts` |
| **Variable Naming** | `z1-z5`, `v1-v7` (analysis vars) | `z_*` prefix for standardized vars |
| **Output Class** | `"gbsg_dgm"` | `"aft_dgm_flex"` |
| **Spline Support** | No | Yes (`spline_spec` parameter) |
| **Censoring Model** | Weibull or uniform | Weibull, lognormal, or uniform (auto-selected) |

### 1.2 Output Structure Comparison

```r
# create_gbsg_dgm() output
list(
  df_super_rand,      # Super-population data
  hr_H_true,          # HR in harm subgroup
  hr_Hc_true,         # HR in complement
  hr_causal,          # Overall causal HR
  AHR, AHR_H_true, AHR_Hc_true,  # AHR metrics (aligned version)
  hazard_ratios,      # Unified list (aligned version)
  model_params,       # mu, sigma, gamma, b_true, b_hr
  cens_params,        # type, mu, sigma
  subgroup_info,      # fs_harm_true, grf_harm_true, definition
  analysis_vars,      # c("v1", "v2", ...)
  model_type,         # "alt" or "null"
  n_super, seed
)

# generate_aft_dgm_flex() output
list(
  df_super,           # Super-population data
  model_params,       # mu, tau, gamma, b0, censoring, spline_info
  subgroup_info,      # vars, cuts, definitions, size, proportion
  hazard_ratios,      # overall, harm_subgroup, no_harm_subgroup, AHR, AHR_harm, AHR_no_harm
  analysis_vars,      # list(continuous, factor, covariates, treatment, outcome, event)
  model_type,         # "alt" or "null"
  n_super, seed
)
```

### 1.3 Key Differences in Implementation

#### a) Subgroup Definition

**GBSG (hardcoded)**:
```r
# In create_gbsg_dgm()
z1_cut <- quantile(gbsg$er, probs = z1_quantile)
df$z1 <- ifelse(df$er <= z1_cut, 1, 0)
df$z3 <- df$meno  # Already 0/1
flag_harm <- (df$z1 == 1) & (df$z3 == 1)
```

**Flex (configurable)**:
```r
# In generate_aft_dgm_flex()
subgroup_vars = c("er", "meno")
subgroup_cuts = list(
  er = list(type = "quantile", value = 0.25),
  meno = 0
)
```

#### b) Design Matrix Construction

Both use similar approaches for building design matrices with potential outcomes:
- `X_treat` (treatment arm, A=1)
- `X_control` (control arm, A=0)
- Both compute `lin_pred_0`, `lin_pred_1`, `theta_0`, `theta_1`, `loghr_po`

#### c) Hazard Ratio Computation

Both use the **stacked potential outcomes** approach:
```r
# Both functions do this:
df_temp <- data.frame(
  time = c(T_1, T_0),
  event = 1,
  treat = c(rep(1, n_super), rep(0, n_super)),
  flag_harm = rep(flag_harm, 2)
)
hr_overall <- exp(coxph(Surv(time, event) ~ treat, data = df_temp)$coef)
```

---

## 2. Simulation: `simulate_from_gbsg_dgm()` vs `simulate_from_dgm()`

### 2.1 Parameter Comparison

| Parameter | `simulate_from_gbsg_dgm()` | `simulate_from_dgm()` |
|-----------|---------------------------|----------------------|
| **Sample size** | `n` | `n` |
| **Randomization** | `rand_ratio` | `rand_ratio` |
| **Admin censoring** | `max_follow` | `analysis_time - entry_time` |
| **Censoring adjust** | `muC_adj` | `cens_adjust` |
| **Treatment redraw** | `draw_treatment` | `draw_treatment` |
| **Entry time** | Not supported | `entry_var`, `max_entry` |
| **Seed** | `sim_id` (offset from base) | `seed` |

### 2.2 Output Column Naming

| Concept | GBSG Version | Flex Version |
|---------|--------------|--------------|
| Observed time | `y.sim` | `y_sim` |
| Event indicator | `event.sim` | `event_sim` |
| True time | `t.sim` | `t_true` |
| Treatment | `treat` | `treat_sim` |
| Subgroup flag | `flag.harm` | `flag_harm` |

### 2.3 Censoring Implementation

Both handle Weibull censoring similarly:
```r
# GBSG version
log_C <- mu_cens + sigma_cens * epsilon_cens + df_sim$linC.conf

# Flex version  
logC_sim <- params$censoring$mu + cens_adjust + params$censoring$tau * epsilon_cens + lin_pred_cens
```

---

## 3. Integration Strategy

### 3.1 Option A: Wrapper Functions (Recommended)

Create `create_gbsg_dgm()` as a **thin wrapper** around `generate_aft_dgm_flex()`:

```r
create_gbsg_dgm <- function(
    model = "alt",
    k_treat = 1,
    k_inter = 1,
    z1_quantile = 0.25,
    n_super = 5000,
    cens_type = "weibull",
    seed = 8316951,
    verbose = TRUE
) {
  
  # Load GBSG data
  data(gbsg, package = "survival")
  
  # Call generate_aft_dgm_flex with GBSG-specific defaults
  dgm_flex <- generate_aft_dgm_flex(
    data = gbsg,
    continuous_vars = c("age", "er", "pgr", "nodes", "size"),
    factor_vars = c("meno", "grade"),
    outcome_var = "rfstime",
    event_var = "status",
    treatment_var = "hormon",
    subgroup_vars = c("er", "meno"),
    subgroup_cuts = list(
      er = list(type = "quantile", value = z1_quantile),
      meno = 0
    ),
    model = model,
    k_treat = k_treat,
    k_inter = k_inter,
    n_super = n_super,
    cens_type = cens_type,
    seed = seed,
    verbose = verbose
  )
  
  # Convert to GBSG format with backward-compatible column names
  dgm_gbsg <- convert_to_gbsg_format(dgm_flex, z1_quantile)
  
  return(dgm_gbsg)
}
```

### 3.2 Option B: Column Mapping Layer

Create a mapping layer for `oc_analyses_gbsg.R` to work with `simulate_from_dgm()`:

```r
simulate_from_gbsg_dgm <- function(dgm, n, sim_id, max_follow, muC_adj, ...) {
  
  # If dgm is gbsg_dgm, convert to aft_dgm_flex format
  if (inherits(dgm, "gbsg_dgm")) {
    dgm_flex <- convert_gbsg_to_flex(dgm)
  } else {
    dgm_flex <- dgm
  }
  
  # Call generic simulate_from_dgm
  sim_data <- simulate_from_dgm(
    dgm = dgm_flex,
    n = n,
    analysis_time = max_follow,
    max_entry = 0,  # No staggered entry
    cens_adjust = muC_adj,
    seed = sim_id
  )
  
  # Rename columns to GBSG convention
  sim_data <- rename_columns_gbsg(sim_data)
  
  return(sim_data)
}
```

---

## 4. Required Changes for Full Integration

### 4.1 Changes to `generate_aft_dgm_flex()`

1. **Add backward-compatible HR accessors**:
```r
# In assemble_results():
results$hr_H_true <- hr_results$harm_subgroup
results$hr_Hc_true <- hr_results$no_harm_subgroup
results$hr_causal <- hr_results$overall
results$AHR_H_true <- hr_results$AHR_harm
results$AHR_Hc_true <- hr_results$AHR_no_harm
results$AHR <- hr_results$AHR
```

2. **Support for analysis variable aliases**:
```r
# Allow mapping like v1 -> z_er, v2 -> z_meno, etc.
analysis_vars$aliases <- list(v1 = "z_er", v2 = "z_meno", ...)
```

### 4.2 Changes to `simulate_from_dgm()`

1. **Add `max_follow` parameter** (equivalent to admin censoring without staggered entry):
```r
simulate_from_dgm <- function(..., max_follow = Inf) {
  # When max_entry = 0 and max_follow is specified:
  # follow_up <- max_follow (constant for all subjects)
}
```

2. **Support for seed offset pattern**:
```r
simulate_from_dgm <- function(..., seed = NULL, sim_id = NULL, seed_base = 8316951L) {
  if (!is.null(sim_id) && is.null(seed)) {
    seed <- seed_base + 1000L * sim_id
  }
  if (!is.null(seed)) set.seed(seed)
}
```

### 4.3 Changes to `oc_analyses_gbsg.R`

1. **Use `get_dgm_hr()` helper consistently** (already implemented):
```r
get_dgm_hr(dgm, "hr_H")  # Works with both gbsg_dgm and aft_dgm_flex
```

2. **Column name adapter**:
```r
standardize_column_names <- function(df, source_format = "flex") {
  if (source_format == "flex") {
    names(df) <- gsub("y_sim", "y.sim", names(df))
    names(df) <- gsub("event_sim", "event.sim", names(df))
    names(df) <- gsub("flag_harm", "flag.harm", names(df))
    # etc.
  }
  return(df)
}
```

---

## 5. Migration Path

### Phase 1: Add Compatibility Layer (Low Risk)

1. Add `get_dgm_hr()` helper to handle both formats ✓ (done)
2. Add column name standardization functions
3. Ensure `oc_analyses_gbsg.R` works with both DGM types

### Phase 2: Create GBSG Wrapper (Medium Risk)

1. Create `create_gbsg_dgm_v2()` that wraps `generate_aft_dgm_flex()`
2. Create `simulate_from_gbsg_dgm_v2()` that wraps `simulate_from_dgm()`
3. Validate output equivalence with original functions
4. Add deprecation warnings to original functions

### Phase 3: Full Migration (Higher Risk)

1. Replace `create_gbsg_dgm()` with wrapper implementation
2. Replace `simulate_from_gbsg_dgm()` with wrapper implementation
3. Update all vignettes and examples
4. Remove redundant code

---

## 6. Code Reuse Analysis

### Functions that can be shared:

| Function | Current Location | Reuse Potential |
|----------|------------------|-----------------|
| `calculate_hazard_ratios()` | generate_aft_dgm_helpers.R | High - identical logic |
| `compute_ahr()` | oc_analyses_gbsg.R | High - can use flex version |
| `calibrate_k_inter()` | sim_aft_gbsg_refactored.R | Medium - wrap `find_k_inter_for_target_hr()` |
| `validate_k_inter_effect()` | sim_aft_gbsg_refactored.R | Medium - GBSG-specific but could generalize |

### Functions that remain GBSG-specific:

| Function | Reason |
|----------|--------|
| `print.gbsg_dgm()` | Display format specific to GBSG output |
| `create_summary_table()` | GBSG variable labels |

---

## 7. Recommended Implementation

### 7.1 New File: `R/gbsg_dgm_wrapper.R`

```r
#' Create GBSG-Based DGM (Wrapper)
#'
#' Wrapper around generate_aft_dgm_flex() with GBSG-specific defaults.
#' Provides backward compatibility with existing GBSG simulation code.
#'
#' @inheritParams create_gbsg_dgm
#' @return Object of class c("gbsg_dgm", "aft_dgm_flex", "list")
#' @export
create_gbsg_dgm_flex <- function(
    model = "alt",
    k_treat = 1,
    k_inter = 1,
    z1_quantile = 0.25,
    n_super = 5000,
    cens_type = "weibull",
    seed = 8316951,
    verbose = TRUE
) {
  
  # Load and prepare GBSG data
  data(gbsg, package = "survival")
  
  # Create DGM using flex infrastructure

  dgm <- generate_aft_dgm_flex(
    data = gbsg,
    continuous_vars = c("age", "er", "pgr", "nodes", "size"),
    factor_vars = c("meno", "grade"),
    outcome_var = "rfstime",
    event_var = "status",
    treatment_var = "hormon",
    subgroup_vars = c("er", "meno"),
    subgroup_cuts = list(
      er = list(type = "quantile", value = z1_quantile),
      meno = 0
    ),
    model = model,
    k_treat = k_treat,
    k_inter = k_inter,
    n_super = n_super,
    cens_type = cens_type,
    select_censoring = TRUE,
    seed = seed,
    verbose = verbose
  )
  
  # Add GBSG-specific accessors for backward compatibility
  dgm$hr_H_true <- dgm$hazard_ratios$harm_subgroup
  dgm$hr_Hc_true <- dgm$hazard_ratios$no_harm_subgroup
  dgm$hr_causal <- dgm$hazard_ratios$overall
  dgm$AHR_H_true <- dgm$hazard_ratios$AHR_harm
  dgm$AHR_Hc_true <- dgm$hazard_ratios$AHR_no_harm
  dgm$AHR <- dgm$hazard_ratios$AHR
  
  # Rename df_super to df_super_rand for backward compatibility
  dgm$df_super_rand <- dgm$df_super
  
  # Add GBSG-style analysis variables (v1-v7 aliases)
  dgm$df_super_rand <- add_gbsg_aliases(dgm$df_super_rand)
  
  # Update class
  class(dgm) <- c("gbsg_dgm", class(dgm))
  
  return(dgm)
}


#' Simulate from GBSG DGM (Wrapper)
#'
#' Wrapper around simulate_from_dgm() with GBSG column naming.
#'
#' @inheritParams simulate_from_gbsg_dgm
#' @return data.frame with GBSG-style column names
#' @export
simulate_from_gbsg_dgm_flex <- function(
    dgm,
    n = NULL,
    rand_ratio = 1,
    sim_id = 1L,
    max_follow = Inf,
    muC_adj = 0,
    draw_treatment = TRUE
) {
  
  # Set seed using GBSG convention
  set.seed(8316951L + 1000L * sim_id)
  
  # Call generic simulation
  sim_data <- simulate_from_dgm(
    dgm = dgm,
    n = n,
    rand_ratio = rand_ratio,
    max_entry = 0,
    analysis_time = max_follow,
    cens_adjust = muC_adj,
    draw_treatment = draw_treatment
  )
  
  # Rename columns to GBSG convention
  sim_data <- rename_to_gbsg_format(sim_data)
  
  return(sim_data)
}


# Helper: Add v1-v7 aliases
add_gbsg_aliases <- function(df) {
  var_map <- c(
    v1 = "z_er", v2 = "z_meno", v3 = "z_age",
    v4 = "z_pgr", v5 = "z_nodes", v6 = "z_size", v7 = "z_grade_1"
  )
  
  for (alias in names(var_map)) {
    source_col <- var_map[alias]
    if (source_col %in% names(df)) {
      df[[alias]] <- df[[source_col]]
    }
  }
  
  return(df)
}


# Helper: Rename columns to GBSG format
rename_to_gbsg_format <- function(df) {
  col_map <- c(
    "y_sim" = "y.sim",
    "event_sim" = "event.sim",
    "t_true" = "t.sim",
    "treat_sim" = "treat",
    "flag_harm" = "flag.harm"
  )
  
  for (old in names(col_map)) {
    if (old %in% names(df)) {
      names(df)[names(df) == old] <- col_map[old]
    }
  }
  
  return(df)
}
```

---

## 8. Summary

| Component | Current State | Recommended Action |
|-----------|--------------|-------------------|
| `create_gbsg_dgm()` | Standalone implementation | Wrap `generate_aft_dgm_flex()` |
| `simulate_from_gbsg_dgm()` | Standalone implementation | Wrap `simulate_from_dgm()` |
| `calibrate_k_inter()` | GBSG-specific | Wrap `find_k_inter_for_target_hr()` |
| `oc_analyses_gbsg.R` | Uses GBSG column names | Add column standardization layer |
| `run_simulation_analysis()` | Calls GBSG-specific functions | Use `get_dgm_hr()` helper (done) |

**Benefits of Integration**:
1. **Code Reuse**: ~60% reduction in DGM-related code
2. **Consistency**: Single source of truth for HR calculations
3. **Features**: GBSG users get access to spline treatment effects, lognormal censoring
4. **Maintenance**: Bug fixes apply to all use cases
5. **Testing**: Unified test suite for DGM infrastructure
