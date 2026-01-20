# ForestSearch Package Changes to Implement

## Summary

This document outlines all changes discussed and implemented during this session. Two files have been modified:

1. `plot_subgroup_results_forestplot.R`
2. `subgroup_consistency_main.R`

---

## File 1: `plot_subgroup_results_forestplot.R`

### 1.1 Add `conf.level` Parameter for Configurable Confidence Intervals

**Location:** Function signature and documentation

**Changes:**
- Add parameter `conf.level = 0.95` to function signature
- Add documentation for the parameter
- Calculate `z_alpha <- qnorm(1 - (1 - conf.level) / 2)` after input validation
- Create dynamic CI column label: `ci_label <- sprintf("HR (%d%% CI)", round(conf.level * 100))`
- Replace all hardcoded `1.96` with `z_alpha`
- Update `summary(fit)` calls to use `summary(fit, conf.int = conf.level)`

**Affected functions:**
- `plot_subgroup_results_forestplot()` - main function
- `create_hr_row()` - nested helper function
- `compute_sg_hr()` - exported helper function
- `create_subgroup_summary_df()` - exported helper function

---

### 1.2 Update Cox Model Estimation to Use Robust Standard Errors

**Rationale:** Align with `cox_summary()` used in `analyze_subgroup()` for consistency across the package.

**Changes to `create_hr_row()`:**

```r
# BEFORE:
fit <- tryCatch(
  survival::coxph(cox.formula, data = dfa),
  error = function(e) { ... }
)

# AFTER:
fit <- tryCatch(
  survival::coxph(
    cox.formula, 
    data = dfa,
    robust = TRUE,
    model = FALSE,
    x = FALSE,
    y = FALSE
  ),
  error = function(e) { ... }
)
```

**Changes to `compute_sg_hr()`:**

```r
# BEFORE:
fit <- survival::coxph(cox.formula, data = df)

# AFTER:
fit <- tryCatch({
  survival::coxph(
    cox.formula, 
    data = df,
    robust = TRUE,
    model = FALSE,
    x = FALSE,
    y = FALSE
  )
}, error = function(e) {
  warning("Subgroup '", sg_name, "': Cox model failed: ", e$message)
  return(NULL)
})
```

---

### 1.3 Add Input Validation to Cox Model Functions

**Changes to `create_hr_row()`:**

Add after sample size check:
```r
# Extract vectors for validation
Y <- dfa[[outcome.name]]
E <- dfa[[event.name]]
Treat <- dfa[[treat.name]]

# Check for sufficient events
n_events <- sum(E)
if (n_events < 2) {
  warning(paste("Subgroup", sg_name, "has fewer than 2 events"))
  return(NULL)
}

# Check treatment variation
if (length(unique(Treat)) < 2) {
  warning(paste("Subgroup", sg_name, "has no variation in treatment"))
  return(NULL)
}
```

**Changes to `compute_sg_hr()`:**

Add input validation section:
```r
# Extract vectors
Y <- df[[outcome.name]]
E <- df[[event.name]]
Treat <- df[[treat.name]]

# Check for sufficient events
n_events <- sum(E)
if (n_events < 2) {
  warning("Subgroup '", sg_name, "': Fewer than 2 events; returning NA")
  ntreat <- sum(Treat)
  ncontrol <- sum(1 - Treat)
  return(data.frame(
    Subgroup = sg_name,
    n_treat = ntreat,
    n_control = ncontrol,
    est = NA_real_,
    low = NA_real_,
    hi = NA_real_,
    se = NA_real_
  ))
}

# Check treatment variation
if (length(unique(Treat)) < 2) {
  warning("Subgroup '", sg_name, "': No variation in treatment; returning NA")
  # ... similar NA return
}
```

---

### 1.4 Add Edge Case Handling for Confidence Interval Extraction

**Changes to both `create_hr_row()` and `compute_sg_hr()`:**

```r
# BEFORE:
hr <- summary(fit, conf.int = conf.level)$conf.int[c(1, 3, 4)]

# AFTER:
fit_summary <- summary(fit, conf.int = conf.level)
conf_int <- fit_summary$conf.int

# Handle edge case where conf.int might not exist
if (is.null(conf_int) || nrow(conf_int) == 0) {
  warning("Subgroup '", sg_name, "': No confidence interval available")
  return(NULL)  # or return data frame with NA values
}

hr <- conf_int[1, c(1, 3, 4)]
```

---

## File 2: `subgroup_consistency_main.R`

### 2.1 Add `seed` Parameter for Reproducible Consistency Splits

**Location:** `subgroup.consistency()` function

**Changes to documentation:**

```r
#' @param seed Integer. Random seed for reproducible consistency splits
#'   (default: 8316951). Set to NULL for non-reproducible random splits.
```

**Changes to function signature:**

```r
# BEFORE:
subgroup.consistency <- function(df,
                                 ...
                                 twostage_args = list()) {

# AFTER:
subgroup.consistency <- function(df,
                                 ...
                                 twostage_args = list(),
                                 seed = 8316951) {
```

**Add set.seed() call after input validation:**

```r
# Set random seed for reproducible splits
if (!is.null(seed)) {
  set.seed(seed)
  if (details) {
    cat("Random seed set to:", seed, "\n")
  }
}
```

**Update future_lapply for parallel execution:**

```r
# BEFORE:
results_list <- future.apply::future_lapply(
  seq_len(n_candidates),
  eval_fun,
  future.seed = TRUE,
  ...
)

# AFTER:
results_list <- future.apply::future_lapply(
  seq_len(n_candidates),
  eval_fun,
  future.seed = if (!is.null(seed)) seed else TRUE,
  ...
)
```

**Add seed to output list:**

```r
output <- list(
  out_sg = out_sg,
  sg_focus = sg_focus,
  df_flag = df_flag,
  sg.harm = sg.harm,
  sg.harm.id = sg.harm.id,
  algorithm = ifelse(use_twostage, "twostage", "fixed"),
  n_candidates_evaluated = n_candidates,
  n_passed = n_passed,
  seed = seed  # NEW
)
```

**Update @return documentation:**

```r
#'     \item{seed}{Integer seed used for reproducibility (NULL if not set)}
```

---

## File 3: `subgroup_consistency_helpers.R` (Discussed but NOT yet implemented)

### 3.1 Differentiate `hrMaxSG`/`hrMinSG` from `maxSG`/`minSG` Sorting

**Location:** `sort_subgroups()` function

**Proposed change:**

```r
# CURRENT:
sort_subgroups <- function(result_new, sg_focus) {
  if (sg_focus == "hr") data.table::setorder(result_new, -Pcons, -hr, K)
  if (sg_focus %in% c("hrMaxSG", "maxSG")) data.table::setorder(result_new, -N, -Pcons, K)
  if (sg_focus %in% c("hrMinSG", "minSG")) data.table::setorder(result_new, N, -Pcons, K)
  result_new
}

# PROPOSED:
sort_subgroups <- function(result_new, sg_focus) {
  if (sg_focus == "hr") data.table::setorder(result_new, -Pcons, -hr, K)
  if (sg_focus == "maxSG") data.table::setorder(result_new, -N, -Pcons, K)
  if (sg_focus == "hrMaxSG") data.table::setorder(result_new, -N, -hr, K)
  if (sg_focus == "minSG") data.table::setorder(result_new, N, -Pcons, K)
  if (sg_focus == "hrMinSG") data.table::setorder(result_new, N, -hr, K)
  result_new
}
```

**Updated sorting criteria:**

| sg_focus | Primary Sort | Secondary Sort | Tertiary Sort |
|----------|--------------|----------------|---------------|
| `"hr"` | Consistency (↓) | Hazard Ratio (↓) | # Factors (↑) |
| `"maxSG"` | Sample Size (↓) | Consistency (↓) | # Factors (↑) |
| `"hrMaxSG"` | Sample Size (↓) | **Hazard Ratio (↓)** | # Factors (↑) |
| `"minSG"` | Sample Size (↑) | Consistency (↓) | # Factors (↑) |
| `"hrMinSG"` | Sample Size (↑) | **Hazard Ratio (↓)** | # Factors (↑) |

---

## Files Ready for Download

The following modified files are available in `/mnt/user-data/outputs/`:

1. **`plot_subgroup_results_forestplot.R`** - All changes from sections 1.1-1.4 implemented
2. **`subgroup_consistency_main.R`** - All changes from section 2.1 implemented

---

## Testing Recommendations

After implementing these changes:

1. **Test conf.level parameter:**
   ```r
   # Should show "HR (90% CI)" in column header
   plot_subgroup_results_forestplot(..., conf.level = 0.90)
   ```

2. **Test reproducibility of consistency splits:**
   ```r
   # These should produce identical results
   result1 <- subgroup.consistency(..., seed = 12345)
   result2 <- subgroup.consistency(..., seed = 12345)
   identical(result1$out_sg, result2$out_sg)  # Should be TRUE
   ```

3. **Test robust SE produces consistent results with analyze_subgroup():**
   ```r
   # Compare HR estimates from forest plot vs SG_tab_estimates
   ```

4. **Test edge cases:**
   - Subgroups with < 2 events
   - Subgroups with no treatment variation
   - Very small subgroups (< 10 observations)
