# ForestSearch Package Architecture

## Overview

ForestSearch is a comprehensive R package for exploratory subgroup identification in clinical trials with survival endpoints. The package implements advanced methods from León et al. (2024) published in *Statistics in Medicine*.

---

## Function Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           MAIN ENTRY POINT                                   │
│                                                                              │
│                         ┌──────────────────────┐                            │
│                         │   forestsearch()     │                            │
│                         │   (Main Algorithm)   │                            │
│                         └──────────┬───────────┘                            │
└────────────────────────────────────┼────────────────────────────────────────┘
                                     │
        ┌────────────────────────────┼────────────────────────────┐
        │                            │                            │
        ▼                            ▼                            ▼
┌───────────────┐          ┌───────────────┐          ┌───────────────────┐
│   VARIABLE    │          │     DATA      │          │    SUBGROUP       │
│   SELECTION   │          │  PREPARATION  │          │     SEARCH        │
│               │          │               │          │                   │
│ • GRF         │          │ • get_FSdata  │          │ • subgroup.search │
│ • LASSO       │          │ • Cuts        │          │ • Combinations    │
└───────┬───────┘          └───────┬───────┘          └─────────┬─────────┘
        │                          │                            │
        └──────────────────────────┼────────────────────────────┘
                                   │
                                   ▼
                    ┌──────────────────────────────┐
                    │    CONSISTENCY EVALUATION    │
                    │                              │
                    │  • subgroup.consistency()    │
                    │  • Split-sample validation   │
                    │  • Two-stage algorithm       │
                    └──────────────┬───────────────┘
                                   │
        ┌──────────────────────────┼──────────────────────────┐
        │                          │                          │
        ▼                          ▼                          ▼
┌───────────────┐          ┌───────────────┐          ┌───────────────┐
│   BOOTSTRAP   │          │    CROSS-     │          │    OUTPUT     │
│   ANALYSIS    │          │  VALIDATION   │          │  & SUMMARY    │
│               │          │               │          │               │
│ • Bias corr.  │          │ • K-fold      │          │ • Tables      │
│ • Inference   │          │ • Ten-fold    │          │ • Plots       │
└───────────────┘          └───────────────┘          └───────────────┘
```

---

## 1. Main Entry Point

### `forestsearch()`
**Location:** `R/forest_search_revised.R`

The main orchestrating function that coordinates all analysis steps.

**Key Parameters:**
| Parameter | Description | Default |
|-----------|-------------|---------|
| `df.analysis` | Analysis dataset | Required |
| `confounders.name` | Candidate variables | NULL (auto-detect) |
| `use_lasso` | Enable LASSO selection | TRUE |
| `use_grf` | Enable GRF for cuts | TRUE |
| `hr.threshold` | Minimum HR for candidates | 1.25 |
| `pconsistency.threshold` | Consistency requirement | 0.90 |
| `fs.splits` | Number of validation splits | 1000 |
| `maxk` | Max factors per subgroup | 2 |
| `use_twostage` | Two-stage algorithm | FALSE |

**Returns:** Object of class `"forestsearch"` containing:
- `grp.consistency` - Consistency evaluation results
- `find.grps` - Subgroup search results
- `sg.harm` - Selected subgroup definition
- `df.est` - Data with treatment recommendations

---

## 2. Variable Selection Module

### `grf.subg.harm.survival()`
**Location:** `R/grf_main.R`

Uses Generalized Random Forests (causal survival forest) to identify variables with treatment effect heterogeneity.

**Helper Functions:**
| Function | Purpose |
|----------|---------|
| `create_grf_config()` | Initialize GRF parameters |
| `fit_causal_forest()` | Wrapper for `grf::causal_survival_forest()` |
| `compute_node_metrics()` | Calculate metrics per tree node |
| `select_best_subgroup()` | Choose optimal subgroup based on criterion |
| `extract_all_tree_cuts()` | Extract cut expressions from trees |
| `fit_policy_trees()` | Fit policy trees at each depth |
| `assign_subgroup_membership()` | Assign observations to subgroups |

### `lasso_selection()`
**Location:** `R/get_FSdata_helpers.R`

Uses Cox-LASSO for dimension reduction.

```r
# Returns:
list(
  selected = character(),   # Variables with non-zero coefficients
  omitted = character(),    # Variables excluded
  coefficients = numeric(), # LASSO coefficients
  lambda_min = numeric(),   # Optimal lambda
  cvfit = object,          # cv.glmnet result
  fit = object             # Final glmnet fit
)
```

---

## 3. Data Preparation Module

### `get_FSdata()`
**Location:** `R/get_FSdata_refactored.R`

Prepares the dataset for ForestSearch by creating binary cut indicators.

**Processing Steps:**
1. Classify variables as continuous/categorical
2. Apply LASSO selection (if enabled)
3. Apply GRF cuts (if provided)
4. Create binary indicator columns
5. Handle forced cuts and exclusions

**Helper Functions:**
| Function | Purpose |
|----------|---------|
| `get_conf_force()` | Generate forced cut expressions |
| `evaluate_cuts_once()` | Evaluate all cuts and cache results |
| `FS_labels()` | Convert q-codes to readable labels |
| `filter_by_lassokeep()` | Filter by LASSO-selected variables |
| `dummy()` | Create dummy variables |
| `is.continuous()` | Check if variable is continuous |
| `process_conf_force_expr()` | Process forced cut expressions |

---

## 4. Subgroup Search Module

### `subgroup.search()`
**Location:** `R/subgroup_search.R`

Exhaustive search over all factor combinations up to `maxk`.

**Algorithm:**
```
1. Generate all combinations of 1 to maxk factors
2. For each combination:
   a. Check minimum prevalence (minp)
   b. Check minimum subgroup size (n.min)
   c. Check minimum events per arm (d0.min, d1.min)
   d. Fit Cox model and compute HR
   e. Store if HR > hr.threshold
3. Return sorted candidates
```

**Helper Functions:**
| Function | Purpose |
|----------|---------|
| `prepare_search_data()` | Clean data, remove missing values |
| `generate_combination_indices()` | Create all k-combinations |
| `search_combinations()` | Main loop over combinations |
| `evaluate_combination()` | Test single combination |
| `get_covs_in()` | Get factor indicators for combination |
| `get_subgroup_membership()` | Compute subgroup membership |
| `calculate_max_combinations()` | Count total combinations |
| `format_search_results()` | Format output as data.table |
| `fit_subgroup_cox()` | Fit Cox model for subgroup |
| `is_time_exceeded()` | Check time limit |

---

## 5. Consistency Evaluation Module

### `subgroup.consistency()`
**Location:** `R/subgroup_consistency_main.R`

Split-sample validation to ensure reproducibility of identified subgroups.

**Two Algorithms:**

**Fixed Algorithm (default):**
- Run exactly `n.splits` random splits
- Each split: 50/50 random assignment
- Consistency = proportion where both halves show HR > threshold

**Two-Stage Algorithm (`use_twostage = TRUE`):**
```
Stage 1: Quick screening (n.splits.screen splits)
         Eliminate clearly non-viable candidates
         
Stage 2: Sequential batched evaluation
         Early stopping when outcome determined
         Uses Wilson confidence intervals
```

**Helper Functions:**
| Function | Purpose |
|----------|---------|
| `evaluate_subgroup_consistency()` | Evaluate single subgroup (fixed) |
| `evaluate_consistency_twostage()` | Evaluate with early stopping |
| `get_split_hr_fast()` | Fast HR calculation for split |
| `wilson_ci()` | Wilson score confidence interval |
| `early_stop_decision()` | Determine if can stop early |
| `sort_subgroups()` | Sort by sg_focus criterion |
| `extract_subgroup()` | Extract subgroup definition |
| `sg_consistency_out()` | Format consistency output |
| `remove_near_duplicate_subgroups()` | Remove redundant subgroups |
| `identify_near_duplicates()` | Find near-identical rows |
| `setup_consistency_parallel()` | Configure parallel processing |

---

## 6. Bootstrap Analysis Module

### `forestsearch_bootstrap_dofuture()`
**Location:** `R/bootstrap_dofuture_main.R`

Bootstrap inference and bias correction using infinitesimal jackknife methods.

**Bias Correction Methods:**
```
Method 1 (Simple Optimism):
  H_adj1 = H_obs - (H*_* - H*_obs)

Method 2 (Double Bootstrap):
  H_adj2 = 2 × H_obs - (H_* + H*_* - H*_obs)
```

**Helper Functions:**
| Function | Location | Purpose |
|----------|----------|---------|
| `bootstrap_results()` | `R/bootstrap_analysis_dofuture.R` | Coordinate bootstrap iterations |
| `run_single_bootstrap()` | `R/bootstrap_analysis_dofuture.R` | Execute one iteration |
| `summarize_bootstrap_results()` | `R/bootstrap_summaries_helpers.R` | Aggregate results |
| `create_bootstrap_timing_plots()` | `R/bootstrap_summaries_helpers.R` | Timing visualizations |
| `summarize_bootstrap_subgroups()` | `R/summarize_bootstrap_subgroups.R` | Subgroup agreement |
| `resolve_bootstrap_parallel_args()` | `R/bootstrap_dofuture_main.R` | Configure parallelization |

---

## 7. Cross-Validation Module

### `forestsearch_Kfold()`
**Location:** `R/forestsearch_cross-validation.R`

K-fold cross-validation for subgroup stability assessment.

### `forestsearch_tenfold()`
**Location:** `R/forestsearch_cross-validation.R`

Repeated K-fold cross-validation (wrapper).

**Helper Functions:**
| Function | Purpose |
|----------|---------|
| `forestsearch_KfoldOut()` | Summarize K-fold results |
| `find_covariate_any_match()` | Match covariates across folds |
| `resolve_cv_parallel_args()` | Configure CV parallelization |

---

## 8. Simulation & DGM Module

### Data Generating Mechanism

| Function | Location | Purpose |
|----------|----------|---------|
| `generate_aft_dgm_flex()` | `R/generate_aft_dgm_flex.R` | Create AFT data generating mechanism |
| `simulate_from_dgm()` | `R/simulate_from_dgm.R` | Generate simulated data from DGM |

### Simulation Functions

| Function | Location | Purpose |
|----------|----------|---------|
| `run_mrct_simulation()` | `R/mrct_simulation.R` | Multi-regional CT simulation |
| `run_single_simulation()` | `R/oc_analyses_gbsg.R` | Single simulation iteration |
| `run_forestsearch_analysis()` | `R/oc_analyses_gbsg.R` | FS analysis within simulation |
| `generate_bootstrap_synthetic()` | `R/synthetic_data_perturbation.R` | Synthetic data generation |

---

## 9. Output & Summary Module

### S3 Methods
| Function | Purpose |
|----------|---------|
| `print.forestsearch()` | Print summary of results |
| `summary.forestsearch()` | Detailed summary statistics |

### Visualization Functions
| Function | Location | Purpose |
|----------|----------|---------|
| `sg_tables()` | `R/summary_utility_functions.R` | gt-formatted summary tables |
| `plot_subgroup_results_forestplot()` | `R/plot_subgroup_results_forestplot.R` | Forest plots with CV metrics |
| `cox_cs_fit()` | `R/cox_spline_fit.R` | Cox model with spline interactions |

### Utility Functions
| Function | Location | Purpose |
|----------|----------|---------|
| `filter_call_args()` | `R/summary_utility_functions.R` | Filter function arguments |
| `get_dfpred()` | `R/subgroup_consistency_helpers.R` | Get prediction dataset |
| `cox_summary()` | `R/summary_utility_functions.R` | Summarize Cox model |

---

## 10. Parallel Processing Infrastructure

ForestSearch supports multiple parallel backends:

| Backend | Use Case |
|---------|----------|
| `"sequential"` | Single-threaded (debugging) |
| `"multisession"` | Cross-platform parallel |
| `"multicore"` | Unix/Mac fork-based |
| `"callr"` | Separate R processes |

**Key Parallel Functions:**
- `setup_consistency_parallel()` - For consistency evaluation
- `resolve_bootstrap_parallel_args()` - For bootstrap analysis
- `resolve_cv_parallel_args()` - For cross-validation

---

## Dependencies

### Core Statistical
- `survival` - Cox models, survival objects
- `grf` - Generalized Random Forests
- `glmnet` - LASSO regularization
- `policytree` - Policy tree algorithms

### Data Manipulation
- `data.table` - Fast data operations

### Parallel Computing
- `future` - Parallel framework
- `doFuture` - foreach integration
- `foreach` - Parallel loops
- `callr` - R subprocess execution

### Visualization
- `ggplot2` - Plotting
- `gt` - Table formatting
- `forestploter` - Forest plots

---

## Typical Workflow

```r
# 1. Run main analysis
fs_result <- forestsearch(
  df.analysis = trial_data,
  confounders.name = c("age", "biomarker", "stage"),
  hr.threshold = 1.25,
  pconsistency.threshold = 0.90,
  use_grf = TRUE,
  use_lasso = TRUE,
  details = TRUE
)

# 2. Bootstrap bias correction
fs_boot <- forestsearch_bootstrap_dofuture(
  fs.est = fs_result,
  nb_boots = 1000,
  parallel_args = list(plan = "callr", workers = 6)
)

# 3. Cross-validation
fs_cv <- forestsearch_Kfold(
  fs.est = fs_result,
  Kfolds = 10
)

# 4. Summary and visualization
print(fs_result)
sg_tables(fs_result)
```

---

## File Organization

```
R/
├── forest_search_revised.R      # Main forestsearch() function
├── get_FSdata_refactored.R      # Data preparation
├── get_FSdata_helpers.R         # Data prep helpers
├── subgroup_search.R            # Exhaustive search
├── subgroup_consistency_main.R  # Consistency evaluation
├── subgroup_consistency_helpers.R # Consistency helpers
├── grf_main.R                   # GRF subgroup identification
├── grf_helpers.R                # GRF helper functions
├── bootstrap_dofuture_main.R    # Bootstrap main function
├── bootstrap_analysis_dofuture.R # Bootstrap iteration logic
├── bootstrap_summaries_helpers.R # Bootstrap summarization
├── forestsearch_cross-validation.R # K-fold CV
├── summary_utility_functions.R  # Output utilities
├── plot_subgroup_results_forestplot.R # Forest plots
├── generate_aft_dgm_flex.R      # DGM generation
├── simulate_from_dgm.R          # Data simulation
└── ...
```

---

*Document generated for ForestSearch R package - CRAN submission preparation*
