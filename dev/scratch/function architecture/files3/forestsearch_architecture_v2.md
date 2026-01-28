# ForestSearch Package Architecture

## Function Relationships and Module Organization

**ForestSearch Development Team**

---

## Overview

ForestSearch is a comprehensive R package for exploratory subgroup identification in clinical trials with survival endpoints. The package implements advanced methods from León et al. (2024) published in *Statistics in Medicine*.

**Key Methods:**

- Generalized Random Forests (GRF) for variable importance
- LASSO regularization for dimension reduction
- Exhaustive combinatorial search for subgroup discovery
- Split-sample validation for consistency assessment
- Bootstrap bias correction using infinitesimal jackknife methods

---

## High-Level Workflow

The ForestSearch algorithm proceeds through five main phases:

```mermaid
flowchart LR
    A[Input Data] --> B[Variable Selection]
    B --> C[Data Preparation]
    C --> D[Subgroup Search]
    D --> E[Consistency Evaluation]
    E --> F[Output]
```

After the core analysis, optional validation steps include:

```mermaid
flowchart LR
    A[ForestSearch Result] --> B[Bootstrap Bias Correction]
    A --> C[K-Fold Cross-Validation]
    B --> D[Final Report]
    C --> D
```

---

## Phase 1: Main Entry Point

### `forestsearch()` - The Orchestrator

**Location:** `R/forest_search_revised.R`

```mermaid
flowchart TB
    FS[[forestsearch]]
    FS --> P1[1. Validate Inputs]
    P1 --> P2[2. GRF Variable Selection]
    P2 --> P3[3. LASSO Selection]
    P3 --> P4[4. Data Preparation]
    P4 --> P5[5. Subgroup Search]
    P5 --> P6[6. Consistency Evaluation]
    P6 --> P7[7. Return Results]
```

**Key Parameters:**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `df.analysis` | Analysis dataset | Required |
| `confounders.name` | Candidate variables | NULL |
| `use_lasso` | Enable LASSO | TRUE |
| `use_grf` | Enable GRF | TRUE |
| `hr.threshold` | Min HR for candidates | 1.25 |
| `pconsistency.threshold` | Consistency requirement | 0.90 |
| `fs.splits` | Validation splits | 1000 |
| `maxk` | Max factors per subgroup | 2 |

---

## Phase 2: Variable Selection

### GRF-Based Selection

**Function:** `grf.subg.harm.survival()`  
**Location:** `R/grf_main.R`

```mermaid
flowchart TB
    A[Input Data] --> B[fit_causal_forest]
    B --> C[Variable Importance]
    C --> D[Fit Policy Trees]
    D --> E[select_best_subgroup]
    E --> F[extract_all_tree_cuts]
    F --> G[GRF Cuts Output]
```

**GRF Helper Functions:**

| Function | Purpose |
|----------|---------|
| `create_grf_config()` | Initialize parameters |
| `fit_causal_forest()` | Fit causal survival forest |
| `compute_node_metrics()` | Metrics per tree node |
| `select_best_subgroup()` | Choose optimal subgroup |
| `extract_all_tree_cuts()` | Get cut expressions |

### LASSO-Based Selection

**Function:** `lasso_selection()`  
**Location:** `R/get_FSdata_helpers.R`

```mermaid
flowchart TB
    A[Confounders] --> B[cv.glmnet Cox-LASSO]
    B --> C[Select lambda.min]
    C --> D[Extract Coefficients]
    D --> E{Coef != 0?}
    E -->|Yes| F[Selected]
    E -->|No| G[Omitted]
```

---

## Phase 3: Data Preparation

### `get_FSdata()` - Creating Binary Indicators

**Location:** `R/get_FSdata_refactored.R`

```mermaid
flowchart TB
    A[Raw Data] --> B[Classify Variables]
    B --> C{Continuous?}
    C -->|Yes| D[Create Cut Indicators]
    C -->|No| E[Create Dummy Variables]
    D --> F[Apply LASSO Filter]
    E --> F
    F --> G[Apply GRF Cuts]
    G --> H[Handle Forced Cuts]
    H --> I[Binary Matrix Z]
```

**Data Preparation Helpers:**

| Function | Purpose |
|----------|---------|
| `is.continuous()` | Check variable type |
| `dummy()` | Create dummy variables |
| `get_conf_force()` | Generate forced cuts |
| `evaluate_cuts_once()` | Evaluate & cache cuts |
| `FS_labels()` | Convert q-codes to labels |
| `filter_by_lassokeep()` | Apply LASSO filter |

---

## Phase 4: Subgroup Search

### `subgroup.search()` - Exhaustive Combinatorial Search

**Location:** `R/subgroup_search.R`

```mermaid
flowchart TB
    A[Generate Combinations] --> B[Loop Each Combination]
    B --> C{Prevalence >= minp?}
    C -->|No| B
    C -->|Yes| D{Size >= n.min?}
    D -->|No| B
    D -->|Yes| E{Events OK?}
    E -->|No| B
    E -->|Yes| F[Fit Cox Model]
    F --> G{HR > threshold?}
    G -->|No| B
    G -->|Yes| H[Store Candidate]
    H --> B
    B -->|Done| I[Return Sorted Candidates]
```

**Search Helper Functions:**

| Function | Purpose |
|----------|---------|
| `prepare_search_data()` | Clean data |
| `generate_combination_indices()` | Create k-combinations |
| `search_combinations()` | Main search loop |
| `evaluate_combination()` | Test one combination |
| `get_covs_in()` | Factor indicators |
| `get_subgroup_membership()` | Membership vector |
| `format_search_results()` | Format as data.table |

---

## Phase 5: Consistency Evaluation

### `subgroup.consistency()` - Split-Sample Validation

**Location:** `R/subgroup_consistency_main.R`

#### Fixed Algorithm (Default)

```mermaid
flowchart TB
    A[Candidate Subgroup] --> B[Repeat n.splits times]
    B --> C[Random 50/50 Split]
    C --> D[HR in Split 1]
    C --> E[HR in Split 2]
    D --> F{Both HR > threshold?}
    E --> F
    F -->|Yes| G[Consistent]
    F -->|No| H[Inconsistent]
    G --> B
    H --> B
    B -->|Done| I[Pcons = Consistent/Total]
```

#### Two-Stage Algorithm (Optional)

```mermaid
flowchart TB
    A[Candidate] --> B[Stage 1: Quick Screen]
    B --> C{Pass threshold?}
    C -->|No| D[REJECT]
    C -->|Yes| E[Stage 2: Batched]
    E --> F[Run 20 splits]
    F --> G[Wilson CI]
    G --> H{CI lower > threshold?}
    H -->|Yes| I[ACCEPT]
    H -->|No| J{CI upper < threshold?}
    J -->|Yes| D
    J -->|No| E
```

**Consistency Helper Functions:**

| Function | Purpose |
|----------|---------|
| `evaluate_subgroup_consistency()` | Fixed algorithm |
| `evaluate_consistency_twostage()` | Two-stage algorithm |
| `get_split_hr_fast()` | Fast HR for split |
| `wilson_ci()` | Wilson confidence interval |
| `early_stop_decision()` | Early stopping logic |
| `sort_subgroups()` | Sort by sg_focus |
| `extract_subgroup()` | Get subgroup definition |

---

## Phase 6: Bootstrap Bias Correction

### `forestsearch_bootstrap_dofuture()`

**Location:** `R/bootstrap_dofuture_main.R`

```mermaid
flowchart TB
    A[Original Result] --> B[Generate B Bootstrap Samples]
    B --> C[For each bootstrap]
    C --> D[Run forestsearch]
    D --> E[Compute HRs]
    E --> C
    C -->|Done| F[Aggregate Results]
    F --> G[Bias Correction]
    G --> H[Adjusted HR]
```

**Bias Correction Formulas:**

| Method | Formula |
|--------|---------|
| Simple Optimism | H_adj1 = H_obs - (H*_* - H*_obs) |
| Double Bootstrap | H_adj2 = 2×H_obs - (H_* + H*_* - H*_obs) |

**Bootstrap Helper Functions:**

| Function | Location |
|----------|----------|
| `bootstrap_results()` | Coordinate iterations |
| `run_single_bootstrap()` | One bootstrap iteration |
| `summarize_bootstrap_results()` | Aggregate statistics |
| `summarize_bootstrap_subgroups()` | Subgroup agreement |

---

## Phase 7: Cross-Validation

### `forestsearch_Kfold()`

**Location:** `R/forestsearch_cross-validation.R`

```mermaid
flowchart TB
    A[Full Dataset] --> B[Split into K Folds]
    B --> C[For fold k = 1 to K]
    C --> D[Train: K-1 folds]
    C --> E[Test: fold k]
    D --> F[forestsearch on training]
    F --> G[Evaluate on test]
    G --> C
    C -->|Done| H[Aggregate Results]
```

**CV Helper Functions:**

| Function | Purpose |
|----------|---------|
| `forestsearch_KfoldOut()` | Summarize K-fold results |
| `forestsearch_tenfold()` | Repeated K-fold wrapper |
| `find_covariate_any_match()` | Match covariates across folds |

---

## Phase 8: Output & Visualization

### S3 Methods and Summary Functions

```mermaid
flowchart TB
    A[forestsearch result] --> B[print.forestsearch]
    A --> C[summary.forestsearch]
    A --> D[sg_tables]
    A --> E[plot_subgroup_results_forestplot]
    B --> F[Console Summary]
    C --> G[Detailed Stats]
    D --> H[gt Tables]
    E --> I[Forest Plot]
```

**Output Functions:**

| Function | Purpose |
|----------|---------|
| `print.forestsearch()` | Basic summary |
| `summary.forestsearch()` | Detailed stats |
| `sg_tables()` | gt-formatted tables |
| `plot_subgroup_results_forestplot()` | Forest plots |
| `cox_cs_fit()` | Spline Cox models |
| `filter_call_args()` | Filter function args |

---

## Simulation & DGM Module

### Data Generation for Simulations

```mermaid
flowchart TB
    A[Define DGM Parameters] --> B[generate_aft_dgm_flex]
    B --> C[simulate_from_dgm]
    C --> D[Simulated Trial Data]
    D --> E[run_forestsearch_analysis]
    E --> F[Evaluate Performance]
```

**Simulation Functions:**

| Function | Location | Purpose |
|----------|----------|---------|
| `generate_aft_dgm_flex()` | `generate_aft_dgm_flex.R` | Create AFT DGM |
| `simulate_from_dgm()` | `simulate_from_dgm.R` | Generate data |
| `run_mrct_simulation()` | `mrct_simulation.R` | MRCT simulation |
| `run_single_simulation()` | `oc_analyses_gbsg.R` | Single iteration |

---

## Parallel Processing

ForestSearch supports multiple parallel backends:

| Backend | Description | Platform |
|---------|-------------|----------|
| `sequential` | Single thread | All |
| `multisession` | Background R sessions | All |
| `multicore` | Fork processes | Unix/Mac |
| `callr` | Separate R processes | All |

---

## Package Dependencies

### Core Statistical
- `survival` - Cox models
- `grf` - Generalized Random Forests
- `glmnet` - LASSO regularization
- `policytree` - Policy trees

### Data Handling
- `data.table` - Fast data operations
- `stringr` - String manipulation

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

## File Organization

```
R/
├── forest_search_revised.R          # Main forestsearch()
├── get_FSdata_refactored.R          # Data preparation
├── get_FSdata_helpers.R             # Data prep helpers
├── subgroup_search.R                # Exhaustive search
├── subgroup_consistency_main.R      # Consistency evaluation
├── subgroup_consistency_helpers.R   # Consistency helpers
├── grf_main.R                       # GRF subgroup ID
├── grf_helpers.R                    # GRF helpers
├── bootstrap_dofuture_main.R        # Bootstrap main
├── bootstrap_analysis_dofuture.R    # Bootstrap iteration
├── bootstrap_summaries_helpers.R    # Bootstrap summary
├── summarize_bootstrap_subgroups.R  # Subgroup summary
├── forestsearch_cross-validation.R  # K-fold CV
├── summary_utility_functions.R      # Output utilities
├── summary_forestsearch.R           # S3 summary
├── plot_subgroup_results_forestplot.R # Forest plots
├── cox_spline_fit.R                 # Spline Cox
├── generate_aft_dgm_flex.R          # DGM generation
├── simulate_from_dgm.R              # Simulation
├── mrct_simulation.R                # MRCT sim
└── oc_analyses_gbsg.R               # OC analyses
```

---

## Typical Workflow Example

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

## References

- León LF, et al. (2024). Exploratory subgroup identification in the heterogeneous Cox model. *Statistics in Medicine*.
- Athey S, Imbens G (2016). Recursive partitioning for heterogeneous causal effects. *PNAS*.
- Wager S, Athey S (2018). Estimation and inference of heterogeneous treatment effects using random forests. *JASA*.

---

*Document generated for ForestSearch R package - CRAN submission preparation*
