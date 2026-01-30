# ForestSearch Package Architecture

## Function Relationships and Module Organization

**ForestSearch Development Team**

---

## Overview

ForestSearch is a comprehensive R package for exploratory subgroup identification in clinical trials with survival endpoints. The package implements advanced methods from León et al. (2024) published in *Statistics in Medicine*.

---

## High-Level Workflow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         FORESTSEARCH PIPELINE                                │
└─────────────────────────────────────────────────────────────────────────────┘

    ┌──────────┐     ┌──────────┐     ┌──────────┐     ┌──────────┐     ┌──────────┐
    │  INPUT   │────▶│ VARIABLE │────▶│   DATA   │────▶│ SUBGROUP │────▶│CONSISTENCY│
    │   DATA   │     │ SELECTION│     │   PREP   │     │  SEARCH  │     │   EVAL   │
    └──────────┘     └──────────┘     └──────────┘     └──────────┘     └─────┬────┘
                                                                              │
                                                                              ▼
                     ┌──────────────────────────────────────────────────────────┐
                     │                      OUTPUT                              │
                     │  • Identified Subgroup (sg.harm)                         │
                     │  • Treatment Recommendations                             │
                     │  • Consistency Metrics                                   │
                     └──────────────────────────────────────────────────────────┘


    ┌─────────────────────────────────────────────────────────────────────────┐
    │                     OPTIONAL VALIDATION STEPS                            │
    └─────────────────────────────────────────────────────────────────────────┘

                              ┌──────────────┐
                         ┌───▶│  BOOTSTRAP   │───┐
                         │    │ BIAS CORRECT │   │
    ┌──────────────┐     │    └──────────────┘   │    ┌──────────────┐
    │  FORESTSEARCH│─────┤                       ├───▶│    FINAL     │
    │    RESULT    │     │    ┌──────────────┐   │    │    REPORT    │
    └──────────────┘     └───▶│   K-FOLD     │───┘    └──────────────┘
                              │     CV       │
                              └──────────────┘
```

---

## Phase 1: Main Entry Point

### `forestsearch()` - The Orchestrator

**Location:** `R/forest_search_revised.R`

```
┌─────────────────────────────────────────────────────────────────┐
│                      forestsearch()                              │
│                   Main Entry Function                            │
└───────────────────────────┬─────────────────────────────────────┘
                            │
        ┌───────────────────┼───────────────────┐
        │                   │                   │
        ▼                   ▼                   ▼
┌───────────────┐   ┌───────────────┐   ┌───────────────┐
│ 1. Validate   │   │ 2. Variable   │   │ 3. Data       │
│    Inputs     │   │    Selection  │   │    Preparation│
└───────┬───────┘   └───────┬───────┘   └───────┬───────┘
        │                   │                   │
        └───────────────────┼───────────────────┘
                            │
        ┌───────────────────┼───────────────────┐
        │                   │                   │
        ▼                   ▼                   ▼
┌───────────────┐   ┌───────────────┐   ┌───────────────┐
│ 4. Subgroup   │   │ 5. Consistency│   │ 6. Return     │
│    Search     │   │    Evaluation │   │    Results    │
└───────────────┘   └───────────────┘   └───────────────┘
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

```
┌─────────────────────────────────────────────────────────────────┐
│                  GRF VARIABLE SELECTION                          │
└─────────────────────────────────────────────────────────────────┘

┌─────────────┐
│ Input Data  │
│ + Confounders│
└──────┬──────┘
       │
       ▼
┌──────────────────────┐
│ fit_causal_forest()  │  ◀── grf::causal_survival_forest()
└──────────┬───────────┘
           │
           ▼
┌──────────────────────┐
│ Compute Variable     │  ◀── grf::variable_importance()
│ Importance           │
└──────────┬───────────┘
           │
           ▼
┌──────────────────────┐
│ Fit Policy Trees     │  ◀── policytree::policy_tree()
│ (depth 1, 2, 3)      │      at each depth level
└──────────┬───────────┘
           │
           ▼
┌──────────────────────┐
│ select_best_subgroup()│  ◀── Choose by criterion (mDiff/Nsg)
└──────────┬───────────┘
           │
           ▼
┌──────────────────────┐
│extract_all_tree_cuts()│  ◀── Get cut expressions
└──────────┬───────────┘
           │
           ▼
┌─────────────────┐
│  GRF Cuts       │
│  (e.g., "age>50")│
└─────────────────┘
```

**GRF Helper Functions:**

| Function | Purpose |
|----------|---------|
| `create_grf_config()` | Initialize parameters |
| `fit_causal_forest()` | Fit causal survival forest |
| `compute_node_metrics()` | Metrics per tree node |
| `select_best_subgroup()` | Choose optimal subgroup |
| `extract_all_tree_cuts()` | Get cut expressions |

---

### LASSO-Based Selection

**Function:** `lasso_selection()`  
**Location:** `R/get_FSdata_helpers.R`

```
┌─────────────────────────────────────────────────────────────────┐
│                  LASSO VARIABLE SELECTION                        │
└─────────────────────────────────────────────────────────────────┘

┌─────────────┐
│ Confounders │
└──────┬──────┘
       │
       ▼
┌──────────────────────┐
│ glmnet::cv.glmnet()  │  ◀── Cox-LASSO with cross-validation
│ family = "cox"       │
└──────────┬───────────┘
           │
           ▼
┌──────────────────────┐
│ Select lambda.min    │  ◀── Optimal regularization
└──────────┬───────────┘
           │
           ▼
┌──────────────────────┐
│ Extract Coefficients │
└──────────┬───────────┘
           │
           ▼
     ┌─────┴─────┐
     │           │
     ▼           ▼
┌─────────┐ ┌─────────┐
│ Coef≠0  │ │ Coef=0  │
│SELECTED │ │ OMITTED │
└─────────┘ └─────────┘
```

---

## Phase 3: Data Preparation

### `get_FSdata()` - Creating Binary Indicators

**Location:** `R/get_FSdata_refactored.R`

```
┌─────────────────────────────────────────────────────────────────┐
│                    DATA PREPARATION                              │
└─────────────────────────────────────────────────────────────────┘

┌─────────────┐
│  Raw Data   │
└──────┬──────┘
       │
       ▼
┌──────────────────────┐
│ Classify Variables   │  ◀── is.continuous()
│ (continuous vs cat)  │
└──────────┬───────────┘
           │
     ┌─────┴─────┐
     │           │
     ▼           ▼
┌─────────┐ ┌─────────┐
│Continuous│ │Categorical│
└────┬────┘ └────┬────┘
     │           │
     ▼           ▼
┌─────────┐ ┌─────────┐
│ Create  │ │ Create  │
│Cut Indic│ │ Dummy   │  ◀── dummy()
└────┬────┘ └────┬────┘
     │           │
     └─────┬─────┘
           │
           ▼
┌──────────────────────┐
│ Apply LASSO Filter   │  ◀── filter_by_lassokeep()
└──────────┬───────────┘
           │
           ▼
┌──────────────────────┐
│ Apply GRF Cuts       │  ◀── Add GRF-identified cuts
└──────────┬───────────┘
           │
           ▼
┌──────────────────────┐
│ Handle Forced Cuts   │  ◀── get_conf_force()
└──────────┬───────────┘
           │
           ▼
┌─────────────────────┐
│ Binary Indicator    │
│ Matrix Z (n × L)    │
└─────────────────────┘
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

```
┌─────────────────────────────────────────────────────────────────┐
│                    SUBGROUP SEARCH ALGORITHM                     │
└─────────────────────────────────────────────────────────────────┘

┌───────────────────────────┐
│ Generate All Combinations │  ◀── generate_combination_indices()
│ (1 to maxk factors)       │      L choose 1, L choose 2, ...
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│   LOOP: Each Combination  │◀─────────────────────────────────┐
└────────────┬──────────────┘                                  │
             │                                                 │
             ▼                                                 │
┌───────────────────────────┐                                  │
│ Check Prevalence ≥ minp   │──── NO ──────────────────────────┤
└────────────┬──────────────┘                                  │
             │ YES                                             │
             ▼                                                 │
┌───────────────────────────┐                                  │
│ Check Size ≥ n.min        │──── NO ──────────────────────────┤
└────────────┬──────────────┘                                  │
             │ YES                                             │
             ▼                                                 │
┌───────────────────────────┐                                  │
│ Check Events:             │                                  │
│   d0 ≥ d0.min (control)   │──── NO ──────────────────────────┤
│   d1 ≥ d1.min (treatment) │                                  │
└────────────┬──────────────┘                                  │
             │ YES                                             │
             ▼                                                 │
┌───────────────────────────┐                                  │
│ Fit Cox Model             │  ◀── survival::coxph()           │
│ Compute HR                │                                  │
└────────────┬──────────────┘                                  │
             │                                                 │
             ▼                                                 │
┌───────────────────────────┐                                  │
│ HR > hr.threshold?        │──── NO ──────────────────────────┤
└────────────┬──────────────┘                                  │
             │ YES                                             │
             ▼                                                 │
┌───────────────────────────┐                                  │
│ *** STORE CANDIDATE ***   │──────────────────────────────────┘
└────────────┬──────────────┘
             │
             ▼ (when loop done)
┌───────────────────────────┐
│ Sort by HR (descending)   │
│ Return Candidate List     │  ◀── format_search_results()
└───────────────────────────┘
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

### Fixed Algorithm (Default)

```
┌─────────────────────────────────────────────────────────────────┐
│              FIXED CONSISTENCY ALGORITHM                         │
└─────────────────────────────────────────────────────────────────┘

┌───────────────────────┐
│  Candidate Subgroup   │
│  (from search phase)  │
└───────────┬───────────┘
            │
            ▼
┌───────────────────────┐
│ REPEAT n.splits times │◀──────────────────────────────┐
└───────────┬───────────┘                               │
            │                                           │
            ▼                                           │
┌───────────────────────┐                               │
│ Random 50/50 Split    │                               │
│ of subgroup members   │                               │
└───────────┬───────────┘                               │
            │                                           │
      ┌─────┴─────┐                                     │
      │           │                                     │
      ▼           ▼                                     │
┌──────────┐ ┌──────────┐                               │
│ Split 1  │ │ Split 2  │                               │
│ Compute  │ │ Compute  │  ◀── get_split_hr_fast()     │
│   HR_1   │ │   HR_2   │                               │
└────┬─────┘ └────┬─────┘                               │
     │            │                                     │
     └─────┬──────┘                                     │
           │                                            │
           ▼                                            │
┌───────────────────────┐                               │
│ Both HR > threshold?  │                               │
└───────────┬───────────┘                               │
            │                                           │
      ┌─────┴─────┐                                     │
      │           │                                     │
      ▼           ▼                                     │
┌──────────┐ ┌──────────┐                               │
│   YES    │ │    NO    │                               │
│Consistent│ │Inconsist.│                               │
│  Count++ │ │          │                               │
└────┬─────┘ └────┬─────┘                               │
     │            │                                     │
     └─────┬──────┘                                     │
           │                                            │
           └────────────────────────────────────────────┘

            │ (when loop done)
            ▼
┌───────────────────────────────┐
│ Pcons = Consistent / n.splits │
│                               │
│ If Pcons ≥ threshold: PASS    │
│ Otherwise: FAIL               │
└───────────────────────────────┘
```

---

### Two-Stage Algorithm (Optional)

```
┌─────────────────────────────────────────────────────────────────┐
│            TWO-STAGE CONSISTENCY ALGORITHM                       │
│            (use_twostage = TRUE)                                 │
└─────────────────────────────────────────────────────────────────┘

┌───────────────────────┐
│  Candidate Subgroup   │
└───────────┬───────────┘
            │
            ▼
┌───────────────────────────────────────────┐
│         STAGE 1: QUICK SCREEN             │
│         (n.splits.screen = 30)            │
└───────────────────────┬───────────────────┘
                        │
                        ▼
              ┌─────────────────┐
              │ Pcons_screen >  │
              │screen.threshold?│
              └────────┬────────┘
                       │
          ┌────────────┴────────────┐
          │                         │
          ▼                         ▼
   ┌────────────┐            ┌────────────┐
   │     NO     │            │    YES     │
   │  ──────▶   │            │            │
   │   REJECT   │            │ Continue   │
   └────────────┘            └─────┬──────┘
                                   │
                                   ▼
┌───────────────────────────────────────────┐
│       STAGE 2: BATCHED EVALUATION         │
│       with Early Stopping                 │
└───────────────────────┬───────────────────┘
                        │
                        ▼
              ┌─────────────────┐
              │ Run batch of    │◀─────────────────┐
              │ 20 splits       │                  │
              └────────┬────────┘                  │
                       │                           │
                       ▼                           │
              ┌─────────────────┐                  │
              │ Compute Wilson  │  ◀── wilson_ci() │
              │ Confidence Int. │                  │
              └────────┬────────┘                  │
                       │                           │
                       ▼                           │
              ┌─────────────────┐                  │
              │ CI_lower >      │                  │
              │ threshold?      │                  │
              └────────┬────────┘                  │
                       │                           │
          ┌────────────┴────────────┐              │
          │                         │              │
          ▼                         ▼              │
   ┌────────────┐            ┌────────────┐       │
   │    YES     │            │     NO     │       │
   │  ──────▶   │            │            │       │
   │   ACCEPT   │            └─────┬──────┘       │
   └────────────┘                  │              │
                                   ▼              │
                         ┌─────────────────┐      │
                         │ CI_upper <      │      │
                         │ threshold?      │      │
                         └────────┬────────┘      │
                                  │               │
                     ┌────────────┴────────────┐  │
                     │                         │  │
                     ▼                         ▼  │
              ┌────────────┐            ┌─────────┴───┐
              │    YES     │            │     NO      │
              │  ──────▶   │            │  CONTINUE   │
              │   REJECT   │            │  (more      │
              └────────────┘            │   batches)  │
                                        └─────────────┘
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

```
┌─────────────────────────────────────────────────────────────────┐
│                 BOOTSTRAP BIAS CORRECTION                        │
└─────────────────────────────────────────────────────────────────┘

┌───────────────────────────┐
│    Original ForestSearch  │
│    Result (H_obs)         │
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ Generate B Bootstrap      │
│ Samples (default B=1000)  │
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ FOR each bootstrap b:     │◀───────────────────────────┐
└────────────┬──────────────┘                            │
             │                                           │
             ▼                                           │
┌───────────────────────────┐                            │
│ Draw bootstrap sample     │                            │
│ with replacement          │                            │
└────────────┬──────────────┘                            │
             │                                           │
             ▼                                           │
┌───────────────────────────┐                            │
│ Run forestsearch() on     │  (may find different      │
│ bootstrap sample          │   subgroup than original) │
└────────────┬──────────────┘                            │
             │                                           │
             ▼                                           │
┌───────────────────────────┐                            │
│ Compute:                  │                            │
│ • H*_b  = orig SG on boot │                            │
│ • H**_b = new SG on boot  │                            │
│ • H*obs = new SG on orig  │                            │
└────────────┬──────────────┘                            │
             │                                           │
             └───────────────────────────────────────────┘

             │ (when loop done)
             ▼
┌───────────────────────────┐
│ Aggregate Results         │
│ Compute Bias Correction   │
└────────────┬──────────────┘
             │
             ▼
┌─────────────────────────────────────────────────────────────────┐
│                    BIAS CORRECTION FORMULAS                      │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Method 1 (Simple Optimism):                                    │
│      H_adj1 = H_obs - (H**_mean - H*obs_mean)                   │
│                                                                 │
│  Method 2 (Double Bootstrap):                                   │
│      H_adj2 = 2 × H_obs - (H*_mean + H**_mean - H*obs_mean)     │
│                                                                 │
│  Where:                                                         │
│      H_obs     = Original subgroup HR on original data          │
│      H*        = Original subgroup HR on bootstrap data         │
│      H*obs     = New subgroup (from bootstrap) HR on orig data  │
│      H**       = New subgroup (from bootstrap) HR on boot data  │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

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

```
┌─────────────────────────────────────────────────────────────────┐
│                    K-FOLD CROSS-VALIDATION                       │
└─────────────────────────────────────────────────────────────────┘

┌───────────────────────────┐
│      Full Dataset         │
│      (n observations)     │
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ Split into K Folds        │
│ (default K = 10)          │
│                           │
│ ┌───┬───┬───┬───┬───┐     │
│ │ 1 │ 2 │ 3 │...│ K │     │
│ └───┴───┴───┴───┴───┘     │
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ FOR fold k = 1 to K:      │◀───────────────────────────┐
└────────────┬──────────────┘                            │
             │                                           │
             ▼                                           │
┌───────────────────────────────────────┐                │
│  Training Set: Folds 1,...,k-1,k+1,...K│               │
│  Test Set:     Fold k                  │               │
│                                        │               │
│  ┌───┬───┬───┬───┬───┐                │               │
│  │ T │ T │ T │TST│ T │  (K=5 example) │               │
│  └───┴───┴───┴───┴───┘                │               │
└────────────┬──────────────────────────┘                │
             │                                           │
             ▼                                           │
┌───────────────────────────┐                            │
│ Run forestsearch() on     │                            │
│ Training Set              │                            │
└────────────┬──────────────┘                            │
             │                                           │
             ▼                                           │
┌───────────────────────────┐                            │
│ Apply identified subgroup │                            │
│ to Test Set               │                            │
│ • Compute HR in test      │                            │
│ • Check if SG matches     │                            │
└────────────┬──────────────┘                            │
             │                                           │
             └───────────────────────────────────────────┘

             │ (when loop done)
             ▼
┌───────────────────────────┐
│ forestsearch_KfoldOut()   │
│ • Agreement across folds  │
│ • Stability metrics       │
│ • Covariate matching      │
└───────────────────────────┘
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

```
┌─────────────────────────────────────────────────────────────────┐
│                    OUTPUT FUNCTIONS                              │
└─────────────────────────────────────────────────────────────────┘

                    ┌───────────────────────┐
                    │   forestsearch        │
                    │   result object       │
                    └───────────┬───────────┘
                                │
        ┌───────────────────────┼───────────────────────┐
        │                       │                       │
        ▼                       ▼                       ▼
┌───────────────┐       ┌───────────────┐       ┌───────────────┐
│    print()    │       │  summary()    │       │  sg_tables()  │
│               │       │               │       │               │
│ Basic console │       │ Detailed      │       │ gt-formatted  │
│ summary       │       │ statistics    │       │ tables        │
└───────────────┘       └───────────────┘       └───────────────┘
        │                       │                       │
        ▼                       ▼                       ▼
┌───────────────┐       ┌───────────────┐       ┌───────────────┐
│ • SG defn     │       │ • Parameters  │       │ • Estimates   │
│ • HR estimate │       │ • Var select  │       │ • Top 10 SGs  │
│ • Pcons       │       │ • Consistency │       │ • Search info │
│ • N, Events   │       │ • By sg_focus │       │               │
└───────────────┘       └───────────────┘       └───────────────┘


┌───────────────────────────────────────────────────────────────┐
│              plot_subgroup_results_forestplot()                │
│                                                               │
│  Creates publication-ready forest plot with:                  │
│  • ITT population estimate                                    │
│  • Reference subgroups                                        │
│  • Post-hoc identified subgroups                              │
│  • Bias-corrected estimates                                   │
│  • Cross-validation metrics                                   │
└───────────────────────────────────────────────────────────────┘
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

```
┌─────────────────────────────────────────────────────────────────┐
│                    SIMULATION WORKFLOW                           │
└─────────────────────────────────────────────────────────────────┘

┌───────────────────────────┐
│ Define DGM Parameters     │
│ • Treatment effect        │
│ • Subgroup definition     │
│ • Censoring distribution  │
│ • Covariate distributions │
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ generate_aft_dgm_flex()   │  ◀── Create AFT data generating
│                           │      mechanism (super population)
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ simulate_from_dgm()       │  ◀── Draw samples from DGM
│ • Sample n observations   │
│ • Apply censoring         │
│ • Randomize treatment     │
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ Simulated Trial Data      │
│ (with known truth)        │
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ run_forestsearch_analysis()│  ◀── Apply ForestSearch
└────────────┬──────────────┘
             │
             ▼
┌───────────────────────────┐
│ Evaluate Performance      │
│ • Sensitivity/Specificity │
│ • PPV/NPV                 │
│ • HR bias                 │
│ • Subgroup agreement      │
└───────────────────────────┘
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

```
┌─────────────────────────────────────────────────────────────────┐
│                   PARALLEL BACKENDS                              │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────┐     ┌─────────────────────────────────────────┐
│   sequential    │────▶│ Single thread (debugging/testing)       │
└─────────────────┘     └─────────────────────────────────────────┘

┌─────────────────┐     ┌─────────────────────────────────────────┐
│  multisession   │────▶│ Background R sessions (all platforms)   │
└─────────────────┘     └─────────────────────────────────────────┘

┌─────────────────┐     ┌─────────────────────────────────────────┐
│   multicore     │────▶│ Fork processes (Unix/Mac only)          │
└─────────────────┘     └─────────────────────────────────────────┘

┌─────────────────┐     ┌─────────────────────────────────────────┐
│     callr       │────▶│ Separate R processes (recommended)      │
└─────────────────┘     └─────────────────────────────────────────┘
```

---

## Package Dependencies

```
┌─────────────────────────────────────────────────────────────────┐
│                   PACKAGE DEPENDENCIES                           │
└─────────────────────────────────────────────────────────────────┘

                         ┌─────────────┐
                         │ forestsearch│
                         └──────┬──────┘
                                │
        ┌───────────────────────┼───────────────────────┐
        │                       │                       │
        ▼                       ▼                       ▼
┌───────────────┐       ┌───────────────┐       ┌───────────────┐
│    CORE       │       │     DATA      │       │   PARALLEL    │
│  STATISTICAL  │       │   HANDLING    │       │  COMPUTING    │
├───────────────┤       ├───────────────┤       ├───────────────┤
│ • survival    │       │ • data.table  │       │ • future      │
│ • grf         │       │ • stringr     │       │ • doFuture    │
│ • glmnet      │       │               │       │ • foreach     │
│ • policytree  │       │               │       │ • callr       │
└───────────────┘       └───────────────┘       └───────────────┘
        │
        ▼
┌───────────────┐
│ VISUALIZATION │
├───────────────┤
│ • ggplot2     │
│ • gt          │
│ • forestploter│
└───────────────┘
```

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
