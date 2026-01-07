
mypackage/
├── R/
│   └── *.R                    # Package functions
├── man/
│   └── *.Rd                   # Generated documentation (roxygen2)
├── src/
│   └── *.cpp                  # C/C++ code if needed
├── data/
│   └── *.rda                  # Package datasets
├── data-raw/
│   └── *.R                    # Scripts to prepare datasets
├── inst/
│   ├── extdata/               # External data files
│   └── scripts/               # Utility scripts not part of package
├── vignettes/
│   └── *.Rmd                  # Package vignettes (installed with package)
├── tests/
│   └── testthat/
│       ├── testthat.R
│       └── test-*.R           # Unit tests
├── dev/
│   ├── scratch/               # Experimental code, notes
│   ├── vignettes-working/     # Vignette drafts before moving to vignettes/
│   └── dev-workflow.R         # Package development helper script
├── .github/
│   └── workflows/
│       └── R-CMD-check.yaml   # CI/CD workflows
├── .Rbuildignore              # Exclude dev/ and other non-package files
├── .gitignore
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── NEWS.md
└── README.md

Key points:
The dev/ folder is your working area—add it to .Rbuildignore so it doesn't get bundled with the package. Use dev/vignettes-working/ to draft vignettes before polishing and moving them to vignettes/. The data-raw/ folder (also in .Rbuildignore) stores the scripts that generate your data/*.rda files, which is useful for reproducibility.
For your .Rbuildignore, include:

^dev$
^data-raw$
^\.github$
^.*\.Rproj$
^README\.Rmd$

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`forestsearch` is an R package for identifying subgroups with heterogeneous treatment effects in survival data. It uses causal forest methods (GRF), LASSO-based dimension reduction, and consistency-based validation to discover subgroups where treatment may be harmful or less effective.

## Core Workflow

The typical analysis flow follows this pattern:

1. **forestsearch()** (R/forest_search.R:137-540) - Main entry point that orchestrates the entire analysis
2. **GRF screening** (R/grf_main.R:38-193) - Optional causal forest for identifying important cut points
3. **get_FSdata()** (R/get_FSdata_refactored.r) - Prepares data, applies LASSO screening, creates candidate factors
4. **subgroup.search()** (R/subgroup_search.R:174-244) - Exhaustive search through factor combinations
5. **subgroup.consistency()** (R/subgroup_consistency_main.r:293-845) - Validates subgroups via random splits
6. **forestsearch_bootstrap_dofuture()** (R/bootstrap_dofuture_main.r:115-417) - Bootstrap analysis with bias correction

## Key Architecture Patterns

### 1. Subgroup Search Strategy (maxk parameter)
The search evaluates combinations of binary factors up to size `maxk` (default: 2):
- maxk=1: Single factors only (e.g., "age>65")
- maxk=2: Up to 2-factor combinations (e.g., "age>65 & male")
- maxk=3: Up to 3-factor combinations (computational cost increases significantly)

Combinations are generated in R/subgroup_search.R:274-302 using combn() and evaluated sequentially.

### 2. Factor Creation Pipeline
Confounders are converted to binary factors via:
- **Continuous variables**: Cut at median or GRF-identified thresholds
- **Categorical variables**: Dummy encoding (reference level dropped)
- **GRF cuts**: Data-driven cut points from causal_survival_forest (optional)
- **Forced cuts**: User-specified cut expressions (conf_force parameter)

The cut selection logic is in get_FSdata() which uses:
- `cut_type="default"`: Reference level for factors, no cuts for continuous unless specified
- `cut_type="median"`: Median cuts for continuous variables
- GRF cuts override median cuts when `replace_med_grf=TRUE`

### 3. Consistency Validation
After initial HR-based screening, candidates undergo consistency testing (R/subgroup_consistency_main.r:28-260):
- Dataset split randomly into halves (n.splits times, default 1000)
- Subgroup must show HR > hr.consistency in BOTH halves
- Proportion meeting threshold must exceed pconsistency.threshold (default 0.90)
- This step uses either sequential or parallel processing via evaluate_subgroup_consistency()

### 4. Parallel Processing
Two separate parallel systems:
- **Consistency evaluation**: Uses doFuture/foreach (configured via parallel_args)
- **Bootstrap analysis**: Uses future.apply (inherits or overrides parallel_args)

Plans: "multisession", "multicore", "callr", "sequential"
Setup in R/bootstrap_dofuture_setup_helpers.R

### 5. Bias Correction Methods
Bootstrap implements two bias correction approaches (R/bootstrap_dofuture_main.r:6-25):
- **Method 1**: Simple optimism correction: H_adj1 = H_obs - (H*_* - H*_obs)
- **Method 2**: Double bootstrap: H_adj2 = 2*H_obs - (H_* + H*_* - H*_obs)

Where:
- H_obs: Original subgroup HR on full data
- H_*: Original subgroup HR on bootstrap sample
- H*_*: New subgroup HR on bootstrap sample
- H*_obs: New subgroup HR evaluated on full data

## Important Function Signatures

### Main Analysis
```r
forestsearch(
  df.analysis,                    # Main dataset
  outcome.name = "tte",           # Time-to-event variable
  event.name = "event",           # Event indicator (0/1)
  treat.name = "treat",           # Treatment (0/1)
  confounders.name,               # Character vector of covariate names
  parallel_args = list(plan = "callr", workers = 6),
  use_grf = TRUE,                 # Use GRF for cut selection
  use_lasso = TRUE,               # Use LASSO for screening
  hr.threshold = 1.25,            # Min HR for subgroup candidate
  hr.consistency = 1.0,           # Min HR for consistency splits
  pconsistency.threshold = 0.90,  # Min proportion of consistent splits
  fs.splits = 1000,               # Number of consistency splits
  maxk = 2,                       # Max subgroup factors
  sg_focus = "hr"                 # "hr", "maxSG", or "minSG"
)
```

### Bootstrap Analysis
```r
forestsearch_bootstrap_dofuture(
  fs.est,                         # Output from forestsearch()
  nb_boots = 1000,                # Number of bootstrap samples
  parallel_args = list(),         # Inherits from forestsearch if empty
  details = FALSE,
  show_three = FALSE              # Verbose output for first 3 iterations
)
```

## Development Commands

### Building and Checking
```r
# Generate documentation from roxygen comments
devtools::document()

# Build package
devtools::build()

# Check package
devtools::check()

# Install locally
devtools::install()

# Load for interactive development
devtools::load_all()
```

### Rendering README
The README.Rmd needs to be knitted to update README.md:
```r
devtools::build_readme()
```

## Code Organization

### Helper Function Structure
Most main functions follow this pattern:
1. Input validation section
2. Data preparation
3. Core algorithm (often parallelized)
4. Post-processing/formatting
5. Return structured list

Helper functions are typically defined in separate files:
- `*_helpers.R`: Utility functions for specific modules
- `*_main.R`: Primary exported functions

### Data.table Usage
The package extensively uses data.table for performance. Key patterns:
- := for in-place modification
- .SD for subsetting operations
- Global variables declared in R/forestsearch-package.R:28-100

### Error Handling
- try() with silent=TRUE for model fitting
- Graceful degradation: returns NULL components rather than failing entirely
- Warning messages for partial failures

## Important Implementation Details

### Subgroup Membership Logic
Subgroups are defined by requiring ALL selected factors = 1:
```r
# In R/subgroup_search.R:64-68
get_subgroup_membership <- function(zz, covs.in) {
  selected_cols <- which(covs.in == 1)
  x <- zz[, selected_cols, drop = FALSE]
  rowSums(x) == length(selected_cols)  # AND logic
}
```

### Redundancy Filtering
Combinations are filtered if adding a factor doesn't reduce sample size by at least `rmin`:
```r
# R/subgroup_search.R:128-142
extract_idx_flagredundancy(x, rmin)
```

### GRF Time Horizon
Restricted Mean Survival Time (RMST) horizon calculated as:
```r
tau.rmst <- frac.tau * min(max(Y[W==1 & D==1]), max(Y[W==0 & D==1]))
```
Default frac.tau=0.6 for screening, 0.9 for final models.

### Treatment Recommendation Flag
After subgroup identification:
- treat.recommend = 0: Harm/questionable subgroup (avoid treatment)
- treat.recommend = 1: Benefit/recommend subgroup (give treatment)

## Testing and Simulation

The dev/ directory contains:
- private/: Internal development scripts
- public/: Example analyses and vignettes
- legacy_code/: Older implementations preserved for reference

Simulation functions in R/simulate_from_dgm.R and R/sim_weibspline_functions_v1.R for validation studies.

## Dependencies

Core dependencies (see DESCRIPTION:11-27):
- **grf**: Causal forest implementation
- **policytree**: Policy tree for interpretable splits
- **survival**: Cox models and survival analysis
- **data.table**: High-performance data manipulation
- **glmnet**: LASSO for dimension reduction
- **doFuture/foreach/future**: Parallel processing
- **weightedsurv**: Weighted survival analysis (custom remote: larry-leon/weightedsurv)

## Common Issues

1. **GRF returns NULL**: Increase frac.tau or decrease dmin.grf if no subgroups meet RMST threshold
2. **No subgroups found**: Lower hr.threshold or increase maxk
3. **Consistency fails**: Decrease pconsistency.threshold or fs.splits if subgroups are small
4. **Parallel errors**: Check that workers < available cores; "callr" plan is most stable
5. **Memory issues**: Reduce nb_boots or workers in bootstrap analysis

## Variable Naming Conventions

- **confounders.name**: Original covariate names in dataset
- **FSconfounders.name**: Factor names after dummy encoding
- **sg.harm**: Selected subgroup factor names (harmful/questionable subgroup)
- **H / Hc**: Harm subgroup / Complement (recommended) subgroup
- **_obs**: Estimate from original data
- **_star**: Estimate from bootstrap sample
- **treat.recommend**: Binary flag (0=avoid treatment, 1=recommend)
- **flag.harm**: True harm status in simulations

## Output Structure

Main forestsearch() returns list with:
- **grp.consistency**: Full consistency analysis results
- **find.grps**: Initial search results (before consistency)
- **df.est**: Dataset with treat.recommend flag
- **df.predict**: Prediction dataset (if provided)
- **sg.harm**: Selected subgroup factor names
- **grf_cuts**: GRF-identified cut expressions
- **args_call_all**: All arguments for reproducibility

Bootstrap returns:
- **results**: Data.table with all bootstrap iterations
- **SG_CIs**: Confidence intervals (raw and bias-corrected)
- **FSsg_tab**: Formatted summary table
- **H_estimates / Hc_estimates**: Detailed estimate objects
