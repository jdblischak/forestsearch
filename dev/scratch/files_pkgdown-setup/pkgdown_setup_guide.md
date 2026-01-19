# pkgdown Setup Guide for ForestSearch

This guide provides step-by-step instructions for setting up a pkgdown website for the ForestSearch package.

---

## 1. Installation and Initial Setup

### Install pkgdown

```r
install.packages("pkgdown")
```

### Initialize pkgdown in your package

```r
# From your package root directory
usethis::use_pkgdown()
```

This creates:
- `_pkgdown.yml` - Configuration file
- `pkgdown/` directory - For custom assets (optional)
- Adds `^_pkgdown\.yml$` and `^docs$` to `.Rbuildignore`

### Build the site locally

```r
# Build entire site
pkgdown::build_site()

# Preview in browser
pkgdown::preview_site()

# Build only the reference (faster for testing)
pkgdown::build_reference()
```

---

## 2. Recommended _pkgdown.yml Configuration

Create or replace `_pkgdown.yml` in your package root:

```yaml
url: https://yourusername.github.io/forestsearch/

template:

  bootstrap: 5
  bootswatch: flatly
  bslib:
    primary: "#0d6efd"
    border-radius: 0.5rem
    btn-border-radius: 0.25rem

home:
  title: ForestSearch
  description: >
    Exploratory subgroup identification for clinical trials with survival endpoints,
    featuring bootstrap bias correction and cross-validation methods.

authors:
  Larry Leon:
    href: https://github.com/yourusername

navbar:
  structure:
    left:  [intro, reference, articles, tutorials, news]
    right: [search, github]
  components:
    github:
      icon: fa-github
      href: https://github.com/yourusername/forestsearch
      aria-label: GitHub

reference:
# ==============================================================================
# CORE ANALYSIS
# ==============================================================================
- title: "Core ForestSearch Analysis"
  desc: >
    Main functions for running subgroup identification analysis
  contents:
  - forestsearch
  - print.forestsearch
  - summary.forestsearch
  - plot.forestsearch

- title: "Data Preparation"
  desc: >
    Functions to prepare data for ForestSearch analysis, including
    LASSO selection, GRF cuts, and factor creation
  contents:
  - get_FSdata
  - add_id_column
  - starts_with("get_conf")
  - starts_with("get_cut")
  - starts_with("dummy")
  - add_unprocessed_vars

# ==============================================================================
# BOOTSTRAP METHODS
# ==============================================================================
- title: "Bootstrap Analysis"
  desc: >
    Bootstrap bias correction for subgroup-specific treatment effects
  contents:
  - forestsearch_bootstrap_dofuture
  - bootstrap_results
  - bootstrap_ystar
  - summarize_bootstrap_results
  - summarize_bootstrap_subgroups

- title: "Bootstrap Helpers"
  desc: >
    Utility functions for bootstrap calculations
  contents:
  - count_boot_id
  - calc_cov
  - ci_est
  - format_ci
  - bc_ci_result
  - target_est_se

- title: "Bootstrap Summaries & Tables"
  desc: >
    Functions for summarizing and formatting bootstrap results
  contents:
  - format_bootstrap_table
  - format_bootstrap_timing_table
  - format_bootstrap_diagnostics_table
  - summarize_bootstrap_timing
  - summarize_event_counts
  - generate_bootstrap_caption

# ==============================================================================
# CROSS-VALIDATION
# ==============================================================================
- title: "Cross-Validation"
  desc: >
    K-fold cross-validation for subgroup stability assessment
  contents:
  - forestsearch_Kfold
  - starts_with("CV_")
  - contains("kfold")
  - sens_text

# ==============================================================================
# SUBGROUP CONSISTENCY
# ==============================================================================
- title: "Subgroup Consistency Analysis"
  desc: >
    Functions for evaluating subgroup consistency across data splits
  contents:
  - evaluate_subgroup_consistency
  - run_single_consistency_split
  - get_split_hr_fast
  - wilson_ci
  - early_stop_decision

- title: "Duplicate Detection"
  desc: >
    Functions for identifying and handling duplicate/near-duplicate subgroups
  contents:
  - identify_near_duplicates
  - remove_near_duplicate_subgroups
  - show_duplicate_subgroups

# ==============================================================================
# COX MODEL UTILITIES
# ==============================================================================
- title: "Cox Model Estimation"
  desc: >
    Functions for fitting Cox models and extracting estimates
  contents:
  - get_Cox_sg
  - build_cox_formula
  - Cox_sg_recommendation
  - cox_summary
  - cox_summary_vectorized
  - cox_summary_batch
  - cox_summary_legacy

- title: "Cox Model Extensions"
  desc: >
    Advanced Cox model functionality including splines and AHR
  contents:
  - cox_ahr_cde_analysis
  - print.cox_ahr_cde
  - summary.cox_ahr_cde
  - cox_cs_fit
  - plot_subgroup_effects

# ==============================================================================
# GRF INTEGRATION
# ==============================================================================
- title: "Generalized Random Forests"
  desc: >
    Integration with the grf package for causal forest-based subgroup identification
  contents:
  - grf_subgroup_search
  - create_grf_config
  - fit_causal_forest
  - fit_policy_trees
  - select_best_subgroup
  - compute_node_metrics
  - extract_tree_cuts
  - extract_all_tree_cuts
  - find_leaf_split
  - assign_subgroup_membership
  - create_success_result
  - create_null_result
  - print_grf_details
  - validate_grf_data

# ==============================================================================
# SIMULATION & DATA GENERATION
# ==============================================================================
- title: "Data Generation"
  desc: >
    Functions for generating simulated clinical trial data with
    known subgroup effects using AFT models
  contents:
  - generate_aft_dgm_flex
  - simulate_from_dgm
  - plot_spline_treatment_effect
  - compare_multiple_survreg

- title: "Synthetic Data"
  desc: >
    Functions for creating perturbed/synthetic versions of datasets
  contents:
  - starts_with("perturb_")
  - starts_with("synthetic_")

# ==============================================================================
# VISUALIZATION
# ==============================================================================
- title: "Forest Plots"
  desc: >
    Publication-ready forest plot visualization
  contents:
  - plot_subgroup_results_forestplot
  - print.fs_forestplot
  - plot.fs_forestplot
  - create_subgroup_summary_df

- title: "Subgroup Visualization"
  desc: >
    Additional plotting functions for subgroup analysis
  contents:
  - plot_subgroup

# ==============================================================================
# SUMMARY TABLES
# ==============================================================================
- title: "Summary Tables"
  desc: >
    Functions for creating publication-ready summary tables
  contents:
  - sg_tables
  - create_summary_table
  - create_summary_table_compact
  - create_summary_table_publication
  - create_summary_table_presentation
  - create_summary_table_minimal
  - create_factor_summary_tables
  - format_subgroup_definition_table
  - format_subgroup_agreement_table

# ==============================================================================
# SUBGROUP SEARCH
# ==============================================================================
- title: "Subgroup Search Algorithms"
  desc: >
    Core search algorithms for identifying subgroups
  contents:
  - find_k_interaction
  - find_markers
  - starts_with("subgroup_")

# ==============================================================================
# INTERNAL/UTILITIES
# ==============================================================================
- title: "Utility Functions"
  desc: >
    Internal helper functions (exported for advanced use)
  contents:
  - setup_parallel_backend
  - safe_expr
  - sort_dt_by_cols
  - FS_labels
  - SG_tab_estimates
  - analyze_subgroup

articles:
- title: Getting Started
  navbar: ~
  contents:
  - forestsearch-introduction
  - data-preparation

- title: Analysis Workflows
  contents:
  - basic-analysis
  - bootstrap-bias-correction
  - cross-validation
  - multi-regional-trials

- title: Advanced Topics
  contents:
  - simulation-studies
  - custom-subgroup-definitions
  - grf-integration

- title: Reference
  contents:
  - methodology
  - leon-2024-replication

news:
  releases:
  - text: "Version 1.0.0"
    href: https://github.com/yourusername/forestsearch/releases/tag/v1.0.0

footer:
  structure:
    left: developed_by
    right: built_with
```

---

## 3. Recommended Vignettes/Articles

Create vignettes in `vignettes/` directory. pkgdown will automatically convert them to articles.

### Essential Vignettes

| File | Title | Content |
|------|-------|---------|
| `forestsearch-introduction.Rmd` | Introduction to ForestSearch | Package overview, methodology, quick start |
| `data-preparation.Rmd` | Data Preparation | How to use `get_FSdata()`, LASSO, GRF cuts |
| `basic-analysis.Rmd` | Basic Subgroup Analysis | Step-by-step analysis with GBSG data |
| `bootstrap-bias-correction.Rmd` | Bootstrap Bias Correction | Using `forestsearch_bootstrap_dofuture()` |
| `cross-validation.Rmd` | Cross-Validation | Using `forestsearch_Kfold()` for stability |

### Advanced Vignettes

| File | Title | Content |
|------|-------|---------|
| `simulation-studies.Rmd` | Simulation Studies | Using `generate_aft_dgm_flex()` for simulations |
| `multi-regional-trials.Rmd` | Multi-Regional Trials | Regional consistency evaluation |
| `methodology.Rmd` | Statistical Methodology | Mathematical details, references to León et al. |
| `leon-2024-replication.Rmd` | Replicating León et al. (2024) | Code to reproduce paper results |

### Vignette Template

```r
---
title: "Introduction to ForestSearch"
output: rmarkdown::html_vignette
vignette: >
%\VignetteIndexEntry{Introduction to ForestSearch}
%\VignetteEngine{knitr::rmarkdown}
%\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>",
fig.width = 7,
fig.height = 5
)
library(forestsearch)
```

## Overview

ForestSearch is an R package for exploratory subgroup identification...

## Quick Start

```{r example, eval = FALSE}
# Load example data
data(gbsg, package = "survival")

# Prepare data
fs_data <- get_FSdata(
df.analysis = gbsg,
confounders.name = c("age", "size", "grade"),
outcome.name = "rfstime",
event.name = "status"
)

# Run analysis
results <- forestsearch(...)
```
```

---

## 4. Adding a Logo

Create a hex sticker logo and add it:

```r
# Create logo (save as man/figures/logo.png)
# Then reference in _pkgdown.yml:
```

Add to `_pkgdown.yml`:

```yaml
template:
logo:
  src: man/figures/logo.png
  width: 120
```

Or use `usethis::use_logo("path/to/logo.png")`.

---
## 5. Custom Styling (Optional)

Create `pkgdown/extra.css` for custom styles:

```css
/* Custom ForestSearch styling */

/* Syntax highlighting for R code */
pre {
background-color: #f8f9fa;
border: 1px solid #e9ecef;
border-radius: 0.25rem;
}

/* Function reference cards */
.card {
border-left: 4px solid #0d6efd;
}

/* Section headers in reference */
.section h2 {
border-bottom: 2px solid #0d6efd;
padding-bottom: 0.5rem;
}
```

Reference it in `_pkgdown.yml`:

```yaml
template:
includes:
  in_header: |
    <!-- Custom meta tags -->
  after_body: |
    <!-- Custom scripts -->
```

---

## 6. GitHub Pages Deployment

### Option A: GitHub Actions (Recommended)

Create `.github/workflows/pkgdown.yaml`:
```yaml
on:
push:
  branches: [main, master]
pull_request:
  branches: [main, master]
release:
  types: [published]
workflow_dispatch:

name: pkgdown

jobs:
pkgdown:
  runs-on: ubuntu-latest
  concurrency:
    group: pkgdown-${{ github.event_name != 'pull_request' || github.run_id }}
  env:
    GITHUB_PAT: ${{ secrets.GITHUB_TOKEN }}
  permissions:
    contents: write
  steps:
    - uses: actions/checkout@v4

    - uses: r-lib/actions/setup-pandoc@v2

    - uses: r-lib/actions/setup-r@v2
      with:
        use-public-rspm: true

    - uses: r-lib/actions/setup-r-dependencies@v2
      with:
        extra-packages: any::pkgdown, local::.
        needs: website

    - name: Build site
      run: pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)
      shell: Rscript {0}

    - name: Deploy to GitHub pages 🚀
      if: github.event_name != 'pull_request'
      uses: JamesIves/github-pages-deploy-action@v4.5.0
      with:
        clean: false
        branch: gh-pages
        folder: docs
```

Or use the helper:
```r
usethis::use_pkgdown_github_pages()
```

### Option B: Manual Deployment

```r
# Build site
pkgdown::build_site()

# The site is in docs/ - push to gh-pages branch
# Or configure GitHub Pages to serve from docs/ on main branch
```

---

## 7. Reference Index Organization Tips

With 119 exported functions, good organization is crucial. The configuration above groups functions by:

1. **Core Analysis** - Main user-facing functions
2. **Data Preparation** - Pre-processing utilities
3. **Bootstrap Methods** - Bias correction
4. **Cross-Validation** - Stability assessment
5. **Subgroup Consistency** - Split-sample validation
6. **Cox Model Utilities** - Statistical estimation
7. **GRF Integration** - Causal forest methods
8. **Simulation** - Data generation for studies
9. **Visualization** - Plotting functions
10. **Summary Tables** - Formatted output
11. **Utilities** - Helper functions

### Selector Patterns

pkgdown supports these selectors in `contents:`:

```yaml
contents:
- function_name          # Exact match
- starts_with("prefix_") # Functions starting with prefix
- ends_with("_suffix")   # Functions ending with suffix
- contains("pattern")    # Functions containing pattern
- matches("regex")       # Regex pattern matching
- has_concept("concept") # Functions with @concept tag
- has_keyword("internal")# Functions with @keywords
```

---

## 8. Enhancing Documentation for pkgdown

### Add @family tags for cross-linking

```r
#' @family bootstrap functions
#' @family subgroup analysis
```

### Add @seealso for related functions

```r
#' @seealso
#' \code{\link{forestsearch}} for main analysis
#' \code{\link{forestsearch_Kfold}} for cross-validation
```

### Add @concept for grouping

```r
#' @concept subgroup-identification
#' @concept bootstrap
```

---

## 9. Build Checklist

Before publishing:

- [ ] Run `pkgdown::build_site()` locally
- [ ] Check all function documentation renders correctly
- [ ] Verify vignettes build without errors
- [ ] Test navigation and search
- [ ] Check responsive design on mobile
- [ ] Verify external links work
- [ ] Review function groupings in reference index
- [ ] Add Google Analytics (optional)

```yaml
# In _pkgdown.yml for analytics:
template:
includes:
  in_header: |
    <!-- Google tag (gtag.js) -->
    <script async src="https://www.googletagmanager.com/gtag/js?id=G-XXXXXXX"></script>
```

---

## 10. Quick Start Commands

```r
# One-time setup
usethis::use_pkgdown()
usethis::use_pkgdown_github_pages()  # If using GitHub

# Development workflow
pkgdown::build_site()                # Full build
pkgdown::build_reference()           # Reference only
pkgdown::build_articles()            # Articles only
pkgdown::build_home()                # Home page only
pkgdown::preview_site()              # Preview in browser

# Check for issues
pkgdown::check_pkgdown()
```

---

## Summary

The pkgdown site will provide:

1. **Professional landing page** with package description
2. **Organized function reference** grouped by topic (119 functions in ~12 categories)
3. **Searchable documentation** with full-text search
4. **Vignettes as articles** for tutorials and methodology
5. **News/changelog** tracking
6. **GitHub integration** with source links

This makes ForestSearch much more accessible to clinical researchers who may not be familiar with navigating R package documentation directly.
