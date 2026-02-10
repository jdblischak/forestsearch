# Export Audit: Reducing the Public API Surface

## Overview

The NAMESPACE currently exports **119 functions** plus **15 S3 methods**
(S3 methods are auto-registered and not counted toward the export total).
This audit categorizes every export into one of three tiers:

| Tier | Description | Recommendation |
|------|-------------|----------------|
| **KEEP** | User-facing entry points | Remain `@export` |
| **INTERNAL** | Algorithm internals, helpers, parallel plumbing | Change to `@keywords internal` (still documented but not in main index) |
| **NORD** | Pure implementation detail, no user value | Change to `@noRd` (no documentation page generated) |

**Important constraint:** Functions called by parallel workers via
`future`/`doFuture`/`callr` backends *must* remain exported, because
workers only see the installed package namespace. These are marked with
⚡ below. This is a legitimate reason to keep otherwise-internal functions
exported, but they should still get `@keywords internal` so they don't
clutter the user-facing reference index.

---

## Tier 1: KEEP Exported (~50 functions)

These are the functions a user would reasonably call directly.

### Core Algorithm (5)
| Function | Rationale |
|----------|-----------|
| `forestsearch` | Primary entry point |
| `subgroup.search` | Core search engine |
| `get_FSdata` | Data preparation entry point |
| `subgroup.consistency` | Consistency evaluation entry point |
| `forestsearch_bootstrap_dofuture` | Bootstrap bias correction entry point |

### Cross-Validation (4)
| Function | Rationale |
|----------|-----------|
| `forestsearch_Kfold` | K-fold CV entry point |
| `forestsearch_tenfold` | Repeated K-fold entry point |
| `forestsearch_KfoldOut` | CV output summarization |
| `CV_sgs` | CV subgroup comparison |

### GRF Integration (2)
| Function | Rationale |
|----------|-----------|
| `grf.subg.harm.survival` | GRF subgroup identification entry point |
| `grf.subg.eval` | GRF evaluation entry point |

### Cox Model Utilities (4)
| Function | Rationale |
|----------|-----------|
| `cox_summary` | Primary Cox model wrapper |
| `cox_ahr_cde_analysis` | AHR/CDE analysis entry point |
| `cox_cs_fit` | Cox spline fitting entry point |
| `rmst_calculation` | RMST calculation |

### Visualization (10)
| Function | Rationale |
|----------|-----------|
| `plot_subgroup_results_forestplot` | Publication forest plot |
| `render_forestplot` | Render forest plot |
| `save_forestplot` | Save forest plot to file |
| `create_forest_theme` | Theme customization |
| `plot_km_band_forestsearch` | KM band plot |
| `quick_km_band_plot` | Convenience KM wrapper |
| `plot_sg_results` | Subgroup results visualization |
| `plot_sg_weighted_km` | Weighted KM curves |
| `plot_spline_treatment_effect` | Spline treatment effect plot |
| `plot_subgroup_effects` | Subgroup effect visualization |

### Detection Curves (3)
| Function | Rationale |
|----------|-----------|
| `plot_detection_curve` | Detection curve visualization |
| `compare_detection_curves` | Multi-curve comparison |
| `generate_detection_curve` | Generate detection data |

### Tables & Summaries (7)
| Function | Rationale |
|----------|-----------|
| `sg_tables` | Primary subgroup tables (gt) |
| `SG_tab_estimates` | Subgroup estimate tables |
| `SGplot_estimates` | Plot-ready estimates |
| `FS_labels` | Subgroup label generation |
| `create_summary_table` | Primary summary table |
| `summarize_bootstrap_results` | Bootstrap summary entry point |
| `format_bootstrap_table` | Bootstrap table formatting |

### Simulation Framework (6)
| Function | Rationale |
|----------|-----------|
| `generate_aft_dgm_flex` | DGM creation entry point |
| `create_gbsg_dgm` | GBSG-based DGM |
| `simulate_from_dgm` | Simulate from DGM |
| `simulate_from_gbsg_dgm` | Simulate from GBSG DGM |
| `run_simulation_analysis` | Run OC simulation |
| `summarize_simulation_results` | Summarize OC results |

### MRCT (3)
| Function | Rationale |
|----------|-----------|
| `mrct_region_sims` | MRCT simulation entry point |
| `summaryout_mrct` | MRCT summary |
| `validate_mrct_data` | MRCT data validation |

### Data Preparation (3)
| Function | Rationale |
|----------|-----------|
| `get_dfpred` | Prediction data preparation |
| `dummy_encode` | Primary dummy encoding function |
| `lasso_selection` | LASSO variable selection |

### Formatting Helpers (3)
| Function | Rationale |
|----------|-----------|
| `format_CI` | CI formatting for reports |
| `hrCI_format` | HR + CI formatting |
| `n_pcnt` | N (%) formatting |

**Tier 1 total: ~50 functions**

---

## Tier 2: INTERNAL — Keep Exported but Add `@keywords internal` (~40 functions)

These serve other package functions or parallel workers. They must remain
in the namespace but should not appear in the main reference index.

### Parallel Worker Requirements ⚡

These are called inside `foreach`/`doFuture` loops and **must** be exported
for `callr`/`multisession` backends to find them:

| Function | Called from | Why exported |
|----------|------------|--------------|
| `bootstrap_results` ⚡ | Bootstrap loop | Worker function |
| `bootstrap_ystar` ⚡ | Bootstrap loop | Bootstrap resampling |
| `count_boot_id` ⚡ | Bootstrap loop | Bootstrap counting |
| `get_dfRes` ⚡ | Bootstrap loop | Result extraction |
| `build_cox_formula` ⚡ | Bootstrap/CV workers | Formula construction |
| `fit_cox_models` ⚡ | Bootstrap/CV workers | Model fitting |
| `get_Cox_sg` ⚡ | Consistency splits | Cox model for subgroup |
| `get_split_hr_fast` ⚡ | Consistency splits | Fast HR calculation |
| `run_single_consistency_split` ⚡ | Consistency parallelism | Single split worker |
| `setup_parallel_SGcons` ⚡ | Consistency | Parallel setup |
| `sg_consistency_out` ⚡ | Consistency | Output formatting |
| `evaluate_consistency_twostage` ⚡ | Two-stage consistency | Worker |
| `wilson_ci` ⚡ | Workers | CI calculation |
| `early_stop_decision` ⚡ | Workers | Stopping rule |
| `evaluate_comparison` ⚡ | Search workers | Cut evaluation |
| `evaluate_cuts_once` ⚡ | Search workers | Batch cut evaluation |
| `sort_subgroups` ⚡ | Search workers | Subgroup sorting |
| `select_best_subgroup` ⚡ | Search workers | Selection |
| `analyze_subgroup` ⚡ | Search workers | Subgroup analysis |
| `assign_subgroup_membership` ⚡ | GRF workers | Membership assignment |
| `extract_subgroup` ⚡ | Workers | Subgroup extraction |
| `get_subgroup_membership` ⚡ | Workers | Membership lookup |
| `prepare_subgroup_data` ⚡ | Workers | Data prep |
| `evaluate_subgroup_consistency` ⚡ | Workers | Consistency check |
| `generate_bootstrap_synthetic` ⚡ | Bootstrap workers | Synthetic data |
| `generate_bootstrap_with_noise` ⚡ | Bootstrap workers | Noise bootstrap |
| `generate_gbsg_bootstrap_general` ⚡ | Bootstrap workers | GBSG bootstrap |

### Algorithm Internals (not parallel, but used by other exported functions)

| Function | Rationale |
|----------|-----------|
| `filter_call_args` | Argument filtering utility |
| `filter_by_lassokeep` | LASSO filter helper |
| `get_combinations_info` | Combination counting |
| `get_conf_force` | Forced confounder handling |
| `get_covs_in` | Covariate extraction |
| `process_conf_force_expr` | Expression processing |
| `create_grf_config` | GRF configuration builder |
| `validate_grf_data` | GRF data validation |
| `fit_causal_forest` | Causal forest fitting |
| `fit_policy_trees` | Policy tree fitting |
| `compute_node_metrics` | Node metric computation |
| `find_leaf_split` | Leaf split finding |
| `print_grf_details` | GRF diagnostic printing |

**Tier 2 total: ~40 functions**

---

## Tier 3: NORD or INTERNAL — Strong Candidates to Un-export (~29 functions)

These have no clear reason to be exported. They are not called by parallel
workers, not user-facing entry points, and serve purely as internal
implementation helpers. **Verify each is not called in parallel before
un-exporting.**

### Bitwise / Encoding Helpers
| Function | Current use | Recommendation |
|----------|-------------|----------------|
| `one.zero` | Bitwise helper for combination indexing | `@noRd` |
| `ztrail` | Trailing zero counter | `@noRd` |
| `acm.disjctif` | French-named disjunctive coding | `@noRd` — superseded by `dummy_encode` |
| `dummy` | Alias for `dummy_encode` | `@noRd` or deprecate |
| `dummy2` | Second alias for `dummy_encode` | `@noRd` or deprecate |

### Data Preparation Internals
| Function | Current use | Recommendation |
|----------|-------------|----------------|
| `get_cut_name` | Extract variable name from cut string | `@keywords internal` |
| `is_flag_continuous` | Check if cut is continuous | `@keywords internal` |
| `is_flag_drop` | Check if cut should be dropped | `@keywords internal` |
| `is.continuous` | Check variable continuity | `@keywords internal` |
| `cut_var` | Cut variable helper | `@keywords internal` |
| `detect_variable_types` | Variable type detection | `@keywords internal` |
| `add_id_column` | Add ID column to data frame | `@keywords internal` |

### Tree Cut Extraction
| Function | Current use | Recommendation |
|----------|-------------|----------------|
| `extract_all_tree_cuts` | Extract cuts from all trees | `@keywords internal` |
| `extract_selected_tree_cuts` | Extract selected cuts | `@keywords internal` |
| `extract_tree_cuts` | Base tree cut extraction | `@keywords internal` |
| `extract_idx_flagredundancy` | Redundancy flagging | `@keywords internal` |

### Subgroup Sorting / Filtering (if NOT called in parallel)
| Function | Current use | Recommendation |
|----------|-------------|----------------|
| `sort_subgroups_preview` | Preview sort results | `@keywords internal` |
| `remove_near_duplicate_subgroups` | Duplicate removal | `@keywords internal` |
| `remove_redundant_subgroups` | Redundancy removal | `@keywords internal` |

### Summary Table Variants
| Function | Current use | Recommendation |
|----------|-------------|----------------|
| `create_summary_table_compact` | Compact variant | `@keywords internal` |
| `create_summary_table_minimal` | Minimal variant | `@keywords internal` |
| `create_summary_table_presentation` | Presentation variant | `@keywords internal` |
| `create_summary_table_publication` | Publication variant | `@keywords internal` |
| `create_sample_size_table` | Sample size table | `@keywords internal` |
| `create_subgroup_summary_df` | Summary data frame | `@keywords internal` |

### Bootstrap / Formatting Internals
| Function | Current use | Recommendation |
|----------|-------------|----------------|
| `summarize_bootstrap_subgroups` | Bootstrap subgroup summary | `@keywords internal` |
| `summarize_bootstrap_events` | Bootstrap event summary | `@keywords internal` |
| `summarize_factor_presence_robust` | Factor presence | `@keywords internal` |
| `format_bootstrap_diagnostics_table` | Diagnostics table | `@keywords internal` |
| `format_bootstrap_timing_table` | Timing table | `@keywords internal` |
| `format_subgroup_summary_tables` | Summary tables | `@keywords internal` |
| `create_factor_summary_tables` | Factor tables | `@keywords internal` |
| `format_results` | General formatting | `@keywords internal` |
| `format_oc_results` | OC formatting | `@keywords internal` |

### Simulation Internals
| Function | Current use | Recommendation |
|----------|-------------|----------------|
| `get_dgm_with_output` | DGM wrapper | `@keywords internal` |
| `calibrate_k_inter` | k_inter calibration | Keep if user-facing, else `@keywords internal` |
| `find_k_inter_for_target_hr` | Target HR finding | Keep if user-facing, else `@keywords internal` |
| `validate_k_inter_effect` | Validation | `@keywords internal` |
| `sensitivity_analysis_k_inter` | Sensitivity analysis | Keep if user-facing |
| `compute_detection_probability` | Detection probability | `@keywords internal` |
| `find_required_sample_size` | Sample size finding | Keep if user-facing |
| `create_null_result` | Null result constructor | `@keywords internal` |
| `create_success_result` | Success result constructor | `@keywords internal` |
| `default_fs_params` | Default FS parameters | `@keywords internal` |
| `default_grf_params` | Default GRF parameters | `@keywords internal` |

### Miscellaneous
| Function | Current use | Recommendation |
|----------|-------------|----------------|
| `ci_est` | CI estimation helper | `@keywords internal` |
| `calc_cov` | Coverage calculation | `@keywords internal` |
| `get_targetEst` | Target estimand extraction | `@keywords internal` |
| `calculate_counts` | Count calculation | `@keywords internal` |
| `calculate_potential_hr` | Potential HR calculation | `@keywords internal` |
| `density_threshold_both` | Density threshold | `@keywords internal` |
| `find_quantile_for_proportion` | Quantile finder | `@keywords internal` |
| `qlow` | Lower quartile | `@keywords internal` |
| `qhigh` | Upper quartile | `@keywords internal` |
| `get_best_survreg` | Best survreg selection | `@keywords internal` |
| `compare_multiple_survreg` | Survreg comparison | Keep if user-facing |
| `cox_summary_batch` | Batch Cox summaries | `@keywords internal` |
| `cox_summary_legacy` | Legacy Cox wrapper | `@keywords internal` — deprecate |
| `cox_summary_vectorized` | Vectorized Cox | `@keywords internal` |
| `sens_text` | Sensitivity text | `@keywords internal` |
| `figure_note` | Figure note text | `@keywords internal` |
| `km_summary` | KM summary extraction | `@keywords internal` |
| `plot_subgroup` | Single subgroup plot | `@keywords internal` |
| `cv_summary_tables` | CV summary tables | Keep if user-facing |
| `cv_metrics_tables` | CV metrics tables | Keep if user-facing |
| `cv_summary_text` | CV text summary | Keep if user-facing |
| `cv_compare_results` | CV comparison | Keep if user-facing |
| `create_dgm_for_mrct` | MRCT DGM creation | Keep if user-facing |

---

## How to Implement

### Step 1: Tag functions with `@keywords internal`

For Tier 2 and Tier 3 `@keywords internal` candidates, change the roxygen:

```r
#' @export
my_internal_function <- function(...) { ... }
```

to:

```r
#' @keywords internal
#' @export
my_internal_function <- function(...) { ... }
```

This keeps the function exported (accessible via `forestsearch::my_function()`
and available to parallel workers), generates a help page, but **excludes it
from the pkgdown reference index** by default.

### Step 2: For `@noRd` candidates, fully un-export

```r
#' Brief description for maintainers
#' @noRd
my_helper <- function(...) { ... }
```

Remove `@export`. The function becomes truly internal — not in NAMESPACE,
no help page. **Only do this after verifying the function is not called
by parallel workers.**

### Step 3: Update pkgdown `_pkgdown.yml`

After tagging, the pkgdown reference sections should be simplified. Functions
with `@keywords internal` won't appear in the default reference index unless
explicitly listed. You can optionally add an "Internal" section:

```yaml
- title: internal
  contents:
    - has_keyword("internal")
```

### Step 4: Rebuild and verify

```r
devtools::document()
devtools::check()
pkgdown::build_reference_index()
```

---

## Summary

| Tier | Count | Action |
|------|-------|--------|
| **KEEP** (user-facing) | ~50 | No change |
| **INTERNAL** (parallel/infrastructure) | ~40 | Add `@keywords internal`, keep `@export` |
| **NORD/un-export** | ~29 | Remove `@export`, add `@noRd` or `@keywords internal` |

This reduces the **visible** public API from 119 to ~50 user-facing functions,
while preserving all functionality for parallel workers and power users who
know to call `forestsearch::internal_function()` directly.
