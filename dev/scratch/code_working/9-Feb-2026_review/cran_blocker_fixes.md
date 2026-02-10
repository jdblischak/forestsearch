# ForestSearch: CRAN Blocker Fixes — Step-by-Step

## Priority Order

| # | Item | File(s) | Risk | Effort |
|---|------|---------|------|--------|
| 1 | ~~Duplicate globalVariables()~~ | ~~globals.R + 4 files~~ | ~~Critical~~ | ~~Done~~ |
| 2 | install.packages() in error messages | 3 files | Critical | 10 min |
| 3 | set.seed() in exported functions | 2 files | Critical | 20 min |
| 4 | Dead code: BOOTSTRAP_REQUIRED_FUNCTIONS | 2 files | Critical | 10 min |
| 5 | bootstrap_results() fragile globals | 1 file | High/Bug | 15 min |

---

## Fix 2: Remove `install.packages()` from Error Messages

**CRAN policy:** Packages must not suggest `install.packages()` in messages,
warnings, or errors. Use `requireNamespace()` checks only.

### 2a. `R/bootstrap_dofuture_setup_helpers.R` — `ensure_packages()`

**Find this** (around line 60):

```r
ensure_packages <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing) > 0) {
    stop(
      "Required package(s) not installed: ",
      paste(missing, collapse = ", "),
      "\nPlease install with: install.packages(c('",
      paste(missing, collapse = "', '"), "'))",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
```

**Replace with:**

```r
ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0) {
    stop(
      "Required package(s) not installed: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
```

**Changes:** Removed the `\nPlease install with:` line. Also switched `sapply`
→ `vapply` (safer, CRAN-preferred).

### 2b. `R/plot_subgroup_results_forestplot.R`

**Find this:**

```r
  if (!requireNamespace("forestploter", quietly = TRUE)) {
    stop("Package 'forestploter' is required. Install with: install.packages('forestploter')")
  }
```

**Replace with:**

```r
  if (!requireNamespace("forestploter", quietly = TRUE)) {
    stop("Package 'forestploter' is required but not installed.", call. = FALSE)
  }
```

### 2c. `R/bootstrap_summaries_helpers.R` — `format_bootstrap_table()`

**Find this:**

```r
 if (!requireNamespace("gt", quietly = TRUE)) {
    stop(
      "Package 'gt' is required for table formatting. ",
      "Install with: install.packages('gt')"
    )
  }
```

**Replace with:**

```r
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required for table formatting.", call. = FALSE)
  }
```

### 2d. Verify no others remain

Run in Terminal:

```bash
grep -rn "install.packages" R/
```

Should return **zero** results. If any others appear, apply the same pattern.

---

## Fix 3: Add `seed` Parameter to Exported Functions with `set.seed()`

**CRAN policy:** Exported functions must not call `set.seed()` without giving
the user control. The fix is to add a `seed` parameter with the current value
as the default, preserving backward compatibility.

### 3a. `R/bootstrap_calculations_helpers.R` — `bootstrap_ystar()`

**Find this:**

```r
#' Bootstrap Ystar Matrix
#'
#' Generates a bootstrap matrix for Ystar using parallel processing.
#'
#' @param df Data frame.
#' @param nb_boots Integer. Number of bootstrap samples.
#' @return Matrix of bootstrap samples.
#' @importFrom foreach foreach
#' @export

bootstrap_ystar <- function(df, nb_boots) {
  NN <- nrow(df)
  # do not modify seed below it need to align with main bootstrap
  # using manual seeding to allow reproducibility when qc-ing
  set.seed(8316951)
```

**Replace with:**

```r
#' Bootstrap Ystar Matrix
#'
#' Generates a bootstrap matrix for Ystar using parallel processing.
#'
#' @param df Data frame.
#' @param nb_boots Integer. Number of bootstrap samples.
#' @param seed Integer. Random seed for reproducibility. Default 8316951L.
#'   Must match the seed used in \code{\link{bootstrap_results}} to ensure
#'   bootstrap index alignment with the Ystar matrix.
#' @return Matrix of bootstrap samples (nb_boots x nrow(df)).
#' @importFrom foreach foreach
#' @export

bootstrap_ystar <- function(df, nb_boots, seed = 8316951L) {
  NN <- nrow(df)
  set.seed(seed)
```

### 3b. `R/bootstrap_analysis_dofuture.R` — `bootstrap_results()`

**Find this** (function signature + first lines):

```r
bootstrap_results <- function(fs.est, df_boot_analysis, cox.formula.boot,
                              nb_boots, show_three, H_obs, Hc_obs) {
  # =========================================================================
  # SECTION: INITIALIZE TIMING
  # =========================================================================
  t_start_bootstrap <- proc.time()[3]

  set.seed(8316951)
```

**Replace with:**

```r
bootstrap_results <- function(fs.est, df_boot_analysis, cox.formula.boot,
                              nb_boots, show_three, H_obs, Hc_obs,
                              seed = 8316951L) {
  # =========================================================================
  # SECTION: INITIALIZE TIMING
  # =========================================================================
  t_start_bootstrap <- proc.time()[3]

  set.seed(seed)
```

Also update the roxygen block for `bootstrap_results()`. Find `@param Hc_obs`
and add after it:

```r
#' @param seed Integer. Random seed for reproducibility. Default 8316951L.
#'   Must match the seed used in \code{\link{bootstrap_ystar}} to ensure
#'   bootstrap index alignment.
```

### 3c. Update the caller in `R/bootstrap_dofuture_main.R`

In `forestsearch_bootstrap_dofuture()`, both functions are called. The caller
doesn't need to change because both default to 8316951L. But for clarity you
can optionally make the alignment explicit:

**Find this:**

```r
  Ystar_mat <- bootstrap_ystar(fs.est$df.est, nb_boots)
```

**Replace with:**

```r
  boot_seed <- 8316951L
  Ystar_mat <- bootstrap_ystar(fs.est$df.est, nb_boots, seed = boot_seed)
```

And find this:

```r
  results <- bootstrap_results(
    fs.est = fs.est,
    df_boot_analysis = fs.est$df.est,
    cox.formula.boot = cox.formula.boot,
    nb_boots = nb_boots,
    show_three = show_three,
    H_obs = H_obs,
    Hc_obs = Hc_obs
  )
```

**Replace with:**

```r
  results <- bootstrap_results(
    fs.est = fs.est,
    df_boot_analysis = fs.est$df.est,
    cox.formula.boot = cox.formula.boot,
    nb_boots = nb_boots,
    show_three = show_three,
    H_obs = H_obs,
    Hc_obs = Hc_obs,
    seed = boot_seed
  )
```

### 3d. Also fix `forestsearch_tenfold()` in `R/forestsearch_cross-validation.R`

There is a `set.seed()` inside the foreach loop body. **Find this:**

```r
    set.seed(8316951 + 1000 * ksim)
```

This one is inside the foreach body (not the function signature), and
`seed = TRUE` in `.options.future` already handles parallel RNG. The manual
`set.seed()` inside the worker is safe per CRAN because the function-level
`.options.future = list(seed = TRUE)` controls reproducibility. However, for
consistency you can add a `seed` parameter to `forestsearch_tenfold()` too:

**Find the function signature:**

```r
forestsearch_tenfold <- function(fs_args, dfnew, Kfolds = 10, sims = 10,
                                  details = FALSE, ...) {
```

**Replace with:**

```r
forestsearch_tenfold <- function(fs_args, dfnew, Kfolds = 10, sims = 10,
                                  details = FALSE, seed = 8316951L, ...) {
```

Then **find:**

```r
    set.seed(8316951 + 1000 * ksim)
```

**Replace with:**

```r
    set.seed(seed + 1000L * ksim)
```

And add roxygen:

```r
#' @param seed Integer. Base random seed for fold shuffling. Default 8316951L.
#'   Each simulation uses seed + 1000 * ksim for reproducibility.
```

### 3e. Verify no uncontrolled seeds remain

```bash
grep -rn "set\.seed" R/ | grep -v "function\|param\|@param\|seed ="
```

Any remaining `set.seed(LITERAL_NUMBER)` without a `seed` parameter should be
addressed the same way.

---

## Fix 4: Remove Dead Code `BOOTSTRAP_REQUIRED_FUNCTIONS`

**Problem:** `BOOTSTRAP_REQUIRED_FUNCTIONS` is a static list that:
- References functions that don't exist: `"clean_data"`, `"prepare_data"`,
  `"run_grf"`, `"evaluate_subgroups"`, `"summarize_results"`, `"run_bootstrap"`,
  `"ci.est"`, `"count_boot_it"`
- Is never called at runtime (`get_bootstrap_exports()` does dynamic discovery)
- Is exported in `_pkgdown.yml` which would confuse users

### 4a. `R/bootstrap_dofuture_setup_helpers.R`

**Delete the entire block** (approximately 60 lines). Find and remove:

```r
#' Functions required in parallel bootstrap environment
#'
#' Organized by functional category for maintainability
BOOTSTRAP_REQUIRED_FUNCTIONS <- list(
  statistics = c(
    ...
  ),
  subgroup_analysis = c(
    ...
  ),
  data_prep = c(
    ...
  ),
  forestsearch_core = c(
    ...
  ),
  grf_policy = c(
    ...
  ),
  search_consistency = c(
    ...
  ),
  bootstrap_parallel = c(
    ...
  )
)
```

Everything from `#' Functions required in parallel bootstrap environment` down
to the closing `)` of the list.

### 4b. `_pkgdown.yml`

**Find and remove this line** from the "Internal Helpers" section:

```yaml
      - BOOTSTRAP_REQUIRED_FUNCTIONS
```

### 4c. Check for any references

```bash
grep -rn "BOOTSTRAP_REQUIRED_FUNCTIONS" R/ man/ _pkgdown.yml
```

Should return **zero** results after cleanup.

---

## Fix 5 (High/Bug): Clean Up `bootstrap_results()` `.options.future` Globals

**Problem:** The manual globals list in `bootstrap_results()` contains:
1. **Typo:** `"args_foresearch_call"` — missing 't', should be
   `"args_forestsearch_call"`
2. **Out-of-scope variables:** `"confounders_candidate"`,
   `"pconsistency.threshold"`, `"pconsistency.digits"`, `"hr.consistency"` —
   these aren't in `bootstrap_results()`'s local scope; they live inside
   `fs.est$args_call_all`

Since `get_bootstrap_exports()` already dynamically discovers all package
functions, and `globals = TRUE` auto-detects data variables, the manual list
is both buggy and redundant.

### 5a. `R/bootstrap_analysis_dofuture.R`

**Find this entire `.options.future` block:**

```r
    .options.future = list(
      seed = TRUE,
      globals = structure(TRUE, add = c(
      # Functions
      get_bootstrap_exports(),
      # Variables
      "fs.est", "df_boot_analysis", "cox.formula.boot","confounders_candidate",
      "H_obs", "Hc_obs", "nb_boots", "show_three", "args_foresearch_call","pconsistency.threshold",
      "pconsistency.digits","hr.consistency"
      ))
  ),
```

**Replace with:**

```r
    .options.future = list(
      seed = TRUE,
      globals = TRUE
    ),
```

**Why this is safe:**
- `globals = TRUE` is the doFuture default — it auto-detects all referenced
  variables through static code analysis
- `get_bootstrap_exports()` is only needed when manually listing globals; with
  automatic detection, all package functions are found via the namespace
- The typo and out-of-scope references are eliminated
- `seed = TRUE` still handles parallel RNG

**If you prefer to keep an explicit list** (belt-and-suspenders), use only
variables that actually exist in `bootstrap_results()`'s scope:

```r
    .options.future = list(
      seed = TRUE,
      globals = structure(TRUE, add = c(
        get_bootstrap_exports(),
        "fs.est", "df_boot_analysis", "cox.formula.boot",
        "H_obs", "Hc_obs", "nb_boots", "show_three"
      ))
    ),
```

---

## Verification Checklist

After all fixes, run these in order:

```r
# 1. Rebuild documentation
devtools::document()

# 2. Full check (the real test)
devtools::check()
# Target: 0 errors, 0 warnings, 0 notes

# 3. Verify no stray references
system('grep -rn "install.packages" R/')
system('grep -rn "BOOTSTRAP_REQUIRED_FUNCTIONS" R/ man/ _pkgdown.yml')
system('grep -rn "globalVariables" R/ | grep -v "globals.R"')

# 4. Spot-check set.seed
system('grep -rn "set\\.seed" R/')
# All hits should show seed = <parameter>, not a hardcoded literal

# 5. Install and verify pkgdown
devtools::install()
pkgdown::check_pkgdown()
```

---

## What NOT to Change

- **`BOOTSTRAP_REQUIRED_PACKAGES`** (the *packages* list) — keep it, it's used
  by `ensure_packages()` at runtime
- **`get_bootstrap_exports()`** — keep it, it's the correct dynamic approach
- **The seed value 8316951** — keep it as the default; it's the published
  reproducibility seed from your methodology
- **`filter_call_args()` pattern** — it's well-designed, don't touch it
- **S3 method definitions in `forestsearch_methods.R`** — verify consolidation
  separately (not a CRAN blocker)
