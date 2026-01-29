# dev-workflow.R
# Development workflow helper script for ForestSearch package
# Source this file to load helper functions, or run sections interactively

#
# Prior devtools-helpers file
library(usethis)
library(devtools)

# Step 1: Create the package structure
# This will create the template directory
# with Rproject setup
#usethis::create_package("../../GitHub/forestsearch")

# Copy functions into R directory
# Step 3: Add README and MIT license
usethis::use_readme_rmd(open = FALSE)
usethis::use_mit_license("Larry Leon")

# Step 4: Add dependencies to DESCRIPTION
# desc::desc_set_dep("survival", file = "DESCRIPTION")
# desc::desc_set_dep("ggplot2", file = "DESCRIPTION")
# desc::desc_set_dep("grf", file = "DESCRIPTION")
# desc::desc_set_dep("policytree", file = "DESCRIPTION")

pkgs <- c("data.table", "foreach", "future", "doFuture", "future.apply", "glmnet", "gt", "randomForest", "stringr", "survival", "grf",
          "ggplot2", "policytree", "future.callr")
for (pkg in pkgs) {
  desc::desc_set_dep(pkg, file = "DESCRIPTION")
}

#desc::desc_set_dep("weightedSurv", type = "Suggests", file = "DESCRIPTION")

#desc::desc_set_dep(package = "weightedsurv", type = "Imports")

desc::desc_set_remotes("larry-leon/weightedsurv")

usethis::use_package("patchwork", type = "Suggests")

# DiagrammeR for grf graph
usethis::use_package("DiagrammeR", type = "Suggests")

usethis::use_package("tidyr", type = "Suggests")

usethis::use_package("rlang")

usethis::use_package("dplyr")

# Required per mrct_simulation()

usethis::use_package("progressr")

usethis::use_package("cubature", type = "Suggests")


# Step 5: Generate documentation
# Also, run this if revising R files such as @importFrom

rm(list=ls())

.rs.restartR()

#source("apply_cran_fixes.R")
# removed from root after applied

# clean up old documentation
unlink("man/*.Rd")

devtools::document()

devtools::load_all()

devtools::check()

devtools::clean_dll()


# Run these before CRAN submission
devtools::check(cran = TRUE)     # Full CRAN validation

# End prior flow




# =============================================================================
# SETUP - Run once when starting development
# =============================================================================

# Install development dependencies
install.packages(c("devtools", "usethis", "roxygen2", "testthat", "pkgdown", "spelling"))

# Initialize package infrastructure (run once)
if (FALSE) {

  usethis::use_testthat()
  usethis::use_vignette("getting-started")
  usethis::use_readme_rmd()
  usethis::use_news_md()
  usethis::use_github_action_check_standard()
  usethis::use_spell_check()
  usethis::use_code_of_conduct()
}

# =============================================================================
# DAILY WORKFLOW - Common development tasks
# =============================================================================

# Load all package functions into environment (like library() but for development)
devtools::load_all()

# Generate/update documentation from roxygen2 comments
devtools::document()

# Run R CMD check (do this frequently!)
devtools::check()

# Quick check - skips some slower tests
devtools::check(args = "--no-examples")

# Run tests only
devtools::test()

# Run tests for a specific file
devtools::test_active_file()  # In RStudio, tests file related to current file
testthat::test_file("tests/testthat/test-bootstrap.R")

# =============================================================================
# NAMESPACE & IMPORTS
# =============================================================================
# Fix common namespace issues

# Check for missing imports
devtools::check() |>
  # Look for: "no visible global function definition for"
  # Look for: "no visible binding for global variable"

  # Add imports to NAMESPACE via roxygen2 tags:
  # @importFrom data.table data.table .N setnames
  # @importFrom stats median quantile sd

  # For pipe operator
  usethis::use_pipe()

# For data.table .SD, .N, etc. (add to R/package-imports.R)
if (FALSE) {
  cat('
#\' @importFrom data.table .N .SD .I .GRP .BY .EACHI
#\' @importFrom rlang .data
NULL
')
}

# =============================================================================
# CRAN CHECKS
# =============================================================================

# Full CRAN-like check
devtools::check(cran = TRUE)

# Check on multiple R versions (via R-hub)
devtools::check_rhub()

# Check on Windows (if you're on Mac/Linux)
devtools::check_win_devel()

# Spell check
spelling::spell_check_package()

# Check URLs in documentation
urlchecker::url_check()

# =============================================================================
# BUILDING & INSTALLING
# =============================================================================
# Build source package
devtools::build()

# Build binary package
devtools::build(binary = TRUE)
# Install package locally
devtools::install()

# Install with vignettes
devtools::install(build_vignettes = TRUE)

# =============================================================================
# VIGNETTES
# =============================================================================

# Build vignettes
devtools::build_vignettes()

# Preview a specific vignette
rmarkdown::render("vignettes/getting-started.Rmd")

# Move working vignette to official vignettes folder
# file.copy("dev/vignettes-working/my-analysis.Rmd", "vignettes/", overwrite = TRUE)

# =============================================================================
# DOCUMENTATION WEBSITE (pkgdown)
# =============================================================================

# Initialize pkgdown (run once)
usethis::use_pkgdown()

# Build documentation website locally
pkgdown::build_site()

# Preview site
pkgdown::preview_site()

# Build single article for quick preview
pkgdown::build_article("getting-started")

# =============================================================================
# VERSION & RELEASES
# =============================================================================

# Increment version number
usethis::use_version("patch")  # 0.1.0 -> 0.1.1
usethis::use_version("minor")  # 0.1.0 -> 0.2.0
usethis::use_version("major")  # 0.1.0 -> 1.0.0

# Development version
usethis::use_dev_version()     # 0.1.0 -> 0.1.0.9000

# Update NEWS.md with changes before release

# Create GitHub release
usethis::use_github_release()

# =============================================================================
# DATA
# =============================================================================

# Save internal data (not exported, available to package functions)
# usethis::use_data(internal_object, internal = TRUE)

# Save exported data (available to users)
# usethis::use_data(exported_dataset, overwrite = TRUE)

# Document datasets in R/data.R

# =============================================================================
# HELPER FUNCTIONS FOR THIS PROJECT
# =============================================================================

#' Quick development cycle
#' @description Document, load, and optionally test
quick_dev <- function(test = FALSE) {
  devtools::document()
  devtools::load_all()
  if (test) devtools::test()
  invisible(TRUE)
}

#' Check for common issues before committing
pre_commit_check <- function() {
  message("=== Documenting ===")
  devtools::document()


  message("\n=== Running tests ===")
  test_results <- devtools::test()

  message("\n=== Quick check ===
")
  devtools::check(args = c("--no-manual", "--no-vignettes"))

  message("\n=== Spell check ===")
  spelling::spell_check_package()

  invisible(TRUE)
}

#' Full pre-release check
pre_release_check <- function() {
  message("=== Full CRAN check ===")
  devtools::check(cran = TRUE)

  message("\n=== URL check ===")
  urlchecker::url_check()

  message("\n=== Reverse dependency check ===")
  # revdepcheck::revdep_check()  # If you have reverse dependencies

  invisible(TRUE)
}

#' Clean up generated files
clean_package <- function() {
  unlink("man", recursive = TRUE)
  unlink("Meta", recursive = TRUE)
  unlink("doc", recursive = TRUE)
  unlink("inst/doc", recursive = TRUE)
  unlink(list.files(pattern = "\\.tar\\.gz$"))
  message("Cleaned generated files. Run devtools::document() to regenerate.")
}

# =============================================================================
# PROJECT-SPECIFIC WORKFLOWS
# =============================================================================

#' Regenerate synthetic data for examples
regenerate_example_data <- function() {
  source("data-raw/generate-example-data.R")
  message("Example data regenerated")
}

#' Run bootstrap simulation for testing
run_test_simulation <- function(n_boots = 10) {
  devtools::load_all()
  # Your simulation code here
  message(sprintf("Completed %d bootstrap iterations", n_boots))
}

# =============================================================================
# TROUBLESHOOTING NOTES
# =============================================================================

# Common errors and fixes:
#
# 1. "no visible global function definition for 'xyz'"
#    -> Add @importFrom pkg xyz to the function's roxygen block
#
# 2. "no visible binding for global variable 'xyz'"
#    -> For data.table: add utils::globalVariables(c("xyz")) in R/globals.R
#    -> Or use .data$xyz with @importFrom rlang .data
#
# 3. "library() call in package code"
#    -> Replace library(pkg) with pkg::function() or @importFrom
#
# 4. Unicode characters causing issues
#    -> Replace with ASCII: >= instead of ≥, use \u2265 in strings
#
# 5. Examples taking too long
#    -> Wrap slow examples in \donttest{} or \dontrun{}

# =============================================================================
# USEFUL REFERENCES
# =============================================================================

# R Packages book: https://r-pkgs.org/
# CRAN policies: https://cran.r-project.org/web/packages/policies.html
# Writing R Extensions: https://cran.r-project.org/doc/manuals/R-exts.html

message("dev-workflow.R loaded. Key functions: quick_dev(), pre_commit_check()")


library(usethis)
library(devtools)

# Step 1: Create the package structure
# This will create the template directory
# with Rproject setup
#usethis::create_package("../../GitHub/forestsearch")

# Copy functions into R directory
# Step 3: Add README and MIT license
usethis::use_readme_rmd(open = FALSE)
usethis::use_mit_license("Larry Leon")

# Step 4: Add dependencies to DESCRIPTION
# desc::desc_set_dep("survival", file = "DESCRIPTION")
# desc::desc_set_dep("ggplot2", file = "DESCRIPTION")
# desc::desc_set_dep("grf", file = "DESCRIPTION")
# desc::desc_set_dep("policytree", file = "DESCRIPTION")

pkgs <- c("data.table", "foreach", "future", "doFuture", "future.apply", "glmnet", "gt", "randomForest", "stringr", "survival", "grf",
          "ggplot2", "policytree", "future.callr")
for (pkg in pkgs) {
  desc::desc_set_dep(pkg, file = "DESCRIPTION")
}

#desc::desc_set_dep("weightedSurv", type = "Suggests", file = "DESCRIPTION")

#desc::desc_set_dep(package = "weightedsurv", type = "Imports")

desc::desc_set_remotes("larry-leon/weightedsurv")

usethis::use_package("patchwork", type = "Suggests")

# DiagrammeR for grf graph
usethis::use_package("DiagrammeR", type = "Suggests")

usethis::use_package("tidyr", type = "Suggests")

usethis::use_package("rlang")

usethis::use_package("dplyr")

# Step 5: Generate documentation
# Also, run this if revising R files such as @importFrom

rm(list=ls())

.rs.restartR()

# clean up old documentation
unlink("man/*.Rd")

devtools::document()

devtools::load_all()

devtools::check()

devtools::clean_dll()


# Run these before CRAN submission
devtools::check(cran = TRUE)     # Full CRAN validation
devtools::check_rhub()            # Check on multiple platforms
devtools::check_win_devel()       # Windows checks


tools::showNonASCIIfile("~/Documents/GitHub/forestsearch/R/summary_utility_functions.R")



# gitignore additions
#.Rproj.user
#vignettes/results/
#  dev/private/
#  *_files/
# vignettes/working/



# Stop tracking files in GitHub

#git rm -r --cached *_files/

# Issue with trying to remove weightedSurv and replace with weightedsurv
# Remove from search path and unload
library(forestsearch)  # load it first if not loaded
detach("package:forestsearch", unload = TRUE, force = TRUE)
# Clear workspace and restart
rm(list = ls())
.rs.restartR()

# Notes
# Incorporate in AI prompts when documenting
# Every \item in a \describe{} block must have both a label and a description.
# Do not leave a lone \item{...} without {...} after it.
# Do not next \item inside another \item
# Use plain text, dashes, or a single paragraph for subpoints.
#
#
# León LF, Jemielita T, Guo Z, Marceau West R, Anderson KM.
# Exploratory subgroup identification in the heterogeneous Cox model: A relatively simple procedure.
# Statistics in Medicine. 2024; 43(20): 3921-3942. doi: 10.1002/sim.10163


#roxygen2::roxygenise()


# Remove and reinstall the package
remove.packages("forestsearch")
# Clear package cache
.libPaths()  # Find your library path
# Then reinstall from your source
# Either from GitHub or local source:
#devtools::install_github("larry-leon/forestsearch")
# OR
devtools::load_all()  # if developing locally

#devtools::install()

# run this in terminal (next to console [go to tools terminal tab])
#git pull --no-rebase



gitcreds::gitcreds_set()

usethis::use_git()

usethis::use_github()

# If functions are in namespace but not directly loaded
# devtools::load_all()

# Or  access hidden files:  mypackage::my_function




#library(codetools)
#codetools::findGlobals(forestsearch_bootstrap_dofuture, merge = FALSE)$functions

codetools::findGlobals(forestsearch, merge = FALSE)$functions

codetools::findGlobals(subgroup.search, merge = FALSE)$functions

codetools::findGlobals(get_FSdata, merge = FALSE)$functions


check <- c("calc_cov",  "calculate_counts", "analyze_subgroups", "calculate_potential_hr","ci.est","count.id","CV_sgs",
           "cox_summary","df_counting","double_robust_scores", "extract_subgroup","format_results", "get_targetEst","getci_Cox",
           "getCIs","grf.estimates.out","hrCI_format","km_summary","n_pcnt","plot_subgroup","plot_weighted_km",
           "prepare_subgroup_data","quiet","rmst_calculation","sg_tables","sort_subgroups","SummaryStat","var_summary",
           "get_FSdata", "dummy","run_bootstrap",
           "forestsearch", "forestsearch_bootstrap_dofuture","get_combinations_info",
           "get_dfpred",
           "grf.subg.harm.survival",
           "subgroup.search",
           "subgroup.consistency",
           "lasso_selection",
           "get_Cox_sg",
           "get_conf_force",
           "filter_by_lassokeep",
           "is.continuous",
           "process_conf_force_expr",
           "is_flag_continuous",
           "is_flag_drop", "acm.disjctif",  "acm.util.df2", "acm.util.df", "dummy2","ztrail","one.zero",
           "get_dfRes", "get_subgroup_membership",
           "SG_tab_estimates",
           "prepare_data",
           "run_grf",
           "evaluate_subgroups",
           "summarize_results",
           "clean_data", "qlow", "qhigh","FS_labels","thiscut","get_cut_name",
           "bootstrap_results", "remove_redundant_subgroups", "sg_consistency_out","get_split_hr","cut_var",
           "bootstrap_ystar", "ensure_packages", "fit_cox_models", "build_cox_formula",
           "format_CI","setup_parallel_SGcons", "get_covs_in", "extract_idx_flagredundancy"
)

duplicated_elements <- check[duplicated(check)]
duplicated_elements


# Per Claude

# In your local forestsearch repo
# git checkout -b fix/must-fix-implementations
#
# # Copy the new files to R/ directory (as shown above)
#
# # Stage changes
# git add R/bootstrap_helpers.R
# git add R/summary_utility_functions.R
# git add R/improved_grf_functions.r
# git add R/input_validation_utils.R

# Commit
# git commit -m "Implement must-fix recommendations
#
# - Add comprehensive input validation to all functions
# - Fix division-by-zero issues throughout
# - Standardize variable naming conventions
# - Add input_validation_utils.R with validation helpers
# - Maintain 100% backward compatibility"
#
# # Push to GitHub
# git push origin fix/must-fix-implementations

# Then create a Pull Request on GitHub


# Generate documentation
devtools::document()

# Load your updated package
devtools::load_all()


# Run tests
devtools::test()

# Check package
devtools::check()

# Search all R files in your package
# grep -rn "enumerate" ~/path/to/ForestSearch/R/ --include="*.R"
#
# # Or from the package root
# cd ~/path/to/ForestSearch
# grep -rn "enumerate" R/
#
#
#   How to Update Project Files
#
# Go to your Claude Project (in the left sidebar, click on the project name)
# Open Project Settings (click the gear icon or "Project settings")
# Remove outdated files:
#
#   Find the files in the "Project knowledge" or "Files" section
# Delete the old versions
#
#
# Upload current files:
#
#   Click "Add content" or "Upload files"
# Select the current versions from your R directory
#
#
# Return to this chat - the new files will be available in /mnt/project
#

# Finding Your Project Knowledge (Files)
# You'll find the project knowledge base on the right side of your project's main page. Click on the "+" button to add content to the project.
# Step-by-Step:
#
#   Go to your Project page:
#
#   Click on "Projects" in the left sidebar (or go to claude.ai/projects)
# Click on your ForestSearch project
#
#
# Look to the RIGHT side of the screen:
#
#   You should see a "Project Knowledge" or "Knowledge" panel
# Your 6 uploaded files should be listed there
#
#
# To remove old files:
#
#   Hover over a file name
# Click the trash/delete icon or "×" that appears
#
#
# To add updated files:
#
#   Click the "+" button in that panel
# Select your current R files from your computer
