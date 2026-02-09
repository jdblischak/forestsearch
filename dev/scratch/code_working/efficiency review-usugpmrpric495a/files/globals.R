# =============================================================================
# ForestSearch Package - Consolidated Global Variables
# =============================================================================
#
# CRAN COMPLIANCE NOTE:
# All globalVariables() declarations are consolidated in this single file.
# This ensures maintainability and avoids duplicate declarations.
#
# DO NOT add globalVariables() in other R files.
# Add new variables to the appropriate section below.
#
# =============================================================================

#' @import data.table
#' @import survival
#' @importFrom stats AIC BIC model.frame vcov sd median quantile
#' @importFrom stats cor coef predict residuals fitted formula terms update
#' @importFrom stats pchisq qnorm qt pnorm pt
#' @importFrom stats rexp runif rnorm rbinom rgamma
#' @importFrom stats na.omit complete.cases aggregate
#' @importFrom stats as.formula setNames
#' @importFrom graphics box plot.new plot lines points abline legend par
#' @importFrom graphics hist barplot mtext axis title rug text
#' @importFrom utils data object.size head tail str
#' @importFrom utils capture.output write.table read.table
#' @importFrom utils packageVersion sessionInfo hasName
#' @importFrom splines ns
NULL

# =============================================================================
# Notes on Base R Functions (No Import Needed)
# =============================================================================
# The following functions are from base R and don't require importing:
# - Math: mean, sum, min, max, exp, log, sqrt, abs, round, floor, ceiling
# - Structure: length, nrow, ncol, dim, names, colnames, rownames
# - Creation: c, list, data.frame, matrix, as.matrix, as.data.frame
# - Manipulation: subset, which, order, sort, unique, duplicated
# - Text: paste, paste0, sprintf, cat, print
# - Control: if, else, for, while, return, stop, warning, message
# - Logic: is.null, is.na, is.numeric, is.character, is.factor, is.data.frame
# - Constants: TRUE, FALSE, NULL, NA, NaN, Inf

# =============================================================================
# Global Variables Declaration
# =============================================================================
# These variables are used in non-standard evaluation (NSE) contexts such as
# data.table operations, formula construction, and tidyverse-style programming.
# Declaring them here suppresses R CMD check NOTEs.

utils::globalVariables(c(

# =============================================================================
# Summary Table Columns
# =============================================================================
  "Category",
  "Metric",
  "Value",
  "Count",
  "Percent",
  "Base_Factor",
  "Factor_Definition",
  "Factor",
  "Consistency Range",
  "Size Range",
  "N_positions",
  "N_total",
  "Positions",
  "Position",
  "Rank",

# =============================================================================
# Cross-Validation Variables
# =============================================================================
  "cv_index",
  "cvindex",
  "fold",
  "fold_id",

# =============================================================================
# ForestSearch Parameters (Used in NSE Contexts)
# =============================================================================
  "est.scale",
  "insplit1",
  "confounders.name",
  "outcome.name",
  "event.name",
  "treat.name",
  "id.name",

# =============================================================================
# Duplicate Detection Variables
# =============================================================================
  "dup_key",
  "dup_group",
  "dup_count",
  "dup_rank",

# =============================================================================
# Survival Analysis & Cox Model Variables
# =============================================================================
  "treat",
  "event",
  "Y",
  "E",
  "Treat",
  "stratum",
  "Strata",
  "z",
  "loghr_po",
  "theta_1",
  "theta_0",
  "treat_harm",
  "time",
  "status",
  "hr",
  "HR",
  "se",
  "lower",
  "upper",
  "logHR",
  "selogHR",
  "pvalue",
  "conf.int",

# =============================================================================
# AFT DGM / Simulation Variables
# =============================================================================
  "y",
  "y_sim",
  "event_sim",
  "t_true",
  "c_time",
  "lin_pred_1",
  "lin_pred_0",
  "lin_pred_obs",
  "lin_pred_cens_1",
  "lin_pred_cens_0",
  "epsilon",
  "mu",
  "sigma",
  "tau",
  "gamma",

# =============================================================================
# GBSG Dataset Column Names
# =============================================================================
  "id",
  "pid",
  "rfstime",
  "age",
  "size",
  "nodes",
  "pgr",
  "er",
  "meno",
  "grade",
  "hormon",
  "desc",

# =============================================================================
# Standardized GBSG Covariate Names (z_* prefix)
# =============================================================================
  "z_age",
  "z_size",
  "z_nodes",
  "z_pgr",
  "z_er",
  "z_meno",
  "z_grade_1",
  "z_grade_2",
  "z_hormon",

# =============================================================================
# GBSG DGM-Specific Variables
# =============================================================================
  "zh",
  "flag.harm",
  "flag_harm",
  "v1", "v2", "v3", "v4", "v5", "v6", "v7",
  "z1", "z2", "z3", "z4", "z5",
  "grade3",
  "lin.conf.true",
  "lin1.conf",
  "lin0.conf",
  "linC1.conf",
  "linC0.conf",
  "hlin.conf.1",
  "hlin.conf.0",
  "hlin.ratio",
  "h1.potential",
  "h0.potential",
  "Ts",
  "es",
  "t.sim",

# =============================================================================
# Operating Characteristics Analysis Variables
# =============================================================================
  "any.H",
  "ppv",
  "npv",
  "sensitivity",
  "specificity",
  "size.H",
  "size.Hc",
  "hr.H.true",
  "hr.H.hat",
  "hr.Hc.true",
  "hr.Hc.hat",
  "hr.itt",
  "hr.adj.itt",
  "p.cens",
  "taumax",
  "analysis",
  "sim",
  "aa",
  "sg_hat",
  "ahr.H.true",
  "ahr.Hc.true",
  "ahr.H.hat",
  "ahr.Hc.hat",

# =============================================================================
# MRCT Simulation Variables
# =============================================================================
  "hr_test",
  "hr_sg",
  "any_found",
  "sg_found",
  "hr_sg_null",
  "regAflag",
  "sg_le85",
  "regAflag2",
  "regAflag3",
  "found",
  "sg_biomarker",
  "sg_age",
  "sg_male",
  "sg_ecog",
  "sg_histology",
  "sg_CTregimen",
  "sg_region",
  "sg_surgery",
  "sg_prior_treat",
  "est",
  "region_var",
  "z_regA",

# =============================================================================
# Subgroup Consistency Variables
# =============================================================================
  "Pcons",
  "N",
  "K",
  "g",
  "m",
  "d0",
  "d1",
  "grp",
  "sg.harm",
  "sg.harm.id",
  "sg.harm_label",
  "treat.recommend",
  "M.1", "M.2", "M.3", "M.4", "M.5", "M.6", "M.7",

# =============================================================================
# Bootstrap Analysis Variables
# =============================================================================
  "boot",
  "iteration",
  "H_obs",
  "Hc_obs",
  "H_biasadj_2",
  "Hc_biasadj_2",
  "H_boot",
  "Hc_boot",
  "events_Hstar_0",
  "events_Hstar_1",
  "events_Hcstar_0",
  "events_Hcstar_1",
  "tmins_search",
  "tmins_iteration",
  "n_found",
  "success_rate",

# =============================================================================
# Model Comparison Variables
# =============================================================================
  "By_AIC",
  "By_BIC",
  "By_LogLik",
  "AIC_value",
  "BIC_value",
  "LogLik_value",
  "model_name",
  "converged",

# =============================================================================
# GRF-Related Variables
# =============================================================================
  "vi.grf",
  "importance",
  "variable",
  "grf_cuts",
  "tau.hat",
  "predictions",

# =============================================================================
# Plotting Variables
# =============================================================================
  "Subgroup",
  "n_treat",
  "n_control",
  "low",
  "hi",
  "x",
  "xmin",
  "xmax",
  "ymin",
  "ymax",

# =============================================================================
# Data.table Special Symbols
# =============================================================================
# Note: .SD, .N, .I, .GRP are data.table special symbols
# They are handled by data.table imports, but including for documentation
  ".",
  ".SD",
  ".SDcols",
  ".N",
  ".I",
  ".GRP",
  ".BY",

# =============================================================================
# Foreach Loop Variables
# =============================================================================
  "i",
  "j",
  "k",
  "split",
  "candidate"

))

# =============================================================================
# Package Constants (Not Variables, But Useful to Define Here)
# =============================================================================

# Default random seed for reproducibility
FORESTSEARCH_DEFAULT_SEED <- 8316951L

# Minimum recommended bootstrap iterations
BOOTSTRAP_MIN_RECOMMENDED <- 100L
BOOTSTRAP_MIN_PRODUCTION <- 500L

# Days per month (for time conversions)
DAYS_PER_MONTH <- 30.4375

# Default super-population size for DGM
DEFAULT_N_SUPER <- 5000L

# Consistency threshold precision
PCONSISTENCY_DIGITS <- 6L
