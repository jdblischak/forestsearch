#' Global Variables Declaration for ForestSearch Package
#'
#' Consolidated declarations to suppress R CMD check NOTEs for
#' non-standard evaluation (NSE) in data.table, ggplot2, and dplyr operations.
#'
#' These variables are used in contexts where column names are referenced
#' directly without quotes. R CMD check cannot detect these are valid
#' column references, so we declare them here to suppress false positive warnings.
#'
#' @keywords internal
#' @noRd

utils::globalVariables(c(


# ============================================================================
# Data.table special symbols
# ============================================================================
  ".",
  ".N",
  ".SD",
  ".I",
  ".GRP",
  ".BY",
  "..col_names",

# ============================================================================
# Future package
# ============================================================================
  ".Last.future.plan",

# ============================================================================
# Bootstrap & bias correction
# ============================================================================
  "boot",
  "boot_id",
  "iteration",
  "converged",
  "error_message",
  "time_elapsed",
  "H",
  "H_star",
  "H_bootstrap",
  "H_biasadj_1",
  "H_biasadj_2",
  "Hc_biasadj_1",
  "Hc_biasadj_2",

# ============================================================================
# Subgroup identification
# ============================================================================
  "Subgroup",
  "match_string",
  "treat.recommend",
  "flag.harm",
  "flag_harm",
  "M.1",
  "M.2",
  "M.3",
  "M.4",
  "M.5",
  "M.6",
  "M.7",
  "grp",
  "g_sg",
  "m_sg",

# ============================================================================
# Metrics: HR, consistency, sample size
# ============================================================================
  "HR",
  "hr",
  "hr_sg",
  "hr_individual",
  "Pcons",
  "Pcons_bin",
  "N",
  "N_sg",
  "N_sg_bin",
  "K",
  "K_sg",
  "E_sg",
  "Percent",
  "Percent_of_successful",
  "Rank",

# ============================================================================
# Summary table columns
# ============================================================================
  "Category",
  "Metric",
  "Value",
  "Count",
  "Base_Factor",
  "Factor_Definition",
  "Factor",
  "Consistency Range",
  "Size Range",
  "N_positions",
  "N_total",
  "Positions",
 "Position",

# ============================================================================
# Cross-validation
# ============================================================================
  "cv_index",
  "cvindex",

# ============================================================================
# ForestSearch parameters (used in NSE contexts)
# ============================================================================
  "est.scale",
  "insplit1",
  "confounders.name",
  "outcome.name",
  "event.name",
  "treat.name",

# ============================================================================
# Duplicate detection
# ============================================================================
  "dup_key",
  "dup_group",
  "dup_count",
  "dup_rank",

# ============================================================================
# Survival analysis & Cox models
# ============================================================================
  "treat",
  "event",
  "Y",
  "E",
  "Treat",
  "stratum",
  "z",
  "loghr_po",
  "theta_1",
  "theta_0",
  "treat_harm",

# ============================================================================
# AFT DGM / Simulation variables
# ============================================================================
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
  "z_age",
  "z_size",
  "z_nodes",
  "z_pgr",
  "z_er",
  "z_meno",
  "z_grade_1",
  "z_grade_2",
  "z_hormon",

# ============================================================================
# GBSG dataset column names (used in examples/tests)
# ============================================================================
  "id",
  "pid",
  "status",
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

# ============================================================================
# Model comparison
# ============================================================================
  "By_AIC",
  "By_BIC",
  "By_LogLik",
  "AIC_value",
  "BIC_value",
  "LogLik_value"
))


# ============================================================================
# Additional package-level imports
# ============================================================================
# These imports supplement those in forestsearch-package.R and ensure
# all commonly used functions are available without explicit namespacing.

#' @import data.table
#' @import survival
#' @importFrom stats AIC BIC model.frame vcov sd median quantile
#' @importFrom stats cor coef predict residuals fitted formula terms update
#' @importFrom stats pchisq qnorm qt pnorm pt
#' @importFrom stats rexp runif rnorm rbinom rgamma
#' @importFrom stats na.omit complete.cases aggregate
#' @importFrom graphics box plot.new plot lines points abline legend par
#' @importFrom graphics hist barplot mtext axis title rug text
#' @importFrom utils data object.size head tail str
#' @importFrom utils capture.output write.table read.table
#' @importFrom utils packageVersion sessionInfo
#' @importFrom splines ns
NULL

# ============================================================================
# Notes on base R functions (no import needed)
# ============================================================================
# The following functions are from base R and don't need importing:
# - Math: mean, sum, min, max, exp, log, sqrt, abs, round, floor, ceiling
# - Structure: length, nrow, ncol, dim, names, colnames, rownames
# - Creation: c, list, data.frame, matrix, as.matrix, as.data.frame
# - Manipulation: subset, which, order, sort, unique, duplicated
# - Text: paste, paste0, sprintf, cat, print
# - Control: if, else, for, while, return, stop, warning, message
# - Constants: TRUE, FALSE, NULL, NA, NaN, Inf
