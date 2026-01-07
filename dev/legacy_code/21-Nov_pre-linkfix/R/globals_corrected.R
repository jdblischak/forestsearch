#' Global Variables Declaration for ForestSearch Package
#' 
#' This file declares global variables to avoid R CMD check warnings.
#' Place this file in your R/ directory.
#' 
#' @keywords internal
#' @noRd

# Suppress R CMD check notes about undefined global variables
# These are primarily used in data.table operations and other contexts
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
  
  # ============================================================================
  # Future package related
  # ============================================================================
  ".Last.future.plan",
  
  # ============================================================================
  # Variables from bootstrap analysis functions
  # ============================================================================
  "flag_harm",
  "Pcons",
  "hr_sg",
  "N_sg",
  "K_sg",
  "Subgroup",
  "M.1",
  "M.2",
  "M.3",
  "M.4",
  "M.5",
  "Percent_of_successful",
  "Rank",
  "Pcons_bin",
  "N_sg_bin",
  
  # ============================================================================
  # Variables from summary table functions
  # ============================================================================
  "Consistency Range",
  "Size Range",
  "Factor",
  "N_positions",
  "N_total",
  "Positions",
  "Count",
  "Percent",
  "Position",
  "Base_Factor",
  "Factor_Definition",
  "Metric",
  "Value",
  
  # ============================================================================
  # Variables from survival analysis functions
  # ============================================================================
  "treat",
  "event",
  "y_sim",
  "event_sim",
  "t_true",
  "c_time",
  "lin_pred_1",
  "lin_pred_0",
  "lin_pred_obs",
  "treat_harm",
  "hr_individual",
  "lin_pred_cens_1",
  "lin_pred_cens_0",
  
  # ============================================================================
  # Variables from AFT DGM functions
  # ============================================================================
  "y",
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
  # Variables from factor processing
  # ============================================================================
  "desc",
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
  
  # ============================================================================
  # Additional variables that might appear in your functions
  # ============================================================================
  "H",
  "H_star",
  "H_bootstrap",
  "boot_id",
  "iteration",
  "converged",
  "error_message",
  "time_elapsed",
  
  # ============================================================================
  # Variables from cox_cs_fit and other functions
  # ============================================================================
  "Y",
  "E",
  "Treat",
  "stratum",
  "z",
  "loghr_po",
  "theta_1",
  "theta_0",
  
  # ============================================================================
  # Variables from data.table operations in summaries
  # ============================================================================
  "By_AIC",
  "By_BIC",
  "By_LogLik",
  "AIC_value",
  "BIC_value",
  "LogLik_value"
))

# Package documentation (NULL is conventional for package-level documentation)
# This also serves to document package-level imports

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

# Note: The following functions are from base R and don't need importing:
# mean, sum, min, max, exp, log, sqrt, abs, round, floor, ceiling
# length, nrow, ncol, dim, names, colnames, rownames
# c, list, data.frame, matrix, as.matrix, as.data.frame
# subset, which, order, sort, unique, duplicated
# paste, paste0, sprintf, cat, print
# if, else, for, while, return, stop, warning, message
# TRUE, FALSE, NULL, NA, NaN, Inf
