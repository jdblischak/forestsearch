# =============================================================================
# ForestSearch Helper Functions
# =============================================================================
#
# General utility functions for ForestSearch package.
# Includes data preparation, subgroup assignment, and validation helpers.
#
# =============================================================================

# =============================================================================
# Required Package Check
# =============================================================================

#' Check and Report Missing Packages
#'
#' Verifies that required packages are installed. Returns missing package names.
#'
#' @param packages Character vector of package names to check.
#' @param quiet Logical. Suppress messages? Default FALSE.
#'
#' @return Character vector of missing package names (empty if all installed).
#'
#' @keywords internal
check_required_packages <- function(packages, quiet = FALSE) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]

  if (length(missing) > 0L && !quiet) {
    warning(
      "Missing required packages: ", paste(missing, collapse = ", "),
      "\nInstall with: install.packages(c('", paste(missing, collapse = "', '"), "'))",
      call. = FALSE
    )
  }

  missing
}


# =============================================================================
# ID Column Management
# =============================================================================

#' Add ID Column to Data Frame
#'
#' Ensures that a data frame has a unique ID column. Creates one if missing.
#'
#' @param df_analysis Data frame to which the ID column will be added.
#' @param id_name Character. Name of the ID column (default: "id").
#'
#' @return Data frame with the ID column added if necessary.
#'
#' @examples
#' df <- data.frame(x = 1:5, y = letters[1:5])
#' df <- add_id_column(df)
#' # df now has column 'id' with values 1:5
#'
#' @export
add_id_column <- function(df_analysis, id_name = "id") {
  if (!id_name %in% names(df_analysis)) {
    df_analysis[[id_name]] <- seq_len(nrow(df_analysis))
  }
  df_analysis
}


# =============================================================================
# Subgroup Assignment Functions
# =============================================================================

#' Generate Prediction Dataset with Subgroup Treatment Recommendation
#'
#' Creates a prediction dataset with a treatment recommendation flag based on
#' subgroup definition. Uses vectorized logic instead of eval(parse()) for
#' safety and clarity.
#'
#' @param df_predict Data frame for prediction (test or validation set).
#' @param sg_harm Character vector of subgroup-defining covariate names.
#'   These should be binary (0/1) indicator variables in the dummy-coded data.
#' @param dummy_version Integer. Version of dummy coding to apply:
#'   1 uses \code{dummy()}, 2 uses \code{dummy2()}. Default: 1.
#'
#' @return Data frame with added \code{treat_recommend} column:
#'   \itemize{
#'     \item 0 = do not recommend treatment (in harm subgroup)
#'     \item 1 = recommend treatment (not in harm subgroup)
#'   }
#'
#' @details
#' The harm subgroup is defined as observations where ALL variables in
#' \code{sg_harm} equal 1. Observations not meeting all criteria are
#' assigned to the no-harm (complement) subgroup.
#'
#' @examples
#' \dontrun{
#' # After ForestSearch identifies sg_harm = c("age_gt50", "stage_high")
#' df_pred <- get_dfpred(test_data, sg_harm = c("age_gt50", "stage_high"))
#'
#' # Patients with age_gt50==1 AND stage_high==1 get treat_recommend=0
#' # All others get treat_recommend=1
#' }
#'
#' @export
get_dfpred <- function(df_predict, sg_harm, dummy_version = 1L) {

  # Validate inputs
  if (!is.data.frame(df_predict)) {
    stop("df_predict must be a data frame")
  }

  if (length(sg_harm) == 0L || is.null(sg_harm)) {
    stop("sg_harm must be a non-empty character vector")
  }

  # Apply dummy coding
  df_pred <- if (dummy_version == 1L) {
    dummy(df_predict)
  } else {
    dummy2(df_predict)
  }

  # Validate that all sg_harm variables exist in dummy-coded data
  missing_vars <- setdiff(sg_harm, names(df_pred))
  if (length(missing_vars) > 0L) {
    stop(
      "sg_harm variables not found in dummy-coded data: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  # Identify harm subgroup using vectorized AND logic
  # All sg_harm columns must equal 1 for membership in harm group
  in_harm <- Reduce(`&`, lapply(sg_harm, function(col) {
    df_pred[[col]] == 1L
  }))

  # Assign treatment recommendations
  # treat_recommend = 0: do not treat (in harm subgroup)
  # treat_recommend = 1: treat (not in harm subgroup)
  df_pred$treat_recommend <- ifelse(in_harm, 0L, 1L)

  df_pred
}


#' Get Subgroup Membership Indicator
#'
#' Determines subgroup membership based on multiple indicator variables.
#'
#' @param df Data frame with indicator variables.
#' @param indicator_vars Character vector of variable names (all must equal 1).
#'
#' @return Logical vector indicating subgroup membership.
#'
#' @export
get_subgroup_membership <- function(df, indicator_vars) {
  if (length(indicator_vars) == 0L) {
    return(rep(TRUE, nrow(df)))
  }

  # Validate all variables exist
  missing <- setdiff(indicator_vars, names(df))
  if (length(missing) > 0L) {
    stop("Variables not found in data: ", paste(missing, collapse = ", "))
  }

  # Vectorized AND across all indicators
  Reduce(`&`, lapply(indicator_vars, function(v) df[[v]] == 1L))
}


# =============================================================================
# Parameter Extraction Utilities
# =============================================================================

#' Get Parameter with Default Fallback
#'
#' Safely extracts a parameter from a list, returning a default if not found.
#'
#' @param args_list List of arguments/parameters.
#' @param param_name Character. Name of parameter to extract.
#' @param default_value Default value if parameter is NULL or missing.
#'
#' @return Parameter value or default.
#'
#' @examples
#' args <- list(a = 1, b = NULL, c = 3)
#' get_param(args, "a", 0)  # Returns 1
#' get_param(args, "b", 0)  # Returns 0 (b is NULL)
#' get_param(args, "d", 0)  # Returns 0 (d doesn't exist)
#'
#' @keywords internal
get_param <- function(args_list, param_name, default_value) {
  if (hasName(args_list, param_name) && !is.null(args_list[[param_name]])) {
    args_list[[param_name]]
  } else {
    default_value
  }
}


#' Filter Call Arguments to Relevant Parameters
#'
#' Extracts only the parameters relevant to a specific function from a
#' larger argument list.
#'
#' @param args_full List of all arguments.
#' @param valid_params Character vector of valid parameter names.
#'
#' @return List containing only valid parameters.
#'
#' @keywords internal
filter_call_args <- function(args_full, valid_params) {
  args_full[intersect(names(args_full), valid_params)]
}


# =============================================================================
# Covariate Processing
# =============================================================================

#' Get Covariates Included in Subgroup Definition
#'
#' Extracts the base covariate names from a subgroup factor definition.
#'
#' @param factor_labels Character vector of factor labels
#'   (e.g., "age_gt50", "stage==2").
#'
#' @return Character vector of unique base covariate names.
#'
#' @keywords internal
get_covs_in <- function(factor_labels) {
  # Extract covariate names by removing comparison operators and values
  covs <- gsub("[<>=!]+.*$", "", factor_labels)
  covs <- gsub("_[^_]+$", "", covs)  # Remove suffixes like _gt50

  unique(covs)
}


#' Get Cut Variable Name from Expression
#'
#' Parses a cut expression to extract the variable name.
#'
#' @param cut_expr Character. Cut expression (e.g., "age <= 50").
#'
#' @return Character. Variable name.
#'
#' @keywords internal
get_cut_name <- function(cut_expr) {
  # Remove whitespace and extract variable name before operator
  expr_clean <- trimws(cut_expr)
  parts <- strsplit(expr_clean, "\\s*[<>=!]+\\s*")[[1L]]
  trimws(parts[1L])
}


# =============================================================================
# Combination Information
# =============================================================================

#' Get Combinations Information
#'
#' Calculates the number of possible k-factor combinations from n factors.
#'
#' @param n_factors Integer. Number of available factors.
#' @param maxk Integer. Maximum combination size.
#'
#' @return List with:
#'   \item{n_combinations}{Total number of combinations}
#'   \item{by_k}{Vector of combinations per k value}
#'
#' @examples
#' # 10 factors, up to 2-way combinations
#' get_combinations_info(10, 2)
#' # Returns list(n_combinations = 55, by_k = c(10, 45))
#'
#' @export
get_combinations_info <- function(n_factors, maxk = 2L) {
  by_k <- vapply(seq_len(maxk), function(k) {
    choose(n_factors, k)
  }, numeric(1L))

  names(by_k) <- paste0("k=", seq_len(maxk))

  list(
    n_combinations = sum(by_k),
    by_k = by_k
  )
}


# =============================================================================
# Formatting Functions
# =============================================================================

#' ForestSearch Factor Labels
#'
#' Creates formatted labels for factors based on cut expressions.
#'
#' @param cut_expr Character. Cut expression.
#' @param label_map Named vector mapping variables to labels (optional).
#'
#' @return Character. Formatted label.
#'
#' @keywords internal
FS_labels <- function(cut_expr, label_map = NULL) {
  var_name <- get_cut_name(cut_expr)

  if (!is.null(label_map) && var_name %in% names(label_map)) {
    # Replace variable name with label
    gsub(var_name, label_map[var_name], cut_expr, fixed = TRUE)
  } else {
    cut_expr
  }
}


#' Format CI as String
#'
#' Formats a confidence interval as "est (lower, upper)".
#'
#' @param est Point estimate.
#' @param lower Lower bound.
#' @param upper Upper bound.
#' @param digits Number of decimal places.
#'
#' @return Formatted string.
#'
#' @export
format_CI <- function(est, lower, upper, digits = 2L) {
  sprintf("%.*f (%.*f, %.*f)", digits, est, digits, lower, digits, upper)
}


#' Format Results for Display
#'
#' Creates formatted strings from numeric results.
#'
#' @param x Numeric vector or single value.
#' @param digits Number of decimal places.
#' @param percent Logical. Format as percentage?
#'
#' @return Formatted character string(s).
#'
#' @export
format_results <- function(x, digits = 2L, percent = FALSE) {
  if (percent) {
    sprintf("%.*f%%", digits, x * 100)
  } else {
    sprintf("%.*f", digits, x)
  }
}


# =============================================================================
# Quantile Functions
# =============================================================================

#' Lower Quantile
#'
#' @param x Numeric vector.
#' @param alpha Significance level (default: 0.05 for 95% CI).
#'
#' @return Lower quantile value.
#'
#' @export
qlow <- function(x, alpha = 0.05) {
  stats::quantile(x, probs = alpha / 2, na.rm = TRUE)
}


#' Upper Quantile
#'
#' @param x Numeric vector.
#' @param alpha Significance level (default: 0.05 for 95% CI).
#'
#' @return Upper quantile value.
#'
#' @export
qhigh <- function(x, alpha = 0.05) {
  stats::quantile(x, probs = 1 - alpha / 2, na.rm = TRUE)
}


# =============================================================================
# Quiet Execution
# =============================================================================
#' Execute Expression Quietly
#'
#' Suppresses messages and warnings from an expression.
#'
#' @param expr Expression to evaluate.
#'
#' @return Result of expression.
#'
#' @keywords internal
quiet <- function(expr) {
  suppressWarnings(suppressMessages(expr))
}
