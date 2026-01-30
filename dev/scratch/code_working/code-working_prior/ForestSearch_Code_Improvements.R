# ForestSearch: Ready-to-Use Code Improvements
# =============================================
# Copy these implementations directly into your package

# =============================================================================
# 1. UTILITY FUNCTIONS (Add to a new file: R/utils.R)
# =============================================================================

#' Ensure Object is a data.table
#'
#' Safely converts data.frames, matrices, or data.tables to data.table,
#' always returning a copy to avoid modifying the original.
#'
#' @param x Object to convert
#' @return A data.table copy
#' @keywords internal
#' @importFrom data.table as.data.table copy
ensure_data_table <- function(x) {

  if (inherits(x, "data.table")) {
    return(data.table::copy(x))
  }
  if (is.data.frame(x)) {
    return(data.table::as.data.table(x))
  }
  if (is.matrix(x)) {
    return(data.table::as.data.table(as.data.frame(x)))
  }
  stop(sprintf("Cannot convert object of class '%s' to data.table", class(x)[1]))
}


#' Validate Common Function Inputs
#'
#' Performs standard validation checks for data and variable specifications.
#'
#' @param data A data.frame
#' @param continuous_vars Character vector of continuous variable names
#' @param cat_vars Character vector of categorical variable names
#' @param outcome_var Optional outcome variable name
#' @param event_var Optional event indicator variable name
#' @param treatment_var Optional treatment variable name
#' @return Invisible TRUE if all checks pass
#' @keywords internal
validate_dgm_inputs <- function(data, continuous_vars = NULL, cat_vars = NULL,
                                outcome_var = NULL, event_var = NULL,
                                treatment_var = NULL) {
  # Check data type
  if (!is.data.frame(data)) {
    stop(sprintf("'data' must be a data.frame, received: %s", class(data)[1]))
  }
  
  # Check data has rows
  if (nrow(data) == 0) {
    stop("'data' contains no rows")
  }
  
  # Collect all specified variables
  all_vars <- c(continuous_vars, cat_vars, outcome_var, event_var, treatment_var)
  all_vars <- unique(all_vars[!is.null(all_vars) & nzchar(all_vars)])
  
  # Check variable existence
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(sprintf(
      "Variable(s) not found in data: %s\nAvailable columns: %s",
      paste(missing_vars, collapse = ", "),
      paste(head(names(data), 15), collapse = ", ")
    ))
  }
  
  # Check continuous variables are numeric
  if (!is.null(continuous_vars)) {
    non_numeric <- continuous_vars[!sapply(data[continuous_vars], is.numeric)]
    if (length(non_numeric) > 0) {
      warning(sprintf(
        "Non-numeric variables specified as continuous: %s",
        paste(non_numeric, collapse = ", ")
      ))
    }
  }
  
  invisible(TRUE)
}


#' Report Progress Message
#'
#' Conditionally prints progress messages based on verbosity setting.
#'
#' @param msg Message to display
#' @param verbose Logical indicating whether to display
#' @param level Indentation level (1 = no indent, 2 = one indent, etc.)
#' @keywords internal
report_progress <- function(msg, verbose = TRUE, level = 1L) {
  if (isTRUE(verbose)) {
    indent <- strrep("  ", max(0L, level - 1L))
    message(paste0(indent, msg))
  }
  invisible(NULL)
}


# =============================================================================
# 2. IMPROVED FACTOR PROCESSING (Replace in generate_aft_dgm*.R)
# =============================================================================

#' Process Factor Variables with Consistent Reference Level
#'
#' Creates dummy variables for factor/categorical variables using the largest
#' value as the reference level, which is standard in clinical research.
#'
#' @param data Original dataset
#' @param factor_vars Character vector of factor variable names
#' @param df_work Working data.frame to add processed variables to
#' @param prefix Prefix for dummy variable names (default: "z_")
#' @param verbose Print processing details
#' @return df_work with added dummy variables
#' @keywords internal
#' @importFrom stats model.matrix
process_factor_variables_v2 <- function(data, factor_vars, df_work,
                                        prefix = "z_", verbose = FALSE) {
  
  if (length(factor_vars) == 0) return(df_work)
  
  for (var_name in factor_vars) {
    if (verbose) report_progress(sprintf("Processing: %s", var_name), TRUE)
    
    var_data <- data[[var_name]]
    
    # Determine all unique levels
    if (is.factor(var_data)) {
      all_levels <- levels(var_data)
    } else {
      all_levels <- sort(unique(var_data[!is.na(var_data)]))
    }
    
    n_levels <- length(all_levels)
    
    if (n_levels <= 1) {
      if (verbose) report_progress("  Skipping - single level", TRUE, 2)
      next
    }
    
    # Determine reference level (largest value)
    if (is.numeric(all_levels)) {
      ref_level <- max(all_levels)
      other_levels <- sort(setdiff(all_levels, ref_level))
    } else {
      # For character/factor: last alphabetically
      sorted_levels <- sort(as.character(all_levels))
      ref_level <- sorted_levels[length(sorted_levels)]
      other_levels <- sorted_levels[-length(sorted_levels)]
    }
    
    if (verbose) {
      report_progress(sprintf("  Reference level: %s", ref_level), TRUE, 2)
    }
    
    # Create dummy variables
    if (n_levels == 2) {
      # Binary: single indicator
      dummy_name <- paste0(prefix, var_name)
      df_work[[dummy_name]] <- as.numeric(var_data == other_levels[1])
    } else {
      # Multi-level: one dummy per non-reference level
      for (lvl in other_levels) {
        dummy_name <- paste0(prefix, var_name, "_", lvl)
        df_work[[dummy_name]] <- as.numeric(var_data == lvl)
      }
    }
  }
  
  return(df_work)
}


# =============================================================================
# 3. VECTORIZED CATEGORICAL PERTURBATION (Replace in generate_bootstrap_synthetic)
# =============================================================================

#' Perturb Categorical Variable (Vectorized)
#'
#' Efficiently perturbs categorical values by flipping a random subset
#' to different values from the same category set.
#'
#' @param values Vector of categorical values
#' @param flip_prob Probability of flipping each value
#' @param is_ordinal If TRUE, only flip to adjacent values
#' @return Perturbed vector
#' @keywords internal
perturb_categorical_vectorized <- function(values, flip_prob, is_ordinal = FALSE) {
  n <- length(values)
  flip_mask <- runif(n) < flip_prob
  n_flip <- sum(flip_mask)
  
  if (n_flip == 0) return(values)
  
  unique_vals <- sort(unique(values))
  n_unique <- length(unique_vals)
  
  if (n_unique < 2) return(values)
  
  result <- values
  flip_idx <- which(flip_mask)
  
  if (is_ordinal && is.numeric(unique_vals)) {
    # Ordinal: move to adjacent value
    for (i in flip_idx) {
      current_pos <- match(values[i], unique_vals)
      candidates <- c()
      if (current_pos > 1) candidates <- c(candidates, unique_vals[current_pos - 1])
      if (current_pos < n_unique) candidates <- c(candidates, unique_vals[current_pos + 1])
      if (length(candidates) > 0) {
        result[i] <- sample(candidates, 1)
      }
    }
  } else if (n_unique == 2) {
    # Binary: simple flip
    other_val <- setdiff(unique_vals, values[flip_idx])
    if (is.numeric(values)) {
      result[flip_idx] <- 1 - values[flip_idx]
    } else {
      result[flip_idx] <- ifelse(
        values[flip_idx] == unique_vals[1],
        unique_vals[2],
        unique_vals[1]
      )
    }
  } else {
    # Multi-category: vectorized replacement
    # Pre-compute replacement options for each unique value
    replacement_map <- lapply(unique_vals, function(v) setdiff(unique_vals, v))
    names(replacement_map) <- as.character(unique_vals)
    
    # Apply replacements
    current_vals <- as.character(values[flip_idx])
    new_vals <- vapply(current_vals, function(cv) {
      candidates <- replacement_map[[cv]]
      if (length(candidates) > 0) sample(candidates, 1) else cv
    }, FUN.VALUE = character(1))
    
    # Restore original type
    if (is.numeric(values)) {
      result[flip_idx] <- as.numeric(new_vals)
    } else if (is.factor(values)) {
      result[flip_idx] <- factor(new_vals, levels = levels(values))
    } else {
      result[flip_idx] <- new_vals
    }
  }
  
  return(result)
}


# =============================================================================
# 4. IMPROVED HAZARD RATIO CALCULATION (Extract common pattern)
# =============================================================================

#' Calculate Empirical Hazard Ratios from Potential Outcomes
#'
#' Computes hazard ratios by fitting Cox models to potential outcome data,
#' optionally stratified by subgroup membership.
#'
#' @param T_treat Potential survival times under treatment
#' @param T_control Potential survival times under control
#' @param subgroup_indicator Optional logical/numeric vector indicating subgroup
#' @param verbose Print results
#' @return Named list with overall HR and optionally subgroup HRs
#' @keywords internal
#' @importFrom survival coxph Surv
calculate_empirical_hr <- function(T_treat, T_control,
                                   subgroup_indicator = NULL,
                                   verbose = FALSE) {
  n <- length(T_treat)
  stopifnot(length(T_control) == n)
  
  # Create stacked dataset for Cox model
  df_pot <- data.frame(
    time = c(T_treat, T_control),
    event = 1L,  # All events for potential outcomes
    treat = rep(c(1L, 0L), each = n)
  )
  
  # Overall HR
  fit_overall <- survival::coxph(
    survival::Surv(time, event) ~ treat,
    data = df_pot
  )
  hr_overall <- exp(coef(fit_overall)["treat"])
  
  result <- list(overall = unname(hr_overall))
  
  # Subgroup-specific HRs
  if (!is.null(subgroup_indicator)) {
    stopifnot(length(subgroup_indicator) == n)
    df_pot$subgroup <- rep(as.numeric(subgroup_indicator), 2)
    
    # Subgroup HR (subgroup_indicator == 1)
    df_sub <- df_pot[df_pot$subgroup == 1, ]
    if (nrow(df_sub) >= 10) {  # Minimum observations
      fit_sub <- survival::coxph(survival::Surv(time, event) ~ treat, data = df_sub)
      result$subgroup <- unname(exp(coef(fit_sub)["treat"]))
    }
    
    # Complement HR (subgroup_indicator == 0)
    df_comp <- df_pot[df_pot$subgroup == 0, ]
    if (nrow(df_comp) >= 10) {
      fit_comp <- survival::coxph(survival::Surv(time, event) ~ treat, data = df_comp)
      result$complement <- unname(exp(coef(fit_comp)["treat"]))
    }
  }
  
  if (verbose) {
    cat("\n=== Empirical Hazard Ratios ===\n")
    for (nm in names(result)) {
      cat(sprintf("  %s: %.3f\n", nm, result[[nm]]))
    }
  }
  
  return(result)
}


# =============================================================================
# 5. FLEXIBLE CUTPOINT PROCESSOR (Improved version)
# =============================================================================

#' Process Flexible Cutpoint Specification
#'
#' Evaluates a cutpoint specification against data to produce a logical
#' indicator vector. Supports numeric cutpoints, quantiles, functions,
#' ranges, and custom functions.
#'
#' @param var_data Numeric vector of variable values
#' @param cut_spec Cutpoint specification (see Details)
#' @param var_name Variable name (for error messages)
#' @return Logical vector indicating condition is met
#' 
#' @details
#' The \code{cut_spec} parameter accepts:
#' \itemize{
#'   \item Numeric scalar: \code{cut_spec = 20} means \code{var_data <= 20}
#'   \item Quantile: \code{list(type = "quantile", value = 0.25)}
#'   \item Function: \code{list(type = "function", fun = median)}
#'   \item Range: \code{list(type = "range", min = 10, max = 50)}
#'   \item Greater than: \code{list(type = "greater", value = 30)}
#'   \item Custom: \code{list(type = "custom", fun = function(x) x < mean(x))}
#' }
#'
#' @keywords internal
evaluate_cutpoint <- function(var_data, cut_spec, var_name = "variable") {
  
  # Handle NULL - use median as default

  if (is.null(cut_spec)) {
    cutpoint <- median(var_data, na.rm = TRUE)
    return(var_data <= cutpoint)
  }
  
  # Handle simple numeric cutpoint
  if (is.numeric(cut_spec) && length(cut_spec) == 1) {
    return(var_data <= cut_spec)
  }
  
  # Handle list specifications
  if (!is.list(cut_spec)) {
    stop(sprintf(
      "Invalid cutpoint specification for '%s': expected numeric or list, got %s",
      var_name, class(cut_spec)[1]
    ))
  }
  
  cut_type <- cut_spec$type
  
  switch(cut_type,
    "quantile" = {
      if (is.null(cut_spec$value)) {
        stop(sprintf("Quantile cutpoint for '%s' requires 'value' parameter", var_name))
      }
      cutpoint <- quantile(var_data, probs = cut_spec$value, na.rm = TRUE)
      var_data <= cutpoint
    },
    
    "function" = {
      if (is.null(cut_spec$fun)) {
        stop(sprintf("Function cutpoint for '%s' requires 'fun' parameter", var_name))
      }
      cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
      var_data <= cutpoint
    },
    
    "range" = {
      if (is.null(cut_spec$min) || is.null(cut_spec$max)) {
        stop(sprintf("Range cutpoint for '%s' requires 'min' and 'max'", var_name))
      }
      var_data >= cut_spec$min & var_data <= cut_spec$max
    },
    
    "greater" = {
      cutpoint <- cut_spec$value %||% 
                  (if (!is.null(cut_spec$quantile)) 
                     quantile(var_data, cut_spec$quantile, na.rm = TRUE)
                   else if (!is.null(cut_spec$fun))
                     cut_spec$fun(var_data, na.rm = TRUE)
                   else stop("Greater cutpoint needs value, quantile, or fun"))
      var_data > cutpoint
    },
    
    "multiple" = {
      if (is.null(cut_spec$values)) {
        stop(sprintf("Multiple cutpoint for '%s' requires 'values'", var_name))
      }
      var_data %in% cut_spec$values
    },
    
    "custom" = {
      if (is.null(cut_spec$fun) || !is.function(cut_spec$fun)) {
        stop(sprintf("Custom cutpoint for '%s' requires 'fun' function", var_name))
      }
      result <- cut_spec$fun(var_data)
      if (!is.logical(result)) {
        stop(sprintf("Custom function for '%s' must return logical vector", var_name))
      }
      result
    },
    
    # Default: unknown type
    stop(sprintf("Unknown cutpoint type '%s' for variable '%s'", cut_type, var_name))
  )
}

# Null coalescing operator if not available
`%||%` <- function(a, b) if (is.null(a)) b else a


# =============================================================================
# 6. IMPROVED BOOTSTRAP SUMMARY (Fix for data.table issues)
# =============================================================================

#' Create Frequency Table Using data.table
#'
#' Efficiently creates a frequency table with counts and percentages.
#'
#' @param values Vector of values to tabulate
#' @param total_n Denominator for percentage calculation
#' @param value_col Name for the value column
#' @return data.table with columns: value_col, N, Percent
#' @keywords internal
#' @importFrom data.table data.table .N
create_freq_table_dt <- function(values, total_n = NULL, value_col = "Value") {
  # Remove NA and empty strings
  values <- values[!is.na(values) & values != ""]
  
  if (length(values) == 0) {
    return(data.table::data.table(
      Value = character(),
      N = integer(),
      Percent = numeric()
    ))
  }
  
  if (is.null(total_n)) total_n <- length(values)
  
  # Use data.table aggregation
  dt <- data.table::data.table(val = values)
  freq <- dt[, .(N = .N), by = val]
  freq[, Percent := 100 * N / total_n]
  freq <- freq[order(-N)]
  
  data.table::setnames(freq, "val", value_col)
  return(freq)
}


# =============================================================================
# 7. PACKAGE-LEVEL DOCUMENTATION TEMPLATE
# =============================================================================

# Add this to R/ForestSearch-package.R:

#' ForestSearch: Exploratory Subgroup Identification for Survival Endpoints
#'
#' @description
#' ForestSearch implements methods for exploratory subgroup identification in
#' clinical trials with survival endpoints. The package combines machine learning
#' approaches (Generalized Random Forests, LASSO regularization) with exhaustive
#' combinatorial search algorithms and bootstrap bias correction using
#' infinitesimal jackknife methods.
#'
#' @section Main Analysis Functions:
#' \itemize{
#'   \item \code{\link{forest_search}}: Main subgroup identification algorithm
#'   \item \code{\link{run_bootstrap_analysis}}: Bootstrap validation
#'   \item \code{\link{summarize_bootstrap_subgroups}}: Summarize bootstrap results
#' }
#'
#' @section Data Generating Mechanisms:
#' \itemize{
#'   \item \code{\link{generate_bootstrap_synthetic}}: Non-parametric bootstrap
#'   \item \code{\link{generate_aft_dgm}}: Parametric AFT model
#'   \item \code{\link{generate_aft_dgm_flex}}: AFT with flexible subgroups
#' }
#'
#' @section References:
#' Leon LF, Marceau-West CT, et al. (2024). Exploratory subgroup identification
#' in the heterogeneous Cox model: A relatively simple procedure.
#' \emph{Statistics in Medicine}. doi:10.1002/sim.XXXX
#'
#' @docType package
#' @name ForestSearch-package
#' @aliases ForestSearch
#'
#' @importFrom stats median quantile sd coef
#' @importFrom survival Surv coxph survreg
#' @importFrom data.table data.table as.data.table rbindlist setnames .N
NULL


# =============================================================================
# 8. TEST EXAMPLES
# =============================================================================

if (FALSE) {
  # Test the utility functions
  
  # Test ensure_data_table
  df <- data.frame(x = 1:5, y = letters[1:5])
  dt <- ensure_data_table(df)
  stopifnot(inherits(dt, "data.table"))
  
  # Test validate_dgm_inputs
  test_data <- data.frame(
    time = rexp(100),
    status = rbinom(100, 1, 0.7),
    age = rnorm(100, 50, 10),
    sex = sample(c("M", "F"), 100, TRUE)
  )
  
  validate_dgm_inputs(
    test_data,
    continuous_vars = "age",
    cat_vars = "sex",
    outcome_var = "time",
    event_var = "status"
  )
  
  # Test evaluate_cutpoint
  x <- rnorm(100, 50, 10)
  
  # Simple numeric
  ind1 <- evaluate_cutpoint(x, 50)
  cat("Numeric cutpoint (<=50):", mean(ind1), "\n")
  

  # Quantile
  ind2 <- evaluate_cutpoint(x, list(type = "quantile", value = 0.25))
  cat("25th percentile:", mean(ind2), "\n")
  
  # Range
  ind3 <- evaluate_cutpoint(x, list(type = "range", min = 40, max = 60))
  cat("Range 40-60:", mean(ind3), "\n")
  
  # Test perturb_categorical_vectorized
  cat_var <- sample(c("A", "B", "C"), 1000, replace = TRUE)
  perturbed <- perturb_categorical_vectorized(cat_var, flip_prob = 0.1)
  cat("Perturbation rate:", mean(cat_var != perturbed), "\n")
  
  cat("\nAll tests passed!\n")
}
