# Improved Factor Variable Processing with Largest Value as Reference
# CRAN-compliant version with proper documentation and exports

#' Process Factor Variables with Largest Value as Reference
#'
#' Creates dummy variables from factor/categorical variables using the largest
#' value as the reference level. This is the standard approach in clinical
#' research where higher values often represent more severe conditions.
#'
#' @param data The input dataset containing the factor variables
#' @param factor_vars Character vector of factor/categorical variable names
#' @param df_work Working data.frame to add processed variables to
#' @param verbose Logical. Print processing details. Default is FALSE.
#'
#' @return The input \code{df_work} data.frame with added dummy variables.
#'   Dummy variable names follow the pattern \code{z_varname} for binary
#'   variables or \code{z_varname_level} for multi-level variables.
#'
#' @details
#' For binary variables (2 levels), the function creates a single indicator
#' variable where 1 indicates the smaller/first value and 0 indicates the
#' larger/last value (reference).
#'
#' For multi-level variables (>2 levels), the function creates k-1 dummy
#' variables where k is the number of levels. The largest value (for numeric)
#' or last alphabetically (for character/factor) is the reference level.
#'
#' Variables with only one level are skipped with no dummy created.
#'
#' @examples
#' # Create example data
#' example_data <- data.frame(
#'   binary_var = sample(0:1, 100, replace = TRUE),
#'   grade = sample(1:3, 100, replace = TRUE),
#'   stage = sample(c("I", "II", "III", "IV"), 100, replace = TRUE)
#' )
#'
#' # Initialize working data.frame
#' df_work <- data.frame(id = 1:100)
#'
#' # Process factors
#' df_work <- process_factor_variables(
#'   data = example_data,
#'   factor_vars = c("binary_var", "grade", "stage"),
#'   df_work = df_work,
#'   verbose = TRUE
#' )
#'
#' # Check results
#' names(df_work)
#' # For grade (1,2,3): z_grade_1, z_grade_2 created, grade=3 is reference
#' # For stage (I,II,III,IV): z_stage_I, z_stage_II, z_stage_III created, IV is reference
#'
#' @seealso \code{\link{process_covariates_corrected}} for processing both
#'   continuous and factor variables together
#'
#' @export
process_factor_variables <- function(data, factor_vars, df_work, verbose = FALSE) {

  if (length(factor_vars) == 0) {
    return(df_work)
  }

  for (var in factor_vars) {
    if (verbose) message("Processing factor variable: ", var)

    # Get unique values and determine the reference level (largest)
    if (is.factor(data[[var]])) {
      all_levels <- levels(data[[var]])
    } else {
      all_levels <- sort(unique(data[[var]][!is.na(data[[var]])]))
    }

    n_levels <- length(all_levels)

    if (verbose) {
      message("  Levels found: ", paste(all_levels, collapse = ", "))
    }

    if (n_levels == 0) {
      if (verbose) message("  Skipping - no valid levels")
      next

    } else if (n_levels == 1) {
      if (verbose) message("  Skipping - only one level")
      next

    } else if (n_levels == 2) {
      # Binary variable - create single indicator
      # Use largest value as reference (indicator = 1 for smaller value)
      if (is.numeric(all_levels)) {
        ref_level <- max(all_levels)
        other_level <- min(all_levels)
      } else {
        # For character, last alphabetically is reference
        sorted_levels <- sort(as.character(all_levels))
        ref_level <- sorted_levels[2]
        other_level <- sorted_levels[1]
      }

      if (verbose) {
        message("  Binary variable - Reference level (omitted): ", ref_level)
        message("  Creating indicator for: ", other_level)
      }

      # Create indicator: 1 if equal to smaller value, 0 if equal to larger value
      df_work[[paste0("z_", var)]] <- as.integer(data[[var]] == other_level)

    } else {
      # Multiple levels - create dummy variables with largest as reference

      if (is.numeric(all_levels)) {
        ref_level <- max(all_levels)
        other_levels <- setdiff(all_levels, ref_level)
        other_levels <- sort(other_levels)
      } else {
        sorted_levels <- sort(as.character(all_levels))
        ref_level <- sorted_levels[length(sorted_levels)]
        other_levels <- sorted_levels[-length(sorted_levels)]
      }

      if (verbose) {
        message("  Multi-level variable - Reference level (omitted): ", ref_level)
        message("  Creating indicators for: ", paste(other_levels, collapse = ", "))
      }

      # Create dummy variables for all levels except reference
      for (level in other_levels) {
        dummy_name <- paste0("z_", var, "_", level)
        df_work[[dummy_name]] <- as.integer(data[[var]] == level)

        if (verbose) {
          message("    Created: ", dummy_name, " for level ", level)
        }
      }
    }
  }

  return(df_work)
}


#' Process Continuous and Factor Variables for Modeling
#'
#' Processes covariates by standardizing continuous variables and creating
#' dummy variables for factors. Uses largest value as reference for factors.
#'
#' @param data Original dataset
#' @param continuous_vars Character vector of continuous variable names
#' @param factor_vars Character vector of factor variable names
#' @param df_work Working data.frame to add processed variables to
#' @param standardize Logical. If TRUE, standardize continuous variables
#'   to mean 0 and SD 1. Default is TRUE.
#' @param verbose Logical. Print processing details. Default is FALSE.
#'
#' @return The input \code{df_work} data.frame with added processed variables.
#'   Continuous variables are named \code{z_varname}.
#'   Factor dummy variables follow naming from \code{\link{process_factor_variables}}.
#'
#' @examples
#' # Create example data
#' example_data <- data.frame(
#'   age = rnorm(100, 50, 10),
#'   biomarker = rgamma(100, 2, 1),
#'   grade = sample(1:3, 100, replace = TRUE),
#'   sex = sample(c("M", "F"), 100, replace = TRUE)
#' )
#'
#' # Initialize working data.frame
#' df_work <- data.frame(id = 1:100)
#'
#' # Process all covariates
#' df_work <- process_covariates_corrected(
#'   data = example_data,
#'   continuous_vars = c("age", "biomarker"),
#'   factor_vars = c("grade", "sex"),
#'   df_work = df_work,
#'   verbose = TRUE
#' )
#'
#' head(df_work)
#'
#' @seealso \code{\link{process_factor_variables}}
#'
#' @importFrom stats scale
#' @export
process_covariates_corrected <- function(data,
                                          continuous_vars,
                                          factor_vars,
                                          df_work,
                                          standardize = TRUE,
                                          verbose = FALSE) {

  # Process continuous variables
  for (var in continuous_vars) {
    if (verbose) message("Processing continuous variable: ", var)

    if (standardize) {
      df_work[[paste0("z_", var)]] <- as.numeric(scale(data[[var]]))
    } else {
      df_work[[paste0("z_", var)]] <- data[[var]]
    }
  }

  # Process factor variables
  df_work <- process_factor_variables(
    data = data,
    factor_vars = factor_vars,
    df_work = df_work,
    verbose = verbose
  )

  return(df_work)
}


#' Demonstrate Factor Processing Behavior
#'
#' Creates an example dataset and demonstrates how factor variables are
#' processed with the largest value as reference.
#'
#' @return A data.frame showing the processed dummy variables
#'
#' @examples
#' result <- demonstrate_factor_processing()
#'
#' @export
demonstrate_factor_processing <- function() {

  message("==================================================")
  message("DEMONSTRATION: Factor Processing with Largest Value as Reference")
  message("==================================================\n")

  # Create example dataset
  set.seed(123)
  example_data <- data.frame(
    binary_var = sample(0:1, 100, replace = TRUE),
    sex = sample(c("F", "M"), 100, replace = TRUE),
    grade = sample(1:3, 100, replace = TRUE),
    stage = sample(c("I", "II", "III", "IV"), 100, replace = TRUE),
    severity = factor(
      sample(c("mild", "moderate", "severe"), 100, replace = TRUE),
      levels = c("mild", "moderate", "severe"),
      ordered = TRUE
    )
  )

  # Initialize working dataframe
  df_work <- data.frame(id = seq_len(100))

  # Define factor variables
  factor_vars <- c("binary_var", "sex", "grade", "stage", "severity")

  # Process with verbose output
  df_work <- process_factor_variables(
    data = example_data,
    factor_vars = factor_vars,
    df_work = df_work,
    verbose = TRUE
  )

  message("\n==================================================")
  message("RESULTS")
  message("==================================================\n")

  # Show the created variables
  message("Created dummy variables:")
  dummy_cols <- grep("^z_", names(df_work), value = TRUE)
  for (col in dummy_cols) {
    message(sprintf("  %-20s Mean = %.3f (proportion where indicator = 1)",
                    col, mean(df_work[[col]])))
  }

  # Summary explanation
  message("\n==================================================")
  message("SUMMARY OF PROCESSING RULES")
  message("==================================================\n")

  message("1. BINARY VARIABLES (2 levels):")
  message("   - LARGEST value is reference (omitted)")
  message("   - Smaller value gets indicator = 1")
  message("   - Example: For binary_var (0,1), z_binary_var = 1 when var = 0\n")

  message("2. MULTI-LEVEL VARIABLES (>2 levels):")
  message("   - LARGEST value is reference (omitted)")
  message("   - Other values each get a dummy variable")
  message("   - Example: For grade (1,2,3):")
  message("     * z_grade_1 = 1 when grade = 1")
  message("     * z_grade_2 = 1 when grade = 2")
  message("     * grade = 3 is reference (no dummy)\n")

  message("3. CHARACTER VARIABLES:")
  message("   - Last alphabetically is reference")
  message("   - Example: For stage (I,II,III,IV), IV is reference\n")

  return(df_work)
}


#' Get Reference Levels for Factor Variables
#'
#' Returns the reference level (omitted category) for each factor variable
#' based on the largest-value-as-reference convention.
#'
#' @param data Dataset containing the factor variables
#' @param factor_vars Character vector of factor variable names
#'
#' @return A named character vector where names are variable names and
#'   values are the reference levels
#'
#' @examples
#' df <- data.frame(
#'   grade = sample(1:3, 50, replace = TRUE),
#'   stage = sample(c("I", "II", "III"), 50, replace = TRUE)
#' )
#' get_reference_levels(df, c("grade", "stage"))
#' # Returns: grade = "3", stage = "III"
#'
#' @export
get_reference_levels <- function(data, factor_vars) {

  ref_levels <- character(length(factor_vars))
  names(ref_levels) <- factor_vars

  for (var in factor_vars) {
    if (is.factor(data[[var]])) {
      all_levels <- levels(data[[var]])
    } else {
      all_levels <- sort(unique(data[[var]][!is.na(data[[var]])]))
    }

    if (length(all_levels) == 0) {
      ref_levels[var] <- NA_character_
    } else if (is.numeric(all_levels)) {
      ref_levels[var] <- as.character(max(all_levels))
    } else {
      sorted_levels <- sort(as.character(all_levels))
      ref_levels[var] <- sorted_levels[length(sorted_levels)]
    }
  }

  return(ref_levels)
}
