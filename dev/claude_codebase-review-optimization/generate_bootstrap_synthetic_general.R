# Generalized Bootstrap with Perturbation Function for Synthetic Data Generation
# This function can work with any dataset, not just GBSG
#
# CRAN-compliant version with proper documentation and exports

#' Generate Synthetic Data using Bootstrap with Perturbation
#'
#' Creates synthetic datasets by bootstrap resampling from original data
#' with controlled perturbation of continuous and categorical variables.
#' Useful for simulation studies, data augmentation, and privacy-preserving
#' data sharing.
#'
#' @param data Original dataset to bootstrap from (data.frame)
#' @param continuous_vars Character vector of continuous variable names
#' @param cat_vars Character vector of categorical variable names
#' @param n Integer. Number of synthetic observations to generate.
#'   Default is \code{nrow(data)}.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#' @param noise_level Numeric between 0 and 1. Controls the standard deviation
#'   of noise added to continuous variables as a fraction of the original SD.
#'   Default is 0.1 (10\% of original SD).
#' @param id_var Optional character. Name of ID variable to regenerate
#'   (will be numbered 1:n).
#' @param cat_flip_prob Numeric between 0 and 1. Probability of flipping
#'   categorical values. Default is \code{noise_level / 2}.
#' @param preserve_bounds Logical. If TRUE, continuous variables are constrained
#'   to stay within the original data range. Default is TRUE.
#' @param ordinal_vars Optional character vector of ordinal categorical variables.
#'   These will be perturbed to adjacent values rather than randomly flipped.
#'
#' @return A data.frame with synthetic data having the same structure as input.
#'
#' @details
#' The function works in three stages:
#' \enumerate{
#'   \item Bootstrap resampling with replacement from original data
#'   \item Perturbation of continuous variables with Gaussian noise
#'   \item Perturbation of categorical variables with random flipping
#' }
#'
#' For continuous variables, noise is drawn from \code{N(0, noise_level * sd(x))}.
#' Integer-valued continuous variables are rounded after perturbation.
#'
#' For categorical variables, each value has probability \code{cat_flip_prob}
#' of being changed. Binary variables are flipped to the other value.
#' Multi-category variables are changed to a random different category.
#' Ordinal variables (if specified) are changed only to adjacent categories.
#'
#' @examples
#' # Example with custom dataset
#' set.seed(42)
#' my_data <- data.frame(
#'   id = 1:100,
#'   age = round(rnorm(100, 50, 15)),
#'   weight = rnorm(100, 70, 12),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   stage = sample(1:4, 100, replace = TRUE)
#' )
#'
#' synth <- generate_bootstrap_synthetic(
#'   data = my_data,
#'   continuous_vars = c("age", "weight"),
#'   cat_vars = c("gender"),
#'   ordinal_vars = c("stage"),
#'   id_var = "id",
#'   n = 150,
#'   seed = 123,
#'   noise_level = 0.15
#' )
#'
#' # Compare distributions
#' summary(my_data$age)
#' summary(synth$age)
#'
#' @seealso \code{\link{detect_variable_types}} for automatic variable classification
#'
#' @importFrom stats rnorm runif sd
#' @export
generate_bootstrap_synthetic <- function(data,
                                         continuous_vars,
                                         cat_vars,
                                         n = NULL,
                                         seed = 123,
                                         noise_level = 0.1,
                                         id_var = NULL,
                                         cat_flip_prob = NULL,
                                         preserve_bounds = TRUE,
                                         ordinal_vars = NULL) {


  # ==========================================================================

# Input validation
  # ==========================================================================

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  if (!all(continuous_vars %in% names(data))) {
    missing <- continuous_vars[!continuous_vars %in% names(data)]
    stop("Continuous variables not found in data: ", paste(missing, collapse = ", "))
  }

  if (!all(cat_vars %in% names(data))) {
    missing <- cat_vars[!cat_vars %in% names(data)]
    stop("Categorical variables not found in data: ", paste(missing, collapse = ", "))
  }

  if (!is.null(ordinal_vars) && !all(ordinal_vars %in% names(data))) {
    missing <- ordinal_vars[!ordinal_vars %in% names(data)]
    stop("Ordinal variables not found in data: ", paste(missing, collapse = ", "))
  }

  if (noise_level < 0 || noise_level > 1) {
    stop("'noise_level' must be between 0 and 1")
  }

  # ==========================================================================
  # Set defaults
  # ==========================================================================

  if (is.null(n)) {
    n <- nrow(data)
  }

  if (is.null(cat_flip_prob)) {
    cat_flip_prob <- noise_level / 2
  }

  # Set seed for reproducibility
  set.seed(seed)

  # ==========================================================================
  # Bootstrap sample from original data
  # ==========================================================================

  boot_indices <- sample(nrow(data), n, replace = TRUE)
  synthetic_data <- data[boot_indices, , drop = FALSE]

  # Reset row names
  rownames(synthetic_data) <- NULL

  # ==========================================================================
  # Pre-compute which continuous variables are integers (efficiency fix)
  # ==========================================================================

  is_integer_var <- vapply(continuous_vars, function(var) {
    if (!is.numeric(data[[var]])) return(FALSE)
    all(data[[var]] == round(data[[var]]), na.rm = TRUE)
  }, FUN.VALUE = logical(1))

  # ==========================================================================
  # Perturb continuous variables
  # ==========================================================================

  for (var in continuous_vars) {
    # Check if variable is numeric
    if (!is.numeric(data[[var]])) {
      warning("Variable '", var, "' is not numeric. Skipping continuous perturbation.")
      next
    }

    # Calculate noise based on variable's standard deviation
    var_sd <- sd(data[[var]], na.rm = TRUE)
    noise_sd <- var_sd * noise_level

    # Generate and add noise
    noise <- rnorm(n, mean = 0, sd = noise_sd)
    synthetic_data[[var]] <- synthetic_data[[var]] + noise

    # Preserve bounds if requested
    if (preserve_bounds) {
      var_min <- min(data[[var]], na.rm = TRUE)
      var_max <- max(data[[var]], na.rm = TRUE)
      synthetic_data[[var]] <- pmax(var_min, pmin(var_max, synthetic_data[[var]]))
    }

    # Round if original variable is integer-valued
    if (is_integer_var[var]) {
      synthetic_data[[var]] <- round(synthetic_data[[var]])
    }
  }

  # ==========================================================================
  # Perturb categorical variables (vectorized where possible)
  # ==========================================================================

  for (var in cat_vars) {
    # Determine which indices to flip
    flip_mask <- runif(n) < cat_flip_prob
    flip_indices <- which(flip_mask)

    if (length(flip_indices) == 0) next

    # Check if this is an ordinal variable
    is_ordinal <- !is.null(ordinal_vars) && var %in% ordinal_vars
    unique_vals <- sort(unique(data[[var]]))
    n_unique <- length(unique_vals)

    if (n_unique < 2) next

    if (is_ordinal) {
      # Ordinal perturbation: change to adjacent values
      for (i in flip_indices) {
        current_val <- synthetic_data[[var]][i]
        current_pos <- which(unique_vals == current_val)

        if (length(current_pos) == 0) next

        # Determine possible new values (adjacent only)
        candidates <- c()
        if (current_pos > 1) {
          candidates <- c(candidates, unique_vals[current_pos - 1])
        }
        if (current_pos < n_unique) {
          candidates <- c(candidates, unique_vals[current_pos + 1])
        }

        if (length(candidates) > 0) {
          synthetic_data[[var]][i] <- sample(candidates, 1)
        }
      }

    } else if (n_unique == 2) {
      # Binary variable: efficient flip
      if (is.numeric(synthetic_data[[var]])) {
        # For 0/1 binary
        synthetic_data[[var]][flip_indices] <- 1 - synthetic_data[[var]][flip_indices]
      } else if (is.logical(synthetic_data[[var]])) {
        # For TRUE/FALSE
        synthetic_data[[var]][flip_indices] <- !synthetic_data[[var]][flip_indices]
      } else {
        # For factor or character binary
        current_vals <- synthetic_data[[var]][flip_indices]
        new_vals <- ifelse(current_vals == unique_vals[1], unique_vals[2], unique_vals[1])
        synthetic_data[[var]][flip_indices] <- new_vals
      }

    } else {
      # Multi-category variable: vectorized replacement
      # Pre-compute replacement options
      replacement_map <- lapply(unique_vals, function(v) setdiff(unique_vals, v))
      names(replacement_map) <- as.character(unique_vals)

      current_vals <- as.character(synthetic_data[[var]][flip_indices])
      new_vals <- vapply(current_vals, function(cv) {
        candidates <- replacement_map[[cv]]
        if (length(candidates) > 0) sample(candidates, 1) else cv
      }, FUN.VALUE = character(1))

      # Restore original type
      if (is.factor(synthetic_data[[var]])) {
        synthetic_data[[var]][flip_indices] <- factor(new_vals,
                                                       levels = levels(synthetic_data[[var]]))
      } else if (is.numeric(data[[var]])) {
        synthetic_data[[var]][flip_indices] <- as.numeric(new_vals)
      } else {
        synthetic_data[[var]][flip_indices] <- new_vals
      }
    }
  }

  # ==========================================================================
  # Regenerate ID variable if specified
  # ==========================================================================

  if (!is.null(id_var) && id_var %in% names(synthetic_data)) {
    synthetic_data[[id_var]] <- seq_len(n)
  }

  # ==========================================================================
  # Report on unspecified variables
  # ==========================================================================

  all_specified_vars <- c(continuous_vars, cat_vars)
  if (!is.null(id_var)) {
    all_specified_vars <- c(all_specified_vars, id_var)
  }

  remaining_vars <- setdiff(names(data), all_specified_vars)

  if (length(remaining_vars) > 0) {
    message("Note: The following variables were not specified as continuous or categorical ",
            "and will be kept as-is from bootstrap sample:\n  ",
            paste(remaining_vars, collapse = ", "))
  }

  return(synthetic_data)
}


#' Automatically Detect Variable Types in a Dataset
#'
#' Classifies variables as continuous or categorical based on their type
#' and number of unique values. Useful for preparing datasets for synthetic
#' data generation or modeling.
#'
#' @param data A data.frame to analyze
#' @param max_unique_for_cat Integer. Maximum number of unique values for a
#'   numeric variable to be classified as categorical. Default is 10.
#' @param exclude_vars Character vector of variable names to exclude from
#'   classification (e.g., ID variables, outcomes). Default is NULL.
#'
#' @return A named list with two components:
#' \describe{
#'   \item{continuous_vars}{Character vector of continuous variable names}
#'   \item{cat_vars}{Character vector of categorical variable names}
#' }
#'
#' @examples
#' df <- data.frame(
#'   id = 1:100,
#'   age = rnorm(100, 50, 10),
#'   stage = sample(1:4, 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE)
#' )
#'
#' var_types <- detect_variable_types(df, exclude_vars = "id")
#' var_types$continuous_vars
#' var_types$cat_vars
#'
#' @seealso \code{\link{generate_bootstrap_synthetic}}
#' @export
detect_variable_types <- function(data, max_unique_for_cat = 10, exclude_vars = NULL) {
  continuous_vars <- c()
  cat_vars <- c()

  for (var in names(data)) {
    if (!is.null(exclude_vars) && var %in% exclude_vars) {
      next
    }

    if (is.numeric(data[[var]])) {
      n_unique <- length(unique(data[[var]]))
      # If numeric with few unique values, treat as categorical
      if (n_unique <= max_unique_for_cat) {
        cat_vars <- c(cat_vars, var)
      } else {
        continuous_vars <- c(continuous_vars, var)
      }
    } else {
      # Factor, character, or logical variables
      cat_vars <- c(cat_vars, var)
    }
  }

  return(list(
    continuous_vars = continuous_vars,
    cat_vars = cat_vars
  ))
}


#' Generate Synthetic GBSG Data using Generalized Bootstrap
#'
#' Convenience wrapper for generating synthetic data based on the German
#' Breast Cancer Study Group (GBSG) dataset structure.
#'
#' @param n Integer. Number of observations to generate. Default is 686.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#' @param noise_level Numeric between 0 and 1. Noise level for perturbation.
#'   Default is 0.1.
#'
#' @return A data.frame with synthetic GBSG-like data
#'
#' @examples
#' \donttest{
#' synth_gbsg <- generate_gbsg_bootstrap_general(n = 500, seed = 42)
#' head(synth_gbsg)
#' }
#'
#' @seealso \code{\link{generate_bootstrap_synthetic}}
#' @export
generate_gbsg_bootstrap_general <- function(n = 686, seed = 123, noise_level = 0.1) {
  # Check for survival package (CRAN-compliant approach)
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }

  # Access gbsg data via namespace
  gbsg <- survival::gbsg

  synthetic_data <- generate_bootstrap_synthetic(
    data = gbsg,
    continuous_vars = c("age", "size", "nodes", "pgr", "er", "rfstime"),
    cat_vars = c("meno", "hormon", "status"),
    ordinal_vars = c("grade"),
    id_var = "pid",
    n = n,
    seed = seed,
    noise_level = noise_level,
    preserve_bounds = TRUE
  )

  return(synthetic_data)
}


#' Compare Original and Synthetic Datasets
#'
#' Prints summary statistics comparing distributions of variables between
#' original and synthetic datasets.
#'
#' @param original Original dataset (data.frame)
#' @param synthetic Synthetic dataset (data.frame)
#' @param continuous_vars Character vector of continuous variable names
#' @param cat_vars Character vector of categorical variable names
#'
#' @return Invisibly returns NULL. Called for side effect of printing.
#'
#' @examples
#' \donttest{
#' original <- data.frame(
#'   age = rnorm(100, 50, 10),
#'   sex = sample(c("M", "F"), 100, replace = TRUE)
#' )
#'
#' synthetic <- generate_bootstrap_synthetic(
#'   original,
#'   continuous_vars = "age",
#'   cat_vars = "sex",
#'   n = 100
#' )
#'
#' compare_datasets_general(original, synthetic, "age", "sex")
#' }
#'
#' @export
compare_datasets_general <- function(original, synthetic, continuous_vars, cat_vars) {
  cat("\n=== Dataset Comparison ===\n")
  cat("Original size:", nrow(original), "x", ncol(original), "\n")
  cat("Synthetic size:", nrow(synthetic), "x", ncol(synthetic), "\n\n")

  # Compare continuous variables
  if (length(continuous_vars) > 0) {
    cat("Continuous Variables:\n")
    cat(sprintf("%-15s %10s %10s %10s %10s %10s\n",
                "Variable", "Orig_Mean", "Synth_Mean", "Orig_SD", "Synth_SD", "Mean_Diff%"))
    cat(strrep("-", 75), "\n")

    for (var in continuous_vars) {
      if (var %in% names(original) && var %in% names(synthetic)) {
        orig_mean <- mean(original[[var]], na.rm = TRUE)
        synth_mean <- mean(synthetic[[var]], na.rm = TRUE)
        orig_sd <- sd(original[[var]], na.rm = TRUE)
        synth_sd <- sd(synthetic[[var]], na.rm = TRUE)
        mean_diff_pct <- abs((synth_mean - orig_mean) / orig_mean * 100)

        cat(sprintf("%-15s %10.2f %10.2f %10.2f %10.2f %10.1f\n",
                    var, orig_mean, synth_mean, orig_sd, synth_sd, mean_diff_pct))
      }
    }
  }

  # Compare categorical variables
  if (length(cat_vars) > 0) {
    cat("\nCategorical Variables:\n")

    for (var in cat_vars) {
      if (var %in% names(original) && var %in% names(synthetic)) {
        cat("\n", var, ":\n", sep = "")
        orig_tab <- table(original[[var]])
        synth_tab <- table(synthetic[[var]])

        # Get all unique levels
        all_levels <- union(names(orig_tab), names(synth_tab))

        cat(sprintf("%-10s %15s %15s\n", "Level", "Original_Prop", "Synthetic_Prop"))
        cat(strrep("-", 45), "\n")

        for (level in all_levels) {
          orig_prop <- if (level %in% names(orig_tab)) prop.table(orig_tab)[level] else 0
          synth_prop <- if (level %in% names(synth_tab)) prop.table(synth_tab)[level] else 0
          cat(sprintf("%-10s %15.3f %15.3f\n", level, orig_prop, synth_prop))
        }
      }
    }
  }

  invisible(NULL)
}
