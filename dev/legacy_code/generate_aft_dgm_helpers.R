# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Validate Input Parameters
#' @keywords internal
validate_inputs <- function(data, model, cens_type, outcome_var, event_var,
                            treatment_var, continuous_vars, factor_vars) {

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  if (!model %in% c("alt", "null")) {
    stop("'model' must be either 'alt' or 'null'")
  }

  if (!cens_type %in% c("weibull", "uniform")) {
    stop("'cens_type' must be either 'weibull' or 'uniform'")
  }

  # Check that required variables exist
  required_vars <- c(outcome_var, event_var, treatment_var)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Required variables not found in data: ",
         paste(missing_vars, collapse = ", "))
  }

  # Check continuous and factor variables
  all_covars <- c(continuous_vars, factor_vars)
  missing_covars <- setdiff(all_covars, names(data))
  if (length(missing_covars) > 0) {
    stop("Covariate variables not found in data: ",
         paste(missing_covars, collapse = ", "))
  }
}


#' Prepare Working Dataset with Processed Covariates
#' @keywords internal
prepare_working_dataset <- function(data, outcome_var, event_var, treatment_var,
                                    continuous_vars, factor_vars, standardize,
                                    verbose) {

  # Create base working dataset
  df_work <- data.frame(
    id = 1:nrow(data),
    y = data[[outcome_var]],
    treat = data[[treatment_var]],
    event = ifelse(data[[event_var]] == 1, 1, 0)
  )

  # Process continuous variables
  df_work <- process_continuous_vars(df_work, data, continuous_vars,
                                     standardize)

  # Process factor variables
  df_work <- process_factor_vars(df_work, data, factor_vars)

  # Add unprocessed variables
  df_work <- add_unprocessed_vars(df_work, data, outcome_var, event_var,
                                  treatment_var, continuous_vars, factor_vars,
                                  verbose)

  return(df_work)
}


#' Process Continuous Variables
#' @keywords internal
process_continuous_vars <- function(df_work, data, continuous_vars,
                                    standardize) {

  for (var in continuous_vars) {
    if (standardize) {
      df_work[[paste0("z_", var)]] <- scale(data[[var]])[, 1]  # Standardize
    } else {
      df_work[[paste0("z_", var)]] <- data[[var]]
    }
  }

  return(df_work)
}


#' Process Factor Variables with Largest Value as Reference
#' @keywords internal
process_factor_vars <- function(df_work, data, factor_vars) {

  for (var in factor_vars) {
    # Get unique values
    if (is.factor(data[[var]])) {
      all_levels <- levels(data[[var]])
    } else {
      all_levels <- sort(unique(data[[var]]))
    }

    n_levels <- length(all_levels)

    if (n_levels == 1) {
      # Skip variables with only one level
      next

    } else if (n_levels == 2) {
      # Binary variable - keep as-is (min as reference)
      if (is.numeric(all_levels)) {
        ref_level <- min(all_levels)
        other_level <- max(all_levels)
      } else {
        # For character, first alphabetically
        ref_level <- sort(all_levels, decreasing = FALSE)[1]
        other_level <- setdiff(all_levels, ref_level)
      }

      # Create single indicator for non-reference level
      df_work[[paste0("z_", var)]] <- as.numeric(data[[var]] == other_level)

    } else {
      # Multiple levels - create dummies with largest as reference
      if (is.numeric(all_levels)) {
        ref_level <- max(all_levels)
      } else {
        ref_level <- sort(all_levels, decreasing = TRUE)[1]
      }

      other_levels <- setdiff(all_levels, ref_level)
      other_levels <- sort(other_levels)

      # Create dummy for each non-reference level
      for (level in other_levels) {
        dummy_name <- paste0("z_", var, "_", level)
        df_work[[dummy_name]] <- as.numeric(data[[var]] == level)
      }
    }
  }

  return(df_work)
}


#' Add Unprocessed Variables from Original Data
#' @keywords internal
add_unprocessed_vars <- function(df_work, data, outcome_var, event_var,
                                 treatment_var, continuous_vars, factor_vars,
                                 verbose) {

  # Identify processed variables
  processed <- c(continuous_vars, factor_vars, outcome_var, event_var,
                 treatment_var)

  # Find unprocessed variables
  unprocessed <- setdiff(names(data), processed)

  if (length(unprocessed) > 0) {
    # Add them to df_work
    for (var in unprocessed) {
      if (!var %in% names(df_work)) {  # Avoid duplicates
        df_work[[var]] <- data[[var]]
      }
    }

    if (verbose) {
      cat("\nAdded", length(unprocessed), "unprocessed variables:",
          paste(unprocessed, collapse = ", "), "\n")
    }
  }

  return(df_work)
}


#' Define Subgroups with Flexible Cutpoints
#' @keywords internal
define_subgroups <- function(df_work, data, subgroup_vars, subgroup_cuts,
                             continuous_vars, model, verbose) {

  flag_harm <- rep(0, nrow(df_work))
  interaction_term <- NULL
  definitions <- list()

  if (model == "alt" && !is.null(subgroup_vars)) {
    # Create subgroup indicators
    subgroup_indicators <- list()

    if (verbose) {
      cat("\n=== Subgroup Definitions ===\n")
    }

    for (var in subgroup_vars) {
      # Get the cutpoint specification for this variable
      cut_spec <- subgroup_cuts[[var]]

      if (var %in% continuous_vars || is.numeric(data[[var]])) {
        # Process continuous variable
        result <- process_continuous_subgroup(data[[var]], cut_spec, var,
                                              verbose)
        subgroup_indicators[[var]] <- result$indicator
        definitions[[var]] <- result$definition

      } else {
        # Process factor variable
        result <- process_factor_subgroup(data[[var]], cut_spec, var,
                                          verbose)
        subgroup_indicators[[var]] <- result$indicator
        definitions[[var]] <- result$definition
      }
    }

    # Validate that all indicators are non-NULL
    null_indicators <- sapply(subgroup_indicators, is.null)
    if (any(null_indicators)) {
      null_vars <- names(subgroup_indicators)[null_indicators]
      stop(paste("Failed to create subgroup indicators for variables:",
                 paste(null_vars, collapse = ", ")))
    }

    # Create harm flag (all subgroup conditions met)
    flag_harm <- as.numeric(Reduce("&", subgroup_indicators))

    # Create interaction term
    if (length(subgroup_indicators) > 0) {
      interaction_term <- df_work$treat * flag_harm
    }

    if (verbose) {
      cat("\nOverall subgroup (all conditions met):\n")
      cat("  Size:", sum(flag_harm), "out of", nrow(df_work), "\n")
      cat("  Proportion:", round(mean(flag_harm), 3), "\n")
    }
  }

  return(list(
    flag_harm = flag_harm,
    interaction_term = interaction_term,
    definitions = definitions
  ))
}


#' Process Continuous Variable for Subgroup Definition
#' @keywords internal
process_continuous_subgroup <- function(var_data, cut_spec, var_name, verbose) {

  indicator <- process_cutpoint(var_data, cut_spec, var_name, verbose)
  definition <- format_continuous_definition(var_data, cut_spec, var_name)

  if (verbose) {
    cat("  ", definition, "\n")
    cat("    Proportion in subgroup:",
        round(mean(indicator, na.rm = TRUE), 3), "\n")
  }

  return(list(indicator = indicator, definition = definition))
}


#' Process Factor Variable for Subgroup Definition
#' @keywords internal
process_factor_subgroup <- function(var_data, cut_spec, var_name, verbose) {

  # Initialize to avoid NULL issues
  indicator <- NULL
  definition <- NULL

  if (!is.null(cut_spec)) {
    # Handle simple numeric or character value (e.g., meno = 0 or meno = "1")
    if ((is.numeric(cut_spec) || is.character(cut_spec) || is.factor(cut_spec)) &&
        length(cut_spec) == 1 && !is.list(cut_spec)) {
      indicator <- var_data == cut_spec
      definition <- paste(var_name, "==", cut_spec)

    } else if (is.list(cut_spec) && !is.null(cut_spec$type)) {
      # Handle list specifications
      if (cut_spec$type == "multiple") {
        indicator <- var_data %in% cut_spec$values
        definition <- paste(var_name, "in",
                            paste(cut_spec$values, collapse = ", "))
      } else {
        # Unknown list type, use default
        first_level <- levels(as.factor(var_data))[1]
        indicator <- var_data == first_level
        definition <- paste(var_name, "==", first_level)
      }

    } else {
      # If cut_spec doesn't match expected patterns, use default
      first_level <- levels(as.factor(var_data))[1]
      indicator <- var_data == first_level
      definition <- paste(var_name, "==", first_level)
    }
  } else {
    # Default: use first level
    first_level <- levels(as.factor(var_data))[1]
    indicator <- var_data == first_level
    definition <- paste(var_name, "==", first_level)
  }

  if (verbose) {
    cat("  ", definition, "\n")
    cat("    Proportion in subgroup:",
        round(mean(indicator, na.rm = TRUE), 3), "\n")
  }

  return(list(indicator = indicator, definition = definition))
}


#' Process Cutpoint Specification for Subgroup Definition
#' @keywords internal
process_cutpoint <- function(var_data, cut_spec, var_name = "",
                             verbose = FALSE) {

  # If cut_spec is a simple numeric value, treat as fixed cutpoint
  if (is.numeric(cut_spec) && length(cut_spec) == 1) {
    return(var_data <= cut_spec)
  }

  # If it's a list, process based on type
  if (is.list(cut_spec)) {
    cut_type <- cut_spec$type

    if (cut_type == "quantile") {
      # Quantile-based cutpoint
      cutpoint <- quantile(var_data, probs = cut_spec$value, na.rm = TRUE)
      return(var_data <= cutpoint)

    } else if (cut_type == "function") {
      # Function-based cutpoint (e.g., median, mean)
      cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
      return(var_data <= cutpoint)

    } else if (cut_type == "range") {
      # Range-based (between min and max)
      return(var_data >= cut_spec$min & var_data <= cut_spec$max)

    } else if (cut_type == "greater") {
      # Greater than cutpoint
      cutpoint <- NULL
      if (!is.null(cut_spec$value)) {
        cutpoint <- cut_spec$value
      } else if (!is.null(cut_spec$quantile)) {
        cutpoint <- quantile(var_data, probs = cut_spec$quantile,
                             na.rm = TRUE)
      } else if (!is.null(cut_spec$fun)) {
        cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
      }

      if (is.null(cutpoint)) {
        stop(paste("For 'greater' type, must provide 'value', 'quantile', or 'fun' for variable:",
                   var_name))
      }

      return(var_data > cutpoint)

    } else if (cut_type == "multiple") {
      # Multiple cutpoints (var in specified values)
      return(var_data %in% cut_spec$values)

    } else if (cut_type == "custom") {
      # Custom function that returns logical vector
      return(cut_spec$fun(var_data))

    } else {
      stop(paste("Unknown cutpoint type:", cut_type, "for variable:",
                 var_name))
    }
  }

  # Default: if no specification, use median
  if (is.null(cut_spec)) {
    if (verbose) cat("  Using median as default cutpoint for", var_name, "\n")
    return(var_data <= median(var_data, na.rm = TRUE))
  }

  stop(paste("Invalid cutpoint specification for variable:", var_name))
}


#' Format Continuous Variable Definition for Display
#' @keywords internal
format_continuous_definition <- function(var_data, cut_spec, var_name) {

  if (is.numeric(cut_spec) && length(cut_spec) == 1) {
    actual_cutpoint <- cut_spec
    return(paste(var_name, "<=", actual_cutpoint))

  } else if (is.list(cut_spec)) {
    if (cut_spec$type == "quantile") {
      actual_cutpoint <- quantile(var_data, probs = cut_spec$value,
                                  na.rm = TRUE)
      return(paste(var_name, "<=", round(actual_cutpoint, 2),
                   "(", cut_spec$value * 100, "th percentile)"))

    } else if (cut_spec$type == "function") {
      actual_cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
      fun_name <- deparse(substitute(cut_spec$fun))
      return(paste(var_name, "<=", round(actual_cutpoint, 2),
                   "(", fun_name, ")"))

    } else if (cut_spec$type == "range") {
      return(paste(cut_spec$min, "<=", var_name, "<=", cut_spec$max))

    } else if (cut_spec$type == "greater") {
      actual_cutpoint <- NULL
      if (!is.null(cut_spec$value)) {
        actual_cutpoint <- cut_spec$value
      } else if (!is.null(cut_spec$quantile)) {
        actual_cutpoint <- quantile(var_data, probs = cut_spec$quantile,
                                    na.rm = TRUE)
      } else if (!is.null(cut_spec$fun)) {
        actual_cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
      }

      if (!is.null(actual_cutpoint)) {
        return(paste(var_name, ">", round(actual_cutpoint, 2)))
      } else {
        return(paste(var_name, "> (unspecified cutpoint)"))
      }

    } else if (cut_spec$type == "custom") {
      return(paste(var_name, "(custom function)"))
    }
  }

  return(paste(var_name, "(unknown specification)"))
}




#' Generate Super Population and Calculate Linear Predictors
#' @keywords internal
generate_super_population <- function(df_work, n_super, draw_treatment,
                                      gamma, b0, mu, tau, verbose,
                                      spline_info = NULL) {  # ADD spline_info

  # Sample with replacement to create super population
  idx_sample <- sample(1:nrow(df_work), size = n_super, replace = TRUE)
  df_super <- df_work[idx_sample, ]

  if (verbose) {
    cat("\nOverall subgroup in super-population (all conditions met):\n")
    cat("  Proportion:", round(mean(df_super$flag_harm), 3), "\n")
  }

  # Handle treatment assignment
  if (draw_treatment) {
    # Assign treatment simple 1/2
    n_treat <- round(n_super / 2)
    n_control <- n_super - n_treat
    df_super$treat[1:n_treat] <- 1
    df_super$treat[(n_treat + 1):n_super] <- 0

    # If spline is present, need to recalculate spline interaction terms
    if (!is.null(spline_info)) {
      spline_var <- spline_info$var
      knot <- spline_info$knot

      # Recreate spline terms with new treatment assignment
      df_super[[paste0(spline_var, "_treat")]] <-
        df_super[[spline_var]] * df_super$treat

      df_super[[paste0(spline_var, "_k_treat")]] <-
        df_super[[paste0(spline_var, "_k")]] * df_super$treat
    }
  }

  # Reset IDs
  df_super$id <- 1:n_super

  # Calculate linear predictors for potential outcomes
  covariate_cols <- grep("^z_", names(df_super), value = TRUE)
  df_super <- calculate_linear_predictors(df_super, covariate_cols, gamma, b0,
                                          spline_info)  # PASS spline_info

  return(df_super)
}



#' Calculate Linear Predictors for Potential Outcomes
#' @keywords internal
calculate_linear_predictors <- function(df_super, covariate_cols, gamma, b0,
                                        spline_info = NULL) {

  # Separate spline interaction terms from regular covariates
  if (!is.null(spline_info)) {
    spline_var <- spline_info$var
    spline_interaction_terms <- c(
      paste0(spline_var, "_treat"),
      paste0(spline_var, "_k"),
      paste0(spline_var, "_k_treat")
    )
    # Base covariates (excluding spline interactions)
    base_covariate_cols <- setdiff(covariate_cols, spline_interaction_terms)
  } else {
    base_covariate_cols <- covariate_cols
    spline_interaction_terms <- character(0)
  }

  # Check if subgroup interaction term exists
  has_interaction <- "treat_harm" %in% names(gamma)

  if (has_interaction) {
    # Recalculate interaction for super population
    df_super$treat_harm <- df_super$treat * df_super$flag_harm
  }

  # ==========================================================================
  # Build Design Matrices for Observed Data
  # ==========================================================================

  # Start with base covariates and treatment
  X_super <- as.matrix(df_super[, c("treat", base_covariate_cols)])

  # Add spline terms if present (these were created in df_super already)
  if (!is.null(spline_info)) {
    for (spline_term in spline_interaction_terms) {
      if (spline_term %in% names(df_super)) {
        X_super <- cbind(X_super, df_super[[spline_term]])
        colnames(X_super)[ncol(X_super)] <- spline_term
      }
    }
  }

  # Add subgroup interaction if present
  if (has_interaction) {
    X_super <- cbind(X_super, treat_harm = df_super$treat_harm)
  }

  # ==========================================================================
  # Build Design Matrices for Potential Outcomes Under Treatment (A=1)
  # ==========================================================================

  X_treat <- as.matrix(df_super[, c("treat", base_covariate_cols)])
  X_treat[, "treat"] <- 1  # Set treatment to 1

  # Add spline terms for treatment arm
  if (!is.null(spline_info)) {
    spline_var <- spline_info$var

    # For A=1: spline_var_treat = spline_var * 1 = spline_var
    if (paste0(spline_var, "_treat") %in% names(df_super)) {
      X_treat <- cbind(X_treat,
                       df_super[[spline_var]])  # spline_var * treat where treat=1
      colnames(X_treat)[ncol(X_treat)] <- paste0(spline_var, "_treat")
    }

    # spline_var_k stays the same (doesn't depend on treatment)
    if (paste0(spline_var, "_k") %in% names(df_super)) {
      X_treat <- cbind(X_treat,
                       df_super[[paste0(spline_var, "_k")]])
      colnames(X_treat)[ncol(X_treat)] <- paste0(spline_var, "_k")
    }

    # For A=1: spline_var_k_treat = spline_var_k * 1 = spline_var_k
    if (paste0(spline_var, "_k_treat") %in% names(df_super)) {
      X_treat <- cbind(X_treat,
                       df_super[[paste0(spline_var, "_k")]])  # spline_var_k * treat where treat=1
      colnames(X_treat)[ncol(X_treat)] <- paste0(spline_var, "_k_treat")
    }
  }

  # Add subgroup interaction for treatment arm
  if (has_interaction) {
    X_treat <- cbind(X_treat, treat_harm = df_super$flag_harm)  # treat=1 * flag_harm
  }

  # ==========================================================================
  # Build Design Matrices for Potential Outcomes Under Control (A=0)
  # ==========================================================================

  X_control <- as.matrix(df_super[, c("treat", base_covariate_cols)])
  X_control[, "treat"] <- 0  # Set treatment to 0

  # Add spline terms for control arm
  if (!is.null(spline_info)) {
    spline_var <- spline_info$var

    # For A=0: spline_var_treat = spline_var * 0 = 0
    if (paste0(spline_var, "_treat") %in% names(df_super)) {
      X_control <- cbind(X_control,
                         rep(0, nrow(df_super)))  # spline_var * treat where treat=0
      colnames(X_control)[ncol(X_control)] <- paste0(spline_var, "_treat")
    }

    # spline_var_k stays the same (doesn't depend on treatment)
    if (paste0(spline_var, "_k") %in% names(df_super)) {
      X_control <- cbind(X_control,
                         df_super[[paste0(spline_var, "_k")]])
      colnames(X_control)[ncol(X_control)] <- paste0(spline_var, "_k")
    }

    # For A=0: spline_var_k_treat = spline_var_k * 0 = 0
    if (paste0(spline_var, "_k_treat") %in% names(df_super)) {
      X_control <- cbind(X_control,
                         rep(0, nrow(df_super)))  # spline_var_k * treat where treat=0
      colnames(X_control)[ncol(X_control)] <- paste0(spline_var, "_k_treat")
    }
  }

  # Add subgroup interaction for control arm
  if (has_interaction) {
    X_control <- cbind(X_control, treat_harm = 0)  # treat=0 * flag_harm = 0
  }

  # ==========================================================================
  # Verify Column Alignment
  # ==========================================================================

  # Ensure all matrices have the same columns in the same order
  expected_cols <- names(gamma)

  if (!all(colnames(X_super) == expected_cols)) {
    stop("Column mismatch in X_super")
  }
  if (!all(colnames(X_treat) == expected_cols)) {
    stop("Column mismatch in X_treat")
  }
  if (!all(colnames(X_control) == expected_cols)) {
    stop("Column mismatch in X_control")
  }

  # ==========================================================================
  # Calculate Linear Predictors
  # ==========================================================================

  # Linear predictors per AFT (log-time scale)
  df_super$lin_pred_1 <- as.vector(X_treat %*% gamma)
  df_super$lin_pred_0 <- as.vector(X_control %*% gamma)
  df_super$lin_pred_obs <- as.vector(X_super %*% gamma)

  # PO log-hazards (excluding baseline hazard)
  # Used for calculating empirical version of CDEs (controlled direct effects)
  df_super$theta_0 <- as.vector(X_control %*% b0)
  df_super$theta_1 <- as.vector(X_treat %*% b0)
  df_super$loghr_po <- df_super$theta_1 - df_super$theta_0

  return(df_super)
}

#' Calculate Hazard Ratios from Potential Outcomes
#' @keywords internal
calculate_hazard_ratios <- function(df_super, n_super, mu, tau, model,
                                    verbose) {

  # Generate potential outcomes for HR calculation
  epsilon <- log(rexp(n_super))  # Extreme value distribution

  # Under treatment
  logT_1 <- mu + tau * epsilon + df_super$lin_pred_1
  T_1 <- exp(logT_1)

  # Under control
  logT_0 <- mu + tau * epsilon + df_super$lin_pred_0
  T_0 <- exp(logT_0)

  # Create data frame for Cox models
  df_temp <- data.frame(
    time = c(T_1, T_0),
    event = 1,
    treat = c(rep(1, n_super), rep(0, n_super)),
    flag_harm = rep(df_super$flag_harm, 2)
  )

  # Calculate overall HR
  hr_overall <- exp(coxph(Surv(time, event) ~ treat,
                          data = df_temp)$coefficients)

  # Calculate average hazard ratios (AHR)
  AHR <- with(df_super, exp(mean(loghr_po)))
  AHR_harm <- with(subset(df_super, flag_harm == 1), exp(mean(loghr_po)))
  AHR_no_harm <- with(subset(df_super, flag_harm == 0), exp(mean(loghr_po)))

  hr_results <- list(
    overall = hr_overall,
    AHR = AHR,
    AHR_harm = AHR_harm,
    AHR_no_harm = AHR_no_harm
  )

  # Calculate subgroup-specific HRs if applicable
  if (model == "alt" && sum(df_super$flag_harm) > 0) {
    hr_harm <- exp(coxph(Surv(time, event) ~ treat,
                         data = subset(df_temp, flag_harm == 1))$coefficients)
    hr_no_harm <- exp(coxph(Surv(time, event) ~ treat,
                            data = subset(df_temp, flag_harm == 0))$coefficients)

    hr_results$harm_subgroup <- hr_harm
    hr_results$no_harm_subgroup <- hr_no_harm

    if (verbose) {
      cat("\n=== Hazard Ratios (super popln) ===\n")
      cat("Overall HR:", round(hr_overall, 3), "\n")
      cat("Causal AHR:", round(AHR, 3), "\n")
      cat("Harm subgroup HR:", round(hr_harm, 3), "\n")
      cat("Harm subgroup AHR:", round(AHR_harm, 3), "\n")
      cat("No-harm subgroup HR:", round(hr_no_harm, 3), "\n")
      cat("No-harm subgroup AHR:", round(AHR_no_harm, 3), "\n")
    }
  } else if (verbose) {
    cat("\n=== Hazard Ratios (super popln) ===\n")
    cat("Overall HR:", round(hr_overall, 3), "\n")
    cat("Causal AHR:", round(AHR, 3), "\n")
  }

  return(hr_results)
}


#' Prepare Censoring Model Parameters
#' @keywords internal
prepare_censoring_model <- function(df_work, cens_type, cens_params,
                                    df_super, gamma, b0, spline_info = NULL, verbose) {

  cens_model <- NULL

  # Get covariate columns, excluding spline interaction terms
  covariate_cols <- grep("^z_", names(df_work), value = TRUE)

  # Exclude spline terms if spline was used
  if (!is.null(spline_info)) {
    spline_var <- spline_info$var
    spline_terms_to_exclude <- c(
      paste0(spline_var, "_treat"),
      paste0(spline_var, "_k"),
      paste0(spline_var, "_k_treat")
    )

    covariate_cols <- setdiff(covariate_cols, spline_terms_to_exclude)


    # # Exclude spline terms if spline was used
    # if (!is.null(spline_info)) {
    #   spline_var <- spline_info$var
    #   # Match any column that contains the spline variable followed by interaction terms
    #   spline_pattern <- paste0("^", spline_var, "_(treat|k|k_treat)$")
    #   spline_cols <- grep(spline_pattern, covariate_cols, value = TRUE)
    #   covariate_cols <- setdiff(covariate_cols, spline_cols)
    # }


    if (verbose) {
      cat("\nExcluded spline interaction terms from censoring model:\n")
      cat("  ", paste(spline_terms_to_exclude, collapse = ", "), "\n")
    }
  }

  if (cens_type == "weibull") {
    # Fit censoring model WITHOUT interaction terms (treat + covariates only)
    X_cens <- as.matrix(df_work[, c("treat", covariate_cols)])

    fit_cens <- survreg(Surv(y, 1 - event) ~ X_cens,
                        data = df_work,
                        dist = "weibull")

    mu_cens <- coef(fit_cens)[1]
    tau_cens <- fit_cens$scale
    gamma_cens <- coef(fit_cens)[-1]

    # Store censoring parameters
    cens_model <- list(
      mu = mu_cens,
      tau = tau_cens,
      gamma = gamma_cens,
      type = "weibull"
    )

    # Calculate censoring linear predictors for super population
    # Build matrices for censoring (no interaction term, no spline interactions)
    X_cens_super_treat <- as.matrix(df_super[, c("treat", covariate_cols)])
    X_cens_super_treat[, "treat"] <- 1

    X_cens_super_control <- as.matrix(df_super[, c("treat", covariate_cols)])
    X_cens_super_control[, "treat"] <- 0

    df_super$lin_pred_cens_1 <- X_cens_super_treat %*% gamma_cens
    df_super$lin_pred_cens_0 <- X_cens_super_control %*% gamma_cens

  } else if (cens_type == "uniform") {
    # Use provided or default uniform censoring parameters
    if (is.null(cens_params$min) || is.null(cens_params$max)) {
      # Default: use range of observed times
      cens_params$min <- min(df_work$y) * 0.5
      cens_params$max <- max(df_work$y) * 1.5
    }

    cens_model <- list(
      min = cens_params$min,
      max = cens_params$max,
      type = "uniform"
    )
  }

  return(list(
    cens_model = cens_model,
    df_super = df_super
  ))
}


#' Assemble Final Results Object
#' @keywords internal
assemble_results <- function(df_super,
                             mu,
                             tau,
                             gamma,
                             b0,
                             cens_model,
                             subgroup_vars,
                             subgroup_cuts,
                             subgroup_definitions,
                             hr_results,
                             continuous_vars,
                             factor_vars,
                             model,
                             n_super,
                             seed,
                             spline_info = NULL) {  # DEFAULT = NULL

  # Model parameters
  model_params <- list(
    mu = mu,
    tau = tau,
    gamma = gamma,
    b0 = b0,
    censoring = cens_model,
    spline_info = spline_info  # Can be NULL
  )

  # Subgroup information
  subgroup_info <- list(
    vars = subgroup_vars,
    cuts = subgroup_cuts,
    definitions = subgroup_definitions,
    size = sum(df_super$flag_harm),
    proportion = mean(df_super$flag_harm)
  )

  # Get covariate column names
  covariate_cols <- grep("^z_", names(df_super), value = TRUE)

  # Analysis variables (for downstream use)
  analysis_vars <- list(
    continuous = continuous_vars,
    factor = factor_vars,
    covariates = covariate_cols,
    treatment = "treat",
    outcome = "y_sim",
    event = "event_sim"
  )

  # Return comprehensive results
  results <- list(
    df_super = df_super,
    model_params = model_params,
    subgroup_info = subgroup_info,
    hazard_ratios = hr_results,
    analysis_vars = analysis_vars,
    model_type = model,
    n_super = n_super,
    seed = seed
  )

  class(results) <- c("aft_dgm_flex", "list")

  return(results)
}






#' Assemble Final Results Object
#' @keywords internal
assemble_results_legacy <- function(df_super, mu, tau, gamma, b0, cens_model,
                             subgroup_vars, subgroup_cuts, subgroup_definitions,
                             hr_results, continuous_vars, factor_vars,
                             model, n_super, seed) {

  # Model parameters
  model_params <- list(
    mu = mu,
    tau = tau,
    gamma = gamma,
    b0 = b0,
    censoring = cens_model
  )

  # Subgroup information
  subgroup_info <- list(
    vars = subgroup_vars,
    cuts = subgroup_cuts,
    definitions = subgroup_definitions,
    size = sum(df_super$flag_harm),
    proportion = mean(df_super$flag_harm)
  )

  # Get covariate column names
  covariate_cols <- grep("^z_", names(df_super), value = TRUE)

  # Analysis variables (for downstream use)
  analysis_vars <- list(
    continuous = continuous_vars,
    factor = factor_vars,
    covariates = covariate_cols,
    treatment = "treat",
    outcome = "y_sim",
    event = "event_sim"
  )

  # Return comprehensive results
  results <- list(
    df_super = df_super,
    model_params = model_params,
    subgroup_info = subgroup_info,
    hazard_ratios = hr_results,
    analysis_vars = analysis_vars,
    model_type = model,
    n_super = n_super,
    seed = seed
  )

  class(results) <- c("aft_dgm_flex", "list")

  return(results)
}
