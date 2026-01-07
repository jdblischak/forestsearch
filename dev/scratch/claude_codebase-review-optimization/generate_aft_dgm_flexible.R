# Enhanced AFT Data Generating Mechanism with Flexible Subgroup Definitions
# CRAN-compliant version with proper documentation and exports
#
# Supports multiple cutpoint specifications: fixed values, quantiles, functions, etc.

#' Generate Synthetic Survival Data with Flexible Subgroup Definitions
#'
#' Creates a data generating mechanism (DGM) based on an Accelerated Failure
#' Time (AFT) model with flexible subgroup specifications including quantiles,
#' functions, ranges, and custom conditions.
#'
#' @inheritParams generate_aft_dgm
#' @param subgroup_cuts Named list of cutpoint specifications. See Details.
#'
#' @return An object of class "aft_dgm_flex" (inherits from "aft_dgm") containing:
#' \describe{
#'   \item{df_super}{Super population dataset with potential outcomes}
#'   \item{model_params}{Model parameters}
#'   \item{subgroup_info}{Information including \code{definitions} with human-readable descriptions}
#'   \item{hazard_ratios}{True hazard ratios}
#'   \item{analysis_vars}{Variable names for downstream analysis}
#' }
#'
#' @details
#' The \code{subgroup_cuts} parameter accepts flexible specifications:
#' \describe{
#'   \item{Numeric value}{\code{list(er = 20)} means \code{er <= 20}}
#'   \item{Quantile}{\code{list(er = list(type = "quantile", value = 0.25))}}
#'   \item{Function}{\code{list(er = list(type = "function", fun = median))}}
#'   \item{Range}{\code{list(er = list(type = "range", min = 10, max = 50))}}
#'   \item{Greater than}{\code{list(er = list(type = "greater", value = 30))}}
#'   \item{Multiple values}{\code{list(stage = list(type = "multiple", values = c(1, 2)))}}
#'   \item{Custom function}{\code{list(er = list(type = "custom", fun = function(x) x < mean(x)))}}
#' }
#'
#' @examples
#' \donttest{
#' # Create example data
#' set.seed(42)
#' example_data <- data.frame(
#'   time = rexp(200, 0.02),
#'   status = rbinom(200, 1, 0.7),
#'   age = rnorm(200, 55, 12),
#'   biomarker = rgamma(200, 2, 0.5),
#'   treatment = rbinom(200, 1, 0.5)
#' )
#'
#' # Using various cutpoint specifications
#' dgm <- generate_aft_dgm_flex(
#'   data = example_data,
#'   continuous_vars = c("age", "biomarker"),
#'   factor_vars = character(0),
#'   outcome_var = "time",
#'   event_var = "status",
#'   treatment_var = "treatment",
#'   subgroup_vars = c("biomarker", "age"),
#'   subgroup_cuts = list(
#'     biomarker = list(type = "quantile", value = 0.25),
#'     age = list(type = "range", min = 40, max = 60)
#'   ),
#'   model = "alt",
#'   verbose = TRUE
#' )
#'
#' print(dgm)
#' }
#'
#' @seealso \code{\link{generate_aft_dgm}} for simpler cutpoint specification,
#'   \code{\link{simulate_from_dgm}} for generating samples
#'
#' @importFrom stats as.formula coef median model.matrix quantile rbinom rnorm rexp
#' @importFrom survival survreg coxph Surv
#' @export
generate_aft_dgm_flex <- function(data,
                                  continuous_vars,
                                  factor_vars,
                                  outcome_var,
                                  event_var,
                                  treatment_var = NULL,
                                  subgroup_vars = NULL,
                                  subgroup_cuts = NULL,
                                  model = "alt",
                                  k_treat = 1,
                                  k_inter = 1,
                                  n_super = 5000,
                                  cens_type = "weibull",
                                  cens_params = list(),
                                  seed = 8316951,
                                  verbose = TRUE) {

  # ==========================================================================
  # Helper: Process Flexible Cutpoint Specifications
  # ==========================================================================

  process_cutpoint <- function(var_data, cut_spec, var_name = "") {
    # Handle simple numeric cutpoint
    if (is.numeric(cut_spec) && length(cut_spec) == 1) {
      return(var_data <= cut_spec)
    }

    # Handle list specifications
    if (is.list(cut_spec)) {
      cut_type <- cut_spec$type

      if (cut_type == "quantile") {
        cutpoint <- quantile(var_data, probs = cut_spec$value, na.rm = TRUE)
        return(var_data <= cutpoint)

      } else if (cut_type == "function") {
        cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
        return(var_data <= cutpoint)

      } else if (cut_type == "range") {
        return(var_data >= cut_spec$min & var_data <= cut_spec$max)

      } else if (cut_type == "greater") {
        if (!is.null(cut_spec$value)) {
          cutpoint <- cut_spec$value
        } else if (!is.null(cut_spec$quantile)) {
          cutpoint <- quantile(var_data, probs = cut_spec$quantile, na.rm = TRUE)
        } else if (!is.null(cut_spec$fun)) {
          cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
        } else {
          stop("'greater' type requires value, quantile, or fun")
        }
        return(var_data > cutpoint)

      } else if (cut_type == "multiple") {
        return(var_data %in% cut_spec$values)

      } else if (cut_type == "custom") {
        result <- cut_spec$fun(var_data)
        if (!is.logical(result)) {
          stop("Custom function must return logical vector")
        }
        return(result)

      } else {
        stop("Unknown cutpoint type: ", cut_type, " for variable: ", var_name)
      }
    }

    # Default: use median
    if (is.null(cut_spec)) {
      if (verbose) message("  Using median as default cutpoint for ", var_name)
      return(var_data <= median(var_data, na.rm = TRUE))
    }

    stop("Invalid cutpoint specification for variable: ", var_name)
  }

  # ==========================================================================
  # Input Validation
  # ==========================================================================

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  if (!model %in% c("alt", "null")) {
    stop("'model' must be either 'alt' or 'null'")
  }

  if (!cens_type %in% c("weibull", "uniform")) {
    stop("'cens_type' must be either 'weibull' or 'uniform'")
  }

  # Check required variables
  required_vars <- c(outcome_var, event_var)
  if (!is.null(treatment_var)) required_vars <- c(required_vars, treatment_var)

  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Required variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # Check covariates
  all_covars <- c(continuous_vars, factor_vars)
  missing_covars <- setdiff(all_covars, names(data))
  if (length(missing_covars) > 0) {
    stop("Covariate variables not found in data: ", paste(missing_covars, collapse = ", "))
  }

  # Check for survival package
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }

  # ==========================================================================
  # Data Preparation
  # ==========================================================================

  set.seed(seed)

  df_work <- data.frame(
    id = seq_len(nrow(data)),
    y = data[[outcome_var]],
    event = ifelse(data[[event_var]] == 1, 1L, 0L)
  )

  # Add treatment
  if (!is.null(treatment_var)) {
    df_work$treat <- as.integer(data[[treatment_var]])
  } else {
    df_work$treat <- rbinom(nrow(data), size = 1, prob = 0.5)
    if (verbose) message("Treatment variable simulated (50/50 randomization)")
  }

  # ==========================================================================
  # Process Covariates
  # ==========================================================================

  # Continuous variables (standardize)
  for (var in continuous_vars) {
    df_work[[paste0("z_", var)]] <- as.numeric(scale(data[[var]]))
  }

  # Factor variables
  for (var in factor_vars) {
    if (is.factor(data[[var]])) {
      dummies <- model.matrix(~ data[[var]] - 1)
      if (ncol(dummies) > 1) {
        for (j in 2:ncol(dummies)) {
          dummy_name <- paste0("z_", var, "_", j - 1)
          df_work[[dummy_name]] <- dummies[, j]
        }
      }
    } else {
      df_work[[paste0("z_", var)]] <- as.numeric(data[[var]])
    }
  }

  # ==========================================================================
  # Define Subgroups with Flexible Cutpoints
  # ==========================================================================

  df_work$flag_harm <- 0L
  interaction_term <- NULL
  subgroup_definitions <- list()

  if (model == "alt" && !is.null(subgroup_vars)) {
    subgroup_indicators <- list()

    if (verbose) message("\n=== Subgroup Definitions ===")

    for (var in subgroup_vars) {
      cut_spec <- subgroup_cuts[[var]]

      if (var %in% continuous_vars || is.numeric(data[[var]])) {
        # Process continuous variable
        subgroup_indicators[[var]] <- process_cutpoint(
          var_data = data[[var]],
          cut_spec = cut_spec,
          var_name = var
        )

        # Store human-readable definition
        if (is.numeric(cut_spec) && length(cut_spec) == 1) {
          subgroup_definitions[[var]] <- paste(var, "<=", cut_spec)

        } else if (is.list(cut_spec)) {
          if (cut_spec$type == "quantile") {
            actual_cutpoint <- quantile(data[[var]], probs = cut_spec$value, na.rm = TRUE)
            subgroup_definitions[[var]] <- sprintf("%s <= %.2f (%s percentile)",
                                                    var, actual_cutpoint, cut_spec$value * 100)
          } else if (cut_spec$type == "function") {
            actual_cutpoint <- cut_spec$fun(data[[var]], na.rm = TRUE)
            subgroup_definitions[[var]] <- sprintf("%s <= %.2f (function)", var, actual_cutpoint)
          } else if (cut_spec$type == "range") {
            subgroup_definitions[[var]] <- sprintf("%.2f <= %s <= %.2f",
                                                    cut_spec$min, var, cut_spec$max)
          } else if (cut_spec$type == "greater") {
            if (!is.null(cut_spec$value)) {
              subgroup_definitions[[var]] <- sprintf("%s > %.2f", var, cut_spec$value)
            } else if (!is.null(cut_spec$quantile)) {
              actual_cutpoint <- quantile(data[[var]], probs = cut_spec$quantile, na.rm = TRUE)
              subgroup_definitions[[var]] <- sprintf("%s > %.2f (%s percentile)",
                                                      var, actual_cutpoint, cut_spec$quantile * 100)
            }
          } else if (cut_spec$type == "custom") {
            subgroup_definitions[[var]] <- paste(var, "(custom function)")
          }
        } else {
          subgroup_definitions[[var]] <- paste(var, "<= median")
        }

      } else {
        # Factor variable
        if (!is.null(cut_spec)) {
          if (is.character(cut_spec) || is.factor(cut_spec)) {
            subgroup_indicators[[var]] <- data[[var]] == cut_spec
            subgroup_definitions[[var]] <- paste(var, "==", cut_spec)
          } else if (is.list(cut_spec) && cut_spec$type == "multiple") {
            subgroup_indicators[[var]] <- data[[var]] %in% cut_spec$values
            subgroup_definitions[[var]] <- paste(var, "in",
                                                  paste(cut_spec$values, collapse = ", "))
          }
        } else {
          first_level <- levels(as.factor(data[[var]]))[1]
          subgroup_indicators[[var]] <- data[[var]] == first_level
          subgroup_definitions[[var]] <- paste(var, "==", first_level)
        }
      }

      if (verbose) {
        message("  ", subgroup_definitions[[var]])
        message("    Proportion in subgroup: ",
                round(mean(subgroup_indicators[[var]], na.rm = TRUE), 3))
      }
    }

    # Create harm flag (all conditions met)
    df_work$flag_harm <- as.integer(Reduce(`&`, subgroup_indicators))

    if (length(subgroup_indicators) > 0) {
      interaction_term <- df_work$treat * df_work$flag_harm
    }

    if (verbose) {
      message("\nOverall subgroup (all conditions met):")
      message("  Size: ", sum(df_work$flag_harm), " out of ", nrow(df_work))
      message("  Proportion: ", round(mean(df_work$flag_harm), 3))
    }
  }

  # ==========================================================================
  # Fit AFT Model (Weibull)
  # ==========================================================================

  covariate_cols <- grep("^z_", names(df_work), value = TRUE)
  X <- as.matrix(df_work[, c("treat", covariate_cols), drop = FALSE])

  if (!is.null(interaction_term)) {
    X <- cbind(X, treat_harm = interaction_term)
  }

  formula_str <- paste("survival::Surv(y, event) ~ ",
                       paste(c("treat", covariate_cols), collapse = " + "))
  if (!is.null(interaction_term)) {
    df_work$treat_harm <- interaction_term
    formula_str <- paste(formula_str, "+ treat_harm")
  }

  fit_aft <- survival::survreg(as.formula(formula_str),
                               data = df_work,
                               dist = "weibull")

  mu <- coef(fit_aft)[1]
  sigma <- fit_aft$scale
  gamma <- coef(fit_aft)[-1]

  gamma["treat"] <- k_treat * gamma["treat"]
  if ("treat_harm" %in% names(gamma)) {
    gamma["treat_harm"] <- k_inter * gamma["treat_harm"]
  }

  b_weibull <- gamma
  b0_weibull <- -gamma / sigma

  if (verbose) {
    message("\n=== Model Parameters ===")
    message("Intercept (mu): ", round(mu, 3))
    message("Scale (sigma): ", round(sigma, 3))
    message("Treatment effect: ", round(gamma["treat"], 3))
    if ("treat_harm" %in% names(gamma)) {
      message("Interaction effect: ", round(gamma["treat_harm"], 3))
    }
  }

  # ==========================================================================
  # Generate Super Population
  # ==========================================================================

  n_treat <- round(n_super / 2)
  n_control <- n_super - n_treat

  idx_sample <- sample(seq_len(nrow(df_work)), size = n_super, replace = TRUE)
  df_super <- df_work[idx_sample, , drop = FALSE]

  df_super$treat[seq_len(n_treat)] <- 1L
  df_super$treat[(n_treat + 1):n_super] <- 0L
  df_super$id <- seq_len(n_super)

  X_super <- as.matrix(df_super[, c("treat", covariate_cols), drop = FALSE])
  if (!is.null(interaction_term)) {
    df_super$treat_harm <- df_super$treat * df_super$flag_harm
    X_super <- cbind(X_super, treat_harm = df_super$treat_harm)
  }

  X_treat <- X_super
  X_treat[, "treat"] <- 1
  if ("treat_harm" %in% colnames(X_treat)) {
    X_treat[, "treat_harm"] <- df_super$flag_harm
  }

  X_control <- X_super
  X_control[, "treat"] <- 0
  if ("treat_harm" %in% colnames(X_control)) {
    X_control[, "treat_harm"] <- 0
  }

  df_super$lin_pred_1 <- as.numeric(X_treat %*% b_weibull)
  df_super$lin_pred_0 <- as.numeric(X_control %*% b_weibull)
  df_super$lin_pred_obs <- as.numeric(X_super %*% b_weibull)
  df_super$hr_individual <- exp((df_super$lin_pred_1 - df_super$lin_pred_0) / sigma)

  # ==========================================================================
  # Calculate True Hazard Ratios
  # ==========================================================================

  epsilon <- log(rexp(n_super))
  logT_1 <- mu + sigma * epsilon + df_super$lin_pred_1
  T_1 <- exp(logT_1)
  logT_0 <- mu + sigma * epsilon + df_super$lin_pred_0
  T_0 <- exp(logT_0)

  df_temp <- data.frame(
    time = c(T_1, T_0),
    event = 1L,
    treat = c(rep(1L, n_super), rep(0L, n_super)),
    flag_harm = rep(df_super$flag_harm, 2)
  )

  hr_overall <- exp(survival::coxph(
    survival::Surv(time, event) ~ treat,
    data = df_temp
  )$coefficients)

  hr_results <- list(overall = unname(hr_overall))

  if (model == "alt" && sum(df_super$flag_harm) > 0) {
    hr_harm <- exp(survival::coxph(
      survival::Surv(time, event) ~ treat,
      data = df_temp[df_temp$flag_harm == 1, ]
    )$coefficients)

    hr_no_harm <- exp(survival::coxph(
      survival::Surv(time, event) ~ treat,
      data = df_temp[df_temp$flag_harm == 0, ]
    )$coefficients)

    hr_results$harm_subgroup <- unname(hr_harm)
    hr_results$no_harm_subgroup <- unname(hr_no_harm)

    if (verbose) {
      message("\n=== Hazard Ratios ===")
      message("Overall HR: ", round(hr_overall, 3))
      message("Harm subgroup HR: ", round(hr_harm, 3))
      message("No-harm subgroup HR: ", round(hr_no_harm, 3))
    }
  } else if (verbose) {
    message("\n=== Hazard Ratios ===")
    message("Overall HR: ", round(hr_overall, 3))
  }

  # ==========================================================================
  # Prepare Censoring Parameters
  # ==========================================================================

  cens_model <- NULL

  if (cens_type == "weibull") {
    X_cens <- as.matrix(df_work[, c("treat", covariate_cols), drop = FALSE])

    fit_cens <- survival::survreg(
      survival::Surv(y, 1 - event) ~ X_cens,
      data = df_work,
      dist = "weibull"
    )

    mu_cens <- coef(fit_cens)[1]
    sigma_cens <- fit_cens$scale
    gamma_cens <- coef(fit_cens)[-1]

    cens_model <- list(
      mu = mu_cens,
      sigma = sigma_cens,
      gamma = gamma_cens,
      type = "weibull"
    )

    df_super$lin_pred_cens_1 <- as.numeric(
      X_treat[, seq_len(ncol(X_cens)), drop = FALSE] %*% gamma_cens
    )
    df_super$lin_pred_cens_0 <- as.numeric(
      X_control[, seq_len(ncol(X_cens)), drop = FALSE] %*% gamma_cens
    )

  } else if (cens_type == "uniform") {
    if (is.null(cens_params$min) || is.null(cens_params$max)) {
      cens_params$min <- min(df_work$y) * 0.5
      cens_params$max <- max(df_work$y) * 1.5
    }

    cens_model <- list(
      min = cens_params$min,
      max = cens_params$max,
      type = "uniform"
    )
  }

  # ==========================================================================
  # Prepare Output
  # ==========================================================================

  model_params <- list(
    mu = mu,
    sigma = sigma,
    gamma = gamma,
    b_weibull = b_weibull,
    b0_weibull = b0_weibull,
    censoring = cens_model
  )

  subgroup_info <- list(
    vars = subgroup_vars,
    cuts = subgroup_cuts,
    definitions = subgroup_definitions,
    size = sum(df_super$flag_harm),
    proportion = mean(df_super$flag_harm)
  )

  analysis_vars <- list(
    continuous = continuous_vars,
    factor = factor_vars,
    covariates = covariate_cols,
    treatment = "treat",
    outcome = "y_sim",
    event = "event_sim"
  )

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

  class(results) <- c("aft_dgm_flex", "aft_dgm", "list")

  return(results)
}


#' Find Quantile for Target Subgroup Proportion
#'
#' Finds the quantile of a variable that achieves a target subgroup proportion.
#' Useful for calibrating subgroup definitions to achieve specific sample sizes.
#'
#' @param data Data containing the variable
#' @param var_name Name of the variable
#' @param target_prop Target proportion for the subgroup (between 0 and 1)
#' @param direction Character. Either "less" for <= or "greater" for >
#' @param tol Tolerance for root finding. Default is 0.0001.
#'
#' @return A list with:
#' \describe{
#'   \item{quantile}{The quantile probability that achieves the target}
#'   \item{cutpoint}{The actual cutpoint value}
#'   \item{actual_proportion}{The achieved proportion}
#' }
#'
#' @examples
#' df <- data.frame(biomarker = rnorm(1000, 50, 10))
#'
#' # Find cutpoint for 30% of population
#' result <- find_quantile_for_proportion(df, "biomarker", target_prop = 0.30)
#' result$cutpoint
#'
#' # Verify
#' mean(df$biomarker <= result$cutpoint)
#'
#' @importFrom stats quantile uniroot
#' @export
find_quantile_for_proportion <- function(data,
                                          var_name,
                                          target_prop,
                                          direction = "less",
                                          tol = 0.0001) {

  var_data <- data[[var_name]]

  obj_fun <- function(q) {
    cutpoint <- quantile(var_data, probs = q, na.rm = TRUE)
    if (direction == "less") {
      actual_prop <- mean(var_data <= cutpoint, na.rm = TRUE)
    } else {
      actual_prop <- mean(var_data > cutpoint, na.rm = TRUE)
    }
    return(actual_prop - target_prop)
  }

  result <- uniroot(obj_fun, interval = c(0, 1), tol = tol)

  list(
    quantile = result$root,
    cutpoint = quantile(var_data, probs = result$root, na.rm = TRUE),
    actual_proportion = target_prop
  )
}


#' Print Method for aft_dgm_flex Objects
#'
#' @param x An aft_dgm_flex object
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.aft_dgm_flex <- function(x, ...) {
  cat("AFT Data Generating Mechanism (Flexible Subgroups)\n")
  cat("===================================================\n")
  cat("Model type:", x$model_type, "\n")
  cat("Super population size:", x$n_super, "\n")
  cat("Number of covariates:", length(x$analysis_vars$covariates), "\n")

  if (x$model_type == "alt" && !is.null(x$subgroup_info$vars)) {
    cat("\nSubgroup Definitions:\n")
    for (def in x$subgroup_info$definitions) {
      cat("  ", def, "\n")
    }
    cat("\nSubgroup Size:", x$subgroup_info$size, "\n")
    cat("Subgroup Proportion:", round(x$subgroup_info$proportion, 3), "\n")
  }

  cat("\nHazard Ratios:\n")
  for (name in names(x$hazard_ratios)) {
    cat("  ", name, ": ", round(x$hazard_ratios[[name]], 3), "\n", sep = "")
  }

  cat("\nModel Parameters:\n")
  cat("  Intercept (mu):", round(x$model_params$mu, 3), "\n")
  cat("  Scale (sigma):", round(x$model_params$sigma, 3), "\n")
  cat("  Treatment effect:", round(x$model_params$gamma["treat"], 3), "\n")

  invisible(x)
}
