# Improved and Generalized AFT Data Generating Mechanism
# This is a refactored version of dgm_aftm4_gbsg that can work with any dataset
#
# CRAN-compliant version with proper documentation and exports
# FIXED: enumerate() syntax error, library() calls, proper imports

#' Generate Synthetic Survival Data using AFT Model
#'
#' Creates a data generating mechanism (DGM) based on an Accelerated Failure
#' Time (AFT) model with Weibull distribution. Supports treatment effects
#' with optional subgroup-specific interactions.
#'
#' @param data The input dataset (data.frame)
#' @param continuous_vars Character vector of continuous variable names
#' @param factor_vars Character vector of factor/categorical variable names
#' @param outcome_var Character. Name of the outcome/time variable
#' @param event_var Character. Name of the event/status variable (1 = event)
#' @param treatment_var Character. Name of the treatment variable. If NULL,
#'   treatment will be simulated with 50/50 randomization.
#' @param subgroup_vars Character vector of variables defining the subgroup.
#'   Used only when \code{model = "alt"}.
#' @param subgroup_cuts Named list of cutpoints for subgroup variables.
#'   For continuous variables, observations <= cutpoint are in subgroup.
#'   For factor variables, specify the level(s) in subgroup.
#' @param model Character. Either "alt" (alternative hypothesis with subgroup
#'   effects) or "null" (no subgroup effects). Default is "alt".
#' @param k_treat Numeric. Treatment effect modifier applied to the estimated
#'   treatment coefficient. Default is 1 (no modification).
#' @param k_inter Numeric. Interaction effect modifier for subgroup-treatment
#'   interaction. Default is 1.
#' @param n_super Integer. Size of super population to generate. Default is 5000.
#' @param cens_type Character. Type of censoring distribution: "weibull" or
#'   "uniform". Default is "weibull".
#' @param cens_params List. Parameters for censoring distribution. For uniform,
#'   provide \code{min} and \code{max}. For Weibull, parameters are estimated.
#' @param seed Integer. Random seed for reproducibility. Default is 8316951.
#' @param verbose Logical. Print diagnostic information. Default is TRUE.
#'
#' @return An object of class "aft_dgm" containing:
#' \describe{
#'   \item{df_super}{Super population dataset with potential outcomes}
#'   \item{model_params}{Model parameters (mu, sigma, gamma coefficients)}
#'   \item{subgroup_info}{Information about subgroup definition and size}
#'   \item{hazard_ratios}{True hazard ratios (overall and by subgroup)}
#'   \item{analysis_vars}{Variable names for downstream analysis}
#' }
#'
#' @details
#' The function fits a Weibull AFT model to the input data, then uses the

#' estimated parameters to generate a super population with known treatment
#' effects. This allows for simulation studies where the true effect is known.
#'
#' The Weibull AFT model has the form:
#' \deqn{\log(T) = \mu + X\gamma + \sigma \epsilon}
#' where \eqn{\epsilon} follows a standard extreme value distribution.
#'
#' @examples
#' \donttest
#' # Create example survival data
#' set.seed(42)
#' example_data <- data.frame(
#'   time = rexp(200, rate = 0.02),
#'   status = rbinom(200, 1, 0.7),
#'   age = rnorm(200, 55, 12),
#'   biomarker = rgamma(200, 2, 0.5),
#'   stage = factor(sample(1:3, 200, replace = TRUE)),
#'   treatment = rbinom(200, 1, 0.5)
#' )
#'
#' # Generate DGM with subgroup effect
#' dgm <- generate_aft_dgm(
#'   data = example_data,
#'   continuous_vars = c("age", "biomarker"),
#'   factor_vars = c("stage"),
#'   outcome_var = "time",
#'   event_var = "status",
#'   treatment_var = "treatment",
#'   subgroup_vars = c("biomarker"),
#'   subgroup_cuts = list(biomarker = 3),
#'   model = "alt",
#'   n_super = 2000,
#'   verbose = TRUE
#' )
#'
#' print(dgm)
#' }
#'
#' @seealso \code{\link{simulate_from_dgm}} for generating samples from the DGM
#'
#' @importFrom stats as.formula coef median model.matrix quantile rbinom rnorm
#' @importFrom survival survreg coxph Surv
#' @export
generate_aft_dgm <- function(data,
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

  # Check that required variables exist
  required_vars <- c(outcome_var, event_var)
  if (!is.null(treatment_var)) required_vars <- c(required_vars, treatment_var)

  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Required variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # Check continuous and factor variables
  all_covars <- c(continuous_vars, factor_vars)
  missing_covars <- setdiff(all_covars, names(data))
  if (length(missing_covars) > 0) {
    stop("Covariate variables not found in data: ", paste(missing_covars, collapse = ", "))
  }

  # ==========================================================================
  # Data Preparation
  # ==========================================================================

  set.seed(seed)

  # Create working dataset
  df_work <- data.frame(
    id = seq_len(nrow(data)),
    y = data[[outcome_var]],
    event = ifelse(data[[event_var]] == 1, 1L, 0L)
  )

  # Add treatment (or simulate if not provided)
  if (!is.null(treatment_var)) {
    df_work$treat <- as.integer(data[[treatment_var]])
  } else {
    df_work$treat <- rbinom(nrow(data), size = 1, prob = 0.5)
    if (verbose) message("Treatment variable simulated (50/50 randomization)")
  }

  # ==========================================================================
  # Process Covariates
  # ==========================================================================

  # Process continuous variables (standardize)
  for (var in continuous_vars) {
    df_work[[paste0("z_", var)]] <- as.numeric(scale(data[[var]]))
  }

  # Process factor variables
  for (var in factor_vars) {
    if (is.factor(data[[var]])) {
      # Create dummy variables for factors
      dummies <- model.matrix(~ data[[var]] - 1)
      # Keep all but first level (reference)
      if (ncol(dummies) > 1) {
        for (j in 2:ncol(dummies)) {
          dummy_name <- paste0("z_", var, "_", j - 1)
          df_work[[dummy_name]] <- dummies[, j]
        }
      }
    } else {
      # Treat as binary
      df_work[[paste0("z_", var)]] <- as.numeric(data[[var]])
    }
  }

  # ==========================================================================
  # Define Subgroups (if specified)
  # ==========================================================================

  df_work$flag_harm <- 0L
  interaction_term <- NULL

  if (model == "alt" && !is.null(subgroup_vars)) {
    # Create subgroup indicators
    subgroup_indicators <- list()

    # FIXED: Use seq_along instead of enumerate (which doesn't exist in R)
    for (var in subgroup_vars) {
      if (var %in% continuous_vars) {
        # Use median or provided cutpoint
        cutpoint <- if (!is.null(subgroup_cuts[[var]])) {
          subgroup_cuts[[var]]
        } else {
          median(data[[var]], na.rm = TRUE)
        }
        subgroup_indicators[[var]] <- data[[var]] <= cutpoint
      } else {
        # For factor variables, use first level or provided specification
        if (!is.null(subgroup_cuts[[var]])) {
          subgroup_indicators[[var]] <- data[[var]] == subgroup_cuts[[var]]
        } else {
          subgroup_indicators[[var]] <- data[[var]] == levels(as.factor(data[[var]]))[1]
        }
      }
    }

    # Create harm flag (all subgroup conditions met)
    df_work$flag_harm <- as.integer(Reduce(`&`, subgroup_indicators))

    # Create interaction term
    if (length(subgroup_indicators) > 0) {
      interaction_term <- df_work$treat * df_work$flag_harm
    }

    if (verbose) {
      message("Subgroup defined by: ", paste(subgroup_vars, collapse = " AND "))
      message("Subgroup size: ", sum(df_work$flag_harm), " out of ", nrow(df_work))
      message("Subgroup proportion: ", round(mean(df_work$flag_harm), 3))
    }
  }

  # ==========================================================================
  # Fit AFT Model (Weibull)
  # ==========================================================================

  # Check for survival package
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: install.packages('survival')")
  }

  # Prepare model matrix
  covariate_cols <- grep("^z_", names(df_work), value = TRUE)
  X <- as.matrix(df_work[, c("treat", covariate_cols), drop = FALSE])

  # Add interaction term if needed
  if (!is.null(interaction_term)) {
    X <- cbind(X, treat_harm = interaction_term)
  }

  # Fit Weibull AFT model
  formula_str <- paste("survival::Surv(y, event) ~ ",
                       paste(c("treat", covariate_cols), collapse = " + "))
  if (!is.null(interaction_term)) {
    df_work$treat_harm <- interaction_term
    formula_str <- paste(formula_str, "+ treat_harm")
  }

  fit_aft <- survival::survreg(as.formula(formula_str),
                               data = df_work,
                               dist = "weibull")

  # Extract parameters
  mu <- coef(fit_aft)[1]  # Intercept
  sigma <- fit_aft$scale
  gamma <- coef(fit_aft)[-1]  # Coefficients (excluding intercept)

  # Apply effect modifiers
  gamma["treat"] <- k_treat * gamma["treat"]
  if ("treat_harm" %in% names(gamma)) {
    gamma["treat_harm"] <- k_inter * gamma["treat_harm"]
  }

  # Weibull parameterization
  b_weibull <- gamma
  b0_weibull <- -gamma / sigma

  if (verbose) {
    message("\n=== Model Parameters ===")
    message("Intercept (mu): ", round(mu, 4))
    message("Scale (sigma): ", round(sigma, 4))
    message("Treatment effect: ", round(gamma["treat"], 4))
    if ("treat_harm" %in% names(gamma)) {
      message("Interaction effect: ", round(gamma["treat_harm"], 4))
    }
  }

  # ==========================================================================
  # Generate Super Population
  # ==========================================================================

  # Sample with replacement to create super population
  n_treat <- round(n_super / 2)
  n_control <- n_super - n_treat

  # Sample indices
  idx_sample <- sample(seq_len(nrow(df_work)), size = n_super, replace = TRUE)
  df_super <- df_work[idx_sample, , drop = FALSE]

  # Assign treatment
  df_super$treat[seq_len(n_treat)] <- 1L
  df_super$treat[(n_treat + 1):n_super] <- 0L

  # Reset IDs
  df_super$id <- seq_len(n_super)

  # Calculate linear predictors for potential outcomes
  X_super <- as.matrix(df_super[, c("treat", covariate_cols), drop = FALSE])
  if (!is.null(interaction_term)) {
    # Recalculate interaction for super population
    df_super$treat_harm <- df_super$treat * df_super$flag_harm
    X_super <- cbind(X_super, treat_harm = df_super$treat_harm)
  }

  # Potential outcomes under treatment
  X_treat <- X_super
  X_treat[, "treat"] <- 1
  if ("treat_harm" %in% colnames(X_treat)) {
    X_treat[, "treat_harm"] <- df_super$flag_harm
  }

  # Potential outcomes under control
  X_control <- X_super
  X_control[, "treat"] <- 0
  if ("treat_harm" %in% colnames(X_control)) {
    X_control[, "treat_harm"] <- 0
  }

  # Linear predictors
  df_super$lin_pred_1 <- as.numeric(X_treat %*% b_weibull)
  df_super$lin_pred_0 <- as.numeric(X_control %*% b_weibull)
  df_super$lin_pred_obs <- as.numeric(X_super %*% b_weibull)

  # Hazard ratios (individual level)
  df_super$hr_individual <- exp((df_super$lin_pred_1 - df_super$lin_pred_0) / sigma)

  # ==========================================================================
  # Calculate True Hazard Ratios
  # ==========================================================================

  # Generate potential outcomes for HR calculation
  epsilon <- log(rexp(n_super))  # Extreme value distribution

  # Under treatment
  logT_1 <- mu + sigma * epsilon + df_super$lin_pred_1
  T_1 <- exp(logT_1)

  # Under control
  logT_0 <- mu + sigma * epsilon + df_super$lin_pred_0
  T_0 <- exp(logT_0)

  # Calculate empirical hazard ratios
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
    # Fit censoring model
    X_cens <- as.matrix(df_work[, c("treat", covariate_cols), drop = FALSE])

    fit_cens <- survival::survreg(
      survival::Surv(y, 1 - event) ~ X_cens,
      data = df_work,
      dist = "weibull"
    )

    mu_cens <- coef(fit_cens)[1]
    sigma_cens <- fit_cens$scale
    gamma_cens <- coef(fit_cens)[-1]

    # Store censoring parameters
    cens_model <- list(
      mu = mu_cens,
      sigma = sigma_cens,
      gamma = gamma_cens,
      type = "weibull"
    )

    # Calculate censoring linear predictors for super population
    df_super$lin_pred_cens_1 <- as.numeric(X_treat[, seq_len(ncol(X_cens)), drop = FALSE] %*% gamma_cens)
    df_super$lin_pred_cens_0 <- as.numeric(X_control[, seq_len(ncol(X_cens)), drop = FALSE] %*% gamma_cens)

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

  # ==========================================================================
  # Prepare Output
  # ==========================================================================

  # Model parameters
  model_params <- list(
    mu = mu,
    sigma = sigma,
    gamma = gamma,
    b_weibull = b_weibull,
    b0_weibull = b0_weibull,
    censoring = cens_model
  )

  # Subgroup information
  subgroup_info <- list(
    vars = subgroup_vars,
    cuts = subgroup_cuts,
    size = sum(df_super$flag_harm),
    proportion = mean(df_super$flag_harm)
  )

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

  class(results) <- c("aft_dgm", "list")

  return(results)
}


#' Simulate Data from AFT DGM
#'
#' Generates a sample from a previously created AFT data generating mechanism.
#'
#' @param dgm Object of class "aft_dgm" returned by \code{\link{generate_aft_dgm}}
#' @param n Integer. Sample size. If NULL, uses the super population.
#' @param rand_ratio Numeric. Randomization ratio (treatment:control). Default is 1.
#' @param max_follow Numeric. Maximum follow-up time for administrative censoring.
#'   Default is Inf (no administrative censoring).
#' @param cens_adjust Numeric. Adjustment to censoring rate (on log scale).
#'   Positive values increase censoring. Default is 0.
#' @param seed Integer. Random seed. If NULL, no seed is set.
#'
#' @return A data.frame with simulated survival data including:
#' \describe{
#'   \item{id}{Subject identifier}
#'   \item{treat}{Treatment indicator (0/1)}
#'   \item{flag_harm}{Subgroup membership indicator}
#'   \item{y_sim}{Observed time (min of event and censoring)}
#'   \item{event_sim}{Event indicator (1 = event, 0 = censored)}
#'   \item{t_true}{True event time}
#' }
#'
#' @examples
#' \donttest{
#' # First create a DGM
#' example_data <- data.frame(
#'   time = rexp(100, 0.02),
#'   status = rbinom(100, 1, 0.7),
#'   age = rnorm(100, 55, 12),
#'   treat = rbinom(100, 1, 0.5)
#' )
#'
#' dgm <- generate_aft_dgm(
#'   data = example_data,
#'   continuous_vars = "age",
#'   factor_vars = character(0),
#'   outcome_var = "time",
#'   event_var = "status",
#'   treatment_var = "treat",
#'   model = "null",
#'   verbose = FALSE
#' )
#'
#' # Simulate a sample
#' sim_data <- simulate_from_dgm(dgm, n = 500, seed = 123)
#' head(sim_data)
#' }
#'
#' @seealso \code{\link{generate_aft_dgm}}
#'
#' @importFrom stats rexp runif
#' @export
simulate_from_dgm <- function(dgm,
                              n = NULL,
                              rand_ratio = 1,
                              max_follow = Inf,
                              cens_adjust = 0,
                              seed = NULL) {

  if (!inherits(dgm, "aft_dgm")) {
    stop("dgm must be an object created by generate_aft_dgm()")
  }

  if (!is.null(seed)) set.seed(seed)

  df_super <- dgm$df_super
  params <- dgm$model_params

  # Determine sample size
  if (is.null(n)) {
    df_sim <- df_super
    n <- nrow(df_sim)
  } else {
    # Sample from super population
    n_treat <- round(n * rand_ratio / (1 + rand_ratio))
    n_control <- n - n_treat

    idx_sample <- sample(seq_len(nrow(df_super)), size = n, replace = TRUE)
    df_sim <- df_super[idx_sample, , drop = FALSE]

    # Reassign treatment
    df_sim$treat[seq_len(n_treat)] <- 1L
    df_sim$treat[(n_treat + 1):n] <- 0L

    # Update linear predictors based on assigned treatment
    df_sim$lin_pred_obs <- ifelse(df_sim$treat == 1,
                                  df_sim$lin_pred_1,
                                  df_sim$lin_pred_0)

    # Reset IDs
    df_sim$id <- seq_len(n)
  }

  # Generate survival times
  epsilon <- log(rexp(n))  # Extreme value distribution
  logT_sim <- params$mu + params$sigma * epsilon + df_sim$lin_pred_obs
  T_sim <- exp(logT_sim)

  # Generate censoring times
  if (params$censoring$type == "weibull") {
    # Weibull censoring
    lin_pred_cens <- ifelse(df_sim$treat == 1,
                            df_sim$lin_pred_cens_1,
                            df_sim$lin_pred_cens_0)

    epsilon_cens <- log(rexp(n))
    logC_sim <- params$censoring$mu + cens_adjust +
      params$censoring$sigma * epsilon_cens + lin_pred_cens
    C_sim <- exp(logC_sim)

  } else if (params$censoring$type == "uniform") {
    # Uniform censoring
    C_sim <- runif(n,
                   min = params$censoring$min,
                   max = params$censoring$max)
  }

  # Apply administrative censoring
  C_sim <- pmin(C_sim, max_follow)

  # Observed times and events
  df_sim$y_sim <- pmin(T_sim, C_sim)
  df_sim$event_sim <- as.integer(T_sim <= C_sim)
  df_sim$t_true <- T_sim
  df_sim$c_time <- C_sim

  # Add analysis variables
  analysis_cols <- c("id", "treat", "flag_harm",
                     dgm$analysis_vars$covariates,
                     "y_sim", "event_sim", "t_true")

  # Keep only necessary columns
  df_sim <- df_sim[, intersect(analysis_cols, names(df_sim)), drop = FALSE]
  rownames(df_sim) <- NULL

  return(df_sim)
}


#' Print Method for aft_dgm Objects
#'
#' @param x An aft_dgm object
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns the input object
#'
#' @export
print.aft_dgm <- function(x, ...) {
  cat("AFT Data Generating Mechanism\n")
  cat("=============================\n")
  cat("Model type:", x$model_type, "\n")
  cat("Super population size:", x$n_super, "\n")
  cat("Number of covariates:", length(x$analysis_vars$covariates), "\n")

  if (x$model_type == "alt" && !is.null(x$subgroup_info$vars)) {
    cat("\nSubgroup Information:\n")
    cat("  Variables:", paste(x$subgroup_info$vars, collapse = ", "), "\n")
    cat("  Size:", x$subgroup_info$size, "\n")
    cat("  Proportion:", round(x$subgroup_info$proportion, 3), "\n")
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


#' Summary of Simulated Data
#'
#' Prints summary statistics for simulated survival data.
#'
#' @param data Simulated dataset from \code{\link{simulate_from_dgm}}
#' @param dgm Optional. Original DGM object for comparison of true vs observed HRs
#'
#' @return Invisibly returns NULL
#'
#' @examples
#' \donttest{
#' # See simulate_from_dgm examples
#' }
#'
#' @importFrom stats median
#' @importFrom survival coxph Surv
#' @export
summarize_simulation <- function(data, dgm = NULL) {
  cat("Simulated Data Summary\n")
  cat("======================\n")
  cat("Sample size:", nrow(data), "\n")
  cat("Treatment allocation:", paste(names(table(data$treat)),
                                      table(data$treat), sep = ":", collapse = ", "), "\n")
  cat("Event rate:", round(mean(data$event_sim), 3), "\n")
  cat("Median follow-up:", round(median(data$y_sim), 2), "\n")

  if ("flag_harm" %in% names(data)) {
    cat("\nSubgroup sizes:\n")
    cat("  Harm subgroup:", sum(data$flag_harm), "\n")
    cat("  No-harm subgroup:", sum(1 - data$flag_harm), "\n")
  }

  # Calculate observed hazard ratio
  hr_obs <- exp(survival::coxph(
    survival::Surv(y_sim, event_sim) ~ treat,
    data = data
  )$coefficients)
  cat("\nObserved HR (overall):", round(hr_obs, 3), "\n")

  if ("flag_harm" %in% names(data) && sum(data$flag_harm) > 0) {
    hr_harm_obs <- exp(survival::coxph(
      survival::Surv(y_sim, event_sim) ~ treat,
      data = data[data$flag_harm == 1, ]
    )$coefficients)

    hr_no_harm_obs <- exp(survival::coxph(
      survival::Surv(y_sim, event_sim) ~ treat,
      data = data[data$flag_harm == 0, ]
    )$coefficients)

    cat("Observed HR (harm subgroup):", round(hr_harm_obs, 3), "\n")
    cat("Observed HR (no-harm subgroup):", round(hr_no_harm_obs, 3), "\n")
  }

  if (!is.null(dgm)) {
    cat("\nTrue hazard ratios from DGM:\n")
    for (name in names(dgm$hazard_ratios)) {
      cat("  ", name, ": ", round(dgm$hazard_ratios[[name]], 3), "\n", sep = "")
    }
  }

  invisible(NULL)
}
