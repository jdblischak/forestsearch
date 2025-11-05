# Simulation function for AFT DGM with flexible subgroups and spline models
# Compatible with both standard (alt/null) and spline model types

#' Simulate Data from AFT DGM
#'
#' Generates simulated survival data from a fitted AFT data generating mechanism,
#' supporting both standard subgroup models and spline-based models.
#'
#' @param dgm An object of class "aft_dgm_flex" created by generate_aft_dgm_flex()
#' @param n_sim Integer specifying the number of observations to simulate.
#'   If NULL, uses the size of the super population
#' @param draw_treatment Logical indicating whether to randomly draw treatment
#'   assignments. If FALSE and original treatment exists, uses those assignments
#' @param treatment_prob Numeric probability of treatment assignment when drawing
#'   new treatments. Default is 0.5
#' @param seed Integer random seed for reproducibility. If NULL, no seed is set
#' @param verbose Logical indicating whether to print progress information
#'
#' @return A data.frame containing:
#' \describe{
#'   \item{id}{Observation ID}
#'   \item{[covariates]}{All covariate columns from the super population}
#'   \item{treat}{Treatment assignment (0 or 1)}
#'   \item{time}{Observed survival time}
#'   \item{event}{Event indicator (1 = event, 0 = censored)}
#'   \item{time_uncensored}{True survival time (before censoring)}
#'   \item{time_censor}{Censoring time}
#'   \item{in_subgroup}{Subgroup membership indicator (for standard models)}
#'   \item{z}{Biomarker value (for spline models)}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate DGM
#' dgm <- generate_aft_dgm_flex(
#'   data = gbsg,
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   model = "alt"
#' )
#'
#' # Simulate data
#' sim_data <- simulate_from_dgm(dgm, n_sim = 1000, seed = 123)
#'
#' # Check censoring rate
#' mean(sim_data$event)
#' }
#'
#' @export

simulate_from_dgm_enhanced <- function(dgm,
                             n_sim = NULL,
                             draw_treatment = TRUE,
                             treatment_prob = 0.5,
                             seed = NULL,
                             verbose = FALSE) {

  # Check input
  if (!inherits(dgm, "aft_dgm_flex")) {
    stop("dgm must be an object of class 'aft_dgm_flex' created by generate_aft_dgm_flex()")
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Get super population
  df_super <- dgm$df_super
  n_super <- nrow(df_super)

  # Determine sample size
  if (is.null(n_sim)) {
    n_sim <- n_super
    if (verbose) cat("Using full super population size:", n_sim, "\n")
  }

  # Sample from super population with replacement
  if (n_sim != n_super) {
    sample_ids <- sample(1:n_super, size = n_sim, replace = TRUE)
    df_sim <- df_super[sample_ids, ]
  } else {
    df_sim <- df_super
  }

  # Reset row names
  rownames(df_sim) <- NULL
  df_sim$id <- 1:n_sim

  # Get model parameters
  mu <- dgm$model_params$mu
  sigma <- dgm$model_params$sigma
  gamma <- dgm$model_params$gamma

  if (verbose) {
    cat("\nSimulating", n_sim, "observations from", dgm$model_type, "model\n")
  }

  # ============================================================================
  # Handle treatment assignment
  # ============================================================================

  if (draw_treatment) {
    # Draw new treatment assignments
    df_sim$treat <- rbinom(n_sim, 1, treatment_prob)
    if (verbose) {
      cat("Drew new treatment assignments with probability", treatment_prob, "\n")
      cat("Treatment allocation: ", mean(df_sim$treat), "\n")
    }
  } else {
    # Use original treatment if available
    if ("treat_orig" %in% names(df_sim)) {
      df_sim$treat <- df_sim$treat_orig
      if (verbose) cat("Using original treatment assignments\n")
    } else if (!is.null(dgm$analysis_vars$treatment) &&
               dgm$analysis_vars$treatment %in% names(df_sim)) {
      df_sim$treat <- df_sim[[dgm$analysis_vars$treatment]]
      if (verbose) cat("Using treatment variable:", dgm$analysis_vars$treatment, "\n")
    } else {
      # Default to random assignment
      df_sim$treat <- rbinom(n_sim, 1, treatment_prob)
      if (verbose) cat("No original treatment found, using random assignment\n")
    }
  }

  # ============================================================================
  # Generate survival times based on model type
  # ============================================================================

  if (dgm$model_type == "spline") {
    # SPLINE MODEL
    if (verbose) cat("Generating survival times using spline model\n")

    # Extract spline-specific terms
    if (!"z" %in% names(df_sim)) {
      stop("Biomarker variable 'z' not found in simulated data")
    }

    z <- df_sim$z
    knot <- dgm$spline_info$knot

    # Recreate spline terms
    df_sim$z.treat <- z * df_sim$treat
    df_sim$z.k <- (z - knot) * ifelse(z > knot, 1, 0)
    df_sim$z.k.treat <- df_sim$z.k * df_sim$treat

    # Build design matrix for spline model
    X_spline <- cbind(
      treat = df_sim$treat,
      z = z,
      z.treat = df_sim$z.treat,
      z.k = df_sim$z.k,
      z.k.treat = df_sim$z.k.treat
    )

    # Calculate linear predictor
    # Check if gamma has names or is just a vector
    if (!is.null(names(gamma))) {
      # Extract relevant coefficients
      gamma_spline <- numeric(5)
      coef_names <- c("treat", "z", "z.treat", "z.k", "z.k.treat")
      for (i in 1:5) {
        idx <- which(names(gamma) == coef_names[i])
        if (length(idx) > 0) {
          gamma_spline[i] <- gamma[idx]
        } else {
          gamma_spline[i] <- 0
        }
      }
      lp <- mu + as.vector(X_spline %*% gamma_spline)
    } else if (length(gamma) == 5) {
      # Assume gamma is in the correct order for spline model
      lp <- mu + as.vector(X_spline %*% gamma)
    } else {
      stop("Gamma coefficients not properly structured for spline model")
    }

  } else {
    # STANDARD MODEL (alt or null)
    if (verbose) cat("Generating survival times using", dgm$model_type, "model\n")

    # Build linear predictor using pre-calculated values if available
    if ("lp_treat_1" %in% names(df_sim) || "lp_treat_1_in" %in% names(df_sim)) {
      # Use pre-calculated linear predictors
      if (dgm$model_type == "alt" && "lp_treat_1_in" %in% names(df_sim)) {
        # For alt model with subgroups
        lp <- ifelse(df_sim$treat == 1,
                    ifelse(df_sim$in_subgroup, df_sim$lp_treat_1_in, df_sim$lp_treat_1_out),
                    df_sim$lp_treat_0)
      } else {
        # For null model or when no subgroup
        lp <- ifelse(df_sim$treat == 1, df_sim$lp_treat_1, df_sim$lp_treat_0)
      }
    } else {
      # Fallback: calculate linear predictor from scratch
      lp <- mu

      # Add treatment effect
      treatment_var_name <- dgm$analysis_vars$treatment
      if (is.null(treatment_var_name)) treatment_var_name <- "treat_sim"

      for (i in seq_along(gamma)) {
        coef_name <- names(gamma)[i]

        if (coef_name == treatment_var_name || coef_name == "treat_sim") {
          lp <- lp + df_sim$treat * gamma[i]
        } else if (coef_name == "treat_x_subgroup" && dgm$model_type == "alt") {
          if ("in_subgroup" %in% names(df_sim)) {
            lp <- lp + df_sim$treat * df_sim$in_subgroup * gamma[i]
          }
        } else if (grepl("_std$", coef_name)) {
          # Standardized continuous variable
          var_base <- sub("_std$", "", coef_name)
          if (var_base %in% names(df_sim)) {
            var_vals <- df_sim[[var_base]]
            if (dgm$analysis_vars$standardize) {
              var_vals <- scale(var_vals)[,1]
            }
            lp <- lp + var_vals * gamma[i]
          }
        } else if (coef_name %in% names(df_sim)) {
          # Direct variable (e.g., factor dummies)
          lp <- lp + df_sim[[coef_name]] * gamma[i]
        }
      }
    }
  }

  # Generate error term (extreme value distribution)
  epsilon <- log(rexp(n_sim))

  # Generate log survival times
  log_time <- lp + sigma * epsilon

  # Transform to survival times
  df_sim$time_uncensored <- exp(log_time)

  # ============================================================================
  # Generate censoring times
  # ============================================================================

  cens_params <- dgm$model_params$censoring

  if (cens_params$type == "weibull") {
    # Weibull censoring
    epsilon_c <- log(rexp(n_sim))
    log_censor <- cens_params$mu + cens_params$sigma * epsilon_c
    df_sim$time_censor <- exp(log_censor)

  } else if (cens_params$type == "uniform") {
    # Uniform censoring
    df_sim$time_censor <- runif(n_sim,
                                min = cens_params$min,
                                max = cens_params$max)
  } else {
    stop("Unknown censoring type:", cens_params$type)
  }

  # Apply censoring
  df_sim$time <- pmin(df_sim$time_uncensored, df_sim$time_censor)
  df_sim$event <- as.numeric(df_sim$time_uncensored <= df_sim$time_censor)

  # ============================================================================
  # Calculate additional quantities for analysis
  # ============================================================================

  # Censoring rate
  cens_rate <- 1 - mean(df_sim$event)

  if (verbose) {
    cat("\n--- Simulation Summary ---\n")
    cat("Sample size:", n_sim, "\n")
    cat("Treatment allocation:", round(mean(df_sim$treat), 3), "\n")
    cat("Event rate:", round(mean(df_sim$event), 3), "\n")
    cat("Censoring rate:", round(cens_rate, 3), "\n")

    if (dgm$model_type == "alt" && "in_subgroup" %in% names(df_sim)) {
      cat("Subgroup proportion:", round(mean(df_sim$in_subgroup), 3), "\n")
    }

    if (dgm$model_type == "spline") {
      cat("Biomarker range: [", round(min(df_sim$z), 2), ",",
          round(max(df_sim$z), 2), "]\n")
    }

    cat("\nMedian survival time:", round(median(df_sim$time), 2), "\n")
    cat("IQR:", round(quantile(df_sim$time, 0.25), 2), "-",
        round(quantile(df_sim$time, 0.75), 2), "\n")
  }

  # Select columns to return
  return_cols <- c("id", dgm$analysis_vars$continuous, dgm$analysis_vars$factor,
                   "treat", "time", "event", "time_uncensored", "time_censor")

  if (dgm$model_type == "alt" && "in_subgroup" %in% names(df_sim)) {
    return_cols <- c(return_cols, "in_subgroup")
  }

  if (dgm$model_type == "spline") {
    return_cols <- c(return_cols, "z", "z.k")
  }

  # Return only relevant columns (ensure data frame is returned)
  valid_cols <- intersect(return_cols, names(df_sim))
  df_return <- df_sim[, valid_cols, drop = FALSE]

  return(df_return)
}

# ============================================================================
# Analyze simulated data function
# ============================================================================

#' Analyze Simulated Data from AFT DGM
#'
#' Performs standard survival analysis on simulated data including Cox regression
#' and summary statistics.
#'
#' @param sim_data Data frame from simulate_from_dgm()
#' @param dgm Optional: the original DGM object for comparison with true values
#' @param verbose Logical indicating whether to print results
#'
#' @return List containing analysis results
#'
#' @export

analyze_simulation <- function(sim_data, dgm = NULL, verbose = TRUE) {

  # Load required package
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for this function")
  }

  # Basic Cox model
  cox_fit <- survival::coxph(survival::Surv(time, event) ~ treat, data = sim_data)
  cox_summary <- summary(cox_fit)

  # Extract HR and CI
  hr_est <- exp(coef(cox_fit))
  hr_ci <- exp(confint(cox_fit))

  results <- list(
    n = nrow(sim_data),
    event_rate = mean(sim_data$event),
    treatment_rate = mean(sim_data$treat),
    hr_estimate = hr_est,
    hr_ci_lower = hr_ci[1],
    hr_ci_upper = hr_ci[2],
    p_value = cox_summary$coefficients[, "Pr(>|z|)"],
    cox_fit = cox_fit
  )

  # If DGM provided, compare with true values
  if (!is.null(dgm)) {
    results$true_hr <- dgm$hazard_ratios$overall
    results$bias <- results$hr_estimate - results$true_hr
    results$relative_bias <- (results$hr_estimate - results$true_hr) / results$true_hr
    results$covers_truth <- (results$true_hr >= results$hr_ci_lower) &
                           (results$true_hr <= results$hr_ci_upper)
  }

  if (verbose) {
    cat("\n=== Analysis Results ===\n")
    cat("Sample size:", results$n, "\n")
    cat("Event rate:", round(results$event_rate, 3), "\n")
    cat("Treatment rate:", round(results$treatment_rate, 3), "\n\n")

    cat("Estimated HR:", round(results$hr_estimate, 3), "\n")
    cat("95% CI: [", round(results$hr_ci_lower, 3), ",",
        round(results$hr_ci_upper, 3), "]\n")
    cat("P-value:", format.pval(results$p_value), "\n")

    if (!is.null(dgm)) {
      cat("\nTrue HR:", round(results$true_hr, 3), "\n")
      cat("Bias:", round(results$bias, 3), "\n")
      cat("Relative bias:", round(results$relative_bias * 100, 1), "%\n")
      cat("CI covers truth:", results$covers_truth, "\n")
    }
  }

  return(results)
}
