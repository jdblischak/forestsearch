# =============================================================================
# Cox Model Estimation Helpers - Optimized
# =============================================================================
#
# Optimized functions for Cox model fitting in ForestSearch package.
# These functions are called thousands of times during bootstrap and
# consistency evaluation, so performance is critical.
#
# Key optimizations:
#   1. Minimal memory allocation (model=FALSE, x=FALSE, y=FALSE)
#   2. Direct coefficient extraction (avoid repeated summary() calls)
#   3. Early exit for edge cases
#   4. Vectorized batch processing
#
# =============================================================================

#' Format Hazard Ratio with Confidence Interval
#'
#' @param hr Numeric vector of length 3: c(HR, lower, upper)
#' @param digits Integer number of decimal places (default: 2)
#'
#' @return Character string formatted as "HR (lower, upper)"
#' @keywords internal
hrCI_format <- function(hr, digits = 2L) {
  sprintf("%.*f (%.*f, %.*f)", digits, hr[1], digits, hr[2], digits, hr[3])
}


#' Cox Model Summary for Subgroup (Optimized)
#'
#' Calculates hazard ratio and confidence interval for a subgroup using
#' Cox regression. Optimized for performance in bootstrap iterations.
#'
#' @param Y Numeric vector of outcome times.
#' @param E Numeric vector of event indicators (0/1).
#' @param Treat Numeric vector of treatment indicators (0/1).
#' @param Strata Vector of strata identifiers (optional).
#' @param use_strata Logical. Whether to include strata in model.
#'   Default: TRUE if Strata provided.
#' @param return_format Character. "formatted" for string output or
#'   "numeric" for named vector.
#' @param conf_level Numeric. Confidence level for intervals (default: 0.95).
#'
#' @return If return_format="formatted": Character string "HR (lower, upper)".
#'   If return_format="numeric": Named numeric vector c(HR, Lower, Upper).
#'
#' @importFrom survival coxph Surv
#' @export
cox_summary <- function(Y,
                        E,
                        Treat,
                        Strata = NULL,
                        use_strata = !is.null(Strata),
                        return_format = c("formatted", "numeric"),
                        conf_level = 0.95) {

  return_format <- match.arg(return_format)

  # ---------------------------------------------------------------------------
  # Input Validation (Fail Fast)
  # ---------------------------------------------------------------------------
  n <- length(Y)

  # Check vector lengths match
  if (length(E) != n || length(Treat) != n) {
    stop("Y, E, and Treat must have the same length")
  }

  if (use_strata && !is.null(Strata) && length(Strata) != n) {
    stop("Strata must have the same length as Y")
  }

  # Check for sufficient events (need at least 2 for Cox model)
  n_events <- sum(E)
  if (n_events < 2L) {
    warning("Fewer than 2 events; returning NA")
    return(.return_na(return_format))
  }

  # Check treatment variation
  if (length(unique(Treat)) < 2L) {
    warning("No variation in treatment; returning NA")
    return(.return_na(return_format))
  }

  # ---------------------------------------------------------------------------
  # Fit Cox Model (Memory Optimized)
  # ---------------------------------------------------------------------------
  fit <- tryCatch({
    if (use_strata && !is.null(Strata)) {
      survival::coxph(
        survival::Surv(Y, E) ~ Treat + strata(Strata),
        robust = TRUE,
        model = FALSE,
        x = FALSE,
        y = FALSE
      )
    } else {
      survival::coxph(
        survival::Surv(Y, E) ~ Treat,
        robust = TRUE,
        model = FALSE,
        x = FALSE,
        y = FALSE
      )
    }
  }, error = function(e) {
    warning("Cox model failed: ", e$message)
    NULL
  })

  if (is.null(fit)) {
    return(.return_na(return_format))
  }

  # ---------------------------------------------------------------------------
  # Extract Results (Single summary() Call)
  # ---------------------------------------------------------------------------
  fit_summary <- summary(fit, conf.int = conf_level)
  conf_int <- fit_summary$conf.int

  # Handle edge case where conf.int might not exist
  if (is.null(conf_int) || nrow(conf_int) == 0L) {
    warning("No confidence interval available from Cox model")
    return(.return_na(return_format))
  }

  # Extract HR and CI bounds (columns 1, 3, 4 of conf.int)
  hr <- conf_int[1L, 1L]
  lower <- conf_int[1L, 3L]
  upper <- conf_int[1L, 4L]

  # ---------------------------------------------------------------------------
  # Return Result
  # ---------------------------------------------------------------------------
  if (return_format == "formatted") {
    hrCI_format(c(hr, lower, upper))
  } else {
    c(HR = hr, Lower = lower, Upper = upper)
  }
}


#' Return NA Result in Appropriate Format
#'
#' @param return_format Character. "formatted" or "numeric"
#' @return NA in the requested format
#' @keywords internal
.return_na <- function(return_format) {
  if (return_format == "formatted") {
    "NA (NA, NA)"
  } else {
    c(HR = NA_real_, Lower = NA_real_, Upper = NA_real_)
  }
}


#' Batch Cox Model Summary for Multiple Subgroups
#'
#' Efficiently processes multiple subgroups in a single call.
#' Pre-validates data once to avoid repeated checks.
#'
#' @param Y Numeric vector of outcome (full dataset).
#' @param E Numeric vector of event indicators (full dataset).
#' @param Treat Numeric vector of treatment indicators (full dataset).
#' @param Strata Vector of strata (optional, full dataset).
#' @param subset_indices List of integer vectors, each defining subset indices.
#' @param return_format Character. "formatted" or "numeric".
#'
#' @return List of results, one per subset.
#'
#' @importFrom survival coxph Surv
#' @export
cox_summary_batch <- function(Y,
                              E,
                              Treat,
                              Strata = NULL,
                              subset_indices,
                              return_format = c("formatted", "numeric")) {

  return_format <- match.arg(return_format)

  # Pre-validate full dataset
  n <- length(Y)
  if (length(E) != n || length(Treat) != n) {
    stop("Y, E, and Treat must have the same length")
  }

  use_strata <- !is.null(Strata)
  if (use_strata && length(Strata) != n) {
    stop("Strata must have the same length as Y")
  }

  # Process each subset
  n_subsets <- length(subset_indices)
  results <- vector("list", n_subsets)

  for (i in seq_len(n_subsets)) {
    idx <- subset_indices[[i]]

    # Extract subset
    Y_sub <- Y[idx]
    E_sub <- E[idx]
    Treat_sub <- Treat[idx]
    Strata_sub <- if (use_strata) Strata[idx] else NULL

    # Call optimized function
    results[[i]] <- cox_summary(
      Y_sub, E_sub, Treat_sub, Strata_sub,
      use_strata = use_strata,
      return_format = return_format
    )
  }

  results
}


#' Fit Cox Model and Extract Estimate with SE
#'
#' Fits a Cox model and returns the log hazard ratio and robust standard error.
#' Used for bias correction in bootstrap procedures.
#'
#' @param df_sg Data frame for subgroup.
#' @param cox_formula Cox model formula.
#' @param est_loghr Logical. Return estimate on log(HR) scale? Default TRUE.
#' @param init Numeric. Initial value for coefficient. Default log(1) = 0.
#'
#' @return List with components:
#'   \item{est_obs}{Point estimate (log HR or HR)}
#'   \item{se_obs}{Standard error (robust)}
#'
#' @importFrom survival coxph
#' @export
get_Cox_sg <- function(df_sg, cox_formula, est_loghr = TRUE, init = 0) {

  # Validate inputs
  names_tocheck <- all.vars(cox_formula)
  missing_vars <- setdiff(names_tocheck, names(df_sg))

  if (length(missing_vars) > 0L) {
    stop("df_sg missing required variables: ", paste(missing_vars, collapse = ", "))
  }

  # Fit model with robust SE
  fit <- suppressWarnings(
    survival::coxph(
      cox_formula,
      data = df_sg,
      model = FALSE,
      x = FALSE,
      y = FALSE,
      robust = TRUE,
      init = init
    )
  )

  # Extract coefficients (single summary call)
  fit_sum <- summary(fit)
  coef_matrix <- fit_sum$coefficients

  # Extract log hazard ratio and robust SE
  bhat <- coef_matrix[1L, "coef"]
  se_robust <- coef_matrix[1L, "robust se"]

  # Return on requested scale
  if (est_loghr) {
    list(est_obs = bhat, se_obs = se_robust)
  } else {
    list(est_obs = exp(bhat), se_obs = exp(bhat) * se_robust)
  }
}


#' Build Cox Model Formula
#'
#' Constructs a Cox proportional hazards formula from variable names.
#'
#' @param outcome_name Character. Name of outcome/time variable.
#' @param event_name Character. Name of event indicator variable.
#' @param treat_name Character. Name of treatment variable.
#' @param covariates Character vector. Optional additional covariates.
#' @param strata_name Character. Optional stratification variable.
#'
#' @return An R formula object for Cox regression.
#'
#' @examples
#' # Simple formula
#' build_cox_formula("time", "status", "treatment")
#' # Surv(time, status) ~ treatment
#'
#' # With covariates
#' build_cox_formula("time", "status", "treatment", c("age", "stage"))
#' # Surv(time, status) ~ treatment + age + stage
#'
#' # With stratification
#' build_cox_formula("time", "status", "treatment", strata_name = "center")
#' # Surv(time, status) ~ treatment + strata(center)
#'
#' @export
build_cox_formula <- function(outcome_name,
                              event_name,
                              treat_name,
                              covariates = NULL,
                              strata_name = NULL) {

  # Build LHS
  lhs <- sprintf("Surv(%s, %s)", outcome_name, event_name)

  # Build RHS
  rhs_terms <- treat_name

  if (!is.null(covariates) && length(covariates) > 0L) {
    rhs_terms <- c(rhs_terms, covariates)
  }

  if (!is.null(strata_name)) {
    rhs_terms <- c(rhs_terms, sprintf("strata(%s)", strata_name))
  }

  rhs <- paste(rhs_terms, collapse = " + ")

  # Combine and return formula
  as.formula(paste(lhs, "~", rhs))
}


#' Fit Cox Models for Both Subgroups
#'
#' Fits Cox models for harm and no-harm subgroups defined by treat.recommend.
#'
#' @param df Data frame with treat.recommend column (0 = harm, 1 = no harm).
#' @param formula Cox model formula.
#'
#' @return List with components:
#'   \item{H_obs}{Log HR estimate for harm subgroup}
#'   \item{seH_obs}{SE for harm subgroup}
#'   \item{Hc_obs}{Log HR estimate for no-harm subgroup}
#'   \item{seHc_obs}{SE for no-harm subgroup}
#'
#' @export
fit_cox_models <- function(df, formula) {

  # Validate treat.recommend exists

if (!"treat.recommend" %in% names(df)) {
    stop("df must contain 'treat.recommend' column")
  }

  # Fit harm subgroup (treat.recommend == 0)
  df_H <- df[df$treat.recommend == 0L, , drop = FALSE]
  fit_H <- get_Cox_sg(df_sg = df_H, cox_formula = formula, est_loghr = TRUE)

  # Fit no-harm subgroup (treat.recommend == 1)
  df_Hc <- df[df$treat.recommend == 1L, , drop = FALSE]
  fit_Hc <- get_Cox_sg(df_sg = df_Hc, cox_formula = formula, est_loghr = TRUE)

  list(
    H_obs = fit_H$est_obs,
    seH_obs = fit_H$se_obs,
    Hc_obs = fit_Hc$est_obs,
    seHc_obs = fit_Hc$se_obs
  )
}


#' Fast HR Calculation for Split Samples
#'
#' Optimized function for calculating hazard ratios in split-sample
#' consistency evaluation. Minimizes overhead for repeated calls.
#'
#' @param Y Numeric vector of outcome times (full data).
#' @param E Numeric vector of event indicators (full data).
#' @param Treat Numeric vector of treatment indicators (full data).
#' @param subgroup_idx Logical or integer vector identifying subgroup.
#' @param split_idx Logical or integer vector identifying split sample.
#'
#' @return Numeric hazard ratio, or NA if model fails.
#'
#' @keywords internal
get_split_hr_fast <- function(Y, E, Treat, subgroup_idx, split_idx) {

  # Combine indices
  in_subset <- subgroup_idx & split_idx

  # Quick exit if insufficient data
  n_subset <- sum(in_subset)
  if (n_subset < 10L) return(NA_real_)

  n_events <- sum(E[in_subset])
  if (n_events < 2L) return(NA_real_)

  n_treat <- sum(Treat[in_subset])
  if (n_treat == 0L || n_treat == n_subset) return(NA_real_)

  # Fit minimal Cox model
  fit <- tryCatch({
    survival::coxph(
      survival::Surv(Y[in_subset], E[in_subset]) ~ Treat[in_subset],
      model = FALSE,
      x = FALSE,
      y = FALSE
    )
  }, error = function(e) NULL)

  if (is.null(fit)) return(NA_real_)

  # Extract HR directly
  exp(coef(fit)[1L])
}


#' Compute HR and Confidence Interval for Subgroup
#'
#' Wrapper function with configurable confidence level. Used in
#' forest plot generation and summary tables.
#'
#' @param df Data frame for subgroup.
#' @param outcome_name Character. Outcome variable name.
#' @param event_name Character. Event variable name.
#' @param treat_name Character. Treatment variable name.
#' @param conf_level Numeric. Confidence level (default: 0.95).
#'
#' @return Named list with HR, lower, upper, se, n_events.
#'
#' @export
compute_subgroup_hr <- function(df,
                                outcome_name,
                                event_name,
                                treat_name,
                                conf_level = 0.95) {

  Y <- df[[outcome_name]]
  E <- df[[event_name]]
  Treat <- df[[treat_name]]

  n <- nrow(df)
  n_events <- sum(E)
  n_treat <- sum(Treat)
  n_control <- n - n_treat

  # Check for sufficient data
  if (n_events < 2L || n_treat == 0L || n_control == 0L) {
    return(list(
      hr = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      se = NA_real_,
      n_events = n_events,
      n_treat = n_treat,
      n_control = n_control
    ))
  }

  # Fit Cox model
  fit <- tryCatch({
    survival::coxph(
      survival::Surv(Y, E) ~ Treat,
      robust = TRUE,
      model = FALSE,
      x = FALSE,
      y = FALSE
    )
  }, error = function(e) NULL)

  if (is.null(fit)) {
    return(list(
      hr = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      se = NA_real_,
      n_events = n_events,
      n_treat = n_treat,
      n_control = n_control
    ))
  }

  # Extract results
  fit_sum <- summary(fit, conf.int = conf_level)
  ci <- fit_sum$conf.int[1L, ]

  list(
    hr = ci[1L],
    lower = ci[3L],
    upper = ci[4L],
    se = fit_sum$coefficients[1L, "robust se"],
    n_events = n_events,
    n_treat = n_treat,
    n_control = n_control
  )
}
