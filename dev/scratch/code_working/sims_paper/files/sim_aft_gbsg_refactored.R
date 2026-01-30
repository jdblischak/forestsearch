# =============================================================================
# GBSG-Based AFT Data Generating Mechanism (Refactored)
# =============================================================================
#
# Refactored simulation functions for survival data based on the German Breast
# Cancer Study Group (GBSG) dataset. Aligned with forestsearch package conventions.
#
# Key functions:
#   - create_gbsg_dgm(): Create DGM with treatment effect heterogeneity
#   - simulate_from_gbsg_dgm(): Generate simulated trial data
#   - calibrate_k_inter(): Find k_inter for target subgroup HR
#   - get_dgm_with_output(): Wrapper for DGM creation with file output
#
# =============================================================================

# =============================================================================
# Global Variables Declaration
# =============================================================================

utils::globalVariables(c(

  # GBSG dataset variables
  "rfstime", "status", "hormon", "er", "age", "pgr", "meno", "nodes",
  "grade", "size",
  
 # Derived variables
  "y", "event", "treat", "z1", "z2", "z3", "z4", "z5", "zh", "flag.harm",
  "v1", "v2", "v3", "v4", "v5", "v6", "v7", "grade3",
  "lin.conf.true", "lin1.conf", "lin0.conf", "linC1.conf", "linC0.conf",
  "hlin.conf.1", "hlin.conf.0", "hlin.ratio", "h1.potential", "h0.potential",
  "Ts", "es", "t.sim", "y.sim", "event.sim"
))


# =============================================================================
# Constants
# =============================================================================

# Default seed for reproducibility
SEED_BASE <- 8316951L

# Days per month conversion
DAYS_PER_MONTH <- 30.4375

# Default super-population size
DEFAULT_N_SUPER <- 5000L


# =============================================================================
# Helper Functions
# =============================================================================

#' Discretize Continuous Variable into Quantile-Based Categories
#'
#' @param x Numeric vector to discretize
#' @param probs Numeric vector of probabilities for quantile breaks.
#'   Default: c(0.25, 0.5, 0.75) creates quartiles coded as 1, 2, 3, 4
#'
#' @return Integer vector with category codes (1 = lowest, max = highest)
#'
#' @keywords internal

cut_numeric <- function(x, probs = c(0.25, 0.5, 0.75)) {
  breaks <- stats::quantile(x, probs = probs, na.rm = TRUE)
  result <- rep(1L, length(x))
  for (i in seq_along(breaks)) {
    result[x > breaks[i]] <- i + 1L
  }
  result
}


#' Discretize Continuous Variable by Size Categories
#'
#' @param x Numeric vector (typically tumor size)
#' @param breaks Numeric vector of breakpoints. Default: c(20, 50)
#'
#' @return Integer vector with category codes
#'
#' @keywords internal

cut_size <- function(x, breaks = c(20, 50)) {
  result <- rep(1L, length(x))
  for (i in seq_along(breaks)) {
    result[x > breaks[i]] <- i + 1L
  }
  result
}


# =============================================================================
# Main DGM Creation Function
# =============================================================================

#' Create GBSG-Based AFT Data Generating Mechanism
#'
#' Creates a data generating mechanism (DGM) for survival simulations based on
#' the German Breast Cancer Study Group (GBSG) dataset. Supports heterogeneous
#' treatment effects via treatment-subgroup interactions.
#'
#' @param model Character. Either "alt" for alternative hypothesis with
#'   heterogeneous treatment effects, or "null" for uniform treatment effect.
#'   Default: "alt"
#' @param k_treat Numeric. Treatment effect multiplier applied to the treatment
#'   coefficient from the fitted AFT model. Values > 1 strengthen the treatment
#'   effect. Default: 1
#' @param k_inter Numeric. Interaction effect multiplier for the
#'   treatment-subgroup interaction (z1 * z3). Only used when model = "alt".
#'   Default: 1
#' @param k_z3 Numeric. Effect multiplier for the z3 (menopausal status)
#'   coefficient. Default: 1
#' @param z1_quantile Numeric. Quantile threshold for z1 (estrogen receptor).
#'   Observations with ER <= quantile are coded as z1 = 1. Default: 0.25
#' @param n_super Integer. Size of super-population for empirical HR estimation.
#'   Default: 5000
#' @param cens_type Character. Censoring distribution type: "weibull" or
#'   "uniform". Default: "weibull"
#' @param use_rand_params Logical. If TRUE, modifies confounder coefficients
#'   using estimates from randomized subset (meno == 0). Default: FALSE
#' @param seed Integer. Random seed for super-population generation.
#'   Default: 8316951
#' @param verbose Logical. Print diagnostic information. Default: FALSE
#'
#' @return A list of class "gbsg_dgm" containing:
#' \describe{
#'   \item{df_super_rand}{Data frame with randomized super-population including
#'     potential outcomes}
#'   \item{hr_H_true}{Empirical hazard ratio in harm subgroup (NA if model = "null")}
#'   \item{hr_Hc_true}{Empirical hazard ratio in complement subgroup}
#'   \item{hr_causal}{Overall causal (ITT) hazard ratio}
#'   \item{model_params}{List with AFT model parameters (mu, sigma, gamma, etc.)}
#'   \item{cens_params}{List with censoring model parameters}
#'   \item{subgroup_info}{List with subgroup definitions and true factor names}
#'   \item{analysis_vars}{Character vector of analysis variable names}
#'   \item{model_type}{Character indicating "alt" or "null"}
#' }
#'
#' @details
#' ## Subgroup Definition
#'
#' The harm subgroup H is defined as: z1 = 1 AND z3 = 1, where:
#' - z1: Low estrogen receptor (ER <= 25th percentile by default)
#' - z3: Premenopausal status (meno == 0)
#'
#' ## Model Specification
#'
#' The AFT model uses covariates: treat, z1, z2, z3, z4, z5, and (for "alt")
#' the interaction zh = treat * z1 * z3.
#'
#' @examples
#' \dontrun{
#' # Alternative hypothesis with default parameters
#' dgm_alt <- create_gbsg_dgm(model = "alt", verbose = TRUE)
#'
#' # Null hypothesis
#' dgm_null <- create_gbsg_dgm(model = "null", verbose = TRUE)
#'
#' # Custom subgroup HR via k_inter
#' dgm_custom <- create_gbsg_dgm(
#'   model = "alt",
#'   k_treat = 1.2,
#'   k_inter = 2.0,
#'   verbose = TRUE
#' )
#' }
#'
#' @seealso
#' \code{\link{simulate_from_gbsg_dgm}} for generating data from the DGM
#' \code{\link{calibrate_k_inter}} for finding k_inter to achieve target HR
#'
#' @importFrom survival survreg coxph Surv gbsg
#' @importFrom stats coef quantile
#' @export

create_gbsg_dgm <- function(
    model = c("alt", "null"),
    k_treat = 1,
    k_inter = 1,
    k_z3 = 1,
    z1_quantile = 0.25,
    n_super = DEFAULT_N_SUPER,
    cens_type = c("weibull", "uniform"),
    use_rand_params = FALSE,
    seed = SEED_BASE,
    verbose = FALSE
) {
  
 # -------------------------------------------------------------------------
  # Input Validation
  # -------------------------------------------------------------------------
  model <- match.arg(model)
  cens_type <- match.arg(cens_type)
  
  stopifnot(
    "k_treat must be positive" = k_treat > 0,
    "k_z3 must be numeric" = is.numeric(k_z3),
    "z1_quantile must be between 0 and 1" = z1_quantile > 0 && z1_quantile < 1,
    "n_super must be positive integer" = n_super > 0
 )
  
  # -------------------------------------------------------------------------
  # Load and Prepare GBSG Data
  # -------------------------------------------------------------------------
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' required for GBSG dataset.")
  }
  
  # Get GBSG data
  dfa <- survival::gbsg
  
  # Create derived variables
  dfa$y <- dfa$rfstime / DAYS_PER_MONTH
  dfa$id <- seq_len(nrow(dfa))
  dfa$event <- ifelse(dfa$status == 1, 1L, 0L)
  dfa$treat <- dfa$hormon
  
  # Create subgroup-defining variables
  er_threshold <- stats::quantile(dfa$er, probs = z1_quantile)
  dfa$z1 <- ifelse(dfa$er <= er_threshold, 1L, 0L)
  dfa$z2 <- cut_numeric(dfa$age)
  dfa$z3 <- ifelse(dfa$meno == 0, 1L, 0L)  # Premenopausal
  dfa$z4 <- cut_numeric(dfa$pgr)
  dfa$z5 <- cut_numeric(dfa$nodes)
  
  # Interaction term (treatment * subgroup)
  dfa$zh <- dfa$treat * dfa$z1 * dfa$z3
  
  # Harm subgroup indicator
  dfa$flag.harm <- ifelse(dfa$z1 == 1 & dfa$z3 == 1, 1L, 0L)
  
  # Analysis factors (observed by analyst)
  dfa$v1 <- as.factor(dfa$z1)
  dfa$v2 <- as.factor(dfa$z2)
  dfa$v3 <- as.factor(dfa$z3)
  dfa$v4 <- as.factor(dfa$z4)
  dfa$v5 <- as.factor(dfa$z5)
  dfa$v6 <- as.factor(cut_size(dfa$size))
  dfa$grade3 <- ifelse(dfa$grade == 3, 1L, 0L)
  dfa$v7 <- as.factor(dfa$grade3)
  
  # -------------------------------------------------------------------------
  # Define Subgroup Identities
  # -------------------------------------------------------------------------
  if (model == "alt") {
    fs_harm_true <- c("v1.1", "v3.1")
    grf_harm_true <- c("v1=1", "v3=1")
  } else {
    fs_harm_true <- NULL
    grf_harm_true <- NULL
    dfa$flag.harm <- 0L
  }
  
  # -------------------------------------------------------------------------
  # Display Kaplan-Meier Plot (if verbose)
  # -------------------------------------------------------------------------
  if (verbose) {
    message("=== GBSG DGM Creation ===")
    message(sprintf("Model type: %s", model))
    message(sprintf("Sample size: %d", nrow(dfa)))
    message(sprintf("Events: %d (%.1f%%)", sum(dfa$event), 
                    100 * mean(dfa$event)))
    message(sprintf("Harm subgroup size: %d (%.1f%%)", 
                    sum(dfa$flag.harm), 100 * mean(dfa$flag.harm)))
  }
  
  # -------------------------------------------------------------------------
  # Fit AFT Model
  # -------------------------------------------------------------------------
  covs_true <- c("treat", "z1", "z2", "z3", "z4", "z5")
  loc_inter <- length(covs_true) + 1
  
  # Design matrices for potential outcomes
  if (model == "alt") {
    z_true <- as.matrix(dfa[, c(covs_true, "zh")])
    
    # Potential outcome under treatment
    z_true_1 <- z_true
    z_true_1[, 1] <- 1  # Set treat = 1
    z_true_1[, loc_inter] <- z_true[, "z1"] * z_true[, "z3"]
    
    # Potential outcome under control
    z_true_0 <- z_true
    z_true_0[, 1] <- 0  # Set treat = 0
    z_true_0[, loc_inter] <- 0
    
  } else {
    z_true <- as.matrix(dfa[, covs_true])
    
    z_true_1 <- z_true
    z_true_1[, 1] <- 1
    
    z_true_0 <- z_true
    z_true_0[, 1] <- 0
  }
  
  # Fit Weibull AFT model
  fit_aft <- survival::survreg(
    survival::Surv(y, event) ~ z_true,
    data = dfa,
    dist = "weibull"
  )
  names(fit_aft$coefficients) <- c("(Intercept)", colnames(z_true))
  
  # Extract parameters
  sigma <- fit_aft$scale
  mu <- stats::coef(fit_aft)[1]
  gamma <- stats::coef(fit_aft)[-1]
  
  # -------------------------------------------------------------------------
  # Optionally Modify Parameters Using Randomized Subset
  # -------------------------------------------------------------------------
  if (use_rand_params) {
    df_rand <- subset(dfa, meno == 0)
    fit_rand <- survival::survreg(
      survival::Surv(y, event) ~ treat + z1 + z2 + z3 + z4 + z5,
      data = df_rand,
      dist = "weibull"
    )
    gamma_rand <- stats::coef(fit_rand)[-1]
    gamma[c("z1", "z2", "z4", "z5")] <- gamma_rand[c("z1", "z2", "z4", "z5")]
  }
  
  # -------------------------------------------------------------------------
  # Apply Effect Modifiers
  # -------------------------------------------------------------------------
  gamma["z3"] <- k_z3 * gamma["z3"]
  gamma["treat"] <- k_treat * gamma["treat"]
  
  if (model == "alt") {
    gamma["zh"] <- k_inter * gamma["zh"]
  }
  
  # Convert to hazard scale
  b_true <- gamma
  b0 <- -gamma / sigma
  
  # Linear predictors
  lin_conf <- z_true %*% b_true
  lin1_conf <- z_true_1 %*% b_true
  lin0_conf <- z_true_0 %*% b_true
  
  # Potential outcome hazard ratios
  dfa$hlin.conf.1 <- exp(z_true_1 %*% b0)
  dfa$hlin.conf.0 <- exp(z_true_0 %*% b0)
  
  # -------------------------------------------------------------------------
  # Generate Randomized Super-Population
  # -------------------------------------------------------------------------
  set.seed(seed)
  
  n_treat <- round(n_super / 2)
  n_control <- n_super - n_treat
  
  # Sample with replacement
  id_sample <- sample(seq_len(nrow(dfa)), size = n_super, replace = TRUE)
  df_samp <- dfa[id_sample, ]
  
  # Add potential outcome linear predictors
  df_samp$lin1.conf <- lin1_conf[id_sample]
  df_samp$lin0.conf <- lin0_conf[id_sample]
  
  # Assign treatment groups
  df_treat <- df_samp[seq_len(n_treat), ]
  df_treat$treat <- 1L
  df_treat$lin.conf.true <- df_treat$lin1.conf
  
  df_control <- df_samp[(n_treat + 1):n_super, ]
  df_control$treat <- 0L
  df_control$lin.conf.true <- df_control$lin0.conf
  
  df_big <- rbind(df_treat, df_control)
  
  # -------------------------------------------------------------------------
  # Compute Empirical Hazard Ratios
  # -------------------------------------------------------------------------
  # Generate survival times
  epsilon <- log(stats::rexp(n_super))
  log_Ts <- mu + sigma * epsilon + df_big$lin.conf.true
  df_big$Ts <- exp(log_Ts)
  df_big$es <- 1L
  
  if (model == "alt") {
    hr_H_true <- exp(survival::coxph(
      survival::Surv(Ts, es) ~ treat,
      data = subset(df_big, flag.harm == 1)
    )$coefficients)
    
    hr_Hc_true <- exp(survival::coxph(
      survival::Surv(Ts, es) ~ treat,
      data = subset(df_big, flag.harm == 0)
    )$coefficients)
    
  } else {
    hr_H_true <- NA_real_
    hr_Hc_true <- exp(survival::coxph(
      survival::Surv(Ts, es) ~ treat,
      data = df_big
    )$coefficients)
  }
  
  hr_causal <- exp(survival::coxph(
    survival::Surv(Ts, es) ~ treat,
    data = df_big
  )$coefficients)
  
  if (verbose) {
    message(sprintf("Empirical HR (harm subgroup): %.3f", hr_H_true))
    message(sprintf("Empirical HR (complement): %.3f", hr_Hc_true))
    message(sprintf("Empirical HR (overall/causal): %.3f", hr_causal))
  }
  
  # -------------------------------------------------------------------------
  # Fit Censoring Model
  # -------------------------------------------------------------------------
  cens_model <- NULL
  mu_cens <- NULL
  sigma_cens <- NULL
  
  if (cens_type == "weibull") {
    z_cens <- as.matrix(df_big[, covs_true])
    
    fit_cens <- survival::survreg(
      survival::Surv(y, 1 - event) ~ z_cens,
      data = df_big,
      dist = "weibull"
    )
    
    sigma_cens <- fit_cens$scale
    mu_cens <- stats::coef(fit_cens)[1]
    gamma_cens <- stats::coef(fit_cens)[-1]
    b0_cens <- -gamma_cens / sigma_cens
    b_cens <- -b0_cens * sigma
    
    # Censoring linear predictors for potential outcomes
    z_cens_1 <- z_cens
    z_cens_1[, 1] <- 1
    z_cens_0 <- z_cens
    z_cens_0[, 1] <- 0
    
    df_big$linC1.conf <- as.vector(z_cens_1 %*% b_cens)
    df_big$linC0.conf <- as.vector(z_cens_0 %*% b_cens)
    
    cens_model <- list(
      type = "weibull",
      mu = mu_cens,
      sigma = sigma_cens,
      gamma = gamma_cens
    )
  }
  
  # -------------------------------------------------------------------------
  # Finalize Super-Population Data
  # -------------------------------------------------------------------------
  df_big$hlin.ratio <- df_big$hlin.conf.1 / df_big$hlin.conf.0
  df_big$h1.potential <- df_big$hlin.conf.1
  df_big$h0.potential <- df_big$hlin.conf.0
  
  # -------------------------------------------------------------------------
  # Assemble Output
  # -------------------------------------------------------------------------
  analysis_vars <- c("v1", "v2", "v3", "v4", "v5", "v6", "v7")
  
  result <- list(
    df_super_rand = df_big,
    hr_H_true = hr_H_true,
    hr_Hc_true = hr_Hc_true,
    hr_causal = hr_causal,
    model_params = list(
      mu = mu,
      sigma = sigma,
      gamma = gamma,
      b_true = b_true,
      b_hr = b0
    ),
    cens_params = list(
      type = cens_type,
      mu = mu_cens,
      sigma = sigma_cens
    ),
    subgroup_info = list(
      fs_harm_true = fs_harm_true,
      grf_harm_true = grf_harm_true,
      z1_quantile = z1_quantile,
      definition = "z1 == 1 & z3 == 1 (low ER & premenopausal)"
    ),
    analysis_vars = analysis_vars,
    model_type = model,
    n_super = n_super,
    seed = seed
  )
  
  class(result) <- c("gbsg_dgm", "list")
  
  return(result)
}


# =============================================================================
# Simulation Function
# =============================================================================

#' Simulate Trial Data from GBSG DGM
#'
#' Generates simulated clinical trial data from a GBSG-based data generating
#' mechanism.
#'
#' @param dgm A "gbsg_dgm" object from \code{\link{create_gbsg_dgm}}
#' @param n Integer. Sample size. If NULL, uses full super-population.
#'   Default: NULL
#' @param rand_ratio Numeric. Randomization ratio (treatment:control).
#'   Default: 1 (1:1 randomization)
#' @param sim_id Integer. Simulation ID used for seed offset. Default: 1
#' @param max_follow Numeric. Administrative censoring time (months).
#'   Default: Inf (no administrative censoring)
#' @param muC_adj Numeric. Adjustment to censoring distribution location
#'   parameter. Positive values increase censoring. Default: 0
#' @param min_cens Numeric. Minimum censoring time for uniform censoring.
#'   Required if cens_type = "uniform"
#' @param max_cens Numeric. Maximum censoring time for uniform censoring.
#'   Required if cens_type = "uniform"
#' @param draw_treatment Logical. If TRUE, randomly assigns treatment.
#'   If FALSE, samples from existing treatment arms. Default: TRUE
#'
#' @return A data.frame with simulated trial data containing:
#' \describe{
#'   \item{id}{Subject ID}
#'   \item{y.sim}{Observed time (min of event time and censoring time)}
#'   \item{event.sim}{Event indicator (1 = event, 0 = censored)}
#'   \item{t.sim}{True (latent) event time}
#'   \item{treat}{Treatment indicator}
#'   \item{flag.harm}{True harm subgroup indicator}
#'   \item{v1-v7}{Analysis covariates}
#'   \item{z1-z5, grade3, size}{Underlying covariates}
#' }
#'
#' @examples
#' \dontrun{
#' # Create DGM
#' dgm <- create_gbsg_dgm(model = "alt")
#'
#' # Simulate a trial
#' sim_data <- simulate_from_gbsg_dgm(
#'   dgm = dgm,
#'   n = 500,
#'   sim_id = 1,
#'   max_follow = 60
#' )
#'
#' # Check event rate
#' mean(sim_data$event.sim)
#' }
#'
#' @seealso \code{\link{create_gbsg_dgm}} for creating the DGM
#'
#' @importFrom stats rexp runif
#' @export

simulate_from_gbsg_dgm <- function(
    dgm,
    n = NULL,
    rand_ratio = 1,
    sim_id = 1L,
    max_follow = Inf,
    muC_adj = 0,
    min_cens = NULL,
    max_cens = NULL,
    draw_treatment = TRUE
) {
  
  # -------------------------------------------------------------------------
  # Input Validation
 # -------------------------------------------------------------------------
  stopifnot(
    "dgm must be a gbsg_dgm object" = inherits(dgm, "gbsg_dgm"),
    "rand_ratio must be positive" = rand_ratio > 0,
    "sim_id must be positive integer" = sim_id >= 1
  )
  
  cens_type <- dgm$cens_params$type
  
  if (cens_type == "uniform" && (is.null(min_cens) || is.null(max_cens))) {
    stop("For uniform censoring, min_cens and max_cens must be specified")
  }
  
  # -------------------------------------------------------------------------
  # Set Seed for Reproducibility
  # -------------------------------------------------------------------------
  set.seed(SEED_BASE + 1000L * sim_id)
  
  # -------------------------------------------------------------------------
  # Sample from Super-Population
  # -------------------------------------------------------------------------
  df_super <- dgm$df_super_rand
  
  if (is.null(n)) {
    df_sim <- df_super
  } else {
    # Calculate group sizes
    n0 <- round(n / (1 + rand_ratio))
    n1 <- n - n0
    
    if (draw_treatment) {
      # Random sampling with treatment assignment
      id_sample <- sample(seq_len(nrow(df_super)), size = n, replace = TRUE)
      df_sample <- df_super[id_sample, ]
      
      # Assign treatment
      df_treat <- df_sample[seq_len(n1), ]
      df_treat$treat <- 1L
      df_treat$lin.conf.true <- df_treat$lin1.conf
      if ("linC1.conf" %in% names(df_treat)) {
        df_treat$linC.conf <- df_treat$linC1.conf
      }
      
      df_control <- df_sample[(n1 + 1):n, ]
      df_control$treat <- 0L
      df_control$lin.conf.true <- df_control$lin0.conf
      if ("linC0.conf" %in% names(df_control)) {
        df_control$linC.conf <- df_control$linC0.conf
      }
      
      df_sim <- rbind(df_control, df_treat)
      
    } else {
      # Sample from existing treatment arms
      df1_super <- subset(df_super, treat == 1)
      df0_super <- subset(df_super, treat == 0)
      
      id1 <- sample(seq_len(nrow(df1_super)), size = n1, replace = TRUE)
      id0 <- sample(seq_len(nrow(df0_super)), size = n0, replace = TRUE)
      
      df_sim <- rbind(df0_super[id0, ], df1_super[id1, ])
    }
  }
  
  n_obs <- nrow(df_sim)
  
  # -------------------------------------------------------------------------
  # Generate Event Times
  # -------------------------------------------------------------------------
  mu <- dgm$model_params$mu
  sigma <- dgm$model_params$sigma
  
  epsilon <- log(stats::rexp(n_obs))
  log_T <- mu + sigma * epsilon + df_sim$lin.conf.true
  T_sim <- exp(log_T)
  
  # -------------------------------------------------------------------------
  # Generate Censoring Times
  # -------------------------------------------------------------------------
  if (cens_type == "weibull") {
    mu_cens <- dgm$cens_params$mu + muC_adj
    sigma_cens <- dgm$cens_params$sigma
    
    epsilon_cens <- log(stats::rexp(n_obs))
    log_C <- mu_cens + sigma_cens * epsilon_cens + df_sim$linC.conf
    C_sim <- exp(log_C)
    
  } else {
    C_sim <- stats::runif(n_obs, min = min_cens, max = max_cens)
  }
  
  # Apply administrative censoring
  C_sim <- pmin(C_sim, max_follow)
  
  # -------------------------------------------------------------------------
  # Create Observed Outcomes
  # -------------------------------------------------------------------------
  event_sim <- ifelse(T_sim <= C_sim, 1L, 0L)
  y_sim <- pmin(T_sim, C_sim)
  
  # -------------------------------------------------------------------------
  # Finalize Output
  # -------------------------------------------------------------------------
  df_sim$y.sim <- y_sim
  df_sim$event.sim <- event_sim
  df_sim$t.sim <- T_sim
  df_sim$id <- seq_len(n_obs)
  
  # Select output columns
  keep_cols <- c(
    "id", "y.sim", "event.sim", "t.sim", "treat", "flag.harm",
    "v1", "v2", "v3", "v4", "v5", "v6", "v7",
    "z1", "z2", "z3", "z4", "z5", "grade3", "size",
    "hlin.ratio", "h1.potential", "h0.potential"
  )
  keep_cols <- intersect(keep_cols, names(df_sim))
  
  as.data.frame(df_sim[, keep_cols])
}


# =============================================================================
# Calibration Functions
# =============================================================================

#' Calibrate k_inter for Target Subgroup Hazard Ratio
#'
#' Finds the interaction effect multiplier (k_inter) that achieves a target
#' hazard ratio in the harm subgroup.
#'
#' @param target_hr_harm Numeric. Target hazard ratio for the harm subgroup
#' @param model Character. Model type ("alt" only). Default: "alt"
#' @param k_treat Numeric. Treatment effect multiplier. Default: 1
#' @param cens_type Character. Censoring type. Default: "weibull"
#' @param k_inter_range Numeric vector of length 2. Search range for k_inter.
#'   Default: c(-100, 100)
#' @param tol Numeric. Tolerance for root finding. Default: 1e-8
#' @param verbose Logical. Print diagnostic information. Default: FALSE
#' @param ... Additional arguments passed to \code{\link{create_gbsg_dgm}}
#'
#' @return Numeric value of k_inter that achieves the target HR
#'
#' @examples
#' \dontrun{
#' # Find k_inter for HR = 1.5 in harm subgroup
#' k <- calibrate_k_inter(target_hr_harm = 1.5, verbose = TRUE)
#'
#' # Verify
#' dgm <- create_gbsg_dgm(model = "alt", k_inter = k, verbose = TRUE)
#' }
#'
#' @importFrom stats uniroot
#' @export

calibrate_k_inter <- function(
    target_hr_harm,
    model = "alt",
    k_treat = 1,
    cens_type = "weibull",
    k_inter_range = c(-100, 100),
    tol = 1e-8,
    verbose = FALSE,
    ...
) {
  
  stopifnot(
    "target_hr_harm must be positive" = target_hr_harm > 0,
    "model must be 'alt' for calibration" = model == "alt"
  )
  
  # Objective function
  objective <- function(k_val) {
    dgm <- create_gbsg_dgm(
      model = model,
      k_treat = k_treat,
      k_inter = k_val,
      cens_type = cens_type,
      verbose = FALSE,
      ...
    )
    dgm$hr_H_true - target_hr_harm
  }
  
  if (verbose) {
    message(sprintf("Searching for k_inter to achieve HR(H) = %.3f", target_hr_harm))
    message(sprintf("Search range: [%.1f, %.1f]", k_inter_range[1], k_inter_range[2]))
  }
  
  result <- tryCatch({
    stats::uniroot(objective, interval = k_inter_range, tol = tol)
  }, error = function(e) {
    warning("Root finding failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(result)) {
    return(NA_real_)
  }
  
  k_inter <- result$root
  
  if (verbose) {
    # Verify result
    dgm_verify <- create_gbsg_dgm(
      model = model,
      k_treat = k_treat,
      k_inter = k_inter,
      cens_type = cens_type,
      verbose = FALSE,
      ...
    )
    message(sprintf("Found k_inter = %.6f", k_inter))
    message(sprintf("Achieved HR(H) = %.4f (target: %.4f)", 
                    dgm_verify$hr_H_true, target_hr_harm))
  }
  
  k_inter
}


#' Create DGM with Output File Path
#'
#' Wrapper function that creates a GBSG DGM and generates a standardized
#' output file path for saving results.
#'
#' @param model_harm Character. Model type ("alt" or "null")
#' @param n Integer. Planned sample size (for filename)
#' @param k_treat Numeric. Treatment effect multiplier
#' @param target_hr_harm Numeric. Target HR for harm subgroup (used for
#'   calibration when model = "alt")
#' @param cens_type Character. Censoring type
#' @param out_dir Character. Output directory path. If NULL, no file path
#'   is generated
#' @param file_prefix Character. Prefix for output filename
#' @param file_suffix Character. Suffix for output filename
#' @param include_hr_in_name Logical. Include achieved HR in filename.
#'   Default: FALSE
#' @param verbose Logical. Print diagnostic information. Default: FALSE
#' @param ... Additional arguments passed to \code{\link{create_gbsg_dgm}}
#'
#' @return List with components:
#' \describe{
#'   \item{dgm}{The gbsg_dgm object}
#'   \item{out_file}{Character path to output file (NULL if out_dir is NULL)}
#' }
#'
#' @export

get_dgm_with_output <- function(
    model_harm,
    n,
    k_treat = 1,
    target_hr_harm = NULL,
    cens_type = "weibull",
    out_dir = NULL,
    file_prefix = "sim",
    file_suffix = "",
    include_hr_in_name = FALSE,
    verbose = FALSE,
    ...
) {
  
  # Calibrate k_inter if needed
  k_inter <- NULL
  if (model_harm != "null" && !is.null(target_hr_harm)) {
    k_inter <- calibrate_k_inter(
      target_hr_harm = target_hr_harm,
      model = model_harm,
      k_treat = k_treat,
      cens_type = cens_type,
      verbose = verbose,
      ...
    )
  }
  
  # Create DGM
  dgm <- create_gbsg_dgm(
    model = model_harm,
    k_treat = k_treat,
    k_inter = if (is.null(k_inter)) 1 else k_inter,
    cens_type = cens_type,
    verbose = verbose,
    ...
  )
  
  # Generate output file path
  out_file <- NULL
  if (!is.null(out_dir)) {
    out_file <- file.path(
      out_dir,
      paste0(
        file_prefix, "_",
        "N=", n, "_",
        model_harm, "_",
        "ktreat=", k_treat,
        if (include_hr_in_name && model_harm == "alt") {
          paste0("_hrH=", round(dgm$hr_H_true, 2))
        } else "",
        file_suffix,
        ".Rdata"
      )
    )
  }
  
  list(dgm = dgm, out_file = out_file)
}


# =============================================================================
# Print Method
# =============================================================================

#' Print Method for gbsg_dgm Objects
#'
#' @param x A gbsg_dgm object
#' @param ... Additional arguments (unused)
#'
#' @export

print.gbsg_dgm <- function(x, ...) {
  cat("GBSG-Based AFT Data Generating Mechanism\n")
  cat("=========================================\n\n")
  
  cat("Model type:", x$model_type, "\n")
  cat("Super-population size:", x$n_super, "\n\n")
  
  cat("Hazard Ratios:\n")
  cat("  Overall (causal):", round(x$hr_causal, 3), "\n")
  if (!is.na(x$hr_H_true)) {
    cat("  Harm subgroup (H):", round(x$hr_H_true, 3), "\n")
  }
  cat("  Complement (H^c):", round(x$hr_Hc_true, 3), "\n\n")
  
  cat("Subgroup definition:", x$subgroup_info$definition, "\n")
  cat("Analysis variables:", paste(x$analysis_vars, collapse = ", "), "\n")
  
  invisible(x)
}
