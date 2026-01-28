# =============================================================================
# GBSG DGM Wrapper Functions
# =============================================================================
#
# Wrapper functions that provide GBSG-specific interfaces around the generic
# generate_aft_dgm_flex() and simulate_from_dgm() infrastructure.
#
# Benefits:
#   - Backward compatibility with existing GBSG simulation code
#   - Access to new features (spline effects, lognormal censoring)
#   - Single source of truth for DGM infrastructure
#   - Reduced code duplication
#
# Key Functions:
#   - create_gbsg_dgm_flex(): Wrapper around generate_aft_dgm_flex()
#   - simulate_from_gbsg_dgm_flex(): Wrapper around simulate_from_dgm()
#   - calibrate_k_inter_flex(): Wrapper around find_k_inter_for_target_hr()
#
# =============================================================================

#' @import survival
#' @import data.table
NULL

# =============================================================================
# Constants
# =============================================================================

GBSG_SEED_BASE <- 8316951L

# GBSG variable specifications
GBSG_CONTINUOUS_VARS <- c("age", "er", "pgr", "nodes", "size")
GBSG_FACTOR_VARS <- c("meno", "grade")
GBSG_OUTCOME_VAR <- "rfstime"
GBSG_EVENT_VAR <- "status"
GBSG_TREATMENT_VAR <- "hormon"

# Variable alias mapping (GBSG analysis vars -> flex z_ vars)
GBSG_VAR_ALIASES <- c(

v1 = "z_er",
  v2 = "z_meno",
  v3 = "z_age",
  v4 = "z_pgr",
  v5 = "z_nodes",
  v6 = "z_size",
  v7 = "z_grade_1"
)

# Column name mapping (flex -> GBSG style)
GBSG_COLUMN_MAP <- c(
  "y_sim" = "y.sim",
  "event_sim" = "event.sim",

  "t_true" = "t.sim",
  "treat_sim" = "treat",
  "flag_harm" = "flag.harm",
  "c_time" = "c.sim"
)


# =============================================================================
# Main DGM Creation Function
# =============================================================================

#' Create GBSG-Based DGM Using Flex Infrastructure
#'
#' Wrapper around \code{\link{generate_aft_dgm_flex}} with GBSG-specific defaults.
#' Provides backward compatibility with existing GBSG simulation code while
#' leveraging the more general DGM infrastructure.
#'
#' @param model Character. Model type: "alt" for alternative hypothesis with
#'   heterogeneous treatment effect, "null" for null hypothesis with uniform
#'   treatment effect. Default: "alt"
#' @param k_treat Numeric. Treatment effect multiplier. Values > 1 increase
#'   treatment effect magnitude, < 1 decrease it. Default: 1
#' @param k_inter Numeric. Interaction effect multiplier for treatment-subgroup
#'   interaction. Controls heterogeneity magnitude. Default: 1
#' @param k_z3 Numeric. Effect multiplier for menopausal status (z3).
#'   Default: 1 (not used in current implementation)
#' @param z1_quantile Numeric. Quantile of estrogen receptor (ER) distribution
#'   to use as subgroup cutpoint. Default: 0.25 (25th percentile)
#' @param n_super Integer. Size of super-population. Default: 5000
#' @param cens_type Character. Censoring distribution type: "weibull" or
#'   "uniform". Default: "weibull"
#' @param seed Integer. Random seed for reproducibility. Default: 8316951
#' @param verbose Logical. Print diagnostic information. Default: TRUE
#'
#' @return An object of class \code{c("gbsg_dgm", "aft_dgm_flex", "list")} containing:
#' \describe{
#'   \item{df_super_rand}{Super-population data with GBSG-style column names}
#'   \item{hr_H_true}{Cox-based HR in harm subgroup}
#'   \item{hr_Hc_true}{Cox-based HR in complement subgroup}
#'   \item{hr_causal}{Overall causal HR}
#'   \item{AHR}{Overall Average Hazard Ratio}
#'   \item{AHR_H_true}{AHR in harm subgroup}
#'   \item{AHR_Hc_true}{AHR in complement subgroup}
#'   \item{hazard_ratios}{Unified list of all HR metrics}
#'   \item{model_params}{Model parameters (mu, sigma/tau, gamma, b0)}
#'   \item{cens_params}{Censoring model parameters}
#'   \item{subgroup_info}{Subgroup definition and proportions
}
#'   \item{analysis_vars}{Character vector of analysis variable names (v1-v7)}
#'   \item{model_type}{Character: "alt" or "null"}
#'   \item{n_super}{Super-population size}
#'   \item{seed}{Random seed used}
#' }
#'
#' @details
#' ## Subgroup Definition
#'
#' The harm subgroup H is defined as: low estrogen receptor (ER <= quantile)
#' AND premenopausal (meno == 0). This matches the original GBSG simulation
#' framework from León et al. (2024).
#'
#' ## Model Structure
#'
#' Uses an AFT Weibull model: log(T) = μ + σε + X'γ
#'
#' The treatment effect is modified by:
#' - k_treat: Scales the overall treatment effect
#' - k_inter: Scales the treatment-subgroup interaction (heterogeneity)
#'
#' @examples
#' \dontrun
{
#' # Create DGM with heterogeneous treatment effect
#' dgm <- create_gbsg_dgm_flex(
#'   model = "alt",
#'   k_inter = 2,
#'   z1_quantile = 0.25,
#'   verbose = TRUE
#' )
#'
#' # Check hazard ratios
#' print(dgm$hr_H_true)
#' print(dgm$hr_Hc_true)
#'
#' # Access via unified list
#' dgm$hazard_ratios$harm_subgroup
#' dgm$hazard_ratios$AHR_harm
#' }
#'
#' @seealso
#' \code{\link{generate_aft_dgm_flex}} for the underlying DGM generator
#' \code{\link{simulate_from_gbsg_dgm_flex}} for simulating data from the DGM
#' \code{\link{calibrate_k_inter_flex}} for calibrating k_inter to achieve target HR
#'
#' @export
create_gbsg_dgm_flex <- function(
    model = "alt",
    k_treat = 1,
    k_inter = 1,
    k_z3 = 1,
    z1_quantile = 0.25,
    n_super = 5000L,
    cens_type = "weibull",
    seed = GBSG_SEED_BASE,
    verbose = TRUE
) {
  
 # -------------------------------------------------------------------------
  # Input Validation
  # -------------------------------------------------------------------------
  stopifnot(
    "model must be 'alt' or 'null'" = model %in% c("alt", "null"),
    "k_treat must be positive" = k_treat > 0,
    "z1_quantile must be between 0 and 1" = z1_quantile > 0 && z1_quantile < 1,
    "n_super must be positive integer" = n_super >= 100
  )
  
  # -------------------------------------------------------------------------
  # Load GBSG Data
  # -------------------------------------------------------------------------
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for GBSG data")
  }
  
  # Get GBSG dataset
  gbsg <- survival::gbsg
  
  if (verbose) {
    message("=== Creating GBSG DGM (Flex Wrapper) ===")
    message(sprintf("Model: %s, k_treat: %.2f, k_inter: %.2f", model, k_treat, k_inter))
    message(sprintf("ER quantile: %.2f, n_super: %d", z1_quantile, n_super))
  }
  
  # -------------------------------------------------------------------------
  # Determine Subgroup Specification
  # -------------------------------------------------------------------------
  # For null model, we don't want subgroup interaction
  if (model == "null") {
    subgroup_vars <- NULL
    subgroup_cuts <- NULL
    k_inter_use <- 0
  } else {
    subgroup_vars <- c("er", "meno")
    subgroup_cuts <- list(
      er = list(type = "quantile", value = z1_quantile),
      meno = 0  # Premenopausal
    )
    k_inter_use <- k_inter
  }
  
  # -------------------------------------------------------------------------
  # Call generate_aft_dgm_flex
  # -------------------------------------------------------------------------
  dgm_flex <- generate_aft_dgm_flex(
    data = gbsg,
    continuous_vars = GBSG_CONTINUOUS_VARS,
    factor_vars = GBSG_FACTOR_VARS,
    outcome_var = GBSG_OUTCOME_VAR,
    event_var = GBSG_EVENT_VAR,
    treatment_var = GBSG_TREATMENT_VAR,
    subgroup_vars = subgroup_vars,
    subgroup_cuts = subgroup_cuts,
    model = model,
    k_treat = k_treat,
    k_inter = k_inter_use,
    n_super = n_super,
    cens_type = cens_type,
    select_censoring = TRUE,
    seed = seed,
    verbose = verbose
  )
  
  # -------------------------------------------------------------------------
  # Add GBSG-Style Accessors (Backward Compatibility)
  # -------------------------------------------------------------------------
  
  # Direct HR accessors
  dgm_flex$hr_H_true <- dgm_flex$hazard_ratios$harm_subgroup
  dgm_flex$hr_Hc_true <- dgm_flex$hazard_ratios$no_harm_subgroup
  dgm_flex$hr_causal <- dgm_flex$hazard_ratios$overall
  
  # AHR accessors
  dgm_flex$AHR <- dgm_flex$hazard_ratios$AHR
  dgm_flex$AHR_H_true <- dgm_flex$hazard_ratios$AHR_harm
  dgm_flex$AHR_Hc_true <- dgm_flex$hazard_ratios$AHR_no_harm
  
  # For null model, set harm HR to overall
  if (model == "null") {
    dgm_flex$hr_H_true <- dgm_flex$hr_causal
    dgm_flex$AHR_H_true <- dgm_flex$AHR
  }
  
  # -------------------------------------------------------------------------
  # Add GBSG-Style Variable Aliases
  # -------------------------------------------------------------------------
  dgm_flex$df_super <- add_gbsg_variable_aliases(dgm_flex$df_super)
  
  # Rename df_super to df_super_rand for backward compatibility
  dgm_flex$df_super_rand <- dgm_flex$df_super
  
  # -------------------------------------------------------------------------
  # Add GBSG-Style Analysis Variables
  # -------------------------------------------------------------------------
  dgm_flex$analysis_vars <- c("v1", "v2", "v3", "v4", "v5", "v6", "v7")
  
  # -------------------------------------------------------------------------
  # Add Subgroup Info in GBSG Style
  # -------------------------------------------------------------------------
  if (model == "alt") {
    dgm_flex$subgroup_info$fs_harm_true <- paste(
      sprintf("v1 <= %.2f", quantile(gbsg$er, z1_quantile)),
      "v2 == 0",
      sep = " & "
    )
    dgm_flex$subgroup_info$grf_harm_true <- dgm_flex$subgroup_info$fs_harm_true
    dgm_flex$subgroup_info$z1_quantile <- z1_quantile
    dgm_flex$subgroup_info$definition <- sprintf(
      "z_er <= %.2f (q=%.2f) & z_meno == 0 (premenopausal)",
      quantile(gbsg$er, z1_quantile), z1_quantile
    )
  }
  
  # -------------------------------------------------------------------------
  # Convert Model Params to GBSG Style
  # -------------------------------------------------------------------------
  dgm_flex$cens_params <- list(
    type = dgm_flex$model_params$censoring$type,
    mu = dgm_flex$model_params$censoring$mu,
    sigma = dgm_flex$model_params$censoring$tau  # tau = sigma in AFT
  )
  
  # Add sigma alias for tau
  dgm_flex$model_params$sigma <- dgm_flex$model_params$tau
  
  # -------------------------------------------------------------------------
  # Update Class
  # -------------------------------------------------------------------------
  class(dgm_flex) <- c("gbsg_dgm", class(dgm_flex))
  
  if (verbose) {
    message("\n=== DGM Created Successfully ===")
    message(sprintf("HR(H): %.4f, HR(Hc): %.4f, HR(overall): %.4f",
                    dgm_flex$hr_H_true, dgm_flex$hr_Hc_true, dgm_flex$hr_causal))
    message(sprintf("AHR(H): %.4f, AHR(Hc): %.4f, AHR(overall): %.4f",
                    dgm_flex$AHR_H_true, dgm_flex$AHR_Hc_true, dgm_flex$AHR))
  }
  
  return(dgm_flex)
}


# =============================================================================
# Simulation Function
# =============================================================================

#' Simulate Trial Data from GBSG DGM (Flex Wrapper)
#'
#' Wrapper around \code{\link{simulate_from_dgm}} with GBSG-style column naming
#' and parameter conventions.
#'
#' @param dgm A DGM object from \code{\link{create_gbsg_dgm_flex}} or
#'   \code{\link{generate_aft_dgm_flex}}
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
#'   Required if cens_type = "uniform" in DGM
#' @param max_cens Numeric. Maximum censoring time for uniform censoring.
#'   Required if cens_type = "uniform" in DGM
#' @param draw_treatment Logical. If TRUE, randomly assigns treatment.
#'   If FALSE, samples from existing treatment arms. Default: TRUE
#'
#' @return A data.frame with simulated trial data containing GBSG-style columns:
#' \describe{
#'   \item{id}{Subject ID}
#'   \item{y.sim}{Observed time (min of event time and censoring time)}
#'   \item{event.sim}{Event indicator (1 = event, 0 = censored)}
#'   \item{t.sim}{True (latent) event time}
#'   \item{treat}{Treatment indicator}
#'   \item{flag.harm}{True harm subgroup indicator}
#'   \item{v1-v7}{Analysis covariates}
#'   \item{loghr_po}{Individual log hazard ratio (potential outcome)}
#' }
#'
#' @examples
#' \dontrun{
#' # Create DGM
#' dgm <- create_gbsg_dgm_flex(model = "alt", k_inter = 2)
#'
#' # Simulate a trial
#' sim_data <- simulate_from_gbsg_dgm_flex(
#'   dgm = dgm,
#'   n = 500,
#'   sim_id = 1,
#'   max_follow = 60
#' )
#'
#' # Check event rate
#' mean(sim_data$event.sim)
#'
#' # Fit Cox model
#' library(survival)
#' coxph(Surv(y.sim, event.sim) ~ treat, data = sim_data)
#' }
#'
#' @seealso
#' \code{\link{create_gbsg_dgm_flex}} for creating the DGM
#' \code{\link{simulate_from_dgm}} for the underlying simulation function
#'
#' @export
simulate_from_gbsg_dgm_flex <- function(
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
    "dgm must be a gbsg_dgm or aft_dgm_flex object" = 
      inherits(dgm, c("gbsg_dgm", "aft_dgm_flex")),
    "rand_ratio must be positive" = rand_ratio > 0,
    "sim_id must be positive integer" = sim_id >= 1
  )
  
  # Check for uniform censoring requirements
  cens_type <- dgm$cens_params$type %||% dgm$model_params$censoring$type
  if (identical(cens_type, "uniform") && (is.null(min_cens) || is.null(max_cens))) {
    stop("For uniform censoring, min_cens and max_cens must be specified")
  }
  
  # -------------------------------------------------------------------------
  # Set Seed Using GBSG Convention
  # -------------------------------------------------------------------------
  set.seed(GBSG_SEED_BASE + 1000L * sim_id)
  
  # -------------------------------------------------------------------------
  # Call Generic simulate_from_dgm
  # -------------------------------------------------------------------------
  sim_data <- simulate_from_dgm(
    dgm = dgm,
    n = n,
    rand_ratio = rand_ratio,
    max_entry = 0,              # No staggered entry (GBSG convention)
    analysis_time = max_follow, # Use max_follow as analysis time
    cens_adjust = muC_adj,
    draw_treatment = draw_treatment
  )
  
  # -------------------------------------------------------------------------
  # Handle Uniform Censoring (if applicable)
  # -------------------------------------------------------------------------
  if (identical(cens_type, "uniform") && !is.null(min_cens) && !is.null(max_cens)) {
    n_obs <- nrow(sim_data)
    c_uniform <- runif(n_obs, min = min_cens, max = max_cens)
    
    # Recalculate observed times with uniform censoring
    sim_data$c_time <- pmin(c_uniform, max_follow)
    sim_data$y_sim <- pmin(sim_data$t_true, sim_data$c_time)
    sim_data$event_sim <- as.integer(sim_data$t_true <= sim_data$c_time)
  }
  
  # -------------------------------------------------------------------------
  # Rename Columns to GBSG Convention
  # -------------------------------------------------------------------------
  sim_data <- rename_to_gbsg_columns(sim_data)
  
  # -------------------------------------------------------------------------
  # Add GBSG Variable Aliases (if not present)
  # -------------------------------------------------------------------------
  sim_data <- add_gbsg_variable_aliases(sim_data)
  
  # -------------------------------------------------------------------------
  # Select and Order Output Columns
  # -------------------------------------------------------------------------
  keep_cols <- c(
    "id", "y.sim", "event.sim", "t.sim", "treat", "flag.harm",
    "v1", "v2", "v3", "v4", "v5", "v6", "v7",
    "loghr_po", "theta_0", "theta_1"
  )
  keep_cols <- intersect(keep_cols, names(sim_data))
  
  # Also keep any z_ columns
  z_cols <- grep("^z_", names(sim_data), value = TRUE)
  keep_cols <- unique(c(keep_cols, z_cols))
  
  sim_data <- sim_data[, keep_cols, drop = FALSE]
  
  return(as.data.frame(sim_data))
}


# =============================================================================
# Calibration Function
# =============================================================================

#' Calibrate k_inter for Target Subgroup Hazard Ratio (Flex Wrapper)
#'
#' Wrapper around \code{\link{find_k_inter_for_target_hr}} with GBSG-specific
#' defaults.
#'
#' @param target_hr_harm Numeric. Target hazard ratio for the harm subgroup.
#' @param model Character. Model type ("alt" only). Default: "alt"
#' @param k_treat Numeric. Treatment effect multiplier. Default: 1
#' @param z1_quantile Numeric. ER quantile for subgroup definition. Default: 0.25
#' @param cens_type Character. Censoring type. Default: "weibull"
#' @param k_inter_range Numeric vector of length 2. Search range for k_inter.
#'   Default: c(-100, 100)
#' @param tol Numeric. Tolerance for root finding. Default: 1e-4
#' @param n_super Integer. Super-population size for each evaluation. Default: 5000
#' @param use_ahr Logical. If TRUE, calibrate to AHR instead of Cox-based HR.
#'   Default: FALSE
#' @param verbose Logical. Print diagnostic information. Default: TRUE
#'
#' @return Numeric value of k_inter that achieves the target HR (or AHR).
#'
#' @examples
#' \dontrun{
#' # Find k_inter for HR = 1.5 in harm subgroup
#' k_inter <- calibrate_k_inter_flex(target_hr_harm = 1.5, verbose = TRUE)
#'
#' # Verify
#' dgm <- create_gbsg_dgm_flex(model = "alt", k_inter = k_inter)
#' print(dgm$hr_H_true)  # Should be ~1.5
#' }
#'
#' @seealso
#' \code{\link{find_k_inter_for_target_hr}} for the generic calibration function
#' \code{\link{create_gbsg_dgm_flex}} for creating DGM with calibrated k_inter
#'
#' @export
calibrate_k_inter_flex <- function(
    target_hr_harm,
    model = "alt",
    k_treat = 1,
    z1_quantile = 0.25,
    cens_type = "weibull",
    k_inter_range = c(-100, 100),
    tol = 1e-4,
    n_super = 5000L,
    use_ahr = FALSE,
    verbose = TRUE
) {
  
  # -------------------------------------------------------------------------
  # Input Validation
  # -------------------------------------------------------------------------
  stopifnot(
    "target_hr_harm must be positive" = target_hr_harm > 0,
    "model must be 'alt'" = model == "alt"
  )
  
  # -------------------------------------------------------------------------
  # Load GBSG Data
  # -------------------------------------------------------------------------
  gbsg <- survival::gbsg
  
  if (verbose) {
    metric <- if (use_ahr) "AHR" else "Cox HR"
    message(sprintf("Calibrating k_inter for %s(H) = %.3f", metric, target_hr_harm))
  }
  
  # -------------------------------------------------------------------------
  # Objective Function
  # -------------------------------------------------------------------------
  objective_fn <- function(k_val) {
    dgm_temp <- create_gbsg_dgm_flex(
      model = "alt",
      k_treat = k_treat,
      k_inter = k_val,
      z1_quantile = z1_quantile,
      n_super = n_super,
      cens_type = cens_type,
      verbose = FALSE
    )
    
    if (use_ahr) {
      achieved <- dgm_temp$AHR_H_true
    } else {
      achieved <- dgm_temp$hr_H_true
    }
    
    return(achieved - target_hr_harm)
  }
  
  # -------------------------------------------------------------------------
  # Root Finding
  # -------------------------------------------------------------------------
  result <- tryCatch({
    uniroot(
      f = objective_fn,
      interval = k_inter_range,
      tol = tol
    )
  }, error = function(e) {
    # If root finding fails, try optimization
    if (verbose) {
      message("Root finding failed, trying optimization...")
    }
    
    opt_result <- optim(
      par = 0,
      fn = function(k) abs(objective_fn(k)),
      method = "Brent",
      lower = k_inter_range[1],
      upper = k_inter_range[2]
    )
    
    list(root = opt_result$par)
  })
  
  k_inter_optimal <- result$root
  
  # -------------------------------------------------------------------------
  # Verification
  # -------------------------------------------------------------------------
  if (verbose) {
    dgm_verify <- create_gbsg_dgm_flex(
      model = "alt",
      k_treat = k_treat,
      k_inter = k_inter_optimal,
      z1_quantile = z1_quantile,
      n_super = n_super,
      verbose = FALSE
    )
    
    achieved_hr <- dgm_verify$hr_H_true
    achieved_ahr <- dgm_verify$AHR_H_true
    
    message(sprintf("\nCalibration Result:"))
    message(sprintf("  k_inter = %.6f", k_inter_optimal))
    message(sprintf("  Achieved HR(H) = %.4f", achieved_hr))
    message(sprintf("  Achieved AHR(H) = %.4f", achieved_ahr))
    message(sprintf("  Target = %.4f", target_hr_harm))
  }
  
  return(k_inter_optimal)
}


# =============================================================================
# Validation Function
# =============================================================================

#' Validate k_inter Effect on Treatment Heterogeneity
#'
#' Tests a range of k_inter values to verify the interaction parameter
#' properly modulates treatment effect heterogeneity.
#'
#' @param k_inter_values Numeric vector. Values of k_inter to test.
#'   Default: c(-2, -1, 0, 1, 2, 3)
#' @param k_treat Numeric. Treatment effect multiplier. Default: 1
#' @param z1_quantile Numeric. ER quantile for subgroup. Default: 0.25
#' @param n_super Integer. Super-population size. Default: 5000
#' @param verbose Logical. Print results table. Default: TRUE
#'
#' @return A data.frame with columns: k_inter, hr_H, hr_Hc, ahr_H, ahr_Hc,
#'   ratio_cox, ratio_ahr
#'
#' @details
#' Key expectation: When k_inter = 0, there should be NO heterogeneity,
#' meaning HR(H) ≈ HR(Hc) and the ratio should be approximately 1.
#'
#' @examples
#' \dontrun{
#' validation <- validate_k_inter_effect_flex(
#'   k_inter_values = c(-2, -1, 0, 1, 2),
#'   verbose = TRUE
#' )
#' }
#'
#' @export
validate_k_inter_effect_flex <- function(
    k_inter_values = c(-2, -1, 0, 1, 2, 3),
    k_treat = 1,
    z1_quantile = 0.25,
    n_super = 5000L,
    verbose = TRUE
) {
  
  results <- data.frame(
    k_inter = numeric(),
    hr_H = numeric(),
    hr_Hc = numeric(),
    ahr_H = numeric(),
    ahr_Hc = numeric(),
    ratio_cox = numeric(),
    ratio_ahr = numeric()
  )
  
  for (k_val in k_inter_values) {
    dgm <- create_gbsg_dgm_flex(
      model = "alt",
      k_treat = k_treat,
      k_inter = k_val,
      z1_quantile = z1_quantile,
      n_super = n_super,
      verbose = FALSE
    )
    
    results <- rbind(results, data.frame(
      k_inter = k_val,
      hr_H = dgm$hr_H_true,
      hr_Hc = dgm$hr_Hc_true,
      ahr_H = dgm$AHR_H_true,
      ahr_Hc = dgm$AHR_Hc_true,
      ratio_cox = dgm$hr_H_true / dgm$hr_Hc_true,
      ratio_ahr = dgm$AHR_H_true / dgm$AHR_Hc_true
    ))
  }
  
  if (verbose) {
    message("\n=== k_inter Validation Results ===")
    message(sprintf("%-8s %-8s %-8s %-8s %-8s %-10s %-10s",
                    "k_inter", "HR(H)", "HR(Hc)", "AHR(H)", "AHR(Hc)",
                    "Ratio(Cox)", "Ratio(AHR)"))
    message(paste(rep("-", 70), collapse = ""))
    
    for (i in seq_len(nrow(results))) {
      row <- results[i, ]
      marker <- if (abs(row$k_inter) < 0.01) " <- k=0" else ""
      message(sprintf("%-8.1f %-8.4f %-8.4f %-8.4f %-8.4f %-10.4f %-10.4f%s",
                      row$k_inter, row$hr_H, row$hr_Hc,
                      row$ahr_H, row$ahr_Hc,
                      row$ratio_cox, row$ratio_ahr, marker))
    }
    
    message("\nNote: When k_inter = 0, ratio should be ~1.0 (no heterogeneity)")
  }
  
  return(results)
}


# =============================================================================
# Print Method
# =============================================================================

#' Print Method for GBSG DGM (Flex Version)
#'
#' @param x A gbsg_dgm object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.gbsg_dgm <- function(x, ...) {
  cat("\n")
  cat("====================================================\n")
  cat("GBSG Data Generating Mechanism (Flex Wrapper)\n")
  cat("====================================================\n\n")
  
  cat("Model Type:", x$model_type, "\n")
  cat("Super-population Size:", x$n_super, "\n")
  cat("Seed:", x$seed, "\n\n")
  
  cat("--- Hazard Ratios (Cox-based, stacked PO) ---\n")
  cat(sprintf("  Overall (causal): %.4f\n", x$hr_causal))
  cat(sprintf("  Harm subgroup (H): %.4f\n", x$hr_H_true))
  cat(sprintf("  Complement (Hc): %.4f\n", x$hr_Hc_true))
  cat(sprintf("  Ratio HR(H)/HR(Hc): %.4f\n", x$hr_H_true / x$hr_Hc_true))
  
  cat("\n--- Average Hazard Ratios (from loghr_po) ---\n")
  cat(sprintf("  AHR (overall): %.4f\n", x$AHR))
  cat(sprintf("  AHR_harm (H): %.4f\n", x$AHR_H_true))
  cat(sprintf("  AHR_no_harm (Hc): %.4f\n", x$AHR_Hc_true))
  cat(sprintf("  Ratio AHR(H)/AHR(Hc): %.4f\n", x$AHR_H_true / x$AHR_Hc_true))
  
  if (!is.null(x$subgroup_info)) {
    cat("\n--- Subgroup Definition ---\n")
    cat("  ", x$subgroup_info$definition %||% "See subgroup_info for details", "\n")
    cat(sprintf("  Size: %d (%.1f%%)\n",
                x$subgroup_info$size,
                100 * x$subgroup_info$proportion))
  }
  
  cat("\n--- Model Parameters ---\n")
  cat(sprintf("  mu (intercept): %.4f\n", x$model_params$mu))
  cat(sprintf("  sigma (scale): %.4f\n", x$model_params$sigma %||% x$model_params$tau))
  
  cat("\n--- Censoring ---\n")
  cat(sprintf("  Type: %s\n", x$cens_params$type))
  
  cat("\n--- Analysis Variables ---\n")
  cat("  ", paste(x$analysis_vars, collapse = ", "), "\n")
  
  cat("\n====================================================\n")
  
  invisible(x)
}


# =============================================================================
# Helper Functions
# =============================================================================

#' Add GBSG Variable Aliases (v1-v7) to Data Frame
#'
#' @param df Data frame with z_ prefixed variables
#' @return Data frame with added v1-v7 aliases
#' @keywords internal
add_gbsg_variable_aliases <- function(df) {
  
  for (alias in names(GBSG_VAR_ALIASES)) {
    source_col <- GBSG_VAR_ALIASES[alias]
    if (source_col %in% names(df) && !(alias %in% names(df))) {
      df[[alias]] <- df[[source_col]]
    }
  }
  
  return(df)
}


#' Rename Columns to GBSG Convention
#'
#' @param df Data frame with flex-style column names
#' @return Data frame with GBSG-style column names
#' @keywords internal
rename_to_gbsg_columns <- function(df) {
  
  for (old_name in names(GBSG_COLUMN_MAP)) {
    new_name <- GBSG_COLUMN_MAP[old_name]
    if (old_name %in% names(df)) {
      names(df)[names(df) == old_name] <- new_name
    }
  }
  
  return(df)
}


#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x
