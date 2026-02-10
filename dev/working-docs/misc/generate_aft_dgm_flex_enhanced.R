# Enhanced AFT Data Generating Mechanism with Flexible Subgroup Definitions and Spline Model Option
# Supports multiple cutpoint specifications: fixed values, quantiles, functions, and spline-based hazard ratio models

#' Generate Synthetic Survival Data using AFT Model with Flexible Subgroups or Spline Model
#'
#' Creates a data generating mechanism (DGM) for survival data using an Accelerated
#' Failure Time (AFT) model with Weibull distribution. Supports flexible subgroup
#' definitions, treatment-subgroup interactions, and spline-based hazard ratio models.
#'
#' @param data A data.frame containing the input dataset to base the simulation on
#' @param continuous_vars Character vector of continuous variable names to be
#'   standardized and included as covariates
#' @param factor_vars Character vector of factor/categorical variable names to be
#'   converted to dummy variables (largest value as reference)
#' @param outcome_var Character string specifying the name of the outcome/time variable
#' @param event_var Character string specifying the name of the event/status variable
#'   (1 = event, 0 = censored)
#' @param treatment_var Character string specifying the name of the treatment variable.
#'   If NULL, treatment will be randomly simulated with 50/50 allocation
#' @param subgroup_vars Character vector of variable names defining the subgroup.
#'   Default is NULL (no subgroups)
#' @param subgroup_cuts Named list of cutpoint specifications for subgroup variables.
#'   See Details section for flexible specification options
#' @param draw_treatment Logical indicating whether to redraw treatment assignment
#'   in simulation. Default is FALSE (use original assignments)
#' @param model Character string: "alt" for alternative model with subgroup effects,
#'   "null" for null model without subgroup effects, "spline" for spline-based
#'   hazard ratio model. Default is "alt"
#' @param spline_params List with spline model parameters when model="spline":
#'   \describe{
#'     \item{knot}{Numeric value for knot location. Default is 5}
#'     \item{zeta}{Numeric value for evaluation point. Default is 10}
#'     \item{log_hrs}{Numeric vector of log hazard ratios at z=0, z=knot, and z=zeta.
#'       Default is log(c(0.75, 0.75, 0.75))}
#'     \item{biomarker_var}{Character string specifying which continuous variable to use
#'       as the biomarker z for the spline model. If NULL, uses first continuous variable}
#'   }
#' @param k_treat Numeric treatment effect modifier. Values >1 increase treatment
#'   effect, <1 decrease it. Default is 1 (no modification)
#' @param k_inter Numeric interaction effect modifier for treatment-subgroup interaction.
#'   Default is 1 (no modification)
#' @param n_super Integer specifying size of super population to generate.
#'   Default is 5000
#' @param cens_type Character string specifying censoring distribution: "weibull"
#'   or "uniform". Default is "weibull"
#' @param cens_params List of parameters for censoring distribution. For uniform:
#'   list(min = value, max = value). For Weibull: fitted from data
#' @param seed Integer random seed for reproducibility. Default is 8316951
#' @param verbose Logical indicating whether to print diagnostic information during
#'   execution. Default is TRUE
#' @param standardize Logical indicating whether to standardize continuous variables.
#'   Default is FALSE
#'
#' @details
#' ## Model Types
#' 
#' ### Standard Subgroup Model (model = "alt" or "null")
#' Uses the flexible subgroup definition system with cutpoints as specified in the
#' original function.
#' 
#' ### Spline Model (model = "spline")
#' Creates a spline-based hazard ratio model where the log hazard ratio varies
#' smoothly with a biomarker z, with a change point at the specified knot.
#' The model creates interaction terms z.k = (z - knot) * I(z > knot) and
#' z.k.treat for the treatment interaction.
#'
#' ## Subgroup Cutpoint Specifications (for standard models)
#'
#' The `subgroup_cuts` parameter accepts multiple flexible specifications:
#'
#' ### Fixed Value
#' ```r
#' subgroup_cuts = list(er = 20)  # er <= 20
#' ```
#'
#' ### Quantile-based
#' ```r
#' subgroup_cuts = list(
#'   er = list(type = "quantile", value = 0.25)  # er <= 25th percentile
#' )
#' ```
#'
#' ### Function-based
#' ```r
#' subgroup_cuts = list(
#'   er = list(type = "function", fun = median)  # er <= median
#' )
#' ```
#'
#' ### Range
#' ```r
#' subgroup_cuts = list(
#'   age = list(type = "range", min = 40, max = 60)  # 40 <= age <= 60
#' )
#' ```
#'
#' ### Greater than
#' ```r
#' subgroup_cuts = list(
#'   nodes = list(type = "greater", quantile = 0.75)  # nodes > 75th percentile
#' )
#' ```
#'
#' ### Multiple values (for categorical)
#' ```r
#' subgroup_cuts = list(
#'   grade = list(type = "multiple", values = c(2, 3))  # grade in (2, 3)
#' )
#' ```
#'
#' ### Custom function
#' ```r
#' subgroup_cuts = list(
#'   er = list(
#'     type = "custom",
#'     fun = function(x) x <= quantile(x, 0.3) | x >= quantile(x, 0.9)
#'   )
#' )
#' ```
#'
#' ## Model Structure
#'
#' The AFT model with Weibull distribution is specified as:
#' \deqn{\log(T) = \mu + \gamma' X + \sigma \epsilon}
#'
#' Where:
#' - T is the survival time
#' - μ is the intercept
#' - γ contains the covariate effects
#' - X includes treatment, covariates, and treatment×subgroup interaction
#' - σ is the scale parameter
#' - ε follows an extreme value distribution
#'
#' ## Interaction Term
#'
#' For standard models: creates a SINGLE interaction term representing the treatment 
#' effect modification when ALL subgroup conditions are simultaneously satisfied.
#' 
#' For spline model: creates smooth interaction with the biomarker z through spline terms.
#'
#' @return An object of class `c("aft_dgm_flex", "list")` containing:
#' \describe{
#'   \item{df_super}{Data frame with the super population including all covariates,
#'     linear predictors, and potential outcomes}
#'   \item{model_params}{List containing model parameters:
#'     \describe{
#'       \item{mu}{Intercept from AFT model}
#'       \item{sigma}{Scale parameter}
#'       \item{gamma}{Vector of regression coefficients}
#'       \item{b_weibull}{Weibull parameterization coefficients}
#'       \item{b0_weibull}{Weibull baseline hazard coefficients}
#'       \item{censoring}{Censoring distribution parameters}
#'     }}
#'   \item{subgroup_info}{List with subgroup information (for standard models):
#'     \describe{
#'       \item{vars}{Variables used to define subgroup}
#'       \item{cuts}{Cutpoint specifications used}
#'       \item{definitions}{Human-readable subgroup definitions}
#'       \item{size}{Number of observations in subgroup}
#'       \item{proportion}{Proportion of observations in subgroup}
#'     }}
#'   \item{spline_info}{List with spline information (for spline model):
#'     \describe{
#'       \item{knot}{Knot location}
#'       \item{zeta}{Evaluation point}
#'       \item{log_hrs}{Log hazard ratios at key points}
#'       \item{biomarker_var}{Variable used as biomarker}
#'     }}
#'   \item{hazard_ratios}{List of true hazard ratios:
#'     \describe{
#'       \item{overall}{Overall treatment HR}
#'       \item{harm_subgroup}{HR within subgroup (if model="alt")}
#'       \item{no_harm_subgroup}{HR outside subgroup (if model="alt")}
#'     }}
#'   \item{analysis_vars}{List of variable classifications for analysis}
#'   \item{model_type}{Character: "alt", "null", or "spline"}
#'   \item{n_super}{Size of super population}
#'   \item{seed}{Random seed used}
#' }
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data(gbsg)
#'
#' # Example 1: Simple fixed cutpoints (standard model)
#' dgm1 <- generate_aft_dgm_flex(
#'   data = gbsg,
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   subgroup_vars = c("er", "meno"),
#'   subgroup_cuts = list(
#'     er = 20,         # Fixed value
#'     meno = 0         # Factor level
#'   ),
#'   model = "alt",
#'   verbose = TRUE
#' )
#'
#' # Example 2: Spline model
#' dgm2 <- generate_aft_dgm_flex(
#'   data = gbsg,
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   model = "spline",
#'   spline_params = list(
#'     knot = 5,
#'     zeta = 10,
#'     log_hrs = log(c(0.5, 0.75, 1.0)),
#'     biomarker_var = "nodes"
#'   ),
#'   verbose = TRUE
#' )
#' }
#'
#' @author Enhanced by Claude
#' @export
#' @importFrom survival survreg coxph Surv
#' @importFrom stats quantile median uniroot rexp runif rnorm rbinom model.matrix coef

generate_aft_dgm_flex <- function(data,
                                 continuous_vars,
                                 factor_vars,
                                 outcome_var,
                                 event_var,
                                 treatment_var = NULL,
                                 subgroup_vars = NULL,
                                 subgroup_cuts = NULL,
                                 draw_treatment = FALSE,
                                 model = "alt",
                                 spline_params = list(
                                   knot = 5,
                                   zeta = 10,
                                   log_hrs = log(c(0.75, 0.75, 0.75)),
                                   biomarker_var = NULL
                                 ),
                                 k_treat = 1,
                                 k_inter = 1,
                                 n_super = 5000,
                                 cens_type = "weibull",
                                 cens_params = list(),
                                 seed = 8316951,
                                 verbose = TRUE, 
                                 standardize = FALSE) {

  # ============================================================================
  # Helper function to process cutpoint specifications
  # ============================================================================

  process_cutpoint <- function(var_data, cut_spec, var_name = "") {
    # If cut_spec is a simple numeric value, treat as fixed cutpoint
    if (is.numeric(cut_spec) && length(cut_spec) == 1) {
      return(list(
        indicator = var_data <= cut_spec,
        definition = paste0(var_name, " <= ", cut_spec),
        type = "fixed",
        value = cut_spec
      ))
    }
    
    # If cut_spec is a list, process based on type
    if (is.list(cut_spec)) {
      cut_type <- cut_spec$type
      
      if (cut_type == "quantile") {
        # Quantile-based cutpoint
        q_val <- quantile(var_data, cut_spec$value, na.rm = TRUE)
        return(list(
          indicator = var_data <= q_val,
          definition = paste0(var_name, " <= ", round(q_val, 2),
                            " (", cut_spec$value * 100, "th percentile)"),
          type = "quantile",
          value = q_val,
          quantile = cut_spec$value
        ))
        
      } else if (cut_type == "function") {
        # Function-based cutpoint (e.g., median, mean)
        fun_val <- cut_spec$fun(var_data, na.rm = TRUE)
        fun_name <- deparse(substitute(cut_spec$fun))
        if (grepl("median", fun_name)) fun_name <- "median"
        if (grepl("mean", fun_name)) fun_name <- "mean"
        
        return(list(
          indicator = var_data <= fun_val,
          definition = paste0(var_name, " <= ", round(fun_val, 2),
                            " (", fun_name, ")"),
          type = "function",
          value = fun_val,
          fun_name = fun_name
        ))
        
      } else if (cut_type == "range") {
        # Range-based: min <= var <= max
        return(list(
          indicator = var_data >= cut_spec$min & var_data <= cut_spec$max,
          definition = paste0(cut_spec$min, " <= ", var_name, " <= ", cut_spec$max),
          type = "range",
          min = cut_spec$min,
          max = cut_spec$max
        ))
        
      } else if (cut_type == "greater") {
        # Greater than cutpoint
        if (!is.null(cut_spec$value)) {
          # Fixed value
          return(list(
            indicator = var_data > cut_spec$value,
            definition = paste0(var_name, " > ", cut_spec$value),
            type = "greater",
            value = cut_spec$value
          ))
        } else if (!is.null(cut_spec$quantile)) {
          # Quantile-based
          q_val <- quantile(var_data, cut_spec$quantile, na.rm = TRUE)
          return(list(
            indicator = var_data > q_val,
            definition = paste0(var_name, " > ", round(q_val, 2),
                              " (", cut_spec$quantile * 100, "th percentile)"),
            type = "greater_quantile",
            value = q_val,
            quantile = cut_spec$quantile
          ))
        }
        
      } else if (cut_type == "multiple") {
        # Multiple values (useful for categorical)
        vals_str <- paste(cut_spec$values, collapse = ", ")
        return(list(
          indicator = var_data %in% cut_spec$values,
          definition = paste0(var_name, " in {", vals_str, "}"),
          type = "multiple",
          values = cut_spec$values
        ))
        
      } else if (cut_type == "custom") {
        # Custom function that returns logical vector
        indicator <- cut_spec$fun(var_data)
        return(list(
          indicator = indicator,
          definition = paste0(var_name, ": custom function"),
          type = "custom"
        ))
      }
    }
    
    # Default: treat as fixed value if not recognized
    return(list(
      indicator = var_data <= cut_spec,
      definition = paste0(var_name, " <= ", cut_spec),
      type = "fixed",
      value = cut_spec
    ))
  }

  # ============================================================================
  # Input validation
  # ============================================================================
  
  if (verbose) cat("\n=== AFT DGM with Flexible Subgroups ===\n\n")
  
  # Check model type
  if (!model %in% c("alt", "null", "spline")) {
    stop("model must be 'alt', 'null', or 'spline'")
  }
  
  # For spline model, ensure we have necessary parameters
  if (model == "spline") {
    if (is.null(spline_params$knot)) spline_params$knot <- 5
    if (is.null(spline_params$zeta)) spline_params$zeta <- 10
    if (is.null(spline_params$log_hrs)) spline_params$log_hrs <- log(c(0.75, 0.75, 0.75))
    if (is.null(spline_params$biomarker_var)) {
      # Use first continuous variable as biomarker if not specified
      spline_params$biomarker_var <- continuous_vars[1]
      if (verbose) {
        cat("No biomarker specified for spline model.\n")
        cat("Using", spline_params$biomarker_var, "as biomarker variable.\n\n")
      }
    }
  }
  
  # Check data
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }
  
  # Check variable existence
  all_vars <- c(continuous_vars, factor_vars, outcome_var, event_var, treatment_var)
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(paste("Variables not found in data:", paste(missing_vars, collapse = ", ")))
  }
  
  # Check subgroup variables (for non-spline models)
  if (model != "spline" && !is.null(subgroup_vars)) {
    missing_subgroup <- setdiff(subgroup_vars, names(data))
    if (length(missing_subgroup) > 0) {
      stop(paste("Subgroup variables not found in data:",
                paste(missing_subgroup, collapse = ", ")))
    }
  }
  
  # Validate model type with subgroup specification
  if (model == "alt" && (is.null(subgroup_vars) || is.null(subgroup_cuts))) {
    if (!model == "spline") {
      stop("model='alt' requires both subgroup_vars and subgroup_cuts to be specified")
    }
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # ============================================================================
  # SPLINE MODEL BRANCH
  # ============================================================================
  
  if (model == "spline") {
    if (verbose) cat("Using SPLINE model with biomarker:", spline_params$biomarker_var, "\n\n")
    
    # Prepare data similar to get_dgm_spline
    df_work <- data.frame(data)
    
    # Ensure biomarker variable exists and is numeric
    if (!spline_params$biomarker_var %in% names(df_work)) {
      stop(paste("Biomarker variable", spline_params$biomarker_var, "not found in data"))
    }
    
    # Rename biomarker variable to 'z' for consistency
    df_work$z <- df_work[[spline_params$biomarker_var]]
    
    # Handle treatment variable
    if (is.null(treatment_var)) {
      if (verbose) cat("No treatment variable specified. Generating random 50/50 allocation.\n\n")
      df_work$treat <- rbinom(nrow(df_work), 1, 0.5)
    } else {
      df_work$treat <- df_work[[treatment_var]]
    }
    
    # Standardize continuous variables if requested (excluding z which is handled separately)
    if (standardize) {
      for (var in setdiff(continuous_vars, spline_params$biomarker_var)) {
        if (var %in% names(df_work)) {
          df_work[[var]] <- scale(df_work[[var]])[,1]
        }
      }
    }
    
    # Create spline terms
    knot <- spline_params$knot
    df_work$z.treat <- df_work$z * df_work$treat
    df_work$z.k <- (df_work$z - knot) * ifelse(df_work$z > knot, 1, 0)
    df_work$z.k.treat <- df_work$z.k * df_work$treat
    
    # Build formula for spline model
    weib_formula <- as.formula(paste0("Surv(", outcome_var, ", ", event_var, 
                                      ") ~ treat + z + z.treat + z.k + z.k.treat"))
    
    # Fit the model
    fit_weibk <- survreg(weib_formula, dist = 'weibull', data = df_work)
    
    if (verbose) {
      cat("Spline model fitted successfully.\n")
      cat("Knot location:", knot, "\n")
      cat("Model coefficients:\n")
      print(summary(fit_weibk))
      cat("\n")
    }
    
    # Censoring model
    cens_formula <- as.formula(paste0("Surv(", outcome_var, ", 1 - ", event_var, ") ~ 1"))
    fitC_weib <- survreg(cens_formula, dist = 'weibull', data = df_work)
    tauC <- fitC_weib$scale
    muC <- coef(fitC_weib)[1]
    
    # Extract parameters
    mu <- coef(fit_weibk)[1]
    gamma <- coef(fit_weibk)[-1]
    tau <- fit_weibk$scale
    
    # Apply log hazard ratio constraints
    log_hrs <- spline_params$log_hrs
    loghr_0 <- log_hrs[1]
    loghr_knot <- log_hrs[2]
    loghr_zeta <- log_hrs[3]
    zeta <- spline_params$zeta
    
    # Calculate Weibull hazard-ratio parameters
    b0 <- -gamma / tau
    
    # Solve for b1, b3, b5 to satisfy log(hrs) pattern
    b0[1] <- loghr_0  # treat effect at z=0
    b0[3] <- (loghr_knot - b0[1]) / knot  # z.treat effect
    b0[5] <- (loghr_zeta - b0[1] - zeta * b0[3]) / (zeta - knot)  # z.k.treat effect
    
    # Re-define gamma on AFT (log(T)) parameterization
    gamma_true <- -b0 * tau
    
    # Create super population with spline structure
    n_rep <- ceiling(n_super / nrow(df_work))
    df_super_list <- list()
    
    for (i in 1:n_rep) {
      df_temp <- df_work
      df_temp$rep_id <- i
      df_super_list[[i]] <- df_temp
    }
    
    df_super <- do.call(rbind, df_super_list)
    df_super <- df_super[1:n_super, ]
    
    # Store original treatment if not drawing new
    if (!draw_treatment) {
      df_super$treat_orig <- df_super$treat
    }
    
    # Sort by biomarker
    df_super <- df_super[order(df_super$z), ]
    
    # Calculate theoretical hazard ratios for visualization
    if (verbose) {
      zz <- seq(0, max(df_super$z), by = 1)
      loghr_zz <- b0[1] + b0[3] * zz + b0[5] * (zz - knot) * ifelse(zz > knot, 1, 0)
      
      cat("\nTheoretical log hazard ratios at key points:\n")
      cat("  At z=0:", loghr_0, "\n")
      cat("  At z=", knot, ":", loghr_knot, "\n")
      cat("  At z=", zeta, ":", loghr_zeta, "\n\n")
    }
    
    # Create output structure similar to standard model
    result <- list(
      df_super = df_super,
      model_params = list(
        mu = mu,
        sigma = tau,
        gamma = gamma_true,
        b_weibull = b0,
        b0_weibull = b0,
        censoring = list(
          type = "weibull",
          mu = muC,
          tau = tauC
        )
      ),
      spline_info = list(
        knot = knot,
        zeta = zeta,
        log_hrs = log_hrs,
        biomarker_var = spline_params$biomarker_var
      ),
      analysis_vars = list(
        continuous = continuous_vars,
        factor = factor_vars,
        outcome = outcome_var,
        event = event_var,
        treatment = treatment_var
      ),
      model_type = "spline",
      n_super = n_super,
      seed = seed
    )
    
    class(result) <- c("aft_dgm_flex", "list")
    
    if (verbose) {
      cat("\n=== DGM Creation Complete (Spline Model) ===\n")
      cat("Super population size:", nrow(df_super), "\n")
      cat("Biomarker range: [", min(df_super$z), ",", max(df_super$z), "]\n")
      cat("Knot location:", knot, "\n")
    }
    
    return(result)
  }
  
  # ============================================================================
  # STANDARD MODEL BRANCH (original logic for alt/null models)
  # ============================================================================
  
  # Data preparation
  df_analysis <- data.frame(data)
  n_obs <- nrow(df_analysis)
  
  if (verbose) {
    cat("Data dimensions:", n_obs, "observations\n")
    cat("Continuous variables:", paste(continuous_vars, collapse = ", "), "\n")
    cat("Factor variables:", paste(factor_vars, collapse = ", "), "\n")
    if (!is.null(subgroup_vars)) {
      cat("Subgroup variables:", paste(subgroup_vars, collapse = ", "), "\n")
    }
    cat("\n")
  }
  
  # ============================================================================
  # Process treatment variable
  # ============================================================================
  
  if (is.null(treatment_var)) {
    if (verbose) cat("No treatment variable specified. Generating random 50/50 allocation.\n\n")
    df_analysis$treat_sim <- rbinom(n_obs, 1, 0.5)
    treatment_var_use <- "treat_sim"
  } else {
    treatment_var_use <- treatment_var
  }
  
  # ============================================================================
  # Process subgroups if model is "alt"
  # ============================================================================
  
  if (model == "alt" && !is.null(subgroup_vars) && !is.null(subgroup_cuts)) {
    
    if (verbose) cat("Processing subgroup definitions...\n")
    
    # Initialize subgroup indicator
    df_analysis$in_subgroup <- rep(TRUE, n_obs)
    
    # Store processed cutpoint information
    cutpoint_info <- list()
    subgroup_definitions <- character()
    
    # Process each subgroup variable
    for (i in seq_along(subgroup_vars)) {
      var_name <- subgroup_vars[i]
      var_data <- df_analysis[[var_name]]
      
      # Get cutpoint specification for this variable
      if (!is.null(subgroup_cuts[[var_name]])) {
        cut_spec <- subgroup_cuts[[var_name]]
      } else {
        warning(paste("No cutpoint specification found for", var_name, 
                     ". Using median as default."))
        cut_spec <- list(type = "function", fun = median)
      }
      
      # Process the cutpoint
      cut_result <- process_cutpoint(var_data, cut_spec, var_name)
      
      # Update subgroup membership (AND operation - must satisfy ALL conditions)
      df_analysis$in_subgroup <- df_analysis$in_subgroup & cut_result$indicator
      
      # Store information
      cutpoint_info[[var_name]] <- cut_result
      subgroup_definitions <- c(subgroup_definitions, cut_result$definition)
    }
    
    # Create interaction term
    df_analysis$treat_x_subgroup <- df_analysis[[treatment_var_use]] * df_analysis$in_subgroup
    
    n_in_subgroup <- sum(df_analysis$in_subgroup)
    prop_in_subgroup <- mean(df_analysis$in_subgroup)
    
    if (verbose) {
      cat("\nSubgroup definitions:\n")
      for (def in subgroup_definitions) {
        cat("  -", def, "\n")
      }
      cat("\nSubgroup size:", n_in_subgroup, "observations")
      cat(" (", round(prop_in_subgroup * 100, 1), "% of data)\n\n")
    }
    
    # Check if subgroup is too small or too large
    if (prop_in_subgroup < 0.05) {
      warning("Subgroup contains less than 5% of observations. Model may be unstable.")
    }
    if (prop_in_subgroup > 0.95) {
      warning("Subgroup contains more than 95% of observations. Limited heterogeneity.")
    }
    
  } else {
    # No subgroup for null model
    df_analysis$in_subgroup <- rep(FALSE, n_obs)
    df_analysis$treat_x_subgroup <- rep(0, n_obs)
    n_in_subgroup <- 0
    prop_in_subgroup <- 0
    cutpoint_info <- NULL
    subgroup_definitions <- NULL
  }
  
  # ============================================================================
  # Prepare covariates
  # ============================================================================
  
  # Standardize continuous variables if requested
  continuous_data <- df_analysis[, continuous_vars, drop = FALSE]
  if (standardize) {
    continuous_data <- scale(continuous_data)
  }
  
  # Convert factors to dummy variables (using largest category as reference)
  factor_dummies <- NULL
  factor_dummy_names <- character()
  
  if (length(factor_vars) > 0) {
    for (var in factor_vars) {
      fvar <- factor(df_analysis[[var]])
      levels_var <- levels(fvar)
      
      if (length(levels_var) > 1) {
        # Create contrasts with largest category as reference
        tab <- table(fvar)
        ref_level <- names(tab)[which.max(tab)]
        fvar <- relevel(fvar, ref = ref_level)
        
        # Create dummy variables
        dummies <- model.matrix(~ fvar - 1)[, -1, drop = FALSE]
        colnames(dummies) <- paste0(var, "_", levels_var[-which(levels_var == ref_level)])
        
        if (is.null(factor_dummies)) {
          factor_dummies <- dummies
        } else {
          factor_dummies <- cbind(factor_dummies, dummies)
        }
        
        factor_dummy_names <- c(factor_dummy_names, colnames(dummies))
      }
    }
  }
  
  # Combine all covariates
  if (!is.null(factor_dummies)) {
    X_covariates <- cbind(continuous_data, factor_dummies)
  } else {
    X_covariates <- continuous_data
  }
  
  # ============================================================================
  # Fit AFT model
  # ============================================================================
  
  if (verbose) cat("Fitting AFT Weibull model...\n")
  
  # Build model matrix
  if (model == "alt") {
    X_model <- cbind(
      treat = df_analysis[[treatment_var_use]],
      X_covariates,
      treat_x_subgroup = df_analysis$treat_x_subgroup
    )
  } else {
    X_model <- cbind(
      treat = df_analysis[[treatment_var_use]],
      X_covariates
    )
  }
  
  # Create survival object
  surv_obj <- Surv(time = df_analysis[[outcome_var]], 
                   event = df_analysis[[event_var]])
  
  # Fit model
  aft_fit <- survreg(surv_obj ~ X_model, dist = "weibull")
  
  # Extract parameters
  mu_aft <- aft_fit$coefficients[1]  # Intercept
  gamma_aft <- aft_fit$coefficients[-1]  # Other coefficients
  sigma_aft <- aft_fit$scale
  
  if (verbose) {
    cat("AFT model fitted successfully.\n")
    cat("Scale parameter (sigma):", round(sigma_aft, 4), "\n\n")
  }
  
  # Apply k_treat and k_inter modifications
  if ("treat" %in% names(gamma_aft)) {
    gamma_aft["treat"] <- gamma_aft["treat"] * k_treat
    if (verbose && k_treat != 1) {
      cat("Treatment effect modified by factor:", k_treat, "\n")
    }
  }
  
  if ("treat_x_subgroup" %in% names(gamma_aft) && model == "alt") {
    gamma_aft["treat_x_subgroup"] <- gamma_aft["treat_x_subgroup"] * k_inter
    if (verbose && k_inter != 1) {
      cat("Interaction effect modified by factor:", k_inter, "\n")
    }
  }
  
  # ============================================================================
  # Fit censoring model
  # ============================================================================
  
  if (verbose) cat("Fitting censoring model...\n")
  
  if (cens_type == "weibull") {
    # Weibull censoring - estimate from observed censoring pattern
    cens_fit <- survreg(Surv(df_analysis[[outcome_var]], 
                            1 - df_analysis[[event_var]]) ~ 1, 
                       dist = "weibull")
    
    cens_params_final <- list(
      type = "weibull",
      mu = cens_fit$coefficients[1],
      sigma = cens_fit$scale
    )
    
  } else if (cens_type == "uniform") {
    # Uniform censoring
    if (is.null(cens_params$min) || is.null(cens_params$max)) {
      # Use data range as default
      cens_params$min <- min(df_analysis[[outcome_var]])
      cens_params$max <- max(df_analysis[[outcome_var]])
    }
    
    cens_params_final <- list(
      type = "uniform",
      min = cens_params$min,
      max = cens_params$max
    )
  }
  
  # ============================================================================
  # Calculate hazard ratios
  # ============================================================================
  
  # Convert AFT coefficients to Weibull parameterization
  # HR = exp(-gamma/sigma) for Weibull AFT
  
  b_weibull <- -gamma_aft / sigma_aft
  
  # Overall treatment HR
  hr_overall <- exp(b_weibull["treat"])
  
  if (model == "alt" && !is.null(subgroup_vars)) {
    # HR in harm subgroup (treatment + interaction)
    hr_harm_subgroup <- exp(b_weibull["treat"] + b_weibull["treat_x_subgroup"])
    
    # HR in no-harm subgroup (treatment only)
    hr_no_harm_subgroup <- exp(b_weibull["treat"])
    
    if (verbose) {
      cat("\n=== Hazard Ratios ===\n")
      cat("Overall treatment HR:", round(hr_overall, 3), "\n")
      cat("HR in specified subgroup:", round(hr_harm_subgroup, 3), "\n")
      cat("HR outside subgroup:", round(hr_no_harm_subgroup, 3), "\n\n")
    }
  } else {
    hr_harm_subgroup <- NULL
    hr_no_harm_subgroup <- NULL
    
    if (verbose) {
      cat("\nOverall treatment HR:", round(hr_overall, 3), "\n\n")
    }
  }
  
  # ============================================================================
  # Create super population
  # ============================================================================
  
  if (verbose) cat("Creating super population of size", n_super, "...\n")
  
  # Replicate data to create super population
  n_reps <- ceiling(n_super / n_obs)
  df_super <- df_analysis[rep(seq_len(n_obs), n_reps), ]
  df_super <- df_super[seq_len(n_super), ]
  rownames(df_super) <- NULL
  
  # Add original IDs
  df_super$super_id <- seq_len(n_super)
  df_super$orig_id <- ((seq_len(n_super) - 1) %% n_obs) + 1
  
  # Recreate model matrix for super population
  X_super <- df_super[, c(continuous_vars, factor_vars), drop = FALSE]
  
  # Standardize if needed
  if (standardize) {
    X_super[, continuous_vars] <- scale(X_super[, continuous_vars])
  }
  
  # Store treatment assignment strategy
  if (!draw_treatment && !is.null(treatment_var)) {
    df_super$treat_orig <- df_super[[treatment_var]]
  }
  
  # ============================================================================
  # Store linear predictors for potential outcomes
  # ============================================================================
  
  # Prepare covariate matrix
  continuous_super <- X_super[, continuous_vars, drop = FALSE]
  
  # Recreate factor dummies for super population
  if (length(factor_vars) > 0 && !is.null(factor_dummies)) {
    factor_dummies_super <- NULL
    
    for (var in factor_vars) {
      fvar <- factor(df_super[[var]])
      levels_var <- levels(fvar)
      
      if (length(levels_var) > 1) {
        tab <- table(df_analysis[[var]])  # Use original data for reference
        ref_level <- names(tab)[which.max(tab)]
        fvar <- relevel(fvar, ref = ref_level)
        
        dummies <- model.matrix(~ fvar - 1)[, -1, drop = FALSE]
        colnames(dummies) <- paste0(var, "_", levels_var[-which(levels_var == ref_level)])
        
        if (is.null(factor_dummies_super)) {
          factor_dummies_super <- dummies
        } else {
          factor_dummies_super <- cbind(factor_dummies_super, dummies)
        }
      }
    }
    
    X_covariates_super <- cbind(continuous_super, factor_dummies_super)
  } else {
    X_covariates_super <- continuous_super
  }
  
  # Store in super population dataframe
  df_super <- cbind(df_super, X_covariates_super)
  
  # Calculate linear predictors for treated (Z=1) and control (Z=0)
  # These will be used in simulation
  
  # For Z=1 (treated)
  lp_treat_1 <- mu_aft + gamma_aft["treat"]
  
  # For Z=0 (control)
  lp_treat_0 <- mu_aft
  
  # Add covariate contributions
  X_contrib <- as.matrix(X_covariates_super) %*% 
                gamma_aft[colnames(X_covariates_super)]
  
  df_super$lp_base <- mu_aft + X_contrib
  
  if (model == "alt" && !is.null(subgroup_vars)) {
    df_super$lp_treat_1_in <- df_super$lp_base + gamma_aft["treat"] + 
                              gamma_aft["treat_x_subgroup"]
    df_super$lp_treat_1_out <- df_super$lp_base + gamma_aft["treat"]
    df_super$lp_treat_0 <- df_super$lp_base
  } else {
    df_super$lp_treat_1 <- df_super$lp_base + gamma_aft["treat"]
    df_super$lp_treat_0 <- df_super$lp_base
  }
  
  # ============================================================================
  # Create output object
  # ============================================================================
  
  result <- list(
    df_super = df_super,
    model_params = list(
      mu = mu_aft,
      sigma = sigma_aft,
      gamma = gamma_aft,
      b_weibull = b_weibull,
      b0_weibull = -gamma_aft / sigma_aft,
      censoring = cens_params_final
    ),
    subgroup_info = if (model == "alt" && !is.null(subgroup_vars)) {
      list(
        vars = subgroup_vars,
        cuts = subgroup_cuts,
        cutpoint_info = cutpoint_info,
        definitions = subgroup_definitions,
        size = n_in_subgroup,
        proportion = prop_in_subgroup
      )
    } else {
      NULL
    },
    hazard_ratios = list(
      overall = hr_overall,
      harm_subgroup = hr_harm_subgroup,
      no_harm_subgroup = hr_no_harm_subgroup
    ),
    analysis_vars = list(
      continuous = continuous_vars,
      factor = factor_vars,
      factor_dummies = factor_dummy_names,
      outcome = outcome_var,
      event = event_var,
      treatment = treatment_var
    ),
    model_type = model,
    n_super = n_super,
    seed = seed
  )
  
  class(result) <- c("aft_dgm_flex", "list")
  
  # ============================================================================
  # Print summary if verbose
  # ============================================================================
  
  if (verbose) {
    cat("\n=== DGM Creation Complete ===\n")
    cat("Model type:", model, "\n")
    cat("Super population size:", nrow(df_super), "\n")
    if (model == "alt" && !is.null(subgroup_vars)) {
      cat("Subgroup proportion:", round(prop_in_subgroup * 100, 1), "%\n")
    }
    cat("Censoring type:", cens_type, "\n")
    cat("\n")
  }
  
  return(result)
}

# ============================================================================
# Print method for aft_dgm_flex objects
# ============================================================================

print.aft_dgm_flex <- function(x, ...) {
  cat("\n=== AFT Data Generating Mechanism ===\n\n")
  
  cat("Model Type:", x$model_type, "\n")
  cat("Super Population Size:", x$n_super, "\n")
  cat("Random Seed:", x$seed, "\n\n")
  
  if (x$model_type == "spline") {
    cat("--- Spline Model Parameters ---\n")
    cat("Biomarker Variable:", x$spline_info$biomarker_var, "\n")
    cat("Knot Location:", x$spline_info$knot, "\n")
    cat("Evaluation Point (zeta):", x$spline_info$zeta, "\n")
    cat("Log Hazard Ratios:\n")
    cat("  At z=0:", x$spline_info$log_hrs[1], "\n")
    cat("  At knot:", x$spline_info$log_hrs[2], "\n")
    cat("  At zeta:", x$spline_info$log_hrs[3], "\n\n")
  } else if (!is.null(x$subgroup_info)) {
    cat("--- Subgroup Information ---\n")
    cat("Variables:", paste(x$subgroup_info$vars, collapse = ", "), "\n")
    cat("Definitions:\n")
    for (def in x$subgroup_info$definitions) {
      cat("  -", def, "\n")
    }
    cat("Size:", x$subgroup_info$size, "observations\n")
    cat("Proportion:", round(x$subgroup_info$proportion * 100, 2), "%\n\n")
  }
  
  cat("--- Analysis Variables ---\n")
  cat("Continuous:", paste(x$analysis_vars$continuous, collapse = ", "), "\n")
  cat("Factor:", paste(x$analysis_vars$factor, collapse = ", "), "\n")
  cat("Outcome:", x$analysis_vars$outcome, "\n")
  cat("Event:", x$analysis_vars$event, "\n")
  cat("Treatment:", ifelse(is.null(x$analysis_vars$treatment), 
                          "Simulated", x$analysis_vars$treatment), "\n\n")
  
  if (!is.null(x$hazard_ratios$overall)) {
    cat("--- Hazard Ratios ---\n")
    cat("Overall Treatment HR:", round(x$hazard_ratios$overall, 3), "\n")
    
    if (!is.null(x$hazard_ratios$harm_subgroup)) {
      cat("HR in Subgroup:", round(x$hazard_ratios$harm_subgroup, 3), "\n")
      cat("HR outside Subgroup:", round(x$hazard_ratios$no_harm_subgroup, 3), "\n")
    }
  }
  
  cat("\n--- Model Parameters ---\n")
  cat("Intercept (mu):", round(x$model_params$mu, 4), "\n")
  cat("Scale (sigma):", round(x$model_params$sigma, 4), "\n")
  cat("Censoring Type:", x$model_params$censoring$type, "\n")
  
  invisible(x)
}