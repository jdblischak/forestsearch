#' =============================================================================
#' Improved MRCT Simulation Code: Simsetup Chunk
#' 
#' This is an improved version of the simulation setup code from mrct_legacy.qmd,
#' utilizing the modern patterns and functions from mrct_recent.qmd including:
#' - forestsearch and weightedsurv packages
#' - generate_aft_dgm_flex() for flexible DGM creation
#' - simulate_from_dgm() for data generation
#' - Improved parallel processing patterns
#' - Better error handling and progress reporting
#' =============================================================================

# Required packages
library(survival)
library(data.table)
library(grf)
library(doFuture)
library(doRNG)
library(progressr)
library(forestsearch)
library(weightedsurv)
library(gt)
library(ggplot2)
library(table1)

#' -----------------------------------------------------------------------------
#' mrct_APregion_sims: Main simulation function for MRCT subgroup analysis
#' 
#' Simulates n_sims trials from a specified data generating mechanism (DGM),
#' identifies subgroups based on non-AP data, and evaluates estimates in AP data.
#' 
#' @param dgm Data generating mechanism object from generate_aft_dgm_flex()
#' @param n_sims Number of simulations to run
#' @param n_sample Sample size per simulation (default uses dgm's super-population size)
#' @param wname Name of the region indicator variable (default: "z_AP")
#' @param bw Log hazard ratio for prognostic effect of region (default: -log(5))
#' @param sg_focus Subgroup selection criterion: "minSG", "hr", or "maxSG"
#' @param maxk Maximum number of factors in subgroup combinations (1 or 2)
#' @param hr_threshold Hazard ratio threshold for subgroup identification
#' @param hr_consistency Consistency threshold for hazard ratio
#' @param pconsistency_threshold Probability threshold for consistency
#' @param confounders_name Vector of confounder variable names
#' @param conf_force Vector of forced cuts to consider
#' @param analysis_time Time of analysis for administrative censoring
#' @param cens_adjust Adjustment factor for censoring rate
#' @param parallel_args List of parallel processing arguments
#' @param details Logical: print detailed progress information
#' @param seed Base random seed for reproducibility
#' 
#' @return data.table with simulation results including HR estimates and subgroup info
#' -----------------------------------------------------------------------------
mrct_APregion_sims <- function(
    dgm,
    n_sims,
    n_sample = NULL,
    wname = "z_AP",
    bw = -log(5),
    sg_focus = "minSG",
    maxk = 1,
    hr_threshold = 0.90,
    hr_consistency = 0.80,
    pconsistency_threshold = 0.90,
    confounders_name = NULL,
    conf_force = NULL,
    analysis_time = 60,
    cens_adjust = 0,
    parallel_args = list(plan = "multisession", workers = NULL, show_message = TRUE),
    details = FALSE,
    seed = NULL
) {
  
  t_start <- proc.time()[3]
  
  # Set default sample size from dgm if not provided
  if (is.null(n_sample)) {
    n_sample <- nrow(dgm$df_super)
  }
  
  # Set default confounders from dgm if not provided
  if (is.null(confounders_name)) {
    # Extract covariate names from dgm, excluding outcome-related variables
    all_vars <- names(dgm$df_super)
    exclude_vars <- c("id", "treat", "tte", "event", "y_sim", "treat_sim", 
                      "event_sim", "t_true", "c_time", "lin_pred_0", "lin_pred_1",
                      "loghr_po", "entrytime", "flag_harm")
    confounders_name <- setdiff(all_vars[grepl("^z_|^ecog|^strat", all_vars)], exclude_vars)
  }
  
  # Set default forced cuts if not provided
  if (is.null(conf_force)) {
    # Create default cuts based on dgm variables
    conf_force <- c("z_age <= 65", "z_bm <= 0", "z_bm <= 1", "z_bm <= 2", "z_bm <= 5")
  }
  
  # Set up parallel processing
  setup_parallel(
    plan = parallel_args$plan %||% "multisession",
    workers = parallel_args$workers,
    show_message = parallel_args$show_message %||% TRUE
  )
  
  # Report backend info
  report_parallel_backend()
  
  # Set up progress reporting
  handlers(global = TRUE)
  handlers("progress")
  
  # Run simulations with doFuture
  results <- with_progress({
    p <- progressor(along = seq_len(n_sims))
    
    foreach(
      sim = seq_len(n_sims),
      .options.future = list(seed = TRUE),
      .combine = "rbind",
      .errorhandling = "pass",
      .packages = c("survival", "data.table", "forestsearch", "weightedsurv")
    ) %dofuture% {
      
      # Update progress
      p(sprintf("Simulation %d/%d", sim, n_sims))
      
      # Simulate data from DGM
      dfs <- simulate_from_dgm(
        dgm = dgm,
        n = n_sample,
        rand_ratio = 1,
        draw_treatment = TRUE,
        entry_var = if ("entrytime" %in% names(dgm$df_super)) "entrytime" else NULL,
        analysis_time = analysis_time,
        cens_adjust = cens_adjust,
        seed = if (!is.null(seed)) seed + sim else NULL
      )
      
      # Ensure data.table format
      dfs <- as.data.table(dfs)
      
      # Split into non-AP (training) and AP (testing) populations
      df_nonAP <- dfs[get(wname) == 0]
      df_AP <- dfs[get(wname) == 1]
      
      # Initialize results
      n_test <- nrow(df_AP)
      n_train <- nrow(df_nonAP)
      n_itt <- nrow(dfs)
      
      # Safe Cox model fitting function
      safe_coxph <- function(formula, data) {
        fit <- tryCatch(
          summary(coxph(formula, data = data))$conf.int,
          error = function(e) NULL
        )
        if (!is.null(fit) && is.numeric(fit[1])) return(fit[1])
        return(NA_real_)
      }
      
      # Estimate HR in training (non-AP)
      hr_train <- safe_coxph(Surv(y_sim, event_sim) ~ treat_sim, df_nonAP)
      
      # Estimate HR in ITT (stratified by randomization if available)
      if ("strat" %in% names(dfs)) {
        hr_itt <- safe_coxph(Surv(y_sim, event_sim) ~ treat_sim + strata(strat), dfs)
      } else {
        hr_itt <- safe_coxph(Surv(y_sim, event_sim) ~ treat_sim, dfs)
      }
      
      # Estimate HR in ITT (stratified by region)
      hr_ittX <- safe_coxph(as.formula(paste0("Surv(y_sim, event_sim) ~ treat_sim + strata(", wname, ")")), dfs)
      
      # Initialize testing/subgroup results
      hr_test <- NA_real_
      any_found <- 0
      sg_found <- "none"
      n_sg <- n_test
      hr_sg <- NA_real_
      prev_sg <- 1.0
      POhr_sg <- NA_real_
      
      # Estimate HR in AP (testing)
      hr_test <- safe_coxph(Surv(y_sim, event_sim) ~ treat_sim, df_AP)
      
      # Only proceed with subgroup identification if AP estimates are valid
      if (!is.na(hr_test)) {
        
        # Potential outcome HR for AP population
        if ("loghr_po" %in% names(df_AP)) {
          POhr_sg <- exp(mean(df_AP$loghr_po, na.rm = TRUE))
        }
        
        # Run forestsearch for subgroup identification
        fs_result <- tryCatch({
          forestsearch(
            df.analysis = as.data.frame(df_nonAP),
            df.test = as.data.frame(df_AP),
            confounders.name = confounders_name,
            outcome.name = "y_sim",
            treat.name = "treat_sim",
            event.name = "event_sim",
            id.name = if ("id" %in% names(df_nonAP)) "id" else NULL,
            potentialOutcome.name = if ("loghr_po" %in% names(df_nonAP)) "loghr_po" else NULL,
            hr.threshold = hr_threshold,
            hr.consistency = hr_consistency,
            pconsistency.threshold = pconsistency_threshold,
            sg_focus = sg_focus,
            conf_force = conf_force,
            maxk = maxk,
            n.min = 60,
            d0.min = 10,
            d1.min = 10,
            showten_subgroups = FALSE,
            details = FALSE,
            plot.sg = FALSE,
            parallel_args = list(plan = "sequential")  # Already in parallel context
          )
        }, error = function(e) NULL)
        
        # Process forestsearch results
        if (!is.null(fs_result)) {
          if (is.null(fs_result$sg.harm)) {
            # No subgroup found - use AP population as "subgroup"
            any_found <- 0
            sg_found <- "none"
            n_sg <- n_test
            hr_sg <- hr_test
            prev_sg <- 1.0
          } else {
            # Subgroup found
            any_found <- 1
            sg_found <- paste(fs_result$sg.harm, collapse = " & ")
            
            if (details && sim <= 10) {
              message(sprintf("Simulation %d: Subgroup found: %s", sim, sg_found))
            }
            
            # Get subgroup data from testing
            df_test <- fs_result$df.test
            if (!is.null(df_test) && "treat.recommend" %in% names(df_test)) {
              df_sg <- df_test[df_test$treat.recommend == 1, ]
              n_sg <- nrow(df_sg)
              prev_sg <- n_sg / n_test
              
              # Estimate HR in subgroup
              hr_sg <- safe_coxph(Surv(y_sim, event_sim) ~ treat_sim, df_sg)
              
              # Potential outcome HR for subgroup
              if ("loghr_po" %in% names(df_sg)) {
                POhr_sg <- exp(mean(df_sg$loghr_po, na.rm = TRUE))
              }
            }
          }
        }
      }
      
      # Return results as data.table row
      data.table(
        sim = sim,
        n_itt = n_itt,
        hr_itt = hr_itt,
        hr_ittX = hr_ittX,
        n_train = n_train,
        hr_train = hr_train,
        n_test = n_test,
        hr_test = hr_test,
        any_found = any_found,
        sg_found = sg_found,
        n_sg = n_sg,
        hr_sg = hr_sg,
        POhr_sg = POhr_sg,
        prev_sg = prev_sg
      )
    }
  })
  
  # Reset to sequential processing
  plan(sequential)
  
  # Convert to data.table if needed
  results <- as.data.table(results)
  
  # Add hr_sg_null: missing if no subgroup found
  results[, hr_sg_null := fifelse(any_found == 0, NA_real_, hr_sg)]
  
  # Report summary statistics
  t_now <- proc.time()[3]
  t_min <- (t_now - t_start) / 60
  
  if (details) {
    message(sprintf("\n=== Simulation Complete ==="))
    message(sprintf("Total simulations: %d", n_sims))
    message(sprintf("Time elapsed: %.2f minutes", t_min))
    message(sprintf("Projection per 1000 sims: %.2f minutes", t_min * (1000 / n_sims)))
    message(sprintf("Proportion subgroups found: %.3f", mean(results$any_found, na.rm = TRUE)))
    message(sprintf("Mean HR (ITT): %.3f", mean(results$hr_itt, na.rm = TRUE)))
    message(sprintf("Mean HR (Test/AP): %.3f", mean(results$hr_test, na.rm = TRUE)))
    message(sprintf("Mean HR (Subgroup): %.3f", mean(results$hr_sg_null, na.rm = TRUE)))
  }
  
  return(results)
}


#' -----------------------------------------------------------------------------
#' Helper function: Setup parallel processing
#' -----------------------------------------------------------------------------
setup_parallel <- function(plan = "multisession", workers = NULL, show_message = TRUE) {
  
  if (plan == "sequential") {
    plan(sequential)
    return(invisible())
  }
  
  if (is.null(workers)) {
    workers <- max(1, parallel::detectCores() - 1)
  }
  
  if (plan == "multisession") {
    plan(multisession, workers = workers)
  } else if (plan == "callr") {
    plan(future.callr::callr, workers = workers)
  } else if (plan == "multicore") {
    # Only works on Unix-like systems
    if (.Platform$OS.type == "unix") {
      plan(multicore, workers = workers)
    } else {
      plan(multisession, workers = workers)
    }
  }
  
  if (show_message) {
    message(sprintf("Parallel processing: %d workers using %s backend", 
                    workers, class(plan())[1]))
  }
}


#' -----------------------------------------------------------------------------
#' Helper function: Report parallel backend information
#' -----------------------------------------------------------------------------
report_parallel_backend <- function() {
  if (!is(plan(), "sequential")) {
    message(sprintf("Using %d cores with backend %s", 
                    nbrOfWorkers(), 
                    attr(plan("list")[[1]], "class")[2]))
  } else if (foreach::getDoParWorkers() > 1) {
    message(sprintf("Using %d cores with backend %s", 
                    foreach::getDoParWorkers(), 
                    foreach::getDoParName()))
  } else {
    message("Backend uses sequential processing.")
  }
}


#' -----------------------------------------------------------------------------
#' Null coalescing operator
#' -----------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x


#' -----------------------------------------------------------------------------
#' summaryout_mrct: Create summary tables from simulation results
#' 
#' @param pop_summary Population summary object from large sample approximation
#' @param mrct_sims Simulation results from mrct_APregion_sims()
#' @param out_sgs Vector of subgroup classification columns for sg_type=1
#' @param out_sgs2 Vector of subgroup classification columns for sg_type=2
#' @param sg_type Type of subgroup summary (1 or 2)
#' @param tab_caption Caption for the output table
#' @param showtable Logical: print the table
#' 
#' @return List with res (summary statistics) and out_table (formatted table)
#' -----------------------------------------------------------------------------
summaryout_mrct <- function(
    pop_summary,
    mrct_sims,
    out_sgs = c("sg_found", "sg_biomarker", "sg_age"),
    out_sgs2 = c("sg_biomarker", "sg_age", "sg_male", "sg_ecog", "sg_histology", 
                 "sg_CTregimen", "sg_region", "sg_surgery", "sg_prior_treat"),
    out_est = c("+ APflag + sg_le85 + APflag2 + APflag3 + hr_sg_null"),
    sg_type = 1,
    tab_caption = "Identified subgroups and estimation summaries",
    showtable = TRUE
) {
  
  # Ensure data.table
  mrct_sims <- as.data.table(copy(mrct_sims))
  
  res <- list()
  
  # Calculate summary statistics from population
  if (!is.null(pop_summary)) {
    res$ahr_true <- round(pop_summary$AHR %||% NA, 4)
    res$ahr_unadj <- round(pop_summary$ITT_unadj %||% NA, 4)
    res$ahr_sR <- round(pop_summary$ITT_sR %||% NA, 4)
    res$ahr_sRw <- round(pop_summary$ITT_sRw %||% NA, 4)
    
    # Relative biases
    if (!is.na(res$ahr_true) && res$ahr_true != 0) {
      res$bias_unadj <- round(100 * (res$ahr_unadj - res$ahr_true) / res$ahr_true, 1)
      res$bias_sR <- round(100 * (res$ahr_sR - res$ahr_true) / res$ahr_true, 1)
      res$bias_sRw <- round(100 * (res$ahr_sRw - res$ahr_true) / res$ahr_true, 1)
    }
    
    res$ahr_w1 <- round(pop_summary$W_1 %||% NA, 4)
    if (!is.na(res$ahr_true) && res$ahr_true != 0 && !is.na(res$ahr_w1)) {
      res$bias_w1 <- round(100 * (res$ahr_w1 - res$ahr_true) / res$ahr_true, 1)
    }
  }
  
  # Create derived variables for analysis
  mrct_sims[, `:=`(
    APflag = fifelse(hr_test > 0.9, "AP > 0.9", "AP <= 0.9"),
    sg_le85 = fifelse(hr_sg <= 0.85, "AP(sg) <= 0.85", "AP(sg) > 0.85"),
    APflag2 = fifelse(hr_test > 0.9 & hr_sg <= 0.85, "AP > 0.9 & AP(sg) <= 0.85", "Not"),
    APflag3 = fifelse(hr_test > 0.9 & hr_sg <= 0.80, "AP > 0.9 & AP(sg) <= 0.80", "Not"),
    found = as.factor(any_found)
  )]
  
  # Classify subgroups by factor type
  mrct_sims[, `:=`(
    sg_biomarker = fifelse(grepl("bm|biomarker", sg_found, ignore.case = TRUE),
                           "biomarker", fifelse(sg_found != "none", "other", "none")),
    sg_age = fifelse(grepl("age", sg_found, ignore.case = TRUE),
                     "age", fifelse(sg_found != "none", "other", "none")),
    sg_male = fifelse(grepl("male|sex|gender", sg_found, ignore.case = TRUE),
                      "sex", fifelse(sg_found != "none", "other", "none")),
    sg_ecog = fifelse(grepl("ecog", sg_found, ignore.case = TRUE),
                      "ecog", fifelse(sg_found != "none", "other", "none")),
    sg_histology = fifelse(grepl("histology", sg_found, ignore.case = TRUE),
                           "histology", fifelse(sg_found != "none", "other", "none")),
    sg_CTregimen = fifelse(grepl("CTregimen|chemo", sg_found, ignore.case = TRUE),
                           "CT_regimen", fifelse(sg_found != "none", "other", "none")),
    sg_region = fifelse(grepl("region|EU|asia", sg_found, ignore.case = TRUE),
                        "region", fifelse(sg_found != "none", "other", "none")),
    sg_surgery = fifelse(grepl("surgery", sg_found, ignore.case = TRUE),
                         "surgery", fifelse(sg_found != "none", "other", "none")),
    sg_prior_treat = fifelse(grepl("prior_treat", sg_found, ignore.case = TRUE),
                             "prior_treat", fifelse(sg_found != "none", "other", "none"))
  )]
  
  # Build formula for table1
  outwhat1 <- "~ hr_itt + hr_ittX + hr_test + found + prev_sg + hr_sg + POhr_sg +"
  if (sg_type == 1) {
    outwhat2 <- paste(out_sgs, collapse = " + ")
  } else {
    outwhat2 <- paste(out_sgs2, collapse = " + ")
  }
  out_what <- as.formula(paste(outwhat1, outwhat2, out_est))
  
  # Create summary table
  out_table <- table1(out_what, data = mrct_sims, caption = tab_caption)
  
  if (showtable) print(out_table)
  
  return(list(res = res, out_table = out_table, data = mrct_sims))
}


#' -----------------------------------------------------------------------------
#' SGplot_estimates: Create violin/boxplot visualization of HR estimates
#' 
#' @param df Data frame with simulation results
#' @param label_training Label for training data estimates
#' @param label_testing Label for testing data estimates
#' @param label_itt Label for ITT estimates
#' @param label_sg Label for subgroup estimates
#' 
#' @return List with data frame for plotting and ggplot object
#' -----------------------------------------------------------------------------
SGplot_estimates <- function(
    df,
    label_training = "Training",
    label_testing = "Testing",
    label_itt = "ITT (stratified)",
    label_sg = "Testing (subgroup)"
) {
  
  # Prepare data for plotting
  df_itt <- data.table(est = df$hr_itt, analysis = label_itt)
  df_training <- data.table(est = df$hr_train, analysis = label_training)
  df_testing <- data.table(est = df$hr_test, analysis = label_testing)
  df_sg <- data.table(est = df$hr_sg, analysis = label_sg)
  
  hr_estimates <- rbind(df_itt, df_training, df_testing, df_sg)
  
  # Set factor order
  est_order <- c(label_itt, label_training, label_testing, label_sg)
  hr_estimates[, analysis := factor(analysis, levels = est_order)]
  
  # Create plot
  p <- ggplot(hr_estimates, aes(x = analysis, y = est, fill = analysis)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.8) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      x = "Analysis Population",
      y = "Hazard Ratio Estimate",
      title = "Distribution of HR Estimates Across Simulations"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 15, hjust = 1),
      panel.grid.minor = element_blank()
    )
  
  return(list(dfPlot_estimates = hr_estimates, plot_estimates = p))
}


#' -----------------------------------------------------------------------------
#' create_dgm_for_simulations: Wrapper to create DGM for simulation scenarios
#' 
#' @param df_case Case study data frame
#' @param model_type Either "alt" (alternative with HTE) or "null" (uniform effect)
#' @param log_hrs Vector of log hazard ratios for spline specification
#' @param ap_effect Log hazard ratio for AP region prognostic effect
#' @param verbose Print detailed output
#' 
#' @return DGM object for use with simulate_from_dgm()
#' -----------------------------------------------------------------------------
create_dgm_for_simulations <- function(
    df_case,
    model_type = c("alt", "null"),
    log_hrs = NULL,
    ap_effect = -log(5),
    verbose = FALSE
) {
  
  model_type <- match.arg(model_type)
  
  # Set default log hazard ratios based on model type
  if (is.null(log_hrs)) {
    if (model_type == "alt") {
      # Strong biomarker effect: HR varies from 2 (harmful) to 0.5 (protective)
      log_hrs <- log(c(2, 1.25, 0.5))
    } else {
      # Null: uniform effect HR = 0.7
      log_hrs <- log(c(0.7, 0.7, 0.7))
    }
  }
  
  dgm <- generate_aft_dgm_flex(
    data = df_case,
    continuous_vars = c("age", "bm"),
    factor_vars = c("male", "histology", "prior_treat", "AP"),
    set_beta_spec = list(
      set_var = c("z_AP"),
      beta_var = ap_effect
    ),
    continuous_vars_cens = c("age"),
    factor_vars_cens = c("prior_treat"),
    cens_type = "weibull",
    outcome_var = "tte",
    event_var = "event",
    treatment_var = "treat",
    subgroup_vars = NULL,
    subgroup_cuts = NULL,
    model = model_type,
    spline_spec = list(
      var = "z_bm",
      knot = 5,
      zeta = 10,
      log_hrs = log_hrs
    ),
    k_inter = 0.0,
    verbose = verbose,
    standardize = FALSE
  )
  
  return(dgm)
}


#' =============================================================================
#' EXAMPLE USAGE
#' =============================================================================

# # Load case study data
# df_case <- read.table("../data/dfsynthetic.csv", header = TRUE, sep = ",")
# df_case$AP <- df_case$region_asia
# 
# # Create DGM for alternative hypothesis (strong biomarker effect)
# dgm_alt <- create_dgm_for_simulations(
#   df_case = df_case,
#   model_type = "alt",
#   log_hrs = log(c(3, 1.25, 0.50)),
#   verbose = TRUE
# )
# 
# # Create DGM for null hypothesis (uniform effect)
# dgm_null <- create_dgm_for_simulations(
#   df_case = df_case,
#   model_type = "null",
#   verbose = TRUE
# )
# 
# # Define confounders
# confounders_name <- c("z_age", "z_bm", "z_male", "ecog", "z_histology", "z_prior_treat", "strat")
# 
# # Define forced cuts
# conf_force <- c("z_age <= 65", "z_bm <= 0", "z_bm <= 1", "z_bm <= 2", "z_bm <= 5")
# 
# # Run simulations under alternative hypothesis
# results_alt <- mrct_APregion_sims(
#   dgm = dgm_alt,
#   n_sims = 1000,
#   wname = "z_AP",
#   bw = -log(5),
#   sg_focus = "minSG",
#   maxk = 1,
#   hr_threshold = 0.90,
#   hr_consistency = 0.80,
#   pconsistency_threshold = 0.90,
#   confounders_name = confounders_name,
#   conf_force = conf_force,
#   parallel_args = list(plan = "multisession", workers = NULL, show_message = TRUE),
#   details = TRUE
# )
# 
# # Run simulations under null hypothesis
# results_null <- mrct_APregion_sims(
#   dgm = dgm_null,
#   n_sims = 1000,
#   wname = "z_AP",
#   bw = -log(5),
#   sg_focus = "minSG",
#   maxk = 1,
#   hr_threshold = 0.90,
#   hr_consistency = 0.80,
#   pconsistency_threshold = 0.90,
#   confounders_name = confounders_name,
#   conf_force = conf_force,
#   parallel_args = list(plan = "multisession", workers = NULL, show_message = TRUE),
#   details = TRUE
# )
# 
# # Create summary outputs
# summary_alt <- summaryout_mrct(
#   pop_summary = NULL,
#   mrct_sims = results_alt,
#   tab_caption = "Subgroups under strong biomarker effect"
# )
# 
# summary_null <- summaryout_mrct(
#   pop_summary = NULL,
#   mrct_sims = results_null,
#   tab_caption = "Subgroups under null (uniform effect)"
# )
# 
# # Create visualizations
# plot_alt <- SGplot_estimates(
#   results_alt,
#   label_training = "Non-AP, ITT",
#   label_itt = "Overall, ITT",
#   label_testing = "AP, ITT",
#   label_sg = "AP, identified subgroup"
# )
# print(plot_alt$plot_estimates)
