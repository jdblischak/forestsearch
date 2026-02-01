#' @title ForestSearch: Exploratory Subgroup Identification
#'
#' @description
#' Identifies subgroups with differential treatment effects in clinical trials
#' using a combination of Generalized Random Forests (GRF), LASSO variable
#' selection, and exhaustive combinatorial search with split-sample validation.
#'
#' @param df.analysis Data frame. Analysis dataset with required columns.
#' @param outcome.name Character. Name of time-to-event outcome variable. Default "tte".
#' @param event.name Character. Name of event indicator (1=event, 0=censored). Default "event".
#' @param treat.name Character. Name of treatment variable (1=treatment, 0=control). Default "treat".
#' @param id.name Character. Name of subject ID variable. Default "id".
#' @param potentialOutcome.name Character. Name of potential outcome variable (optional).
#' @param flag_harm.name Character. Name of true harm flag for simulations (optional).
#' @param confounders.name Character vector. Names of candidate subgroup-defining variables.
#' @param parallel_args List. Parallel processing configuration:
#'   \describe{
#'     \item{plan}{Character. One of "multisession", "multicore", "callr", "sequential"}
#'     \item{workers}{Integer. Number of parallel workers}
#'     \item{show_message}{Logical. Show parallel setup messages}
#'   }
#' @param df.predict Data frame. Prediction dataset (optional).
#' @param df.test Data frame. Test dataset (optional).
#' @param is.RCT Logical. Is this a randomized controlled trial? Default TRUE.
#' @param seedit Integer. Random seed. Default 8316951.
#' @param est.scale Character. Estimation scale ("hr" or "rmst"). Default "hr".
#' @param use_lasso Logical. Use LASSO for variable selection. Default TRUE.
#' @param use_grf Logical. Use GRF for variable importance. Default TRUE.
#' @param grf_res GRF results object (optional, for reuse).
#' @param grf_cuts List. Custom GRF cut points (optional).
#' @param max_n_confounders Integer. Maximum confounders to consider. Default 1000.
#' @param grf_depth Integer. GRF tree depth. Default 2.
#' @param dmin.grf Integer. Minimum events for GRF. Default 12.
#' @param frac.tau Numeric. Fraction of tau for RMST. Default 0.6.
#' @param return_selected_cuts_only Logical. If TRUE, GRF returns only cuts from the
#'   tree depth that identified the selected subgroup meeting `dmin.grf`. If FALSE
#'   (default), returns all cuts from all fitted trees (depths 1 through `grf_depth`).
#'   See \code{\link{grf.subg.harm.survival}} for details.
#' @param conf_force Character vector. Variables to force include (optional).
#' @param defaultcut_names Character vector. Default cut variable names (optional).
#' @param cut_type Character. Cut type ("default" or "custom"). Default "default".
#' @param exclude_cuts Character vector. Variables to exclude from cutting (optional).
#' @param replace_med_grf Logical. Replace median with GRF cuts. Default FALSE.
#' @param cont.cutoff Integer. Cutoff for continuous vs categorical. Default 4.
#' @param conf.cont_medians Named numeric vector. Median values for continuous variables (optional).
#' @param conf.cont_medians_force Named numeric vector. Forced median values (optional).
#' @param n.min Integer. Minimum subgroup size. Default 60.
#' @param hr.threshold Numeric. Minimum HR for candidate subgroups. Default 1.25.
#' @param hr.consistency Numeric. Minimum HR for consistency validation. Default 1.0.
#' @param sg_focus Character. Subgroup selection focus. One of "hr", "hrMaxSG", "maxSG",
#'   "hrMinSG", "minSG". Default "hr".
#' @param stop_threshold Numeric. Early stopping threshold for consistency
#'   evaluation. When a candidate subgroup's estimated consistency probability
#'   exceeds this threshold, evaluation stops early. Default 0.95.
#'   \strong{Note:} Automatically reset to NULL when \code{sg_focus} is
#'   "hrMaxSG" or "hrMinSG", as these criteria prioritize hazard ratio in
#'   selection and require full evaluation of all candidates.
#' @param fs.splits Integer. Number of splits for consistency evaluation (or maximum
#'   splits when \code{use_twostage = TRUE}). Default 1000.
#' @param m1.threshold Numeric. Maximum median survival threshold. Default Inf.
#' @param pconsistency.threshold Numeric. Minimum consistency proportion. Default 0.90.
#' @param showten_subgroups Logical. Show top 10 subgroups. Default FALSE.
#' @param d0.min Integer. Minimum control arm events. Default 12.
#' @param d1.min Integer. Minimum treatment arm events. Default 12.
#' @param max.minutes Numeric. Maximum search time in minutes. Default 3.
#' @param minp Numeric. Minimum prevalence threshold. Default 0.025.
#' @param details Logical. Print progress details. Default FALSE.
#' @param maxk Integer. Maximum number of factors per subgroup. Default 2.
#' @param by.risk Integer. Risk table interval. Default 12.
#' @param plot.sg Logical. Plot subgroup survival curves. Default FALSE.
#' @param plot.grf Logical. Plot GRF results. Default FALSE.
#' @param max_subgroups_search Integer. Maximum subgroups to evaluate. Default 10.
#' @param vi.grf.min Numeric. Minimum GRF variable importance. Default -0.2.
#' @param use_twostage Logical. Use two-stage sequential consistency algorithm for
#'   improved performance. Default TRUE. When TRUE, \code{fs.splits} becomes the
#'   maximum number of splits, and early stopping is enabled. See Details.
#' @param twostage_args List. Additional arguments for two-stage consistency:
#'   \describe{
#'     \item{n.splits.screen}{Integer. Splits for screening stage (default: 50)}
#'     \item{batch.size}{Integer. Batch size for sequential updates (default: 25)}
#'     \item{conf.level}{Numeric. Confidence level for early stopping (default: 0.95)}
#'   }
#'
#' @return A list of class "forestsearch" containing:
#'   \item{sg.harm}{Character vector. Selected subgroup definition}
#'   \item{grp.consistency}{List. Consistency evaluation results}
#'   \item{find.grps}{List. All candidate subgroups found}
#'   \item{confounders.candidate}{Character. Candidate confounders}
#'   \item{confounders.evaluated}{Character. Evaluated confounders}
#'   \item{df.est}{Data frame. Estimation dataset with subgroup flags}
#'   \item{df.predict}{Data frame. Prediction dataset with recommendations}
#'   \item{df.test}{Data frame. Test dataset with recommendations}
#'   \item{minutes_all}{Numeric. Total computation time}
#'   \item{grf_res}{List. GRF analysis results}
#'   \item{grf_cuts}{Character. GRF-derived cut points}
#'   \item{args_call_all}{List. All function arguments for reproducibility}
#'
#' @details
#' ForestSearch implements a multi-stage algorithm for exploratory subgroup
#' identification:
#'
#' \enumerate{
#'   \item \strong{Variable Selection}: Uses LASSO and/or GRF to identify
#'     promising candidate variables
#'   \item \strong{Cut Point Generation}: Creates binary indicators at
#'     data-driven or pre-specified cut points
#'   \item \strong{Combinatorial Search}: Exhaustively searches combinations
#'     of up to \code{maxk} factors
#'   \item \strong{Consistency Validation}: Uses split-sample validation to
#'     assess reproducibility
#' }
#'
#' \subsection{Two-Stage Consistency Algorithm}{
#' When \code{use_twostage = TRUE}, ForestSearch uses an adaptive sequential
#' algorithm that can dramatically reduce computation time:
#'
#' \enumerate{
#'   \item \strong{Screening Stage}: Evaluates candidates with a small number
#'     of splits (\code{n.splits.screen}) to quickly identify promising subgroups
#'   \item \strong{Sequential Refinement}: For candidates passing screening,
#'     adds splits in batches until either the consistency threshold is met
#'     with high confidence, or the maximum splits is reached
#' }
#'
#' This approach typically achieves 3-10x speedup while maintaining statistical
#' validity through proper confidence interval construction.
#' }
#'
#' \subsection{GRF Cut Selection}{
#' When \code{use_grf = TRUE}, the function calls \code{grf.subg.harm.survival()}
#' to identify data-driven cut points. The \code{return_selected_cuts_only}
#' parameter controls which cuts are used:
#' \itemize{
#'   \item \code{FALSE} (default): Uses all cuts from all tree depths for maximum
#'     exploration
#'   \item \code{TRUE}: Uses only cuts from the tree that identified the winning
#'     subgroup, providing a more focused set of cuts
#' }
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Standard analysis (backward compatible)
#' result <- forestsearch(
#'   df.analysis = trial_data,
#'   sg_focus = "hr",
#'   hr.threshold = 1.25,
#'   pconsistency.threshold = 0.90,
#'   fs.splits = 400,
#'   details = TRUE
#' )
#'
#' # Example 2: Fast exploratory analysis with two-stage
#' result_fast <- forestsearch(
#'   df.analysis = trial_data,
#'   sg_focus = "maxSG",
#'   hr.threshold = 1.15,
#'   pconsistency.threshold = 0.85,
#'   fs.splits = 500,
#'   use_twostage = TRUE,
#'   details = TRUE
#' )
#'
#' # Example 3: Two-stage with custom parameters
#' result_custom <- forestsearch(
#'   df.analysis = trial_data,
#'   sg_focus = "hr",
#'   hr.threshold = 1.3,
#'   pconsistency.threshold = 0.95,
#'   fs.splits = 600,
#'   use_twostage = TRUE,
#'   twostage_args = list(
#'     n.splits.screen = 50,
#'     batch.size = 25,
#'     conf.level = 0.99
#'   ),
#'   parallel_args = list(plan = "multisession", workers = 4),
#'   details = TRUE
#' )
#'
#' # Example 4: GRF with selected cuts only
#' result_focused <- forestsearch(
#'   df.analysis = trial_data,
#'   use_grf = TRUE,
#'   return_selected_cuts_only = TRUE,
#'   dmin.grf = 0.1,
#'   grf_depth = 2,
#'   details = TRUE
#' )
#' }
#'
#' @references
#' \itemize{
#'   \item FDA Guidance for Industry: Enrichment Strategies for Clinical Trials
#'   \item Athey & Imbens (2016). Recursive partitioning for heterogeneous
#'     causal effects. PNAS.
#'   \item Wager & Athey (2018). Estimation and inference of heterogeneous
#'     treatment effects using random forests. JASA.
#' }
#'
#' @seealso
#' \code{\link{subgroup.consistency}} for consistency evaluation details
#' \code{\link{forestsearch_bootstrap_dofuture}} for bootstrap inference
#' \code{\link{forestsearch_Kfold}} for cross-validation
#' \code{\link{grf.subg.harm.survival}} for GRF subgroup identification
#'
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
                         is.RCT = TRUE,
                         seedit = 8316951,
                         est.scale = "hr",
                         use_lasso = TRUE,
                         use_grf = TRUE,
                         grf_res = NULL,
                         grf_cuts = NULL,
                         max_n_confounders = 1000,
                         grf_depth = 2,
                         dmin.grf = 12,
                         frac.tau = 0.6,
                         return_selected_cuts_only = FALSE,
                         conf_force = NULL,
                         defaultcut_names = NULL,
                         cut_type = "default",
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
                         stop_threshold = 0.95,
                         showten_subgroups = FALSE,
                         d0.min = 12,
                         d1.min = 12,
                         max.minutes = 3,
                         minp = 0.025,
                         details = FALSE,
                         maxk = 2,
                         by.risk = 12,
                         plot.sg = FALSE,
                         plot.grf = FALSE,
                         max_subgroups_search = 10,
                         vi.grf.min = -0.2,
                         # Two-stage consistency parameters
                         use_twostage = TRUE,
                         twostage_args = list()) {

  # ===========================================================================
  # SECTION 1: CAPTURE ALL ARGUMENTS FOR REPRODUCIBILITY
  # ===========================================================================

  args_names <- names(formals())
  args_call_all <- mget(args_names, envir = environment())

  # ===========================================================================
  # SECTION 2: VALIDATE INPUTS
  # ===========================================================================

  # Validate parallel arguments
  if (length(parallel_args) > 0) {
    allowed_plans <- c("multisession", "multicore", "callr", "sequential")
    plan_type <- parallel_args$plan
    n_workers <- parallel_args$workers
    max_cores <- parallel::detectCores()

    if (is.null(plan_type)) {
      stop("parallel_args$plan must be specified.")
    }
    if (!plan_type %in% allowed_plans) {
      stop("parallel_args$plan must be one of: ", paste(allowed_plans, collapse = ", "))
    }
    if (!is.null(n_workers) && n_workers > max_cores) {
      warning(sprintf("Requested workers (%d) exceeds available cores (%d). Using %d.",
                      n_workers, max_cores, max_cores))
      parallel_args$workers <- max_cores
    }
  }

  # Validate sg_focus
  valid_sg_focus <- c("hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG")
  if (!sg_focus %in% valid_sg_focus) {
    stop("sg_focus must be one of: ", paste(valid_sg_focus, collapse = ", "))
  }

  # Reset stop_threshold for sg_focus values that require full evaluation
  if (sg_focus %in% c("hrMaxSG", "hrMinSG")) {
    if (!is.null(stop_threshold) && stop_threshold < 1.0) {
      if (details) {
        cat("Note: stop_threshold reset to NULL for sg_focus =", sg_focus, "\n")
      }
      stop_threshold <- NULL
    }
  }

  # Validate return_selected_cuts_only

  if (!is.logical(return_selected_cuts_only) || length(return_selected_cuts_only) != 1) {
    stop("return_selected_cuts_only must be a single logical value (TRUE or FALSE)")
  }

  # ===========================================================================
  # SECTION 3: INITIALIZATION
  # ===========================================================================

  t0 <- Sys.time()
  set.seed(seedit)

  # Initialize output variables
  grf_plot <- NULL
  df.est_out <- NULL
  df.predict_out <- NULL
  df.test_out <- NULL
  sg.harm <- NULL
  prop_maxk <- NULL
  max_sg_est <- NULL

  # ===========================================================================
  # SECTION 3A: GRF CUT GENERATION (if use_grf = TRUE)
  # ===========================================================================

  # If using grf and cuts not already populated, run grf.subg.harm.survival
  if (use_grf && (is.null(grf_res) || is.null(grf_res$tree.cuts))) {
    if (details) {
      cat("GRF stage for cut selection with dmin, tau =", c(dmin.grf, frac.tau), "\n")
      if (return_selected_cuts_only) {
        cat("  return_selected_cuts_only = TRUE: using cuts from selected tree only\n")
      }
    }

    grf_res <- tryCatch(
      grf.subg.harm.survival(
        data = df.analysis,
        confounders.name = confounders.name,
        outcome.name = outcome.name,
        event.name = event.name,
        id.name = id.name,
        treat.name = treat.name,
        frac.tau = frac.tau,
        n.min = n.min,
        dmin.grf = dmin.grf,
        RCT = is.RCT,
        details = details,
        maxdepth = grf_depth,
        seedit = seedit,
        return_selected_cuts_only = return_selected_cuts_only
      ),
      error = function(e) {
        warning("GRF analysis failed: ", e$message)
        return(NULL)
      }
    )

    # Extract GRF cuts if successful
    if (!is.null(grf_res) && !inherits(grf_res, "try-error")) {
      grf_cuts <- grf_res$tree.cuts

      if (details) {
        cat("GRF cuts identified:", length(grf_cuts), "\n")
        if (length(grf_cuts) > 0) {
          cat("  Cuts:", paste(grf_cuts, collapse = ", "), "\n")
        }
        if (!is.null(grf_res$selected_depth)) {
          cat("  Selected tree depth:", grf_res$selected_depth, "\n")
        }
      }

      # Generate GRF plot if requested
      if (plot.grf && !is.null(grf_res$tree)) {
        grf_plot <- tryCatch({
          plot(grf_res$tree)
        }, error = function(e) {
          warning("GRF plot failed: ", e$message)
          NULL
        })
      }
    } else {
      grf_cuts <- NULL
      if (details) {
        cat("GRF analysis did not produce cuts\n")
      }
    }
  }

  # ===========================================================================
  # SECTION 4: DATA PREPARATION (get_FSdata)
  # ===========================================================================

  FSdata <- tryCatch(
    do.call(
      get_FSdata,
      filter_call_args(
        args_call_all,
        get_FSdata,
        list(df.analysis = df.analysis, grf_cuts = grf_cuts)
      )
    ),
    error = function(e) {
      warning("Error in get_FSdata: ", e$message)
      return(NULL)
    }
  )

  if (inherits(FSdata, "try-error") || is.null(FSdata)) {
    warning("FSdata failure - returning NULL result")
    return(list(sg.harm = NULL))
  }

  # Extract FSdata components
  lassoomit <- FSdata$lassoomit
  lassokeep <- FSdata$lassokeep
  df <- FSdata$df

  Y <- df[, outcome.name]
  Event <- df[, event.name]
  Treat <- df[, treat.name]

  FSconfounders.name <- FSdata$confs_names
  confs_labels <- FSdata$confs

  if (is.null(df.predict)) df.predict <- df

  # ===========================================================================
  # SECTION 5: GRF VARIABLE IMPORTANCE SCREENING
  # ===========================================================================

  if (!is.null(vi.grf.min)) {
    X <- as.matrix(df[, FSconfounders.name])
    X <- apply(X, 2, as.numeric)

    tau.rmst <- min(c(max(Y[Treat == 1 & Event == 1]), max(Y[Treat == 0 & Event == 1])))

    if (!is.RCT) {
      cs.forest <- try(suppressWarnings(
        grf::causal_survival_forest(X, Y, Treat, Event,
                                    horizon = 0.9 * tau.rmst, seed = 8316951)
      ), TRUE)
    } else {
      cs.forest <- try(suppressWarnings(
        grf::causal_survival_forest(X, Y, Treat, Event, W.hat = 0.5,
                                    horizon = 0.9 * tau.rmst, seed = 8316951)
      ), TRUE)
    }

    vi.cs <- round(grf::variable_importance(cs.forest), 4)
    vi.cs2 <- data.frame(confs_labels, FSconfounders.name, vi.cs)
    vi.order <- order(vi.cs, decreasing = TRUE)
    vi.cs2 <- vi.cs2[vi.order, ]

    conf.screen <- vi.cs2[, 2]
    vi_ratio <- vi.cs2[, 3] / max(vi.cs2[, 3])
    selected.vars <- which(vi_ratio > vi.grf.min)
    conf.screen <- conf.screen[selected.vars]
    conf.screen <- conf.screen[seq_len(min(length(conf.screen), max_n_confounders))]
  } else {
    conf.screen <- FSconfounders.name
  }

  # ===========================================================================
  # SECTION 6: SUBGROUP SEARCH
  # ===========================================================================

  # Prepare data for subgroup search
  df.confounders <- df[, conf.screen]
  df.confounders <- dummy(df.confounders)

  id <- df[, c(id.name)]
  df.fs <- data.frame(Y, Event, Treat, id, df.confounders)
  Z <- as.matrix(df.confounders)
  colnames(Z) <- names(df.confounders)

  # Run subgroup search
  find.grps <- tryCatch(
    do.call(
      find_subgroups,
      filter_call_args(
        args_call_all,
        find_subgroups,
        list(
          df.fs = df.fs,
          Z = Z,
          Y = Y,
          Event = Event,
          Treat = Treat
        )
      )
    ),
    error = function(e) {
      warning("Error in find_subgroups: ", e$message)
      return(NULL)
    }
  )

  if (inherits(find.grps, "try-error") || is.null(find.grps)) {
    warning("find_subgroups failure - returning NULL result")
    return(list(sg.harm = NULL, args_call_all = args_call_all))
  }

  # ===========================================================================
  # SECTION 7: CONSISTENCY EVALUATION
  # ===========================================================================

  has_subgroups <- !is.null(find.grps$out_hr) && nrow(find.grps$out_hr) > 0

  if (has_subgroups) {
    # Build consistency arguments
    consistency_args <- filter_call_args(
      args_call_all,
      subgroup.consistency,
      list(
        df = df.fs,
        allgroups.old = find.grps$out_hr,
        Z = Z
      )
    )

    # Add two-stage parameters if enabled
    if (use_twostage) {
      consistency_args$use_twostage <- TRUE
      consistency_args$twostage_args <- twostage_args
    }

    # Run consistency evaluation
    grp.consistency <- tryCatch(
      do.call(subgroup.consistency, consistency_args),
      error = function(e) {
        warning("Error in subgroup.consistency: ", e$message)
        return(NULL)
      }
    )

    # Extract selected subgroup
    if (!is.null(grp.consistency) && !is.null(grp.consistency$out_sg)) {
      sg.harm <- grp.consistency$sg.harm

      if (details && !is.null(grp.consistency$algorithm)) {
        cat("Consistency algorithm:", grp.consistency$algorithm, "\n")
        if (!is.null(grp.consistency$n_candidates_evaluated)) {
          cat("  Candidates evaluated:", grp.consistency$n_candidates_evaluated, "\n")
          cat("  Candidates passed:", grp.consistency$n_passed, "\n")
        }
      }
    }
  } else {
    grp.consistency <- NULL
    if (details) {
      cat("No candidate subgroups found meeting criteria\n")
    }
  }

  # ===========================================================================
  # SECTION 8: ESTIMATION DATASET PREPARATION
  # ===========================================================================

  if (has_subgroups && !is.null(sg.harm)) {
    df.est_out <- tryCatch({
      get_dfpred(
        df.predict = df,
        sg.harm = sg.harm,
        version = 2
      )
    }, error = function(e) {
      warning("Error in get_dfpred for df.est: ", e$message)
      NULL
    })

    # Compute max_sg_est if available
    if (!is.null(df.est_out) && "treat.recommend" %in% names(df.est_out)) {
      max_sg_est <- sum(df.est_out$treat.recommend == 0, na.rm = TRUE)
    }

    # Compute prop_maxk
    if (!is.null(find.grps$out_hr)) {
      n_maxk <- sum(find.grps$out_hr$nfactors == maxk, na.rm = TRUE)
      prop_maxk <- n_maxk / nrow(find.grps$out_hr)
    }
  }

  # ===========================================================================
  # SECTION 9: PREDICTION AND TEST DATASETS
  # ===========================================================================

  if (has_subgroups && !is.null(sg.harm)) {
    # Prediction dataset
    if (!is.null(df.predict) && !identical(df.predict, df)) {
      df.predict_out <- tryCatch({
        get_dfpred(
          df.predict = df.predict,
          sg.harm = grp.consistency$sg.harm,
          version = 2
        )
      }, error = function(e) {
        warning("Error in get_dfpred for df.predict: ", e$message)
        NULL
      })
    }

    # Test dataset
    if (!is.null(df.test)) {
      df.test_out <- tryCatch({
        get_dfpred(
          df.predict = df.test,
          sg.harm = grp.consistency$sg.harm,
          version = 2
        )
      }, error = function(e) {
        warning("Error in get_dfpred for df.test: ", e$message)
        NULL
      })
    }
  }

  # ===========================================================================
  # SECTION 10: COMPILE AND RETURN OUTPUT
  # ===========================================================================

  t.min_all <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  out <- list(
    grp.consistency = grp.consistency,
    find.grps = find.grps,
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
    args_call_all = args_call_all,
    # Include algorithm information
    consistency_algorithm = if (!is.null(grp.consistency)) grp.consistency$algorithm else NA
  )

  # Return early if FSdata or find.grps failed
  if (inherits(FSdata, "try-error") || inherits(find.grps, "try-error")) {
    out <- list(sg.harm = NULL, args_call_all = args_call_all)
  }

  class(out) <- c("forestsearch", "list")
  return(out)
}


#' Print Method for forestsearch Objects
#'
#' @param x A forestsearch object
#' @param ... Additional arguments (unused)
#'
#' @export
print.forestsearch <- function(x, ...) {
  cat("ForestSearch Results\n")
  cat("====================\n\n")

  if (is.null(x$sg.harm)) {
    cat("No subgroup identified.\n")
    return(invisible(x))
  }

  cat("Selected Subgroup:\n")
  cat("  Definition:", paste(x$sg.harm, collapse = " & "), "\n")

  if (!is.null(x$grp.consistency)) {
    gc <- x$grp.consistency
    cat("  sg_focus:", gc$sg_focus, "\n")

    if (!is.null(gc$out_sg) && !is.null(gc$out_sg$result)) {
      top <- gc$out_sg$result[1, ]
      cat("  N:", top$N, "\n")
      cat("  HR:", round(as.numeric(top$hr), 3), "\n")
      cat("  Pcons:", round(as.numeric(top$Pcons), 3), "\n")
    }

    if (!is.null(gc$algorithm)) {
      cat("  Algorithm:", gc$algorithm, "\n")
    }

    if (!is.null(gc$n_candidates_evaluated)) {
      cat("  Candidates evaluated:", gc$n_candidates_evaluated, "\n")
      cat("  Candidates passed:", gc$n_passed, "\n")
    }
  }

  cat("\nComputation time:", round(x$minutes_all, 2), "minutes\n")

  invisible(x)
}


#' Summary Method for forestsearch Objects
#'
#' @param object A forestsearch object to summarize
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the input object
#' @export
summary.forestsearch <- function(object, ...) {
  cat("ForestSearch Summary\n")
  cat("====================\n\n")

  # Parameters section
  cat("Analysis Parameters:\n")
  if (!is.null(object$args_call_all)) {
    params <- object$args_call_all
    cat("  sg_focus:", params$sg_focus, "\n")
    cat("  hr.threshold:", params$hr.threshold, "\n")
    cat("  hr.consistency:", params$hr.consistency, "\n")
    cat("  pconsistency.threshold:", params$pconsistency.threshold, "\n")
    cat("  n.min:", params$n.min, "\n")
    cat("  fs.splits:", params$fs.splits, "\n")
    cat("  maxk:", params$maxk, "\n")
    cat("  use_twostage:", params$use_twostage, "\n")
    cat("  return_selected_cuts_only:", params$return_selected_cuts_only, "\n")
  }
  cat("\n")

  # Variable selection
  cat("Variable Selection:\n")
  cat("  Candidate confounders:", length(object$confounders.candidate), "\n")
  cat("  Confounders evaluated:", length(object$confounders.evaluated), "\n\n")

  # GRF information
  if (!is.null(object$grf_res)) {
    cat("GRF Results:\n")
    cat("  Cuts identified:", length(object$grf_cuts), "\n")
    if (!is.null(object$grf_res$selected_depth)) {
      cat("  Selected tree depth:", object$grf_res$selected_depth, "\n")
    }
    if (!is.null(object$grf_res$return_selected_cuts_only)) {
      cat("  Using selected cuts only:", object$grf_res$return_selected_cuts_only, "\n")
    }
    cat("\n")
  }

  # Consistency results
  if (!is.null(object$grp.consistency)) {
    gc <- object$grp.consistency
    cat("Consistency Evaluation:\n")
    cat("  Algorithm:", ifelse(!is.null(gc$algorithm), gc$algorithm, "fixed"), "\n")

    if (!is.null(gc$n_candidates_evaluated)) {
      cat("  Candidates evaluated:", gc$n_candidates_evaluated, "\n")
      cat("  Candidates passed:", gc$n_passed, "\n")
    }
    cat("\n")
  }

  # Selected subgroup
  if (!is.null(object$sg.harm)) {
    cat("Selected Subgroup:\n")
    cat("  Definition:", paste(object$sg.harm, collapse = " & "), "\n")

    if (!is.null(object$grp.consistency$out_sg$result)) {
      top <- object$grp.consistency$out_sg$result[1, ]
      cat("  Sample size:", top$N, "\n")
      cat("  Hazard ratio:", round(as.numeric(top$hr), 3), "\n")
      cat("  Consistency:", round(as.numeric(top$Pcons) * 100, 1), "%\n")
    }
  } else {
    cat("No subgroup identified.\n")
  }

  invisible(object)
}
