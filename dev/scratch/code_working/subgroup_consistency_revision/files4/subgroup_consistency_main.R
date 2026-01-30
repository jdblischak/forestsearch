#' @title Subgroup Consistency Evaluation
#'
#' @description
#' Evaluates consistency of subgroups found in a survival analysis, using random
#' splits and hazard ratio criteria. Supports both fixed-sample and two-stage
#' sequential evaluation algorithms.
#'
#' @param df The original data.frame with columns Y, Event, Treat, id.
#' @param hr.subgroups Data.table of subgroup hazard ratio results from subgroup search.
#' @param hr.threshold Numeric. Minimum hazard ratio for subgroup inclusion. Default 1.0.
#' @param hr.consistency Numeric. Minimum hazard ratio for consistency in splits. Default 1.0.
#' @param pconsistency.threshold Numeric. Minimum proportion of splits meeting consistency. Default 0.9.
#' @param pconsistency.digits Integer. Significant digits for pconsistency.threshold. Default 2.
#' @param m1.threshold Numeric. Maximum median survival for treatment arm. Default Inf.
#' @param n.splits Integer. Number of random splits for consistency evaluation (or maximum
#'   splits when \code{use_twostage = TRUE}). Default 100.
#' @param details Logical. Print progress details. Default FALSE.
#' @param by.risk Numeric. Risk interval for plotting. Default 12.
#' @param plot.sg Logical. Plot subgroup survival curves. Default FALSE.
#' @param maxk Integer. Maximum number of covariates in a subgroup. Default 7.
#' @param Lsg Integer. Number of covariates (required).
#' @param confs_labels Character vector. Covariate label mapping (required).
#' @param sg_focus Character. Sorting focus for subgroup selection. One of
#'   "hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG". Default "hr".
#' @param stop_Kgroups Integer. Maximum number of subgroups to evaluate. Default 10.
#' @param checking Logical. Enable debugging output. Default FALSE.
#' @param parallel_args List. Parallel processing configuration with elements:
#'   \describe{
#'     \item{plan}{Character. One of "multisession", "multicore", "callr", "sequential"}
#'     \item{workers}{Integer. Number of parallel workers}
#'     \item{show_message}{Logical. Show parallel setup messages}
#'   }
#' @param use_twostage Logical. Use two-stage sequential algorithm for improved
#'   performance. Default FALSE for backward compatibility. When TRUE, enables
#'   early stopping and screening to reduce computation time by 3-10x for most
#'   analyses.
#' @param twostage_args List. Parameters for two-stage algorithm (only used when
#'   \code{use_twostage = TRUE}):
#'   \describe{
#'     \item{n.splits.screen}{Integer. Splits for Stage 1 screening. Default 30.}
#'     \item{screen.threshold}{Numeric. Consistency threshold for Stage 1. Default
#'       is automatically calculated to provide ~2.5 SE margin below pconsistency.threshold.}
#'     \item{batch.size}{Integer. Splits per batch in Stage 2. Default 20.}
#'     \item{conf.level}{Numeric. Confidence level for early stopping decisions. Default 0.95.}
#'     \item{min.valid.screen}{Integer. Minimum valid splits required in Stage 1. Default 10.}
#'   }
#'
#' @return A list containing:
#'   \describe{
#'     \item{out_sg}{Results object from sg_consistency_out()}
#'     \item{sg_focus}{The sg_focus value used}
#'     \item{df_flag}{Data.frame with id and treat.recommend columns}
#'     \item{sg.harm}{Character vector of selected subgroup labels}
#'     \item{sg.harm.id}{Numeric ID of selected subgroup}
#'     \item{algorithm}{Character indicating which algorithm was used}
#'     \item{n_candidates_evaluated}{Number of candidates that underwent consistency evaluation}
#'     \item{n_passed}{Number of candidates meeting consistency threshold}
#'   }
#'
#' @details
#' The function supports two evaluation algorithms:
#'
#' \strong{Fixed-Sample Algorithm} (\code{use_twostage = FALSE}):
#' \itemize{
#'   \item Runs exactly \code{n.splits} random splits for each candidate
#'   \item Consistent and predictable runtime
#'   \item Recommended for final analyses requiring exact reproducibility
#' }
#'
#' \strong{Two-Stage Sequential Algorithm} (\code{use_twostage = TRUE}):
#' \itemize{
#'   \item Stage 1: Quick screening with \code{n.splits.screen} splits
#'   \item Stage 2: Sequential batched evaluation with early stopping
#'   \item Candidates clearly passing/failing stop early
#'   \item 3-10x faster for typical analyses
#'   \item Recommended for exploratory analyses and large candidate sets
#' }
#'
#' @section Performance Considerations:
#' The two-stage algorithm provides significant speedups when:
#' \itemize{
#'   \item Many candidates clearly fail consistency (screened at Stage 1)
#'   \item Many candidates have consistency well above/below threshold
#'   \item \code{n.splits} is large (>200)
#' }
#'
#' Speedup is minimal when most candidates have true consistency near the threshold.
#'
#' @examples
#' \dontrun{
#' # Standard fixed-sample evaluation
#' result <- subgroup.consistency(
#'   df = trial_data,
#'   hr.subgroups = candidates,
#'   n.splits = 400,
#'   pconsistency.threshold = 0.90,
#'   use_twostage = FALSE
#' )
#'
#' # Two-stage sequential evaluation (faster)
#' result_fast <- subgroup.consistency(
#'   df = trial_data,
#'   hr.subgroups = candidates,
#'   n.splits = 400,  # Maximum splits
#'   pconsistency.threshold = 0.90,
#'   use_twostage = TRUE,
#'   twostage_args = list(
#'     n.splits.screen = 30,
#'     batch.size = 20
#'   )
#' )
#' }
#'
#' @seealso
#' \code{\link{evaluate_subgroup_consistency}} for fixed-sample evaluation
#' \code{\link{evaluate_consistency_twostage}} for two-stage evaluation
#' \code{\link{forestsearch}} for the main analysis function
#'
#' @importFrom data.table copy as.data.table is.data.table
#' @importFrom survival coxph Surv
#' @importFrom future.apply future_lapply
#' @export

subgroup.consistency <- function(df,
                                 hr.subgroups,
                                 hr.threshold = 1.0,
                                 hr.consistency = 1.0,
                                 pconsistency.threshold = 0.9,
                                 m1.threshold = Inf,
                                 n.splits = 100,
                                 details = FALSE,
                                 by.risk = 12,
                                 plot.sg = FALSE,
                                 maxk = 7,
                                 Lsg,
                                 confs_labels,
                                 sg_focus = "hr",
                                 stop_Kgroups = 10,
                                 pconsistency.digits = 2,
                                 checking = FALSE,
                                 parallel_args = list(NULL),
                                 # New two-stage parameters
                                 use_twostage = FALSE,
                                 twostage_args = list()) {

  # ===========================================================================
  # SECTION 1: INPUT VALIDATION
  # ===========================================================================

  # Check required inputs exist
  if (missing(df) || !is.data.frame(df)) {
    stop("'df' must be provided and must be a data.frame")
  }

  if (missing(hr.subgroups) || is.null(hr.subgroups)) {
    stop("'hr.subgroups' must be provided")
  }

  if (missing(Lsg) || !is.numeric(Lsg) || Lsg < 1) {
    stop("'Lsg' must be a positive integer representing number of covariates")
  }

  if (missing(confs_labels) || !is.character(confs_labels)) {
    stop("'confs_labels' must be a character vector of covariate labels")
  }

  # Validate use_twostage
 if (!is.logical(use_twostage) || length(use_twostage) != 1) {
    stop("'use_twostage' must be a single logical value (TRUE or FALSE)")
  }

  # Validate twostage_args if provided
  if (use_twostage && length(twostage_args) > 0) {
    valid_ts_args <- c("n.splits.screen", "screen.threshold", "batch.size",
                       "conf.level", "min.valid.screen")
    invalid_args <- setdiff(names(twostage_args), valid_ts_args)
    if (length(invalid_args) > 0) {
      warning("Unknown twostage_args parameters ignored: ",
              paste(invalid_args, collapse = ", "))
    }
  }

  # Check data.table format
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required")
  }

  if (!data.table::is.data.table(hr.subgroups)) {
    hr.subgroups <- data.table::as.data.table(hr.subgroups)
  }

  # ===========================================================================
  # SECTION 2: SET UP TWO-STAGE PARAMETERS (WITH DEFAULTS)
  # ===========================================================================

  # Default two-stage parameters
  ts_defaults <- list(
    n.splits.screen = 30,
    screen.threshold = NULL,  # Will be calculated if NULL
    batch.size = 20,
    conf.level = 0.95,
    min.valid.screen = 10
  )

  # Merge user-provided args with defaults
  ts_params <- modifyList(ts_defaults, twostage_args)

  # Calculate default screen threshold if not provided
  if (is.null(ts_params$screen.threshold)) {
    se_estimate <- sqrt(pconsistency.threshold * (1 - pconsistency.threshold) /
                          ts_params$n.splits.screen)
    ts_params$screen.threshold <- max(0.5, pconsistency.threshold - 2.5 * se_estimate)
  }

  # ===========================================================================
  # SECTION 3: EXTRACT COLUMN NAMES AND VALIDATE STRUCTURE
  # ===========================================================================

  names.Z <- colnames(hr.subgroups)[(ncol(hr.subgroups) - Lsg + 1):ncol(hr.subgroups)]

  required_cols <- c("HR", "n", "E", "grp")
  missing_cols <- setdiff(required_cols, names(hr.subgroups))
  if (length(missing_cols) > 0) {
    stop("hr.subgroups missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ===========================================================================
  # SECTION 4: FILTER SUBGROUPS BY CRITERIA
  # ===========================================================================

  if (nrow(hr.subgroups) == 0) {
    warning("No subgroups provided in hr.subgroups")
    return(list(
      out_sg = NULL, sg_focus = sg_focus,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL,
      algorithm = ifelse(use_twostage, "twostage", "fixed"),
      n_candidates_evaluated = 0, n_passed = 0
    ))
  }

  # Apply m1.threshold filter if specified
  if (is.finite(m1.threshold)) {
    hr.subgroups <- hr.subgroups[!is.na(hr.subgroups$m1), ]
    if (nrow(hr.subgroups) == 0) {
      warning("All subgroups removed after filtering NA m1 values")
      return(list(
        out_sg = NULL, sg_focus = sg_focus,
        df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL,
        algorithm = ifelse(use_twostage, "twostage", "fixed"),
        n_candidates_evaluated = 0, n_passed = 0
      ))
    }
    found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold &
                                hr.subgroups$m1 <= m1.threshold, ]
  } else {
    found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold, ]
  }

  if (nrow(found.hrs) == 0) {
    if (details) {
      cat("No subgroups meet criteria (HR >=", hr.threshold)
      if (is.finite(m1.threshold)) cat(" and m1 <=", m1.threshold)
      cat(")\n")
    }
    return(list(
      out_sg = NULL, sg_focus = sg_focus,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL,
      algorithm = ifelse(use_twostage, "twostage", "fixed"),
      n_candidates_evaluated = 0, n_passed = 0
    ))
  }

  # ===========================================================================
  # SECTION 5: REMOVE DUPLICATES AND SORT BY SG_FOCUS
  # ===========================================================================

  if (nrow(found.hrs) > 1) {
    n_before <- nrow(found.hrs)

    tryCatch({
      found.hrs <- remove_near_duplicate_subgroups(found.hrs, details = details)
    }, error = function(e) {
      warning("Error removing duplicates: ", e$message, ". Proceeding with original subgroups.")
    })

    if (nrow(found.hrs) == 0) {
      stop("All subgroups removed during duplicate removal. This should not happen.")
    }
  }

  # Sort based on sg_focus to prioritize candidates
  if (sg_focus == "maxSG") {
    found.hrs <- found.hrs[order(found.hrs$n, decreasing = TRUE), ]
  } else if (sg_focus == "minSG") {
    found.hrs <- found.hrs[order(found.hrs$n, decreasing = FALSE), ]
  }

  # Extract index matrix
  index.Z <- found.hrs[, names.Z, with = FALSE]

  if (details) {
    cat("# of unique initial candidates:", nrow(found.hrs), "\n")
  }

  # Limit to top stop_Kgroups candidates
  maxsgs <- min(nrow(found.hrs), stop_Kgroups)
  found.hrs <- found.hrs[seq_len(maxsgs), ]
  index.Z <- index.Z[seq_len(maxsgs), ]

  n_candidates <- nrow(found.hrs)

  if (details) {
    cat("# Restricting to top stop_Kgroups =", stop_Kgroups, "\n")
    cat("# of candidates to evaluate:", n_candidates, "\n")
  }

  # ===========================================================================
  # SECTION 6: VALIDATE PARALLEL CONFIGURATION
  # ===========================================================================

  use_parallel <- length(parallel_args) > 0 && !is.null(parallel_args[[1]])

  if (use_parallel) {
    required_parallel <- c("plan", "workers")
    if (!all(required_parallel %in% names(parallel_args))) {
      warning("parallel_args missing required elements. Using sequential processing.")
      use_parallel <- FALSE
    }

    valid_plans <- c("multisession", "multicore", "callr", "sequential")
    if (use_parallel && !parallel_args$plan %in% valid_plans) {
      warning("Invalid parallel plan '", parallel_args$plan,
              "'. Using sequential processing.")
      use_parallel <- FALSE
    }

    if (use_parallel && (!is.numeric(parallel_args$workers) || parallel_args$workers < 1)) {
      warning("Invalid workers value. Using sequential processing.")
      use_parallel <- FALSE
    }
  }

  if (details) {
    algorithm_name <- ifelse(use_twostage, "Two-stage sequential", "Fixed-sample")
    cat("Algorithm:", algorithm_name, "\n")
    if (use_twostage) {
      cat("  Stage 1 splits:", ts_params$n.splits.screen, "\n")
      cat("  Screen threshold:", round(ts_params$screen.threshold, 3), "\n")
      cat("  Max total splits:", n.splits, "\n")
      cat("  Batch size:", ts_params$batch.size, "\n")
    } else {
      cat("  Splits per candidate:", n.splits, "\n")
    }
    if (use_parallel) {
      cat("Parallel processing:", parallel_args$plan,
          "with", parallel_args$workers, "workers\n")
    } else {
      cat("Sequential processing\n")
    }
  }

  # ===========================================================================
  # SECTION 7: INITIALIZE TIMING
  # ===========================================================================

  t.start <- proc.time()[3]

  # ===========================================================================
  # SECTION 8: EVALUATE EACH SUBGROUP
  # ===========================================================================

  # ----- EXECUTE EVALUATION -----

  if (!use_parallel) {
    # SEQUENTIAL EXECUTION
    if (use_twostage) {
      results_list <- lapply(seq_len(n_candidates), function(m) {
        evaluate_consistency_twostage(
          m = m,
          index.Z = index.Z,
          names.Z = names.Z,
          df = df,
          found.hrs = found.hrs,
          hr.consistency = hr.consistency,
          pconsistency.threshold = pconsistency.threshold,
          pconsistency.digits = pconsistency.digits,
          maxk = maxk,
          confs_labels = confs_labels,
          details = details,
          n.splits.screen = ts_params$n.splits.screen,
          screen.threshold = ts_params$screen.threshold,
          n.splits.max = n.splits,
          batch.size = ts_params$batch.size,
          conf.level = ts_params$conf.level,
          min.valid.screen = ts_params$min.valid.screen
        )
      })
    } else {
      results_list <- lapply(seq_len(n_candidates), function(m) {
        evaluate_subgroup_consistency(
          m = m,
          index.Z = index.Z,
          names.Z = names.Z,
          df = df,
          found.hrs = found.hrs,
          n.splits = n.splits,
          hr.consistency = hr.consistency,
          pconsistency.threshold = pconsistency.threshold,
          pconsistency.digits = pconsistency.digits,
          maxk = maxk,
          confs_labels = confs_labels,
          details = details
        )
      })
    }

  } else {
    # PARALLEL EXECUTION
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    setup_parallel_SGcons(parallel_args)

    # Capture function references explicitly for parallel workers
    # These must be captured BEFORE creating the eval_fun closure
    .evaluate_consistency_twostage <- evaluate_consistency_twostage
    .evaluate_subgroup_consistency <- evaluate_subgroup_consistency
    .run_single_consistency_split <- run_single_consistency_split
    .get_split_hr_fast <- get_split_hr_fast
    .wilson_ci <- wilson_ci
    .early_stop_decision <- early_stop_decision
    .FS_labels <- FS_labels

    if (use_twostage) {
      # Create closure that captures all necessary functions and data
      eval_fun <- function(m) {
        .evaluate_consistency_twostage(
          m = m,
          index.Z = index.Z,
          names.Z = names.Z,
          df = df,
          found.hrs = found.hrs,
          hr.consistency = hr.consistency,
          pconsistency.threshold = pconsistency.threshold,
          pconsistency.digits = pconsistency.digits,
          maxk = maxk,
          confs_labels = confs_labels,
          details = FALSE,  # Suppress details in parallel workers
          n.splits.screen = ts_params$n.splits.screen,
          screen.threshold = ts_params$screen.threshold,
          n.splits.max = n.splits,
          batch.size = ts_params$batch.size,
          conf.level = ts_params$conf.level,
          min.valid.screen = ts_params$min.valid.screen
        )
      }
    } else {
      eval_fun <- function(m) {
        .evaluate_subgroup_consistency(
          m = m,
          index.Z = index.Z,
          names.Z = names.Z,
          df = df,
          found.hrs = found.hrs,
          n.splits = n.splits,
          hr.consistency = hr.consistency,
          pconsistency.threshold = pconsistency.threshold,
          pconsistency.digits = pconsistency.digits,
          maxk = maxk,
          confs_labels = confs_labels,
          details = FALSE  # Suppress details in parallel workers
        )
      }
    }

    # Build globals list with all captured functions
    globals_list <- list(
      # Data objects
      index.Z = index.Z,
      names.Z = names.Z,
      df = df,
      found.hrs = found.hrs,
      n.splits = n.splits,
      hr.consistency = hr.consistency,
      pconsistency.threshold = pconsistency.threshold,
      pconsistency.digits = pconsistency.digits,
      maxk = maxk,
      confs_labels = confs_labels,
      ts_params = ts_params,
      # Captured functions with dot prefix
      .evaluate_consistency_twostage = .evaluate_consistency_twostage,
      .evaluate_subgroup_consistency = .evaluate_subgroup_consistency,
      .run_single_consistency_split = .run_single_consistency_split,
      .get_split_hr_fast = .get_split_hr_fast,
      .wilson_ci = .wilson_ci,
      .early_stop_decision = .early_stop_decision,
      .FS_labels = .FS_labels,
      # Also include original names for internal calls
      evaluate_consistency_twostage = .evaluate_consistency_twostage,
      evaluate_subgroup_consistency = .evaluate_subgroup_consistency,
      run_single_consistency_split = .run_single_consistency_split,
      get_split_hr_fast = .get_split_hr_fast,
      wilson_ci = .wilson_ci,
      early_stop_decision = .early_stop_decision,
      FS_labels = .FS_labels
    )

    results_list <- future.apply::future_lapply(
      seq_len(n_candidates),
      eval_fun,
      future.seed = TRUE,
      future.packages = c("survival", "data.table"),
      future.globals = globals_list
    )
  }

  # ===========================================================================
  # SECTION 9: COMPILE RESULTS
  # ===========================================================================

  if (length(results_list) == 0) {
    if (details) cat("No subgroups met consistency criteria\n")
    return(list(
      out_sg = NULL, sg_focus = sg_focus,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL,
      algorithm = ifelse(use_twostage, "twostage", "fixed"),
      n_candidates_evaluated = n_candidates, n_passed = 0
    ))
  }

  # Filter out NULL results
  results_list <- Filter(Negate(is.null), results_list)
  n_passed <- length(results_list)

  if (n_passed == 0) {
    if (details) cat("All subgroup evaluations returned NULL\n")
    return(list(
      out_sg = NULL, sg_focus = sg_focus,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL,
      algorithm = ifelse(use_twostage, "twostage", "fixed"),
      n_candidates_evaluated = n_candidates, n_passed = 0
    ))
  }

  # Combine results into data.table
  res <- tryCatch({
    data.table::as.data.table(do.call(rbind, results_list))
  }, error = function(e) {
    stop("Error combining results: ", e$message,
         "\nThis may indicate inconsistent result structure across subgroups.")
  })

  any.found <- nrow(res)

  if (any.found == 0) {
    if (details) cat("No subgroups found meeting consistency threshold\n")
    return(list(
      out_sg = NULL, sg_focus = sg_focus,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL,
      algorithm = ifelse(use_twostage, "twostage", "fixed"),
      n_candidates_evaluated = n_candidates, n_passed = 0
    ))
  }

  # Convert columns to numeric
  cols_to_numeric <- c("Pcons", "hr", "N", "E", "K")
  missing_cols <- setdiff(cols_to_numeric, names(res))
  if (length(missing_cols) > 0) {
    stop("Result data.table missing expected columns: ",
         paste(missing_cols, collapse = ", "))
  }

  tryCatch({
    res[, (cols_to_numeric) := lapply(.SD, as.numeric), .SDcols = cols_to_numeric]
  }, error = function(e) {
    stop("Error converting result columns to numeric: ", e$message)
  })

  # ===========================================================================
  # SECTION 10: GENERATE OUTPUT
  # ===========================================================================

  out_sg <- NULL
  df_flag <- sg.harm <- sg.harm.id <- NULL

  if (any.found > 0) {
    result_new <- data.table::copy(res)

    # Map sg_focus to output type
    sg_focus_map <- list(
      hr = "hr",
      hrMaxSG = "maxSG",
      maxSG = "maxSG",
      hrMinSG = "minSG",
      minSG = "minSG"
    )

    if (!sg_focus %in% names(sg_focus_map)) {
      stop(sprintf("Unknown sg_focus value: %s", sg_focus))
    }

    focus_type <- sg_focus_map[[sg_focus]]
    sgdetails <- ifelse(plot.sg, TRUE, FALSE)

    out_sg <- tryCatch({
      sg_consistency_out(
        df = df,
        result_new = result_new,
        sg_focus = focus_type,
        details = sgdetails,
        plot.sg = sgdetails,
        index.Z = index.Z,
        names.Z = names.Z,
        by.risk = by.risk,
        confs_labels = confs_labels
      )
    }, error = function(e) {
      warning("Error in sg_consistency_out for '", focus_type, "': ", e$message)
      NULL
    })

    # Extract results if successful
    if (is.null(out_sg)) {
      warning("No valid output for sg_focus='", sg_focus, "'")
    } else {
      required_fields <- c("df_flag", "sg.harm_label", "sg.harm.id")
      missing_fields <- setdiff(required_fields, names(out_sg))
      if (length(missing_fields) > 0) {
        stop("sg_consistency_out result missing fields: ",
             paste(missing_fields, collapse = ", "))
      }

      df_flag <- out_sg$df_flag
      sg.harm <- out_sg$sg.harm_label
      sg.harm.id <- out_sg$sg.harm.id
    }

    if (details) cat("SG focus=", sg_focus, "\n")
  }

  # ===========================================================================
  # SECTION 11: FINAL TIMING AND OUTPUT
  # ===========================================================================

  if (details) {
    t.end <- proc.time()[3]
    t.min <- (t.end - t.start) / 60
    cat("Subgroup Consistency Minutes=", round(t.min, 3), "\n")
    cat("Algorithm used:", ifelse(use_twostage, "Two-stage sequential", "Fixed-sample"), "\n")
    cat("Candidates evaluated:", n_candidates, "\n")
    cat("Candidates passed:", n_passed, "\n")

    if (any.found > 0) {
      cat("Subgroup found (FS) with sg_focus='", sg_focus, "'\n", sep = "")
      if (!is.null(sg.harm)) {
        cat("Selected subgroup:", paste(sg.harm, collapse = " & "), "\n")
      }
    } else {
      cat("NO subgroup found (FS)\n")
    }
  }

  output <- list(
    out_sg = out_sg,
    sg_focus = sg_focus,
    df_flag = df_flag,
    sg.harm = sg.harm,
    sg.harm.id = sg.harm.id,
    algorithm = ifelse(use_twostage, "twostage", "fixed"),
    n_candidates_evaluated = n_candidates,
    n_passed = n_passed
  )

  return(output)
}


# =============================================================================
# FIXED-SAMPLE EVALUATION FUNCTION
# =============================================================================

#' Evaluate Single Subgroup for Consistency (Fixed-Sample)
#'
#' Helper function that evaluates a single subgroup (indexed by m) for consistency
#' across random splits using a fixed number of splits. This function contains
#' the core logic used in both sequential and parallel execution modes.
#'
#' @param m Integer. Index of the subgroup to evaluate (1 to nrow(found.hrs))
#' @param index.Z Data.table or matrix. Factor indicators for all subgroups
#' @param names.Z Character vector. Names of factor columns
#' @param df Data.frame. Original data with Y, Event, Treat, id columns
#' @param found.hrs Data.table. Subgroup hazard ratio results
#' @param n.splits Integer. Number of random splits for consistency evaluation
#' @param hr.consistency Numeric. Minimum HR threshold for consistency
#' @param pconsistency.threshold Numeric. Minimum proportion of splits meeting consistency
#' @param pconsistency.digits Integer. Rounding digits for consistency proportion
#' @param maxk Integer. Maximum number of factors in a subgroup
#' @param confs_labels Character vector. Labels for confounders
#' @param details Logical. Print details during execution
#'
#' @return Named numeric vector with consistency results, or NULL if criteria not met.
#'   Vector contains: Pcons, hr, N, E, g, m, K, and factor labels (M.1, M.2, etc.)
#'
#' @importFrom data.table data.table
#' @importFrom survival coxph Surv
#' @export
evaluate_subgroup_consistency <- function(m, index.Z, names.Z, df, found.hrs,
                                          n.splits, hr.consistency,
                                          pconsistency.threshold, pconsistency.digits,
                                          maxk, confs_labels, details = FALSE) {

  # Internal helper for getting HR from split
  get_split_hr <- function(df, cox_initial = NULL) {
    if (nrow(df) < 2 || sum(df$Event) < 2) {
      return(NA_real_)
    }

    fit <- tryCatch(
      suppressWarnings(
        survival::coxph(
          survival::Surv(Y, Event) ~ Treat,
          data = df,
          init = cox_initial,
          robust = FALSE,
          model = FALSE,
          x = FALSE,
          y = FALSE
        )
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) return(NA_real_)
    return(exp(fit$coefficients[1]))
  }

  # -------------------------------------------------------------------------
  # SECTION 1: VALIDATE SUBGROUP EXTRACTION
  # -------------------------------------------------------------------------

  if (m < 1 || m > nrow(index.Z)) {
    warning("Invalid subgroup index m=", m, ". Skipping.")
    return(NULL)
  }

  indexm <- as.numeric(unlist(index.Z[m, ]))

  if (length(indexm) != length(names.Z)) {
    warning("Subgroup ", m, ": index length mismatch. Skipping.")
    return(NULL)
  }

  if (!all(indexm %in% c(0, 1, NA))) {
    warning("Subgroup ", m, ": invalid index values. Skipping.")
    return(NULL)
  }

  this.m <- names.Z[indexm == 1]

  if (length(this.m) == 0) {
    warning("Subgroup ", m, ": no factors selected. Skipping.")
    return(NULL)
  }

  # -------------------------------------------------------------------------
  # SECTION 2: VALIDATE LABEL CONVERSION
  # -------------------------------------------------------------------------

  this.m_label <- tryCatch({
    unlist(lapply(this.m, FS_labels, confs_labels = confs_labels))
  }, error = function(e) {
    warning("Subgroup ", m, ": error in FS_labels: ", e$message, ". Skipping.")
    return(NULL)
  })

  if (is.null(this.m_label)) return(NULL)

  if (length(this.m_label) != length(this.m)) {
    warning("Subgroup ", m, ": label length mismatch. Skipping.")
    return(NULL)
  }

  # -------------------------------------------------------------------------
  # SECTION 3: CREATE SUBGROUP DEFINITION AND EXTRACT DATA
  # -------------------------------------------------------------------------

  id.m <- paste(paste(this.m, collapse = "==1 & "), "==1")

  df.sub <- tryCatch({
    subset(df, eval(parse(text = id.m)))
  }, error = function(e) {
    warning("Subgroup ", m, ": error extracting data: ", e$message, ". Skipping.")
    return(NULL)
  })

  if (is.null(df.sub)) return(NULL)

  if (nrow(df.sub) == 0) {
    warning("Subgroup ", m, ": no observations match criteria. Skipping.")
    return(NULL)
  }

  df.x <- data.table::data.table(df.sub)
  N.x <- nrow(df.x)

  # -------------------------------------------------------------------------
  # SECTION 4: GET INITIAL COX ESTIMATE
  # -------------------------------------------------------------------------

  cox_init <- log(found.hrs$HR[m])
  if (is.na(cox_init) || is.infinite(cox_init)) {
    cox_init <- 0
  }

  # -------------------------------------------------------------------------
  # SECTION 5: PERFORM CONSISTENCY SPLITS
  # -------------------------------------------------------------------------

  flag.consistency <- sapply(seq_len(n.splits), function(bb) {
    in.split1 <- tryCatch({
      sample(c(TRUE, FALSE), N.x, replace = TRUE, prob = c(0.5, 0.5))
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(in.split1)) return(NA_real_)

    df.x$insplit1 <- in.split1

    df.x.split1 <- subset(df.x, insplit1 == 1)
    df.x.split2 <- subset(df.x, insplit1 == 0)

    if (nrow(df.x.split1) < 5 || nrow(df.x.split2) < 5) {
      return(NA_real_)
    }

    if (sum(df.x.split1$Event) < 2 || sum(df.x.split2$Event) < 2) {
      return(NA_real_)
    }

    hr.split1 <- get_split_hr(df = df.x.split1, cox_initial = cox_init)
    hr.split2 <- get_split_hr(df = df.x.split2, cox_initial = cox_init)

    if (!is.na(hr.split1) && !is.na(hr.split2)) {
      as.numeric(hr.split1 > hr.consistency && hr.split2 > hr.consistency)
    } else {
      NA_real_
    }
  })

  # -------------------------------------------------------------------------
  # SECTION 6: CHECK VALIDITY AND CALCULATE CONSISTENCY
  # -------------------------------------------------------------------------

  n_valid_splits <- sum(!is.na(flag.consistency))

  if (n_valid_splits == 0) {
    if (details) {
      cat("Subgroup ", m, ": No valid consistency splits\n")
    }
    return(NULL)
  }

  if (n_valid_splits < 10) {
    warning("Subgroup ", m, ": only ", n_valid_splits, " valid splits out of ",
            n.splits, ". Results may be unreliable.")
  }

  p.consistency <- tryCatch({
    round(mean(flag.consistency, na.rm = TRUE), pconsistency.digits)
  }, error = function(e) {
    return(NA_real_)
  })

  if (is.na(p.consistency)) {
    return(NULL)
  }

  # -------------------------------------------------------------------------
  # SECTION 7: CHECK CONSISTENCY THRESHOLD
  # -------------------------------------------------------------------------

  if (isTRUE(p.consistency < pconsistency.threshold)) {
    if (details) {
      cat("*** Not met: Subgroup, % Consistency =",
          c(this.m_label, p.consistency), "\n")
    }
    return(NULL)
  }

  # -------------------------------------------------------------------------
  # SECTION 8: FORMAT AND RETURN RESULT
  # -------------------------------------------------------------------------

  k <- length(this.m)
  covsm <- rep("M", maxk)
  mindex <- seq_len(maxk)
  Mnames <- paste(covsm, mindex, sep = ".")

  mfound <- matrix(rep("", maxk))
  mfound[seq_len(k)] <- this.m_label

  resultk <- c(
    p.consistency,
    found.hrs$HR[m],
    found.hrs$n[m],
    found.hrs$E[m],
    found.hrs$grp[m],
    m,
    k,
    mfound
  )

  names(resultk) <- c("Pcons", "hr", "N", "E", "g", "m", "K", Mnames)

  if (details) {
    cat("*** Met: Subgroup, % Consistency =",
        c(this.m_label, p.consistency), "\n")
  }

  return(resultk)
}


# =============================================================================
# TWO-STAGE HELPER FUNCTIONS
# =============================================================================
# These should be placed in subgroup_consistency_twostage.R or
# subgroup_consistency_helpers.R

#' Wilson Score Confidence Interval
#'
#' Computes Wilson score confidence interval for a proportion, which has
#' better coverage properties than the normal approximation for small samples
#' and proportions near 0 or 1.
#'
#' @param x Integer. Number of successes.
#' @param n Integer. Number of trials.
#' @param conf.level Numeric. Confidence level (default 0.95).
#'
#' @return Named numeric vector with elements: estimate, lower, upper.
#'
#' @keywords internal
#' @export
wilson_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) {
    return(c(estimate = NA_real_, lower = 0, upper = 1))
  }

  z <- qnorm(1 - (1 - conf.level) / 2)
  p_hat <- x / n

  denominator <- 1 + z^2 / n
  center <- (p_hat + z^2 / (2 * n)) / denominator
  margin <- (z / denominator) * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))

  lower <- max(0, center - margin)
  upper <- min(1, center + margin)

  c(estimate = p_hat, lower = lower, upper = upper)
}


#' Early Stopping Decision
#'
#' Evaluates whether enough evidence exists to stop early based on
#' confidence interval for consistency proportion.
#'
#' @param n_success Integer. Number of splits meeting consistency.
#' @param n_total Integer. Total number of valid splits.
#' @param threshold Numeric. Target consistency threshold.
#' @param conf.level Numeric. Confidence level for decision (default 0.95).
#' @param min_samples Integer. Minimum samples before allowing early stop.
#'
#' @return Character. One of "continue", "pass", or "fail".
#'
#' @keywords internal
#' @export
early_stop_decision <- function(n_success, n_total, threshold,
                                 conf.level = 0.95, min_samples = 20) {
  if (n_total < min_samples) {
    return("continue")
  }

  ci <- wilson_ci(n_success, n_total, conf.level)

  if (ci["lower"] >= threshold) {
    return("pass")
  }

  if (ci["upper"] < threshold) {
    return("fail")
  }

  return("continue")
}


#' Run Single Consistency Split
#'
#' Performs one random 50/50 split and evaluates whether both halves
#' meet the HR consistency threshold.
#'
#' @param df.x data.table. Subgroup data with columns Y, Event, Treat.
#' @param N.x Integer. Number of observations in subgroup.
#' @param hr.consistency Numeric. Minimum HR threshold for consistency.
#' @param cox_init Numeric. Initial value for Cox model (log HR).
#'
#' @return Numeric. 1 if both splits meet threshold, 0 if not, NA if error.
#'
#' @keywords internal
#' @export
run_single_consistency_split <- function(df.x, N.x, hr.consistency, cox_init = 0) {

  in.split1 <- tryCatch({
    sample(c(TRUE, FALSE), N.x, replace = TRUE, prob = c(0.5, 0.5))
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(in.split1)) return(NA_real_)

  df.x$insplit1 <- in.split1
  df.x.split1 <- df.x[insplit1 == TRUE]
  df.x.split2 <- df.x[insplit1 == FALSE]

  if (nrow(df.x.split1) < 5 || nrow(df.x.split2) < 5) {
    return(NA_real_)
  }

  if (sum(df.x.split1$Event) < 2 || sum(df.x.split2$Event) < 2) {
    return(NA_real_)
  }

  hr.split1 <- get_split_hr_fast(df.x.split1, cox_init)
  hr.split2 <- get_split_hr_fast(df.x.split2, cox_init)

  if (!is.na(hr.split1) && !is.na(hr.split2)) {
    as.numeric(hr.split1 > hr.consistency && hr.split2 > hr.consistency)
  } else {
    NA_real_
  }
}


#' Fast Cox Model HR Estimation
#'
#' Fits a minimal Cox model to estimate hazard ratio with reduced overhead.
#'
#' @param df data.frame or data.table with Y, Event, Treat columns.
#' @param cox_init Numeric. Initial value for coefficient (default 0).
#'
#' @return Numeric. Estimated hazard ratio, or NA if model fails.
#'
#' @importFrom survival coxph Surv
#' @keywords internal
#' @export
get_split_hr_fast <- function(df, cox_init = 0) {
  if (nrow(df) < 2 || sum(df$Event) < 2) {
    return(NA_real_)
  }

  fit <- tryCatch(
    suppressWarnings(
      survival::coxph(
        survival::Surv(Y, Event) ~ Treat,
        data = df,
        init = cox_init,
        robust = FALSE,
        model = FALSE,
        x = FALSE,
        y = FALSE
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NA_real_)
  return(exp(fit$coefficients[1]))
}


#' Two-Stage Sequential Consistency Evaluation
#'
#' Evaluates a single subgroup for consistency using a two-stage approach
#' with sequential early stopping.
#'
#' @param m Integer. Index of subgroup to evaluate.
#' @param index.Z data.table or matrix. Factor indicators for all subgroups.
#' @param names.Z Character vector. Names of factor columns.
#' @param df data.frame. Original data with Y, Event, Treat, id columns.
#' @param found.hrs data.table. Subgroup hazard ratio results.
#' @param hr.consistency Numeric. Minimum HR threshold for consistency.
#' @param pconsistency.threshold Numeric. Final consistency threshold.
#' @param pconsistency.digits Integer. Rounding digits for output.
#' @param maxk Integer. Maximum number of factors in a subgroup.
#' @param confs_labels Character vector. Labels for confounders.
#' @param details Logical. Print progress details.
#' @param n.splits.screen Integer. Number of splits for Stage 1.
#' @param screen.threshold Numeric. Screening threshold for Stage 1.
#' @param n.splits.max Integer. Maximum total splits.
#' @param batch.size Integer. Splits per batch in Stage 2.
#' @param conf.level Numeric. Confidence level for early stopping.
#' @param min.valid.screen Integer. Minimum valid splits in Stage 1.
#'
#' @return Named numeric vector with consistency results, or NULL if not met.
#'
#' @importFrom data.table data.table
#' @importFrom survival coxph Surv
#' @export
evaluate_consistency_twostage <- function(
    m,
    index.Z,
    names.Z,
    df,
    found.hrs,
    hr.consistency,
    pconsistency.threshold,
    pconsistency.digits = 2,
    maxk,
    confs_labels,
    details = FALSE,
    n.splits.screen = 30,
    screen.threshold = NULL,
    n.splits.max = 400,
    batch.size = 20,
    conf.level = 0.95,
    min.valid.screen = 10
) {

  # ===========================================================================
  # EMBEDDED HELPER FUNCTIONS (for parallel execution compatibility)
  # These are defined locally to ensure availability in callr workers
  # ===========================================================================

  # Wilson score confidence interval
  .wilson_ci <- function(x, n, conf.level = 0.95) {
    if (n == 0) {
      return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
    }
    z <- qnorm(1 - (1 - conf.level) / 2)
    p_hat <- x / n
    denom <- 1 + z^2 / n
    center <- (p_hat + z^2 / (2 * n)) / denom
    margin <- (z / denom) * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
    c(estimate = p_hat, lower = max(0, center - margin), upper = min(1, center + margin))
  }

  # Early stopping decision
  .early_stop_decision <- function(n_success, n_total, threshold, conf.level = 0.95, min_samples = 20) {
    if (n_total < min_samples) return("continue")
    ci <- .wilson_ci(n_success, n_total, conf.level)
    if (is.na(ci["lower"]) || is.na(ci["upper"])) return("continue")
    if (ci["lower"] >= threshold) return("pass")
    if (ci["upper"] < threshold) return("fail")
    return("continue")
  }

  # Fast HR calculation from split
  .get_split_hr_fast <- function(df_split, cox_initial = NULL) {
    if (nrow(df_split) < 2 || sum(df_split$Event) < 2) return(NA_real_)
    fit <- tryCatch(
      suppressWarnings(
        survival::coxph(
          survival::Surv(Y, Event) ~ Treat,
          data = df_split,
          init = cox_initial,
          robust = FALSE,
          model = FALSE,
          x = FALSE,
          y = FALSE
        )
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    return(exp(fit$coefficients[1]))
  }

  # Run single consistency split
  .run_single_consistency_split <- function(df.x, N.x, hr.cons, cox_init) {
    in.split1 <- tryCatch({
      sample(c(TRUE, FALSE), N.x, replace = TRUE, prob = c(0.5, 0.5))
    }, error = function(e) NULL)
    if (is.null(in.split1)) return(NA_real_)

    df.x$insplit1 <- in.split1
    df.x.split1 <- df.x[df.x$insplit1 == TRUE, ]
    df.x.split2 <- df.x[df.x$insplit1 == FALSE, ]

    if (nrow(df.x.split1) < 5 || nrow(df.x.split2) < 5) return(NA_real_)
    if (sum(df.x.split1$Event) < 2 || sum(df.x.split2$Event) < 2) return(NA_real_)

    hr.split1 <- .get_split_hr_fast(df_split = df.x.split1, cox_initial = cox_init)
    hr.split2 <- .get_split_hr_fast(df_split = df.x.split2, cox_initial = cox_init)

    if (!is.na(hr.split1) && !is.na(hr.split2)) {
      as.numeric(hr.split1 > hr.cons && hr.split2 > hr.cons)
    } else {
      NA_real_
    }
  }

  # FS_labels helper (embedded version)
  # Converts q-indexed codes to human-readable labels
  # Pattern: q<index>.<action> where action 0 = NOT, action 1 = IN
  .FS_labels <- function(Qsg, confs_labels) {
    pattern <- "^q(\\d+)\\.(\\d)$"
    matches <- regmatches(Qsg, regexec(pattern, Qsg))[[1]]
    if (length(matches) < 3) return(Qsg)
    idx <- as.integer(matches[2])
    action <- matches[3]
    if (idx < 1 || idx > length(confs_labels)) return(Qsg)
    base_label <- confs_labels[idx]
    # Match original FS_labels behavior: wrap with {..} or !{..}
    if (action == "0") {
      paste0("!{", base_label, "}")
    } else {
      paste0("{", base_label, "}")
    }
  }

  # ---------------------------------------------------------------------------
  # Parameter initialization
  # ---------------------------------------------------------------------------

  if (is.null(screen.threshold)) {
    se_estimate <- sqrt(pconsistency.threshold * (1 - pconsistency.threshold) / n.splits.screen)
    screen.threshold <- max(0.5, pconsistency.threshold - 2.5 * se_estimate)
  }

  # ---------------------------------------------------------------------------
  # Validate and extract subgroup
  # ---------------------------------------------------------------------------

  if (m < 1 || m > nrow(index.Z)) {
    warning("Invalid subgroup index m=", m, ". Skipping.")
    return(NULL)
  }

  indexm <- as.numeric(unlist(index.Z[m, ]))

  if (length(indexm) != length(names.Z)) {
    warning("Subgroup ", m, ": index length mismatch. Skipping.")
    return(NULL)
  }

  if (!all(indexm %in% c(0, 1, NA))) {
    warning("Subgroup ", m, ": invalid index values. Skipping.")
    return(NULL)
  }

  this.m <- names.Z[indexm == 1]

  if (length(this.m) == 0) {
    warning("Subgroup ", m, ": no factors selected. Skipping.")
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Convert labels
  # ---------------------------------------------------------------------------

  this.m_label <- tryCatch({
    unlist(lapply(this.m, .FS_labels, confs_labels = confs_labels))
  }, error = function(e) {
    warning("Subgroup ", m, ": error in FS_labels: ", e$message, ". Skipping.")
    return(NULL)
  })

  if (is.null(this.m_label) || length(this.m_label) != length(this.m)) {
    warning("Subgroup ", m, ": label conversion failed. Skipping.")
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Extract subgroup data
  # ---------------------------------------------------------------------------

  id.m <- paste(paste(this.m, collapse = "==1 & "), "==1")

  df.sub <- tryCatch({
    subset(df, eval(parse(text = id.m)))
  }, error = function(e) {
    warning("Subgroup ", m, ": error extracting data: ", e$message, ". Skipping.")
    return(NULL)
  })

  if (is.null(df.sub) || nrow(df.sub) == 0) {
    warning("Subgroup ", m, ": no observations match criteria. Skipping.")
    return(NULL)
  }

  df.x <- data.table::data.table(df.sub)
  N.x <- nrow(df.x)

  cox_init <- log(found.hrs$HR[m])
  if (is.na(cox_init) || is.infinite(cox_init)) {
    cox_init <- 0
  }

  # ---------------------------------------------------------------------------
  # Stage 1: Quick screening
  # ---------------------------------------------------------------------------

  if (details) {
    cat("Subgroup ", m, ": Stage 1 (", n.splits.screen, " splits)\n", sep = "")
  }

  stage1_flags <- numeric(n.splits.screen)
  for (i in seq_len(n.splits.screen)) {
    stage1_flags[i] <- .run_single_consistency_split(df.x, N.x, hr.consistency, cox_init)
  }

  n_valid_stage1 <- sum(!is.na(stage1_flags))
  n_success_stage1 <- sum(stage1_flags == 1, na.rm = TRUE)

  # Screen out clearly non-viable candidates
  if (n_valid_stage1 >= min.valid.screen) {
    p_screen <- n_success_stage1 / n_valid_stage1

    if (p_screen < screen.threshold) {
      if (details) {
        cat("Subgroup ", m, ": SCREENED OUT (Pcons=", round(p_screen, 3),
            " < ", round(screen.threshold, 3), ")\n", sep = "")
      }
      return(NULL)
    }
  }

  # ---------------------------------------------------------------------------
  # Stage 2: Sequential evaluation with early stopping
  # ---------------------------------------------------------------------------

  if (details) {
    cat("Subgroup ", m, ": Stage 2 sequential\n", sep = "")
  }

  all_flags <- stage1_flags
  n_total_valid <- n_valid_stage1
  n_total_success <- n_success_stage1

  n_remaining <- n.splits.max - n.splits.screen
  n_batches <- ceiling(n_remaining / batch.size)

  final_decision <- "continue"

  for (batch_num in seq_len(n_batches)) {

    # Check for early stopping
    decision <- .early_stop_decision(
      n_success = n_total_success,
      n_total = n_total_valid,
      threshold = pconsistency.threshold,
      conf.level = conf.level,
      min_samples = max(20, min.valid.screen)
    )

    if (decision != "continue") {
      final_decision <- decision
      if (details) {
        cat("Subgroup ", m, ": EARLY ", toupper(decision),
            " at n=", n_total_valid, "\n", sep = "")
      }
      break
    }

    # Run next batch
    batch_flags <- numeric(batch.size)
    for (i in seq_len(batch.size)) {
      batch_flags[i] <- .run_single_consistency_split(df.x, N.x, hr.consistency, cox_init)
    }

    n_batch_valid <- sum(!is.na(batch_flags))
    n_batch_success <- sum(batch_flags == 1, na.rm = TRUE)

    n_total_valid <- n_total_valid + n_batch_valid
    n_total_success <- n_total_success + n_batch_success
  }

  # ---------------------------------------------------------------------------
  # Final evaluation
  # ---------------------------------------------------------------------------

  if (final_decision == "continue") {
    final_decision <- .early_stop_decision(
      n_success = n_total_success,
      n_total = n_total_valid,
      threshold = pconsistency.threshold,
      conf.level = conf.level,
      min_samples = 20
    )
  }

  p.consistency <- tryCatch({
    round(n_total_success / n_total_valid, pconsistency.digits)
  }, error = function(e) {
    return(NA_real_)
  })

  if (is.na(p.consistency)) {
    return(NULL)
  }

  if (final_decision == "fail" || p.consistency < pconsistency.threshold) {
    if (details) {
      cat("Subgroup ", m, ": FAILED (Pcons=", p.consistency, ")\n", sep = "")
    }
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Format and return result
  # ---------------------------------------------------------------------------

  k <- length(this.m)
  covsm <- rep("M", maxk)
  mindex <- seq_len(maxk)
  Mnames <- paste(covsm, mindex, sep = ".")

  mfound <- matrix(rep("", maxk))
  mfound[seq_len(k)] <- this.m_label

  resultk <- c(
    p.consistency,
    found.hrs$HR[m],
    found.hrs$n[m],
    found.hrs$E[m],
    found.hrs$grp[m],
    m,
    k,
    mfound
  )

  names(resultk) <- c("Pcons", "hr", "N", "E", "g", "m", "K", Mnames)

  if (details) {
    cat("Subgroup ", m, ": PASSED (Pcons=", p.consistency,
        ", n_splits=", n_total_valid, ")\n", sep = "")
  }

  return(resultk)
}
