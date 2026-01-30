#' ForestSearch K-Fold Cross-Validation
#'
#' Performs K-fold cross-validation for ForestSearch, evaluating subgroup
#' identification and agreement between training and test sets.
#'
#' @description
#' This function assesses the stability and reproducibility of ForestSearch
#' subgroup identification through cross-validation. For each fold:
#' \enumerate{
#'   \item Train ForestSearch on (K-1) folds
#'   \item Apply the identified subgroup to the held-out fold
#'   \item Compare predictions to the original full-data analysis
#' }
#'
#' @section Cross-Validation Types:
#' \itemize{
#'   \item \strong{Leave-One-Out (LOO)}: When \code{Kfolds = nrow(df)}, each
#'     observation is held out once. Most thorough but computationally intensive.
#'   \item \strong{K-Fold}: When \code{Kfolds < nrow(df)}, data is split into K
#'     roughly equal folds. Good balance of bias-variance tradeoff.
#' }
#'
#' @section Output Metrics:
#' The returned \code{resCV} data frame contains:
#' \itemize{
#'   \item \code{treat.recommend}: Prediction from CV model
#'   \item \code{treat.recommend.original}: Prediction from full-data model
#'   \item \code{cvindex}: Fold assignment
#'   \item \code{sg1}, \code{sg2}: Subgroup definitions found in each fold
#' }
#'
#' @param fs.est List. ForestSearch results object from \code{\link{forestsearch}}.
#'   Must contain \code{df.est} (data frame) and \code{args_call_all} (list of arguments).
#' @param Kfolds Integer. Number of folds (default: \code{nrow(fs.est$df.est)} for LOO).
#' @param seedit Integer. Random seed for fold assignment (default: 8316951).
#' @param parallel_args List. Parallelization configuration with elements:
#'   \itemize{
#'     \item \code{plan}: Character. One of "multisession", "multicore", "sequential"
#'     \item \code{workers}: Integer. Number of parallel workers
#'     \item \code{show_message}: Logical. Show parallel setup messages
#'   }
#' @param sg0.name Character. Label for subgroup 0 (default: "Not recommend").
#' @param sg1.name Character. Label for subgroup 1 (default: "Recommend").
#' @param details Logical. Print progress details (default: FALSE).
#'
#' @return List with components:
#' \describe{
#'   \item{resCV}{Data frame with CV predictions for each observation}
#'   \item{cv_args}{Arguments used for CV ForestSearch calls}
#'   \item{timing_minutes}{Execution time in minutes}
#'   \item{prop_SG_found}{Percentage of folds where a subgroup was found}
#'   \item{sg_analysis}{Original subgroup definition from full-data analysis}
#'   \item{sg0.name, sg1.name}{Subgroup labels}
#'   \item{Kfolds}{Number of folds used}
#' }
#'
#' @examples
#' \dontrun
#' # Run initial ForestSearch
#' fs_result <- forestsearch(
#'   df.analysis = trial_data,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment",
#'   confounders.name = c("age", "biomarker")
#' )
#'
#' # Run 10-fold cross-validation
#' cv_results <- forestsearch_Kfold(
#'   fs.est = fs_result,
#'   Kfolds = 10,
#'   parallel_args = list(plan = "multisession", workers = 4),
#'   details = TRUE
#' )
#'
#' # Summarize results
#' cv_summary <- forestsearch_KfoldOut(cv_results, outall = TRUE)
#' }
#'
#' @seealso
#' \code{\link{forestsearch}} for initial subgroup identification
#' \code{\link{forestsearch_KfoldOut}} for summarizing CV results
#' \code{\link{forestsearch_tenfold}} for repeated K-fold simulations
#'
#' @importFrom future plan
#' @importFrom data.table data.table copy
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @export

forestsearch_Kfold <- function(
    fs.est,
    Kfolds = nrow(fs.est$df.est),
    seedit = 8316951,
    parallel_args = list(plan = "multisession", workers = 6, show_message = TRUE),
    sg0.name = "Not recommend",
    sg1.name = "Recommend",
    details = FALSE
) {


  # ===========================================================================
 # SECTION 1: INPUT VALIDATION
 # ===========================================================================

  if (!is.list(fs.est)) {
    stop("fs.est must be a list (ForestSearch results object)")
  }

  required_components <- c("df.est", "args_call_all", "sg.harm")
  missing <- setdiff(required_components, names(fs.est))
  if (length(missing) > 0) {
    stop("fs.est missing required components: ", paste(missing, collapse = ", "))
  }

  if (is.null(fs.est$df.est) || !is.data.frame(fs.est$df.est) || nrow(fs.est$df.est) == 0) {
    stop("fs.est$df.est must be a non-empty data frame")
  }

  if (!is.numeric(Kfolds) || length(Kfolds) != 1 || Kfolds < 2) {
    stop("Kfolds must be an integer >= 2")
  }

  if (Kfolds > nrow(fs.est$df.est)) {
    stop("Kfolds cannot exceed number of observations (", nrow(fs.est$df.est), ")")
  }

  # ===========================================================================
  # SECTION 2: PARALLEL SETUP
  # ===========================================================================

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  # Resolve parallel arguments (similar to bootstrap function)
  resolved_args <- resolve_cv_parallel_args(parallel_args, fs.est$args_call_all, details)
  setup_parallel_SGcons(resolved_args)

  t_start <- proc.time()[3]

  # ===========================================================================
  # SECTION 3: DATA PREPARATION
  # ===========================================================================

  fs_args <- fs.est$args_call_all

  # Extract variable names
  outcome.name <- fs_args$outcome.name
  event.name <- fs_args$event.name
  treat.name <- fs_args$treat.name
  id.name <- fs_args$id.name
  est.scale <- fs_args$est.scale

  get_names <- c(fs_args$confounders.name, outcome.name, event.name, id.name, treat.name)

  # Prepare analysis data
  dfa <- fs.est$df.est[, c(get_names, "treat.recommend")]
  names(dfa)[names(dfa) == "treat.recommend"] <- "treat.recommend.original"
  dfnew <- as.data.frame(dfa)

  # ===========================================================================
  # SECTION 4: FOLD ASSIGNMENT
  # ===========================================================================

  df_scrambled <- data.table::copy(dfnew)

  # Scramble only if not doing LOO
 if (Kfolds < nrow(fs.est$df.est)) {
    set.seed(seedit)
    id_sample <- sample(seq_len(nrow(df_scrambled)), replace = FALSE)
    df_scrambled <- df_scrambled[id_sample, ]
  }

  folds <- cut(seq_len(nrow(df_scrambled)), breaks = Kfolds, labels = FALSE)

  if (details) {
    cat("Cross-validation setup:\n")
    cat("  - Observations:", nrow(df_scrambled), "\n")
    cat("  - Folds:", Kfolds, "\n")
    cat("  - Fold sizes (range):", paste(range(table(folds)), collapse = "-"), "\n")
  }

  # ===========================================================================
  # SECTION 5: CONFIGURE CV FORESTSEARCH ARGUMENTS
  # ===========================================================================

  cv_args <- fs_args
  cv_args$parallel_args <- list(plan = "sequential", workers = 1, show_message = FALSE)
  cv_args$details <- FALSE
  cv_args$plot.sg <- FALSE

  # ===========================================================================
  # SECTION 6: RUN CROSS-VALIDATION (PARALLEL)
  # ===========================================================================

  resCV <- foreach::foreach(
    cv_index = seq_len(Kfolds),
    .options.future = list(seed = TRUE),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {

    # Split data
    testIndexes <- which(folds == cv_index, arr.ind = TRUE)
    x.test <- df_scrambled[testIndexes, ]
    x.train <- df_scrambled[-testIndexes, ]

    # Update training data
    cv_args$df.analysis <- x.train

    # Run ForestSearch on training fold
    fs.train <- suppressWarnings(try(do.call(forestsearch, cv_args), TRUE))

    # Process results
    if (!inherits(fs.train, "try-error") && !is.null(fs.train$sg.harm)) {
      # Subgroup found - apply to test data
      sg1 <- fs.train$sg.harm[1]
      sg2 <- if (length(fs.train$sg.harm) > 1) fs.train$sg.harm[2] else NA

      df.test <- get_dfpred(df.predict = x.test, sg.harm = fs.train$sg.harm, version = 2)
      df.test$cvindex <- cv_index
      df.test$sg1 <- sg1
      df.test$sg2 <- sg2

    } else {
      # No subgroup found - default to ITT (recommend all)
      df.test <- x.test
      df.test$cvindex <- cv_index
      df.test$sg1 <- NA
      df.test$sg2 <- NA
      df.test$treat.recommend <- 1.0
    }

    # Return standardized columns
    data.table::data.table(
      df.test[, c(id.name, outcome.name, event.name, treat.name,
                  "treat.recommend", "treat.recommend.original",
                  "cvindex", "sg1", "sg2")]
    )
  }

  # ===========================================================================
  # SECTION 7: POST-PROCESSING AND VALIDATION
  # ===========================================================================

  t_end <- proc.time()[3]
  t_min <- (t_end - t_start) / 60

  # Validate results
  if (length(unique(resCV[[id.name]])) != nrow(df_scrambled)) {
    stop("K-fold cross-validation sample size does not equal full analysis dataset")
  }

  chk1 <- range(table(folds))
  chk2 <- range(table(resCV$cvindex))
  if (!all(chk1 == chk2)) {
    stop("Mismatch on range of fold sample sizes")
  }

  # Convert to data.frame
  resCV <- as.data.frame(resCV)

  # Handle est.scale for treatment direction
  if (est.scale == "1/hr") {
    resCV$treat2 <- 1 - resCV[, treat.name]
    treat.name <- "treat2"
    sg0.name <- "Recommend"
    sg1.name <- "Not recommend"
  }

  # Calculate proportion of folds where subgroup was found
  sg_found_count <- sum(!is.na(resCV$sg1) | !is.na(resCV$sg2))
  # Adjust for observations per fold
  unique_folds_with_sg <- length(unique(resCV$cvindex[!is.na(resCV$sg1) | !is.na(resCV$sg2)]))
  prop_SG_found <- 100 * round(unique_folds_with_sg / Kfolds, 3)

  if (details) {
    cat("\nCross-validation complete:\n")
    cat("  - Time:", round(t_min, 2), "minutes\n")
    cat("  - Subgroup found in", prop_SG_found, "% of folds\n")
  }

  # ===========================================================================
  # SECTION 8: RETURN RESULTS
  # ===========================================================================

  result <- list(
    resCV = resCV,
    cv_args = cv_args,
    timing_minutes = t_min,
    prop_SG_found = prop_SG_found,
    sg_analysis = fs.est$sg.harm,
    sg0.name = sg0.name,
    sg1.name = sg1.name,
    Kfolds = Kfolds
  )

  class(result) <- c("fs_kfold", "list")
  return(result)
}


#' ForestSearch Repeated K-Fold Cross-Validation
#'
#' Runs repeated K-fold cross-validation simulations for ForestSearch and
#' summarizes subgroup identification stability across repetitions.
#'
#' @description
#' This function performs multiple independent K-fold cross-validations to
#' assess the variability in subgroup identification. Each simulation:
#' \enumerate{
#'   \item Randomly shuffles the data
#'   \item Performs K-fold CV
#'   \item Records sensitivity and agreement metrics
#' }
#' Results are summarized across all simulations.
#'
#' @section Parallelization Strategy:
#' Unlike the single K-fold function which parallelizes across folds, this
#' function parallelizes across simulations for better efficiency when
#' running many repetitions. Each simulation runs its K-fold CV sequentially.
#'
#' @param fs.est List. ForestSearch results object from \code{\link{forestsearch}}.
#' @param sims Integer. Number of simulation repetitions.
#' @param Kfolds Integer. Number of folds per simulation (default: 10).
#' @param details Logical. Print progress details (default: TRUE).
#' @param parallel_args List. Parallelization configuration.
#'
#' @return List with components:
#' \describe{
#'   \item{sens_summary}{Named vector of median sensitivity metrics across simulations}
#'   \item{find_summary}{Named vector of median subgroup-finding metrics}
#'   \item{sens_out}{Matrix of sensitivity metrics (sims x metrics)}
#'   \item{find_out}{Matrix of finding metrics (sims x metrics)}
#'   \item{timing_minutes}{Total execution time}
#'   \item{sims}{Number of simulations run}
#'   \item{Kfolds}{Number of folds per simulation}
#' }
#'
#' @examples
#' \dontrun{
#' # Run 100 repetitions of 10-fold CV
#' tenfold_results <- forestsearch_tenfold(
#'   fs.est = fs_result,
#'   sims = 100,
#'   Kfolds = 10,
#'   parallel_args = list(plan = "multisession", workers = 6),
#'   details = TRUE
#' )
#'
#' # View summary
#' print(tenfold_results$sens_summary)
#' print(tenfold_results$find_summary)
#' }
#'
#' @seealso
#' \code{\link{forestsearch_Kfold}} for single K-fold CV
#' \code{\link{forestsearch_KfoldOut}} for summarizing CV results
#'
#' @importFrom future plan
#' @importFrom data.table data.table copy
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @export

forestsearch_tenfold <- function(
    fs.est,
    sims,
    Kfolds = 10,
    details = TRUE,
    parallel_args = list(plan = "multisession", workers = 6, show_message = TRUE)
) {

  # ===========================================================================
  # SECTION 1: INPUT VALIDATION
  # ===========================================================================

  if (!is.list(fs.est)) {
    stop("fs.est must be a list (ForestSearch results object)")
  }

  required_components <- c("df.est", "args_call_all", "sg.harm")
  missing <- setdiff(required_components, names(fs.est))
  if (length(missing) > 0) {
    stop("fs.est missing required components: ", paste(missing, collapse = ", "))
  }

  if (!is.numeric(sims) || length(sims) != 1 || sims < 1) {
    stop("sims must be a positive integer")
  }

  if (!is.numeric(Kfolds) || length(Kfolds) != 1 || Kfolds < 2) {
    stop("Kfolds must be an integer >= 2")
  }

  # ===========================================================================
  # SECTION 2: PARALLEL SETUP
  # ===========================================================================

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  resolved_args <- resolve_cv_parallel_args(parallel_args, fs.est$args_call_all, details)
  setup_parallel_SGcons(resolved_args)

  t_start <- proc.time()[3]

  if (details) {
    cat("Starting repeated K-fold cross-validation:\n")
    cat("  - Simulations:", sims, "\n")
    cat("  - Folds per simulation:", Kfolds, "\n")
    cat("  - Workers:", resolved_args$workers, "\n")
  }

  # ===========================================================================
  # SECTION 3: DATA PREPARATION
  # ===========================================================================

  fs_args <- fs.est$args_call_all

  # Extract variable names
  outcome.name <- fs_args$outcome.name
  event.name <- fs_args$event.name
  treat.name <- fs_args$treat.name
  id.name <- fs_args$id.name
  confounders.name <- fs_args$confounders.name

  get_names <- c(confounders.name, outcome.name, event.name, id.name, treat.name)

  # Prepare base data
  dfa <- fs.est$df.est[, c(get_names, "treat.recommend")]
  names(dfa)[names(dfa) == "treat.recommend"] <- "treat.recommend.original"
  dfnew <- as.data.frame(dfa)

  # Configure CV arguments (sequential within each simulation)
  cv_args <- fs_args
  cv_args$parallel_args <- list(plan = "sequential", workers = 1, show_message = FALSE)
  cv_args$details <- FALSE
  cv_args$plot.sg <- FALSE

  # ===========================================================================
  # SECTION 4: RUN SIMULATIONS (PARALLEL ACROSS SIMS)
  # ===========================================================================

  simulation_results <- foreach::foreach(
    ksim = seq_len(sims),
    .options.future = list(seed = TRUE),
    .errorhandling = "pass"
  ) %dofuture% {

    t_sim_start <- proc.time()[3]

    # Shuffle data for this simulation
    df_scrambled <- data.table::copy(dfnew)
    set.seed(8316951 + 1000 * ksim)
    id_sample <- sample(seq_len(nrow(df_scrambled)), replace = FALSE)
    df_scrambled <- df_scrambled[id_sample, ]

    # Create fold assignments
    folds <- cut(seq_len(nrow(df_scrambled)), breaks = Kfolds, labels = FALSE)

    # Run K-fold CV sequentially within this simulation
    resCV_list <- vector("list", Kfolds)

    for (cv_index in seq_len(Kfolds)) {
      testIndexes <- which(folds == cv_index, arr.ind = TRUE)
      x.test <- df_scrambled[testIndexes, ]
      x.train <- df_scrambled[-testIndexes, ]

      cv_args$df.analysis <- x.train

      fs.train <- suppressWarnings(try(do.call(forestsearch, cv_args), TRUE))

      if (!inherits(fs.train, "try-error") && !is.null(fs.train$sg.harm)) {
        sg1 <- fs.train$sg.harm[1]
        sg2 <- if (length(fs.train$sg.harm) > 1) fs.train$sg.harm[2] else NA

        df.test <- get_dfpred(df.predict = x.test, sg.harm = fs.train$sg.harm, version = 2)
        df.test$cvindex <- cv_index
        df.test$sg1 <- sg1
        df.test$sg2 <- sg2

      } else {
        df.test <- x.test
        df.test$cvindex <- cv_index
        df.test$sg1 <- NA
        df.test$sg2 <- NA
        df.test$treat.recommend <- 1.0
      }

      resCV_list[[cv_index]] <- df.test[, c(id.name, outcome.name, event.name, treat.name,
                                             "treat.recommend", "treat.recommend.original",
                                             "cvindex", "sg1", "sg2")]
    }

    resCV <- do.call(rbind, resCV_list)
    resCV <- as.data.frame(resCV)

    # Summarize this simulation
    res <- list(
      resCV = resCV,
      Kfolds = Kfolds,
      confounders.name = confounders.name,
      cv_args = cv_args,
      outcome.name = outcome.name,
      event.name = event.name,
      id.name = id.name,
      treat.name = treat.name,
      sg_analysis = fs.est$sg.harm,
      sg0.name = "Not recommend",
      sg1.name = "Recommend"
    )

    out <- forestsearch_KfoldOut(res = res, outall = FALSE, details = FALSE)

    list(
      sens_metrics_original = out$sens_metrics_original,
      find_metrics = out$find_metrics,
      sim_id = ksim
    )
  }

  # ===========================================================================
  # SECTION 5: COMBINE AND SUMMARIZE RESULTS
  # ===========================================================================

  # Extract metrics from successful simulations
  valid_results <- Filter(function(x) !inherits(x, "error") && is.list(x), simulation_results)

  if (length(valid_results) == 0) {
    stop("All simulations failed. Check ForestSearch configuration.")
  }

  if (length(valid_results) < sims && details) {
    warning(sims - length(valid_results), " simulations failed and were excluded.")
  }

  sens_out <- do.call(rbind, lapply(valid_results, `[[`, "sens_metrics_original"))
  find_out <- do.call(rbind, lapply(valid_results, `[[`, "find_metrics"))

  # Compute summaries
  sens_summary <- apply(sens_out, 2, median, na.rm = TRUE)
  find_summary <- apply(find_out, 2, median, na.rm = TRUE)

  t_end <- proc.time()[3]
  t_min <- (t_end - t_start) / 60

  if (details) {
    cat("\nRepeated K-fold CV complete:\n")
    cat("  - Time:", round(t_min, 2), "minutes\n")
    cat("  - Successful simulations:", length(valid_results), "/", sims, "\n")
    cat("  - Projected hours per 100 sims:", round((t_min / 60) * (100 / sims), 2), "\n")
  }

  # ===========================================================================
  # SECTION 6: RETURN RESULTS
  # ===========================================================================

  result <- list(
    sens_summary = sens_summary,
    find_summary = find_summary,
    sens_out = sens_out,
    find_out = find_out,
    timing_minutes = t_min,
    sims = length(valid_results),
    Kfolds = Kfolds
  )

  class(result) <- c("fs_tenfold", "list")
  return(result)
}


#' Resolve Parallel Arguments for Cross-Validation
#'
#' Helper function to resolve and validate parallel processing arguments,
#' similar to bootstrap's \code{resolve_bootstrap_parallel_args}.
#'
#' @param parallel_args List. User-provided parallel arguments.
#' @param fs_args List. Original ForestSearch call arguments.
#' @param details Logical. Print configuration messages.
#'
#' @return List with resolved plan, workers, show_message.
#'
#' @keywords internal

resolve_cv_parallel_args <- function(parallel_args, fs_args, details = FALSE) {

  # Default values
  default_args <- list(
    plan = "multisession",
    workers = max(1, parallel::detectCores() - 1),
    show_message = TRUE
  )

  # Merge user args with defaults
  resolved_args <- default_args

  if (length(parallel_args) > 0) {
    for (nm in names(parallel_args)) {
      resolved_args[[nm]] <- parallel_args[[nm]]
    }
  } else if (!is.null(fs_args$parallel_args) && length(fs_args$parallel_args) > 0) {
    # Fall back to ForestSearch settings
    for (nm in names(fs_args$parallel_args)) {
      resolved_args[[nm]] <- fs_args$parallel_args[[nm]]
    }
    if (details) {
      message("CV parallel config: Using ForestSearch analysis settings")
    }
  }

  # Validate workers
  max_cores <- parallel::detectCores()
  if (resolved_args$workers > max_cores) {
    warning("Requested workers (", resolved_args$workers, ") exceeds available cores (",
            max_cores, "). Using ", max_cores - 1, " workers.")
    resolved_args$workers <- max(1, max_cores - 1)
  }

  if (details && resolved_args$show_message) {
    message("CV will use: ", resolved_args$workers, " workers with '",
            resolved_args$plan, "' plan")
  }

  resolved_args
}


#' Print Method for K-Fold CV Results
#'
#' @param x An fs_kfold object
#' @param ... Additional arguments (ignored)
#' @export

print.fs_kfold <- function(x, ...) {
  cat("ForestSearch K-Fold Cross-Validation Results\n")
  cat("=============================================\n")
  cat("Folds:", x$Kfolds, "\n")
  cat("Observations:", nrow(x$resCV), "\n")
  cat("Subgroup found in:", x$prop_SG_found, "% of folds\n")
  cat("Original subgroup:", paste(x$sg_analysis, collapse = " & "), "\n")
  cat("Time:", round(x$timing_minutes, 2), "minutes\n")
  cat("\nUse forestsearch_KfoldOut() for detailed metrics.\n")
  invisible(x)
}


#' Print Method for Repeated K-Fold CV Results
#'
#' @param x An fs_tenfold object
#' @param ... Additional arguments (ignored)
#' @export

print.fs_tenfold <- function(x, ...) {
  cat("ForestSearch Repeated K-Fold Cross-Validation Results\n")
  cat("======================================================\n")
  cat("Simulations:", x$sims, "\n")
  cat("Folds per simulation:", x$Kfolds, "\n")
  cat("Time:", round(x$timing_minutes, 2), "minutes\n")
  cat("\nSensitivity Summary (median across simulations):\n")
  print(round(x$sens_summary, 3))
  cat("\nSubgroup Finding Summary (median across simulations):\n")
  print(round(x$find_summary, 3))
  invisible(x)
}
