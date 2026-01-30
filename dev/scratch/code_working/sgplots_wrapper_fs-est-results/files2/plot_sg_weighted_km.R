# =============================================================================
# plot_sg_weighted_km.R - Weighted Kaplan-Meier Plots for ForestSearch Subgroups
# =============================================================================
#
# This file provides visualization functions for ForestSearch subgroup results
# using the weightedsurv package for weighted Kaplan-Meier survival curves.
#
# Author: ForestSearch Development Team
# License: GPL-3
# =============================================================================


#' Plot Weighted Kaplan-Meier Curves for ForestSearch Subgroups
#'
#' Creates publication-quality weighted Kaplan-Meier survival curves for the
#' identified subgroups (H and H^c) using the \pkg{weightedsurv} package.
#' Optionally displays bootstrap bias-corrected hazard ratio estimates from
#' \code{\link{forestsearch_bootstrap_dofuture}} results.
#'
#' @param fs.est A forestsearch object or list containing:
#'   \describe{
#'     \item{df.est}{Data frame with analysis data including \code{treat.recommend}}
#'     \item{sg.harm}{Character vector of subgroup-defining variable names}
#'     \item{grp.consistency}{Optional. Consistency results with \code{out_sg}}
#'   }
#' @param fs_bc Optional. Bootstrap results from \code{\link{forestsearch_bootstrap_dofuture}}.
#'   If provided, displays bias-corrected hazard ratio estimates on the plots.
#'   Expected structure includes:
#'   \describe{
#'     \item{SG_CIs}{List with \code{H_bc}, \code{Hc_bc} (formatted CI strings)}
#'     \item{H_estimates}{List with \code{H0}, \code{H2}, and their CIs (numeric)}
#'     \item{Hc_estimates}{List with \code{H0}, \code{H2}, and their CIs (numeric)}
#'   }
#' @param outcome.name Character. Name of time-to-event outcome column.
#'   Default: \code{"Y"}
#' @param event.name Character. Name of event indicator column (1=event, 0=censored).
#'   Default: \code{"Event"}
#' @param treat.name Character. Name of treatment column (1=treatment, 0=control).
#'   Default: \code{"Treat"}
#' @param by.risk Numeric. Time interval for risk table display. If \code{NULL},
#'   automatically calculated as \code{max(outcome) / 12}. Default: \code{NULL}
#' @param est.scale Character. Effect scale for HR interpretation:
#'   \describe{
#'     \item{"hr"}{Standard hazard ratio (HR > 1 favors control)}
#'     \item{"1/hr"}{Inverse hazard ratio (treatment labels reversed)}
#'   }
#'   Default: \code{"hr"}
#' @param sg0_name Character. Label for H subgroup (harm/questionable).
#'   Default: \code{"Questionable (H)"}
#' @param sg1_name Character. Label for H^c subgroup (complement/recommend).
#'   Default: \code{"Recommend (Hc)"}
#' @param treat_labels Character vector of length 2. Labels for treatment arms
#'   in order \code{c(control_label, treatment_label)}.
#'   Default: \code{c("control", "treat")}
#' @param conf.int Logical. Show confidence intervals on survival curves.
#'   Default: \code{TRUE}
#' @param show.logrank Logical. Show log-rank test p-value. Default: \code{TRUE}
#' @param show_hr Logical. Show hazard ratio estimates on plots. Default: \code{TRUE}
#' @param show_bc_hr Logical. If \code{TRUE} and \code{fs_bc} is provided, show
#'   bias-corrected HR instead of (or in addition to) observed HR.
#'   Default: \code{TRUE}
#' @param hr_position Character. Position for HR annotation: \code{"topleft"},
#'   \code{"topright"}, \code{"bottomleft"}, or \code{"bottomright"}.
#'   Default: \code{"topright"}
#' @param put.legend.lr Character. Position for survival curve legend.
#'   Default: \code{"topleft"}
#' @param ymax Numeric. Maximum y-axis value for survival probability.
#'   Default: \code{1.05}
#' @param xmed.fraction Numeric. Fraction of x-axis for median survival lines.
#'   Default: \code{0.65}
#' @param title Character. Overall plot title. Default: \code{NULL} (auto-generated)
#' @param colors Named list. Colors for plot elements:
#'   \describe{
#'     \item{treat}{Color for treatment arm. Default: "#E41A1C"}
#'     \item{control}{Color for control arm. Default: "#377EB8"}
#'   }
#' @param layout Character. Plot layout: \code{"horizontal"} (side-by-side) or
#'   \code{"vertical"} (stacked). Default: \code{"horizontal"}
#' @param use_weightedsurv Logical. If \code{TRUE} (default), use
#'   \code{weightedsurv::plot_weighted_km()} for plotting. If \code{FALSE} or
#'   if weightedsurv package is unavailable, fall back to base
#'   \code{survival::survfit()} plotting. Default: \code{TRUE}
#' @param verbose Logical. Print diagnostic messages. Default: \code{FALSE}
#' @param ... Additional arguments passed to \code{weightedsurv::plot_weighted_km}.
#'
#' @return An object of class \code{fs_weighted_km} (invisibly) containing:
#'   \describe{
#'     \item{dfcount_H}{Counting process data for H subgroup}
#'     \item{dfcount_Hc}{Counting process data for Hc subgroup}
#'     \item{hr_observed}{Data frame with observed HR estimates}
#'     \item{hr_bc}{Data frame with bias-corrected HR estimates (if fs_bc provided)}
#'     \item{summary}{Summary statistics for both subgroups}
#'     \item{call}{The matched call}
#'   }
#'
#' @details
#' This function creates weighted Kaplan-Meier survival curves using the
#' \pkg{weightedsurv} package, which provides proper handling of tied event
#' times and supports inverse probability weighting.
#'
#' @section Subgroup Definition:
#' The function extracts subgroup membership from \code{fs.est$df.est$treat.recommend}:
#' \itemize{
#'   \item \code{treat.recommend == 0}: H subgroup (harm/questionable) -
#'     patients for whom treatment benefit is questionable
#'   \item \code{treat.recommend == 1}: H^c subgroup (complement/recommend) -
#'     patients for whom treatment is recommended
#' }
#'
#' @section Bootstrap Integration:
#' When \code{fs_bc} is provided from \code{\link{forestsearch_bootstrap_dofuture}},
#' the function can display bias-corrected hazard ratio estimates. Two types of
#' estimates are available:
#' \itemize{
#'   \item \strong{Observed (H0)}: Unadjusted HR from original analysis
#'   \item \strong{Bias-corrected (H2)}: HR adjusted for selection optimism using
#'     the double bootstrap method
#' }
#'
#' @section Dependencies:
#' This function requires the \pkg{weightedsurv} package for weighted KM estimation.
#' If not available, the function will stop with an informative error message.
#'
#' @examples
#' \dontrun{
#' library(ForestSearch)
#' library(survival)
#'
#' # Run ForestSearch analysis
#' fs <- forestsearch(
#'   df.analysis = my_data,
#'   outcome.name = "os_time",
#'   event.name = "os_event",
#'   treat.name = "treatment",
#'   confounders.name = c("age_cat", "stage", "biomarker")
#' )
#'
#' # Basic weighted KM plot
#' plot_sg_weighted_km(
#'   fs.est = fs,
#'   outcome.name = "os_time",
#'   event.name = "os_event",
#'   treat.name = "treatment"
#' )
#'
#' # With bootstrap bias-corrected estimates
#' fs_bc <- forestsearch_bootstrap_dofuture(
#'   fs.est = fs,
#'   nb_boots = 500,
#'   parallel_args = list(plan = "multisession", workers = 4)
#' )
#'
#' plot_sg_weighted_km(
#'   fs.est = fs,
#'   fs_bc = fs_bc,
#'   outcome.name = "os_time",
#'   event.name = "os_event",
#'   treat.name = "treatment",
#'   show_bc_hr = TRUE,
#'   sg0_name = "High Risk",
#'   sg1_name = "Standard Risk"
#' )
#'
#' # Customize appearance
#' result <- plot_sg_weighted_km(
#'   fs.est = fs,
#'   fs_bc = fs_bc,
#'   outcome.name = "os_time",
#'   event.name = "os_event",
#'   treat.name = "treatment",
#'   treat_labels = c("Placebo", "Active"),
#'   colors = list(treat = "darkred", control = "darkblue"),
#'   layout = "vertical",
#'   title = "Overall Survival by Subgroup"
#' )
#'
#' # Access returned data
#' result$hr_observed
#' result$hr_bc
#' }
#'
#' @seealso
#' \code{\link{forestsearch}} for running the subgroup analysis
#' \code{\link{forestsearch_bootstrap_dofuture}} for bootstrap bias correction
#' \code{\link{plot_sg_results}} for base R visualization
#' \code{\link{plot_subgroup_results_forestplot}} for publication forest plots
#'
#' @importFrom survival coxph Surv survdiff
#' @importFrom stats confint pchisq as.formula
#' @importFrom graphics par mtext legend text
#' @export
plot_sg_weighted_km <- function(
    fs.est,
    fs_bc = NULL,
    outcome.name = "Y",
    event.name = "Event",
    treat.name = "Treat",
    by.risk = NULL,
    est.scale = c("hr", "1/hr"),
    sg0_name = "Questionable (H)",
    sg1_name = "Recommend (Hc)",
    treat_labels = c("control", "treat"),
    conf.int = TRUE,
    show.logrank = TRUE,
    show_hr = TRUE,
    show_bc_hr = TRUE,
    hr_position = c("topright", "topleft", "bottomleft", "bottomright"),
    put.legend.lr = "topleft",
    ymax = 1.05,
    xmed.fraction = 0.65,
    title = NULL,
    colors = NULL,
    layout = c("horizontal", "vertical"),
    use_weightedsurv = TRUE,
    verbose = FALSE,
    ...
) {

 # ===========================================================================
 # SECTION 1: INPUT VALIDATION
 # ===========================================================================

  est.scale <- match.arg(est.scale)
  hr_position <- match.arg(hr_position)
  layout <- match.arg(layout)

  # Check for weightedsurv package availability
  weightedsurv_available <- requireNamespace("weightedsurv", quietly = TRUE)

  if (use_weightedsurv && !weightedsurv_available) {
    message(
      "Package 'weightedsurv' not available. ",
      "Using base survival::survfit() for plotting."
    )
    use_weightedsurv <- FALSE
  }

  # Validate fs.est structure
  if (!inherits(fs.est, "forestsearch") && !is.list(fs.est)) {
    stop("fs.est must be a forestsearch object or list with df.est component")
  }

  if (is.null(fs.est$df.est)) {
    stop("fs.est$df.est is NULL. No subgroup was identified.")
  }

  df <- fs.est$df.est

  # Check required columns
  required_cols <- c(outcome.name, event.name, treat.name, "treat.recommend")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in df.est: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Validate treat.recommend values
  if (!all(df$treat.recommend %in% c(0, 1, NA))) {
    stop("treat.recommend must contain only 0, 1, or NA values")
  }

  # Validate fs_bc structure if provided
  if (!is.null(fs_bc)) {
    if (!is.list(fs_bc)) {
      warning("fs_bc should be a list from forestsearch_bootstrap_dofuture()")
      fs_bc <- NULL
    } else if (is.null(fs_bc$H_estimates) && is.null(fs_bc$SG_CIs)) {
      warning("fs_bc does not contain expected H_estimates or SG_CIs components")
      fs_bc <- NULL
    }
  }

  if (verbose) {
    message("Input validation passed")
    message(sprintf("  N total: %d", nrow(df)))
    message(sprintf("  H (treat.recommend == 0): %d",
                    sum(df$treat.recommend == 0, na.rm = TRUE)))
    message(sprintf("  Hc (treat.recommend == 1): %d",
                    sum(df$treat.recommend == 1, na.rm = TRUE)))
    if (!is.null(fs_bc)) {
      message("  Bootstrap results provided: YES")
    }
  }

  # ===========================================================================
  # SECTION 2: SET UP PARAMETERS
  # ===========================================================================

  # Default colors
  if (is.null(colors)) {
    colors <- list(
      treat = "#E41A1C",
      control = "#377EB8"
    )
  }

  # Auto-calculate by.risk if not provided
  if (is.null(by.risk)) {
    by.risk <- round(max(df[[outcome.name]], na.rm = TRUE) / 12, 0)
    if (by.risk < 1) by.risk <- 1
    if (verbose) message(sprintf("  Auto-calculated by.risk: %d", by.risk))
  }

  # Handle est.scale for treatment reversal
  if (est.scale == "1/hr") {
    df$treat_plot <- 1 - df[[treat.name]]
    treat_labels <- rev(treat_labels)
    # Swap subgroup names for interpretation
    temp <- sg0_name
    sg0_name <- sg1_name
    sg1_name <- temp
    if (verbose) message("  est.scale = '1/hr': treatment labels reversed")
  } else {
    df$treat_plot <- df[[treat.name]]
  }

  # Map treatment values to labels for weightedsurv
  df$treat_label <- ifelse(df$treat_plot == 1, treat_labels[2], treat_labels[1])

  # ===========================================================================
  # SECTION 3: CREATE SUBGROUP DATA FRAMES
  # ===========================================================================

  # H subgroup (treat.recommend == 0): harm/questionable
  df_H <- df[df$treat.recommend == 0 & !is.na(df$treat.recommend), ]

  # Hc subgroup (treat.recommend == 1): complement/recommend
  df_Hc <- df[df$treat.recommend == 1 & !is.na(df$treat.recommend), ]

  if (nrow(df_H) < 5 || nrow(df_Hc) < 5) {
    stop(
      "Insufficient data in subgroups. H: ", nrow(df_H),
      ", Hc: ", nrow(df_Hc), ". Need at least 5 per group."
    )
  }

  if (verbose) {
    message(sprintf("  H subgroup: n = %d, events = %d",
                    nrow(df_H), sum(df_H[[event.name]])))
    message(sprintf("  Hc subgroup: n = %d, events = %d",
                    nrow(df_Hc), sum(df_Hc[[event.name]])))
  }

  # ===========================================================================
  # SECTION 4: COMPUTE HAZARD RATIO ESTIMATES
  # ===========================================================================

  # Compute observed HRs using Cox model
  hr_observed <- compute_observed_hr(
    df_H = df_H,
    df_Hc = df_Hc,
    outcome.name = outcome.name,
    event.name = event.name,
    treat.name = treat.name,
    verbose = verbose
  )

  # Extract bootstrap bias-corrected HRs if available
  hr_bc <- NULL
  if (!is.null(fs_bc) && show_bc_hr) {
    hr_bc <- extract_bootstrap_hr(
      fs_bc = fs_bc,
      est.scale = est.scale,
      verbose = verbose
    )
  }

  # ===========================================================================
  # SECTION 5: CREATE COUNTING PROCESS DATA FOR WEIGHTEDSURV
  # ===========================================================================

  # Initialize variables for weightedsurv data
  dfcount_H <- NULL
  dfcount_Hc <- NULL
  use_base_survival <- !use_weightedsurv

  if (use_weightedsurv) {
    # Prepare data for weightedsurv - needs treatment column name and arm labels
    # weightedsurv expects: treat.name = column with 0/1, arms = c(exp_label, con_label)
    # where exp_label corresponds to treat==1 and con_label corresponds to treat==0

    # Determine treatment column to use (handle est.scale)
    if (est.scale == "1/hr") {
      # For 1/hr scale, we already created treat_plot with reversed coding
      treat_col_for_ws <- "treat_plot"
    } else {
      treat_col_for_ws <- treat.name
    }

    # arm labels: first is experimental (treat==1), second is control (treat==0)
    ws_arms <- c(treat_labels[2], treat_labels[1])

    if (verbose) {
      message(sprintf("  weightedsurv params: treat_col='%s', arms=c('%s','%s'), by.risk=%d",
                      treat_col_for_ws, ws_arms[1], ws_arms[2], by.risk))
    }

    dfcount_H <- tryCatch({
      weightedsurv::df_counting(
        df_H,
        tte.name = outcome.name,
        event.name = event.name,
        treat.name = treat_col_for_ws,
        arms = ws_arms,
        by.risk = by.risk
      )
    }, error = function(e) {
      if (verbose) message("  H subgroup df_counting error: ", e$message)
      NULL
    })

    dfcount_Hc <- tryCatch({
      weightedsurv::df_counting(
        df_Hc,
        tte.name = outcome.name,
        event.name = event.name,
        treat.name = treat_col_for_ws,
        arms = ws_arms,
        by.risk = by.risk
      )
    }, error = function(e) {
      if (verbose) message("  Hc subgroup df_counting error: ", e$message)
      NULL
    })

    # If weightedsurv fails, fall back to base survival plotting
    if (is.null(dfcount_H) && is.null(dfcount_Hc)) {
      warning(
        "weightedsurv::df_counting() failed for both subgroups. ",
        "Falling back to base survival::survfit() plotting.\n",
        "For more details, run with verbose = TRUE."
      )
      use_base_survival <- TRUE
    }
  }

  # ===========================================================================
  # SECTION 6: SET UP PLOT LAYOUT
  # ===========================================================================

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  if (layout == "horizontal") {
    par(mfrow = c(1, 2), mar = c(5, 4, 4, 2) + 0.1)
  } else {
    par(mfrow = c(2, 1), mar = c(5, 4, 4, 2) + 0.1)
  }

  # ===========================================================================
  # SECTION 7: PLOT H SUBGROUP (QUESTIONABLE)
  # ===========================================================================

  # Track whether plot was successfully created
 plot_success_H <- FALSE

  tryCatch({
    if (use_base_survival) {
      # Fallback: use base survival::survfit plotting
      plot_success_H <- plot_base_km(
        data = df_H,
        outcome.name = outcome.name,
        event.name = event.name,
        treat.name = treat.name,
        title = sprintf("%s (n=%d)", sg0_name, nrow(df_H)),
        treat_labels = treat_labels,
        colors = colors,
        conf.int = conf.int,
        show.logrank = show.logrank,
        put.legend = put.legend.lr
      )
    } else if (!is.null(dfcount_H)) {
      # Build title with sample size
      h_title <- sprintf("%s (n=%d)", sg0_name, nrow(df_H))

      # Plot weighted KM with error handling
      weightedsurv::plot_weighted_km(
        dfcount_H,
        conf.int = conf.int,
        show.logrank = show.logrank,
        put.legend.lr = put.legend.lr,
        ymax = ymax,
        xmed.fraction = xmed.fraction,
        main = h_title,
        ...
      )
      plot_success_H <- TRUE
    } else {
      # No data available - create placeholder
      plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = paste("Could not create plot for", sg0_name))
      text(0, 0, "No data available", cex = 1.2, col = "gray50")
    }
  }, error = function(e) {
    if (verbose) message("  H subgroup plot error: ", e$message)
    # Create placeholder plot on error
    tryCatch({
      plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
           axes = FALSE, xlab = "", ylab = "",
           main = sprintf("%s - Error", sg0_name))
      text(0.5, 0.5, "Plot failed", cex = 1.2, col = "gray50")
    }, error = function(e2) NULL)
  })

  # Add HR annotation only if plot was successful
  if (show_hr && isTRUE(plot_success_H)) {
    tryCatch({
      add_hr_annotation(
        hr_obs = hr_observed$H,
        hr_bc = if (!is.null(hr_bc)) hr_bc$H else NULL,
        position = hr_position,
        show_bc = show_bc_hr && !is.null(hr_bc)
      )
    }, error = function(e) {
      if (verbose) message("  H subgroup annotation error: ", e$message)
    })
  }

  # ===========================================================================
  # SECTION 8: PLOT Hc SUBGROUP (RECOMMEND)
  # ===========================================================================

  # Track whether plot was successfully created
  plot_success_Hc <- FALSE

  tryCatch({
    if (use_base_survival) {
      # Fallback: use base survival::survfit plotting
      plot_success_Hc <- plot_base_km(
        data = df_Hc,
        outcome.name = outcome.name,
        event.name = event.name,
        treat.name = treat.name,
        title = sprintf("%s (n=%d)", sg1_name, nrow(df_Hc)),
        treat_labels = treat_labels,
        colors = colors,
        conf.int = conf.int,
        show.logrank = show.logrank,
        put.legend = put.legend.lr
      )
    } else if (!is.null(dfcount_Hc)) {
      # Build title with sample size
      hc_title <- sprintf("%s (n=%d)", sg1_name, nrow(df_Hc))

      # Plot weighted KM with error handling
      weightedsurv::plot_weighted_km(
        dfcount_Hc,
        conf.int = conf.int,
        show.logrank = show.logrank,
        put.legend.lr = put.legend.lr,
        ymax = ymax,
        xmed.fraction = xmed.fraction,
        main = hc_title,
        ...
      )
      plot_success_Hc <- TRUE
    } else {
      # No data available - create placeholder
      plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = paste("Could not create plot for", sg1_name))
      text(0, 0, "No data available", cex = 1.2, col = "gray50")
    }
  }, error = function(e) {
    if (verbose) message("  Hc subgroup plot error: ", e$message)
    # Create placeholder plot on error
    tryCatch({
      plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
           axes = FALSE, xlab = "", ylab = "",
           main = sprintf("%s - Error", sg1_name))
      text(0.5, 0.5, "Plot failed", cex = 1.2, col = "gray50")
    }, error = function(e2) NULL)
  })

  # Add HR annotation only if plot was successful
  if (show_hr && isTRUE(plot_success_Hc)) {
    tryCatch({
      add_hr_annotation(
        hr_obs = hr_observed$Hc,
        hr_bc = if (!is.null(hr_bc)) hr_bc$Hc else NULL,
        position = hr_position,
        show_bc = show_bc_hr && !is.null(hr_bc)
      )
    }, error = function(e) {
      if (verbose) message("  Hc subgroup annotation error: ", e$message)
    })
  }

  # ===========================================================================
  # SECTION 9: ADD OVERALL TITLE
  # ===========================================================================

  if (!is.null(title)) {
    mtext(title, side = 3, outer = TRUE, line = -1.5, cex = 1.2, font = 2)
  }

  # ===========================================================================
  # SECTION 10: CREATE SUMMARY STATISTICS
  # ===========================================================================

  summary_stats <- data.frame(
    Subgroup = c(sg0_name, sg1_name),
    n = c(nrow(df_H), nrow(df_Hc)),
    n_percent = round(100 * c(nrow(df_H), nrow(df_Hc)) / nrow(df), 1),
    events = c(sum(df_H[[event.name]]), sum(df_Hc[[event.name]])),
    event_rate = round(100 * c(
      sum(df_H[[event.name]]) / nrow(df_H),
      sum(df_Hc[[event.name]]) / nrow(df_Hc)
    ), 1),
    stringsAsFactors = FALSE
  )

  # ===========================================================================
  # SECTION 11: COMPILE AND RETURN OUTPUT
  # ===========================================================================

  result <- list(
    dfcount_H = dfcount_H,
    dfcount_Hc = dfcount_Hc,
    hr_observed = hr_observed,
    hr_bc = hr_bc,
    summary = summary_stats,
    sg_definition = fs.est$sg.harm,
    call = match.call()
  )

  class(result) <- c("fs_weighted_km", "list")

  invisible(result)
}


# =============================================================================
# HELPER FUNCTION: Compute Observed Hazard Ratios
# =============================================================================

#' Compute Observed Hazard Ratios for Subgroups
#'
#' Internal function to compute Cox model hazard ratio estimates with
#' confidence intervals for H and Hc subgroups.
#'
#' @param df_H Data frame for H subgroup
#' @param df_Hc Data frame for Hc subgroup
#' @param outcome.name Character. Outcome variable name
#' @param event.name Character. Event indicator name
#' @param treat.name Character. Treatment variable name
#' @param conf.level Numeric. Confidence level. Default: 0.95
#' @param verbose Logical. Print messages
#'
#' @return List with HR estimates for H and Hc
#'
#' @importFrom survival coxph Surv
#' @importFrom stats confint
#' @keywords internal
compute_observed_hr <- function(
    df_H,
    df_Hc,
    outcome.name,
    event.name,
    treat.name,
    conf.level = 0.95,
    verbose = FALSE
) {

  # Helper function for safe Cox model fitting
  safe_cox <- function(data, label) {
    if (nrow(data) < 10 || sum(data[[event.name]]) < 5) {
      if (verbose) message(sprintf("  %s: Insufficient data for Cox model", label))
      return(list(hr = NA, lower = NA, upper = NA, n = nrow(data),
                  events = sum(data[[event.name]]), pval = NA))
    }

    tryCatch({
      formula_str <- sprintf("survival::Surv(%s, %s) ~ %s",
                             outcome.name, event.name, treat.name)
      fit <- survival::coxph(as.formula(formula_str), data = data)
      ci <- confint(fit, level = conf.level)

      # Get p-value from summary
      fit_summary <- summary(fit)
      pval <- fit_summary$coefficients[1, "Pr(>|z|)"]

      result <- list(
        hr = exp(coef(fit)[1]),
        lower = exp(ci[1]),
        upper = exp(ci[2]),
        n = nrow(data),
        events = sum(data[[event.name]]),
        pval = pval
      )

      if (verbose) {
        message(sprintf("  %s: HR = %.3f (%.3f, %.3f), p = %.4f",
                        label, result$hr, result$lower, result$upper, result$pval))
      }

      result
    }, error = function(e) {
      if (verbose) message(sprintf("  %s: Cox model error: %s", label, e$message))
      list(hr = NA, lower = NA, upper = NA, n = nrow(data),
           events = sum(data[[event.name]]), pval = NA)
    })
  }

  list(
    H = safe_cox(df_H, "H (Questionable)"),
    Hc = safe_cox(df_Hc, "Hc (Recommend)")
  )
}


# =============================================================================
# HELPER FUNCTION: Base R Kaplan-Meier Plotting (Fallback)
# =============================================================================

#' Plot Kaplan-Meier Curves Using Base survival Package
#'
#' Internal fallback function when weightedsurv is unavailable or fails.
#' Creates standard Kaplan-Meier survival curves using survival::survfit.
#'
#' @param data Data frame with survival data
#' @param outcome.name Character. Outcome variable name
#' @param event.name Character. Event indicator name
#' @param treat.name Character. Treatment variable name
#' @param title Character. Plot title
#' @param treat_labels Character vector. Treatment arm labels
#' @param colors List. Colors for treatment arms
#' @param conf.int Logical. Show confidence intervals
#' @param show.logrank Logical. Show log-rank p-value
#' @param put.legend Character. Legend position
#'
#' @importFrom survival survfit Surv survdiff
#' @importFrom stats pchisq as.formula
#' @importFrom graphics plot lines legend mtext
#' @keywords internal
plot_base_km <- function(
    data,
    outcome.name,
    event.name,
    treat.name,
    title,
    treat_labels = c("control", "treat"),
    colors = list(control = "#377EB8", treat = "#E41A1C"),
    conf.int = TRUE,
    show.logrank = TRUE,
    put.legend = "topleft"
) {

  # Always create a plot panel (even if empty) so downstream code doesn't fail
  if (nrow(data) < 5) {
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", main = title)
    text(0, 0, "Insufficient data", cex = 1.2, col = "gray50")
    return(FALSE)  # Return FALSE to indicate plot failed
  }

  # Create survival objects directly (not via formula string)
  # This avoids issues with survival::Surv in formula strings
  Y <- data[[outcome.name]]
  E <- data[[event.name]]
  Treat <- data[[treat.name]]

  # Validate data
  if (any(is.na(Y)) || any(is.na(E)) || any(is.na(Treat))) {
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", main = title)
    text(0, 0, "Missing values in data", cex = 1.2, col = "gray50")
    return(FALSE)
  }

  # Fit survival curves using direct Surv object
  fit <- tryCatch({
    survival::survfit(survival::Surv(Y, E) ~ Treat, data = data)
  }, error = function(e) {
    message("survfit error: ", e$message)
    NULL
  })

  if (is.null(fit)) {
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", main = title)
    text(0, 0, "Could not fit survival curve", cex = 1.2, col = "gray50")
    return(FALSE)
  }

  # Define colors for each arm (order depends on factor levels)
  plot_colors <- c(colors$control, colors$treat)

  # Plot and return success status
  result <- tryCatch({
    plot(fit,
         col = plot_colors,
         lwd = 2,
         xlab = "Time",
         ylab = "Survival Probability",
         main = title,
         conf.int = conf.int,
         mark.time = TRUE)

    # Add legend
    legend(put.legend,
           legend = treat_labels,
           col = plot_colors,
           lwd = 2,
           bty = "n",
           cex = 0.9)

    # Add log-rank test
    if (show.logrank) {
      lr_test <- tryCatch({
        survival::survdiff(survival::Surv(Y, E) ~ Treat, data = data)
      }, error = function(e) NULL)

      if (!is.null(lr_test)) {
        pval <- 1 - pchisq(lr_test$chisq, df = length(lr_test$n) - 1)
        pval_text <- if (pval < 0.001) "p < 0.001" else sprintf("p = %.3f", pval)
        mtext(paste("Log-rank:", pval_text), side = 3, line = 0, adj = 1, cex = 0.8)
      }
    }

    TRUE  # Return TRUE to indicate success

  }, error = function(e) {
    message("Plotting error: ", e$message)
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", main = title)
    text(0, 0, "Plotting failed", cex = 1.2, col = "gray50")
    FALSE
  })

  return(result)
}


# =============================================================================
# HELPER FUNCTION: Extract Bootstrap Hazard Ratios
# =============================================================================

#' Extract Bootstrap Bias-Corrected Hazard Ratios
#'
#' Internal function to extract bias-corrected hazard ratio estimates from
#' \code{forestsearch_bootstrap_dofuture} results.
#'
#' @param fs_bc List. Bootstrap results from \code{forestsearch_bootstrap_dofuture}
#' @param est.scale Character. Effect scale ("hr" or "1/hr")
#' @param verbose Logical. Print messages
#'
#' @return List with bias-corrected HR estimates for H and Hc, or NULL if unavailable
#' @keywords internal
extract_bootstrap_hr <- function(fs_bc, est.scale = "hr", verbose = FALSE) {

  if (is.null(fs_bc)) return(NULL)

  # Try to extract from H_estimates and Hc_estimates (numeric values)
  H_bc <- NULL
  Hc_bc <- NULL

  # Extract H estimates
  if (!is.null(fs_bc$H_estimates)) {
    h_est <- fs_bc$H_estimates
    H_bc <- list(
      hr_obs = if (!is.null(h_est$H0)) h_est$H0 else NA,
      hr_obs_lower = if (!is.null(h_est$H0_lower)) h_est$H0_lower else NA,
      hr_obs_upper = if (!is.null(h_est$H0_upper)) h_est$H0_upper else NA,
      hr_bc = if (!is.null(h_est$H2)) h_est$H2 else NA,
      hr_bc_lower = if (!is.null(h_est$H2_lower)) h_est$H2_lower else NA,
      hr_bc_upper = if (!is.null(h_est$H2_upper)) h_est$H2_upper else NA
    )
    if (verbose) {
      message(sprintf("  H bootstrap: Obs = %.3f, BC = %.3f (%.3f, %.3f)",
                      H_bc$hr_obs, H_bc$hr_bc, H_bc$hr_bc_lower, H_bc$hr_bc_upper))
    }
  }

  # Extract Hc estimates
  if (!is.null(fs_bc$Hc_estimates)) {
    hc_est <- fs_bc$Hc_estimates
    Hc_bc <- list(
      hr_obs = if (!is.null(hc_est$H0)) hc_est$H0 else NA,
      hr_obs_lower = if (!is.null(hc_est$H0_lower)) hc_est$H0_lower else NA,
      hr_obs_upper = if (!is.null(hc_est$H0_upper)) hc_est$H0_upper else NA,
      hr_bc = if (!is.null(hc_est$H2)) hc_est$H2 else NA,
      hr_bc_lower = if (!is.null(hc_est$H2_lower)) hc_est$H2_lower else NA,
      hr_bc_upper = if (!is.null(hc_est$H2_upper)) hc_est$H2_upper else NA
    )
    if (verbose) {
      message(sprintf("  Hc bootstrap: Obs = %.3f, BC = %.3f (%.3f, %.3f)",
                      Hc_bc$hr_obs, Hc_bc$hr_bc, Hc_bc$hr_bc_lower, Hc_bc$hr_bc_upper))
    }
  }

  # Also try to extract formatted strings from SG_CIs
  if (!is.null(fs_bc$SG_CIs)) {
    if (is.null(H_bc)) {
      H_bc <- list(hr_bc_string = fs_bc$SG_CIs$H_bc,
                   hr_obs_string = fs_bc$SG_CIs$H_raw)
    } else {
      H_bc$hr_bc_string <- fs_bc$SG_CIs$H_bc
      H_bc$hr_obs_string <- fs_bc$SG_CIs$H_raw
    }

    if (is.null(Hc_bc)) {
      Hc_bc <- list(hr_bc_string = fs_bc$SG_CIs$Hc_bc,
                    hr_obs_string = fs_bc$SG_CIs$Hc_raw)
    } else {
      Hc_bc$hr_bc_string <- fs_bc$SG_CIs$Hc_bc
      Hc_bc$hr_obs_string <- fs_bc$SG_CIs$Hc_raw
    }
  }

  # Handle est.scale = "1/hr" (swap H and Hc)
  if (est.scale == "1/hr") {
    temp <- H_bc
    H_bc <- Hc_bc
    Hc_bc <- temp
  }

  if (is.null(H_bc) && is.null(Hc_bc)) {
    if (verbose) message("  No bootstrap HR estimates could be extracted")
    return(NULL)
  }

  list(H = H_bc, Hc = Hc_bc)
}


# =============================================================================
# HELPER FUNCTION: Add HR Annotation to Plot
# =============================================================================

#' Add Hazard Ratio Annotation to Survival Plot
#'
#' Internal function to add hazard ratio annotation to the current plot.
#'
#' @param hr_obs List. Observed HR estimates (hr, lower, upper)
#' @param hr_bc List. Bias-corrected HR estimates (hr_bc, hr_bc_lower, hr_bc_upper)
#' @param position Character. Position for annotation
#' @param show_bc Logical. Whether to show bias-corrected estimate
#'
#' @keywords internal
add_hr_annotation <- function(hr_obs, hr_bc = NULL, position = "topright",
                              show_bc = FALSE) {

  # Safety check: verify a plot exists before trying to annotate
  usr <- tryCatch({
    par("usr")
  }, error = function(e) {
    return(NULL)
  })

  # If no plot exists or par("usr") failed, return silently
  if (is.null(usr) || all(usr == c(0, 1, 0, 1))) {
    # Default usr values indicate no real plot exists
    return(invisible(NULL))
  }

  x_range <- usr[2] - usr[1]
  y_range <- usr[4] - usr[3]

  # Additional check: if ranges are zero or invalid, skip annotation

  if (x_range <= 0 || y_range <= 0 || !is.finite(x_range) || !is.finite(y_range)) {
    return(invisible(NULL))
  }

  # Set coordinates based on position
  if (position == "topright") {
    x_pos <- usr[2] - 0.02 * x_range
    y_pos <- usr[4] - 0.05 * y_range
    adj <- c(1, 1)  # right-aligned
  } else if (position == "topleft") {
    x_pos <- usr[1] + 0.02 * x_range
    y_pos <- usr[4] - 0.05 * y_range
    adj <- c(0, 1)  # left-aligned
  } else if (position == "bottomright") {
    x_pos <- usr[2] - 0.02 * x_range
    y_pos <- usr[3] + 0.15 * y_range
    adj <- c(1, 0)
  } else {  # bottomleft
    x_pos <- usr[1] + 0.02 * x_range
    y_pos <- usr[3] + 0.15 * y_range
    adj <- c(0, 0)
  }

  # Build text lines
  lines_text <- character(0)
  line_cols <- character(0)

  # Observed HR
  if (!is.null(hr_obs) && !is.na(hr_obs$hr)) {
    obs_text <- sprintf("HR: %.2f (%.2f, %.2f)",
                        hr_obs$hr, hr_obs$lower, hr_obs$upper)
    lines_text <- c(lines_text, obs_text)
    line_cols <- c(line_cols, "black")
  }

  # Bias-corrected HR (if requested and available)
  if (show_bc && !is.null(hr_bc)) {
    if (!is.null(hr_bc$hr_bc) && !is.na(hr_bc$hr_bc)) {
      bc_text <- sprintf("BC: %.2f (%.2f, %.2f)",
                         hr_bc$hr_bc, hr_bc$hr_bc_lower, hr_bc$hr_bc_upper)
      lines_text <- c(lines_text, bc_text)
      line_cols <- c(line_cols, "darkgreen")
    } else if (!is.null(hr_bc$hr_bc_string)) {
      bc_text <- paste("BC:", hr_bc$hr_bc_string)
      lines_text <- c(lines_text, bc_text)
      line_cols <- c(line_cols, "darkgreen")
    }
  }

  # Add text to plot (with error handling)
  if (length(lines_text) > 0) {
    y_offset <- 0.04 * y_range

    tryCatch({
      for (i in seq_along(lines_text)) {
        y_current <- y_pos - (i - 1) * y_offset
        text(x_pos, y_current, lines_text[i],
             adj = adj, cex = 0.75, col = line_cols[i], font = 1)
      }
    }, error = function(e) {
      # Silently ignore annotation errors
      NULL
    })
  }

  invisible(NULL)
}


# =============================================================================
# PRINT METHOD
# =============================================================================

#' Print Method for fs_weighted_km Objects
#'
#' @param x An fs_weighted_km object
#' @param ... Additional arguments (unused)
#'
#' @export
print.fs_weighted_km <- function(x, ...) {
  cat("ForestSearch Weighted Kaplan-Meier Visualization\n")
  cat("=================================================\n\n")

  if (!is.null(x$sg_definition)) {
    cat("Subgroup Definition:\n")
    cat("  ", paste(x$sg_definition, collapse = " & "), "\n\n")
  }

  cat("Summary Statistics:\n")
  print(x$summary, row.names = FALSE)
  cat("\n")

  cat("Observed Hazard Ratios:\n")
  if (!is.null(x$hr_observed$H) && !is.na(x$hr_observed$H$hr)) {
    cat(sprintf("  H (Questionable): HR = %.3f (%.3f, %.3f)\n",
                x$hr_observed$H$hr, x$hr_observed$H$lower, x$hr_observed$H$upper))
  }
  if (!is.null(x$hr_observed$Hc) && !is.na(x$hr_observed$Hc$hr)) {
    cat(sprintf("  Hc (Recommend):   HR = %.3f (%.3f, %.3f)\n",
                x$hr_observed$Hc$hr, x$hr_observed$Hc$lower, x$hr_observed$Hc$upper))
  }

  if (!is.null(x$hr_bc)) {
    cat("\nBias-Corrected Hazard Ratios:\n")
    if (!is.null(x$hr_bc$H) && !is.null(x$hr_bc$H$hr_bc) && !is.na(x$hr_bc$H$hr_bc)) {
      cat(sprintf("  H (Questionable): HR = %.3f (%.3f, %.3f)\n",
                  x$hr_bc$H$hr_bc, x$hr_bc$H$hr_bc_lower, x$hr_bc$H$hr_bc_upper))
    } else if (!is.null(x$hr_bc$H$hr_bc_string)) {
      cat(sprintf("  H (Questionable): %s\n", x$hr_bc$H$hr_bc_string))
    }
    if (!is.null(x$hr_bc$Hc) && !is.null(x$hr_bc$Hc$hr_bc) && !is.na(x$hr_bc$Hc$hr_bc)) {
      cat(sprintf("  Hc (Recommend):   HR = %.3f (%.3f, %.3f)\n",
                  x$hr_bc$Hc$hr_bc, x$hr_bc$Hc$hr_bc_lower, x$hr_bc$Hc$hr_bc_upper))
    } else if (!is.null(x$hr_bc$Hc$hr_bc_string)) {
      cat(sprintf("  Hc (Recommend):   %s\n", x$hr_bc$Hc$hr_bc_string))
    }
  }

  invisible(x)
}
