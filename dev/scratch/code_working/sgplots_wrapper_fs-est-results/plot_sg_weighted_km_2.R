# =============================================================================
# plot_sg_weighted_km.R - Weighted Kaplan-Meier Plots for ForestSearch Subgroups
# =============================================================================
#
# Matches plot_subgroup() pattern from sg_consistency_out() with flexible
# column name mapping and readable subgroup definitions.
#
# =============================================================================

#' Plot Weighted Kaplan-Meier Curves for ForestSearch Subgroups
#'
#' Creates weighted Kaplan-Meier survival curves for the identified subgroups
#' (H and Hc) using the weightedsurv package, matching the pattern used in
#' \code{sg_consistency_out()}.
#'
#' @param fs.est A forestsearch object containing \code{df.est} with
#'   \code{treat.recommend} column, or a data frame directly.
#' @param fs_bc Optional. Bootstrap results from \code{forestsearch_bootstrap_dofuture()}
#'   containing bias-corrected HR estimates. If provided, bias-corrected HRs
#'   will be annotated on the plots.
#' @param outcome.name Character. Name of time-to-event column.
#'   Default: \code{"Y"}
#' @param event.name Character. Name of event indicator column.
#'   Default: \code{"Event"}
#' @param treat.name Character. Name of treatment column.
#'   Default: \code{"Treat"}
#' @param by.risk Numeric. Risk interval for plotting. Default: \code{NULL}
#'   (auto-calculated as max(outcome)/12)
#' @param sg0_name Character. Label for H subgroup (treat.recommend == 0).
#'   Default: \code{NULL} (auto-extracted from forestsearch object as
#'   "H: {definition}" or "Questionable (H)" if not available)
#' @param sg1_name Character. Label for Hc subgroup (treat.recommend == 1).
#'   Default: \code{NULL} (auto-generated as "Hc: NOT {definition}" or
#'   "Recommend (Hc)" if not available)
#' @param conf.int Logical. Show confidence intervals. Default: \code{TRUE}
#' @param show.logrank Logical. Show log-rank test. Default: \code{TRUE}
#' @param put.legend.lr Character. Legend position. Default: "topleft"
#' @param ymax Numeric. Max y-axis value. Default: 1.05
#' @param xmed.fraction Numeric. Fraction for median lines. Default: 0.65
#' @param hr_bc_position Character. Position for bias-corrected HR annotation.
#'   One of "bottomright", "bottomleft", "topright", "topleft".
#'   Default: "bottomright"
#' @param title Character. Overall plot title. Default: \code{NULL}
#' @param verbose Logical. Print diagnostic messages. Default: \code{FALSE}
#'
#' @details
#' This function uses the exact same calling pattern as \code{plot_subgroup()}
#' in the ForestSearch package. Column names are mapped internally to the
#' standard names (Y, Event, Treat) expected by weightedsurv.
#'
#' Subgroup definitions are automatically extracted from the forestsearch
#' object if available:
#' \itemize{
#'   \item \code{fs$grp.consistency$out_sg$sg.harm_label} - Human-readable labels
#'   \item \code{fs$sg.harm} - Technical factor names (fallback)
#' }
#'
#' If \code{fs_bc} is provided, bias-corrected hazard ratios from bootstrap
#' analysis are displayed in the specified corner of each plot panel.
#'
#' @return Invisibly returns a list with subgroup data frames and counting data
#'
#' @examples
#' \dontrun{
#' # After running forestsearch - auto-extracts subgroup definition
#' plot_sg_weighted_km(fs.est = fs)
#'
#' # With bootstrap bias-corrected estimates
#' plot_sg_weighted_km(fs.est = fs, fs_bc = fs_bootstrap)
#'
#' # With custom column names
#' plot_sg_weighted_km(
#'   fs.est = fs,
#'   outcome.name = "time_months",
#'   event.name = "status",
#'   treat.name = "hormon"
#' )
#' }
#'
#' @importFrom graphics par title mtext text
#' @export
plot_sg_weighted_km <- function(
    fs.est,
    fs_bc = NULL,
    outcome.name = "Y",
    event.name = "Event",
    treat.name = "Treat",
    by.risk = NULL,
    sg0_name = NULL,
    sg1_name = NULL,
    conf.int = TRUE,
    show.logrank = TRUE,
    put.legend.lr = "topleft",
    ymax = 1.05,
    xmed.fraction = 0.65,
    hr_bc_position = "bottomright",
    title = NULL,
    verbose = FALSE
) {

  # ===========================================================================
  # SECTION 1: EXTRACT SUBGROUP DEFINITION
  # ===========================================================================

  # Helper function to add bias-corrected HR annotation
  add_bc_annotation <- function(bc_estimates, position = "bottomright") {
    if (is.null(bc_estimates)) return(invisible(NULL))

    # Extract bias-corrected values (H2 = bias-corrected)
    hr_bc <- bc_estimates$H2
    hr_lower <- bc_estimates$H2_lower
    hr_upper <- bc_estimates$H2_upper

    if (is.na(hr_bc) || is.null(hr_bc)) return(invisible(NULL))

    # Format the annotation text
    hr_text <- sprintf("HR(bc): %.2f (%.2f, %.2f)", hr_bc, hr_lower, hr_upper)

    # Get plot coordinates
    usr <- par("usr")
    x_range <- usr[2] - usr[1]
    y_range <- usr[4] - usr[3]

    # Calculate position based on corner
    if (position == "bottomright") {
      x_pos <- usr[2] - 0.02 * x_range
      y_pos <- usr[3] + 0.08 * y_range
      adj <- c(1, 0)
    } else if (position == "bottomleft") {
      x_pos <- usr[1] + 0.02 * x_range
      y_pos <- usr[3] + 0.08 * y_range
      adj <- c(0, 0)
    } else if (position == "topright") {
      x_pos <- usr[2] - 0.02 * x_range
      y_pos <- usr[4] - 0.08 * y_range
      adj <- c(1, 1)
    } else {
      # topleft
      x_pos <- usr[1] + 0.02 * x_range
      y_pos <- usr[4] - 0.08 * y_range
      adj <- c(0, 1)
    }

    # Add text annotation
    text(x_pos, y_pos, hr_text, adj = adj, cex = 0.85, font = 2, col = "darkblue")

    invisible(NULL)
  }

  # ===========================================================================
  # SECTION 1b: EXTRACT BIAS-CORRECTED ESTIMATES
  # ===========================================================================

  H_bc_estimates <- NULL
  Hc_bc_estimates <- NULL

  if (!is.null(fs_bc)) {
    if (!is.null(fs_bc$H_estimates)) {
      H_bc_estimates <- fs_bc$H_estimates
      if (verbose) {
        cat("Extracted H bias-corrected HR:", 
            sprintf("%.3f (%.3f, %.3f)\n", 
                    H_bc_estimates$H2, H_bc_estimates$H2_lower, H_bc_estimates$H2_upper))
      }
    }
    if (!is.null(fs_bc$Hc_estimates)) {
      Hc_bc_estimates <- fs_bc$Hc_estimates
      if (verbose) {
        cat("Extracted Hc bias-corrected HR:", 
            sprintf("%.3f (%.3f, %.3f)\n", 
                    Hc_bc_estimates$H2, Hc_bc_estimates$H2_lower, Hc_bc_estimates$H2_upper))
      }
    }
  }

  # ===========================================================================
  # SECTION 2: EXTRACT SUBGROUP DEFINITION FROM FORESTSEARCH
  # ===========================================================================

  sg_definition <- NULL
  sg_definition_label <- NULL

  # Try to extract human-readable subgroup definition from forestsearch object
 if (inherits(fs.est, "forestsearch") || is.list(fs.est)) {

    # Priority 1: grp.consistency$out_sg$sg.harm_label (human-readable)
    if (!is.null(fs.est$grp.consistency$out_sg$sg.harm_label)) {
      sg_labels <- fs.est$grp.consistency$out_sg$sg.harm_label
      sg_labels <- sg_labels[!is.na(sg_labels) & sg_labels != ""]
      if (length(sg_labels) > 0) {
        sg_definition_label <- paste(sg_labels, collapse = " & ")
        if (verbose) cat("Extracted sg.harm_label:", sg_definition_label, "\n")
      }
    }

    # Priority 2: sg.harm (may be technical codes or labels)
    if (is.null(sg_definition_label) && !is.null(fs.est$sg.harm)) {
      sg_harm <- fs.est$sg.harm
      if (is.character(sg_harm)) {
        sg_harm <- sg_harm[!is.na(sg_harm) & sg_harm != ""]
        if (length(sg_harm) > 0) {
          sg_definition_label <- paste(sg_harm, collapse = " & ")
          if (verbose) cat("Extracted sg.harm:", sg_definition_label, "\n")
        }
      }
    }

    # Priority 3: Check for sg.harm as list with label element
    if (is.null(sg_definition_label) && is.list(fs.est$sg.harm)) {
      if (!is.null(fs.est$sg.harm$sg.harm_label)) {
        sg_labels <- fs.est$sg.harm$sg.harm_label
        sg_labels <- sg_labels[!is.na(sg_labels) & sg_labels != ""]
        if (length(sg_labels) > 0) {
          sg_definition_label <- paste(sg_labels, collapse = " & ")
        }
      }
    }
  }

  # ===========================================================================
  # SECTION 2: SET SUBGROUP NAMES
  # ===========================================================================

  # Use extracted definition or fallback to defaults
  if (is.null(sg0_name)) {
    if (!is.null(sg_definition_label)) {
      sg0_name <- sg_definition_label
    } else {
      sg0_name <- "Questionable (H)"
    }
  }

  if (is.null(sg1_name)) {
    if (!is.null(sg_definition_label)) {
      sg1_name <- paste0("NOT ", sg_definition_label)
    } else {
      sg1_name <- "Recommend (Hc)"
    }
  }

  # ===========================================================================
  # SECTION 3: EXTRACT AND VALIDATE DATA
  # ===========================================================================

  # Extract data frame
  if (inherits(fs.est, "forestsearch") && !is.null(fs.est$df.est)) {
    df <- fs.est$df.est
  } else if (is.data.frame(fs.est)) {
    df <- fs.est
  } else if (is.list(fs.est) && !is.null(fs.est$df.est)) {
    df <- fs.est$df.est
  } else {
    stop("Cannot extract data frame from fs.est")
  }

  # Check required columns
  required_cols <- c(outcome.name, event.name, treat.name, "treat.recommend")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # ===========================================================================
  # SECTION 4: MAP COLUMN NAMES TO STANDARD NAMES
  # ===========================================================================

  # Create working copy with standard names (Y, Event, Treat) for weightedsurv
  df_work <- df
  df_work$Y <- df[[outcome.name]]
  df_work$Event <- df[[event.name]]
  df_work$Treat <- df[[treat.name]]

  # ===========================================================================
  # SECTION 5: CREATE SUBGROUP DATA FRAMES
  # ===========================================================================

  # Split by treat.recommend (same logic as sg_consistency_out)
  df.sub <- df_work[df_work$treat.recommend == 0, ]   # H subgroup
  df.subC <- df_work[df_work$treat.recommend == 1, ]  # Hc subgroup

  n_H <- nrow(df.sub)
  n_Hc <- nrow(df.subC)

  if (verbose) {
    cat("H subgroup (treat.recommend == 0):", n_H, "\n")
    cat("Hc subgroup (treat.recommend == 1):", n_Hc, "\n")
  }

  if (n_H == 0 || n_Hc == 0) {
    stop("One or both subgroups have no observations")
  }

  # ===========================================================================
  # SECTION 6: CALCULATE by.risk IF NOT PROVIDED
  # ===========================================================================

  if (is.null(by.risk)) {
    by.risk <- round(max(df_work$Y, na.rm = TRUE) / 12, 0)
    if (by.risk < 1) by.risk <- 1
    if (verbose) cat("Auto-calculated by.risk:", by.risk, "\n")
  }

  # ===========================================================================
  # SECTION 7: CHECK FOR WEIGHTEDSURV AND PLOT
  # ===========================================================================

  if (!requireNamespace("weightedsurv", quietly = TRUE)) {
    stop("Package 'weightedsurv' is required. Install with: ",
         "devtools::install_github('larry-leon/weightedsurv')")
  }

  # Use exact same variable names as plot_subgroup()
  tte.name <- "Y"
  event.name.ws <- "Event"
  treat.name.ws <- "Treat"
  con.lab <- "control"
  exp.lab <- "treat"

  # Create counting process data - exact same call as plot_subgroup()
  dfcount <- weightedsurv::df_counting(
    df.sub,
    tte.name = tte.name,
    event.name = event.name.ws,
    treat.name = treat.name.ws,
    arms = c(exp.lab, con.lab),
    by.risk = by.risk
  )

  dfcountC <- weightedsurv::df_counting(
    df.subC,
    tte.name = tte.name,
    event.name = event.name.ws,
    treat.name = treat.name.ws,
    arms = c(exp.lab, con.lab),
    by.risk = by.risk
  )

  # Set up plot layout
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 2))

  # Plot H subgroup - exact same call as plot_subgroup()
  weightedsurv::plot_weighted_km(
    dfcount,
    conf.int = conf.int,
    show.logrank = show.logrank,
    put.legend.lr = put.legend.lr,
    ymax = ymax,
    xmed.fraction = xmed.fraction
  )
  title(main = sprintf("%s (n=%d)", sg0_name, n_H))

  # Add bias-corrected HR annotation for H subgroup
  if (!is.null(H_bc_estimates)) {
    add_bc_annotation(H_bc_estimates, position = hr_bc_position)
  }

  # Plot Hc subgroup - exact same call as plot_subgroup()
  weightedsurv::plot_weighted_km(
    dfcountC,
    conf.int = conf.int,
    show.logrank = show.logrank,
    put.legend.lr = put.legend.lr,
    ymax = ymax,
    xmed.fraction = xmed.fraction
  )
  title(main = sprintf("%s (n=%d)", sg1_name, n_Hc))

  # Add bias-corrected HR annotation for Hc subgroup
  if (!is.null(Hc_bc_estimates)) {
    add_bc_annotation(Hc_bc_estimates, position = hr_bc_position)
  }

  # Add overall title if provided
  if (!is.null(title)) {
    mtext(title, side = 3, outer = TRUE, line = -1.5, cex = 1.2, font = 2)
  }

  # Print summary
  if (verbose) {
    cat("*** Subgroups plotted\n")
    if (!is.null(sg_definition_label)) {
      cat("Subgroup definition:", sg_definition_label, "\n")
    }
    cat("H:", sg0_name, "- n =", n_H, ", events =", sum(df.sub$Event), "\n")
    cat("Hc:", sg1_name, "- n =", n_Hc, ", events =", sum(df.subC$Event), "\n")
  }

  invisible(list(
    df.sub = df.sub,
    df.subC = df.subC,
    dfcount = dfcount,
    dfcountC = dfcountC,
    by.risk = by.risk,
    sg_definition = sg_definition_label,
    sg0_name = sg0_name,
    sg1_name = sg1_name,
    H_bc_estimates = H_bc_estimates,
    Hc_bc_estimates = Hc_bc_estimates
  ))
}
