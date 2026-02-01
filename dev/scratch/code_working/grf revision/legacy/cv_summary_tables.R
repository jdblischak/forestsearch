# =============================================================================
# cv_summary_tables.R - Publication-Ready Tables for CV Results
# =============================================================================
#
# gt-formatted tables for forestsearch_KfoldOut() output
#
# =============================================================================

#' Cross-Validation Summary Tables
#'
#' Creates publication-ready gt tables summarizing ForestSearch cross-validation
#' results, including ITT estimates, subgroup estimates from full analysis vs
#' K-fold analysis, and CV agreement metrics.
#'
#' @param cv_out List. Output from \code{forestsearch_KfoldOut(outall = TRUE)}.
#' @param est.scale Character. Estimate scale: "hr" or "1/hr" (default: "hr").
#' @param sg_definition Character. Subgroup definition text for footnote.
#'   If NULL, attempts to extract from cv_out.
#' @param caption Character. Table caption (default: auto-generated).
#' @param font_size Numeric. Base font size in pixels (default: 12).
#' @param show_cv_metrics Logical. Include CV agreement metrics table
#'   (default: TRUE).
#' @param show_sg_finding Logical. Include subgroup finding summary table
#'   (default: TRUE).
#' @param decimals Integer. Number of decimal places for estimates (default: 2).
#'
#' @return A list containing gt table objects:
#'   \describe{
#'     \item{estimates_table}{Combined ITT/Full-Analysis/K-fold estimates}
#'     \item{cv_metrics_table}{CV agreement metrics (sensitivity, PPV)}
#'     \item{sg_finding_table}{Subgroup finding summary across folds}
#'     \item{tab_all}{The underlying data frame}
#'   }
#'
#' @details
#' The estimates table shows:
#' \itemize{
#'   \item \strong{Overall (ITT)}: Intent-to-treat population
#'   \item \strong{Full-Analysis (FA)}: Using full-data subgroup assignment
#'   \item \strong{K-fold Analysis (KfA)}: Using CV-predicted subgroup assignment
#' }
#'
#' Row naming convention:
#' \itemize{
#'   \item FA_0 / KfA_0: Questionable/Harm subgroup (H)
#'   \item FA_1 / KfA_1: Recommend/Benefit subgroup (Hc)
#' }
#'
#' @examples
#' \dontrun{
#' # Run K-fold CV
#' cv_results <- forestsearch_Kfold(fs.est = fs_result, Kfolds = 10)
#'
#' # Get detailed output
#' cv_out <- forestsearch_KfoldOut(cv_results, outall = TRUE)
#'
#' # Create gt tables
#' tables <- cv_summary_tables(cv_out)
#'
#' # View estimates table
#' tables$estimates_table
#'
#' # View CV metrics
#' tables$cv_metrics_table
#' }
#'
#' @seealso
#' \code{\link{forestsearch_Kfold}} for running cross-validation
#' \code{\link{forestsearch_KfoldOut}} for summarizing CV results
#' \code{\link{sg_tables}} for similar formatting of main ForestSearch results
#'
#' @importFrom gt gt tab_header tab_spanner tab_footnote tab_source_note
#'   tab_options tab_style cells_body cells_column_labels md px fmt_number
#'   cols_label cols_align cell_fill cell_text cell_borders
#' @export
cv_summary_tables <- function(
    cv_out,
    est.scale = "hr",
    sg_definition = NULL,
    caption = NULL,
    font_size = 12,
    show_cv_metrics = TRUE,
    show_sg_finding = TRUE,
    decimals = 2
) {

  # ===========================================================================
 # Input Validation
  # ===========================================================================

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' required. Install with: install.packages('gt')")
  }

  if (is.null(cv_out$tab_all)) {
    stop("cv_out must be from forestsearch_KfoldOut(outall = TRUE)")
  }

  # ===========================================================================
  # Prepare Data
  # ===========================================================================

  tab_all <- cv_out$tab_all

  # Create display-friendly row labels
  row_labels <- rownames(tab_all)
  display_labels <- sapply(row_labels, function(x) {
    switch(x,
           "Overall" = "Overall (ITT)",
           "FA_0" = "Full-Analysis: Questionable (H)",
           "KfA_0" = "K-fold CV: Questionable (H)",
           "FA_1" = "Full-Analysis: Recommend (Hc)",
           "KfA_1" = "K-fold CV: Recommend (Hc)",
           x  # Default: return as-is
    )
  })

  # Add Analysis column
  tab_display <- data.frame(
    Analysis = display_labels,
    tab_all,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(tab_display) <- NULL

  # Determine row types for styling
  row_types <- sapply(row_labels, function(x) {
    if (x == "Overall") return("itt")
    if (grepl("^FA_", x)) return("full_analysis")
    if (grepl("^KfA_", x)) return("kfold")
    return("other")
  })

  # ===========================================================================
  # Create Estimates Table
  # ===========================================================================

  if (is.null(caption)) {
    caption <- "Cross-Validation Subgroup Estimates"
  }

  # Build gt table
  estimates_gt <- gt::gt(tab_display) |>
    gt::tab_header(
      title = gt::md(paste0("**", caption, "**")),
      subtitle = "Comparison of Full-Analysis vs K-fold CV Predictions"
    ) |>
    gt::cols_label(
      Analysis = "Population / Subgroup",
      Subgroup = "Definition",
      n = "N",
      n1 = gt::md("n<sub>E</sub>"),
      m1 = gt::md("Events<sub>E</sub>"),
      m0 = gt::md("Events<sub>C</sub>"),
      RMST = "RMST Diff",
      `Hazard ratio` = "HR (95% CI)"
    ) |>
    gt::cols_align(
      align = "left",
      columns = c("Analysis", "Subgroup")
    ) |>
    gt::cols_align(
      align = "center",
      columns = c("n", "n1", "m1", "m0", "RMST", "Hazard ratio")
    )

  # Style ITT row
  itt_rows <- which(row_types == "itt")
  if (length(itt_rows) > 0) {
    estimates_gt <- estimates_gt |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#fff9c4"),  # Light yellow
          gt::cell_text(weight = "bold")
        ),
        locations = gt::cells_body(rows = itt_rows)
      )
  }

  # Style Full-Analysis rows
  fa_rows <- which(row_types == "full_analysis")
  if (length(fa_rows) > 0) {
    estimates_gt <- estimates_gt |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#e3f2fd")  # Light blue
        ),
        locations = gt::cells_body(rows = fa_rows)
      )
  }

  # Style K-fold rows
  kf_rows <- which(row_types == "kfold")
  if (length(kf_rows) > 0) {
    estimates_gt <- estimates_gt |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#f3e5f5")  # Light purple
        ),
        locations = gt::cells_body(rows = kf_rows)
      )
  }

  # Add column header border
  estimates_gt <- estimates_gt |>
    gt::tab_style(
      style = gt::cell_borders(
        sides = "bottom",
        color = "black",
        weight = gt::px(2)
      ),
      locations = gt::cells_column_labels()
    )

  # Font sizes
  estimates_gt <- estimates_gt |>
    gt::tab_options(
      table.font.size = gt::px(font_size),
      heading.title.font.size = gt::px(font_size + 2),
      heading.subtitle.font.size = gt::px(font_size),
      column_labels.font.size = gt::px(font_size),
      footnotes.font.size = gt::px(font_size - 1),
      source_notes.font.size = gt::px(font_size - 1)
    )

  # Add footnotes
  estimates_gt <- estimates_gt |>
    gt::tab_footnote(
      footnote = gt::md("**Full-Analysis**: Subgroup assignment from full-data model"),
      locations = gt::cells_body(columns = "Analysis", rows = fa_rows[1])
    ) |>
    gt::tab_footnote(
      footnote = gt::md("**K-fold CV**: Subgroup assignment predicted from training folds"),
      locations = gt::cells_body(columns = "Analysis", rows = kf_rows[1])
    )

  # Add subgroup definition if available
  if (!is.null(sg_definition)) {
    estimates_gt <- estimates_gt |>
      gt::tab_source_note(
        source_note = gt::md(paste0("**Identified subgroup (H)**: ", sg_definition))
      )
  }

  # Add legend
  estimates_gt <- estimates_gt |>
    gt::tab_source_note(
      source_note = gt::md(
        paste0(
          "n<sub>E</sub> = sample size experimental arm; ",
          "Events<sub>E</sub>/Events<sub>C</sub> = events in experimental/control; ",
          "RMST = restricted mean survival time difference"
        )
      )
    )

  # ===========================================================================
  # Create CV Metrics Table
  # ===========================================================================

  cv_metrics_gt <- NULL

  if (show_cv_metrics && !is.null(cv_out$sens_metrics_original)) {
    sens <- cv_out$sens_metrics_original

    # Create metrics data frame
    metrics_df <- data.frame(
      Metric = c("Sensitivity (H)", "Sensitivity (Hc)",
                 "PPV (H)", "PPV (Hc)"),
      Description = c(
        "Proportion of FA H patients also predicted H in CV",
        "Proportion of FA Hc patients also predicted Hc in CV",
        "Proportion of CV H predictions that match FA H",
        "Proportion of CV Hc predictions that match FA Hc"
      ),
      Value = sprintf("%.1f%%", 100 * sens),
      stringsAsFactors = FALSE
    )

    cv_metrics_gt <- gt::gt(metrics_df) |>
      gt::tab_header(
        title = gt::md("**Cross-Validation Agreement Metrics**"),
        subtitle = "Agreement between Full-Analysis and K-fold CV Predictions"
      ) |>
      gt::cols_label(
        Metric = "Metric",
        Description = "Description",
        Value = "Value"
      ) |>
      gt::cols_align(align = "left", columns = c("Metric", "Description")) |>
      gt::cols_align(align = "center", columns = "Value") |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_body(columns = "Metric")
      ) |>
      gt::tab_style(
        style = gt::cell_borders(
          sides = "bottom",
          color = "black",
          weight = gt::px(2)
        ),
        locations = gt::cells_column_labels()
      ) |>
      gt::tab_options(
        table.font.size = gt::px(font_size),
        heading.title.font.size = gt::px(font_size + 2),
        heading.subtitle.font.size = gt::px(font_size),
        column_labels.font.size = gt::px(font_size)
      ) |>
      gt::tab_footnote(
        footnote = gt::md("**H** = Harm/Questionable subgroup; **Hc** = Complement/Recommend subgroup"),
        locations = gt::cells_column_labels(columns = "Metric")
      ) |>
      gt::tab_footnote(
        footnote = gt::md("**Sensitivity** = True positive rate; **PPV** = Positive predictive value"),
        locations = gt::cells_column_labels(columns = "Value")
      )

    # Highlight high agreement (>= 70%)
    high_agree_rows <- which(sens >= 0.70)
    if (length(high_agree_rows) > 0) {
      cv_metrics_gt <- cv_metrics_gt |>
        gt::tab_style(
          style = gt::cell_fill(color = "#c8e6c9"),  # Light green
          locations = gt::cells_body(columns = "Value", rows = high_agree_rows)
        )
    }
  }

  # ===========================================================================
  # Create Subgroup Finding Table
  # ===========================================================================

  sg_finding_gt <- NULL

  if (show_sg_finding && !is.null(cv_out$find_metrics)) {
    find <- cv_out$find_metrics

    # Create finding data frame
    finding_df <- data.frame(
      Metric = names(find),
      Description = c(
        "Any subgroup found in fold",
        "Exact match (both factors)",
        "At least 1 factor matches",
        "Covariate 1 found (any threshold)",
        "Covariate 2 found (any threshold)",
        "Both covariates found",
        "Covariate 1 exact match",
        "Covariate 2 exact match"
      )[seq_along(find)],
      Value = sprintf("%.1f%%", 100 * find),
      stringsAsFactors = FALSE
    )

    sg_finding_gt <- gt::gt(finding_df) |>
      gt::tab_header(
        title = gt::md("**Subgroup Finding Summary**"),
        subtitle = "How often the identified subgroup was recovered in CV folds"
      ) |>
      gt::cols_label(
        Metric = "Metric",
        Description = "Description",
        Value = "% of Folds"
      ) |>
      gt::cols_align(align = "left", columns = c("Metric", "Description")) |>
      gt::cols_align(align = "center", columns = "Value") |>
      gt::tab_style(
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_body(columns = "Metric")
      ) |>
      gt::tab_style(
        style = gt::cell_borders(
          sides = "bottom",
          color = "black",
          weight = gt::px(2)
        ),
        locations = gt::cells_column_labels()
      ) |>
      gt::tab_options(
        table.font.size = gt::px(font_size),
        heading.title.font.size = gt::px(font_size + 2),
        heading.subtitle.font.size = gt::px(font_size),
        column_labels.font.size = gt::px(font_size)
      )

    # Highlight key metrics
    any_row <- which(finding_df$Metric == "Any")
    exact_row <- which(finding_df$Metric == "Exact")

    if (length(any_row) > 0) {
      sg_finding_gt <- sg_finding_gt |>
        gt::tab_style(
          style = gt::cell_fill(color = "#fff9c4"),  # Light yellow
          locations = gt::cells_body(rows = any_row)
        )
    }

    if (length(exact_row) > 0) {
      sg_finding_gt <- sg_finding_gt |>
        gt::tab_style(
          style = gt::cell_fill(color = "#c8e6c9"),  # Light green
          locations = gt::cells_body(rows = exact_row)
        )
    }
  }

  # ===========================================================================
  # Return Results
  # ===========================================================================

  result <- list(
    estimates_table = estimates_gt,
    cv_metrics_table = cv_metrics_gt,
    sg_finding_table = sg_finding_gt,
    tab_all = tab_all
  )

  class(result) <- c("cv_summary_tables", "list")

  return(result)
}


#' Print Method for CV Summary Tables
#'
#' @param x A cv_summary_tables object
#' @param which Character. Which table to print: "estimates", "metrics",
#'   "finding", or "all" (default).
#' @param ... Additional arguments (ignored)
#' @export
print.cv_summary_tables <- function(x, which = "all", ...) {

  if (which == "all" || which == "estimates") {
    cat("\n=== Estimates Table ===\n")
    print(x$estimates_table)
  }

  if ((which == "all" || which == "metrics") && !is.null(x$cv_metrics_table)) {
    cat("\n=== CV Agreement Metrics ===\n")
    print(x$cv_metrics_table)
  }

  if ((which == "all" || which == "finding") && !is.null(x$sg_finding_table)) {
    cat("\n=== Subgroup Finding Summary ===\n")
    print(x$sg_finding_table)
  }

  invisible(x)
}


#' Format CV Metrics as Inline Text
#'
#' Creates formatted text summarizing CV metrics for use in reports.
#'
#' @param cv_out List. Output from forestsearch_KfoldOut().
#' @param est.scale Character. "hr" or "1/hr".
#'
#' @return Character string with formatted CV summary.
#'
#' @examples
#' \dontrun{
#' cv_out <- forestsearch_KfoldOut(cv_results, outall = TRUE)
#' cat(format_cv_metrics_text(cv_out))
#' }
#'
#' @export
format_cv_metrics_text <- function(cv_out, est.scale = "hr") {

  if (is.null(cv_out$find_metrics) || is.null(cv_out$sens_metrics_original)) {
    return("CV metrics not available")
  }

  find <- cv_out$find_metrics
  sens <- cv_out$sens_metrics_original

  # Format based on est.scale
  if (est.scale == "hr") {
    sg_label_H <- "Questionable (H)"
    sg_label_Hc <- "Recommend (Hc)"
  } else {
    sg_label_H <- "Recommend (H)"
    sg_label_Hc <- "Questionable (Hc)"
  }

  text <- sprintf(
    paste0(
      "Cross-validation found the subgroup in %.0f%% of folds ",
      "(exact match: %.0f%%). ",
      "Agreement between full-analysis and CV predictions: ",
      "%s sensitivity = %.0f%%, %s sensitivity = %.0f%%."
    ),
    100 * find["Any"],
    100 * find["Exact"],
    sg_label_H, 100 * sens["sens_H"],
    sg_label_Hc, 100 * sens["sens_Hc"]
  )

  return(text)
}


#' Create Combined CV Summary (Estimates + Metrics)
#'
#' Creates a single gt table combining estimates with CV metrics summary.
#' Useful for compact reporting.
#'
#' @param cv_out List. Output from forestsearch_KfoldOut(outall = TRUE).
#' @param sg_definition Character. Subgroup definition for footnote.
#' @param font_size Numeric. Base font size (default: 11).
#'
#' @return A gt table object.
#'
#' @examples
#' \dontrun{
#' cv_out <- forestsearch_KfoldOut(cv_results, outall = TRUE)
#' cv_summary_combined(cv_out, sg_definition = "BM > 1 & Size > 20")
#' }
#'
#' @importFrom gt gt tab_header tab_spanner tab_source_note md px
#' @export
cv_summary_combined <- function(
    cv_out,
    sg_definition = NULL,
    font_size = 11
) {

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' required")
  }

  # Get estimates table
  tables <- cv_summary_tables(
    cv_out,
    font_size = font_size,
    show_cv_metrics = FALSE,
    show_sg_finding = FALSE
  )

  estimates_gt <- tables$estimates_table

  # Add CV metrics as source note
  if (!is.null(cv_out$find_metrics) && !is.null(cv_out$sens_metrics_original)) {
    find <- cv_out$find_metrics
    sens <- cv_out$sens_metrics_original

    cv_text <- sprintf(
      "**CV Summary:** Found in %.0f%% of folds (exact: %.0f%%) | Agreement: Sens(H)=%.0f%%, Sens(Hc)=%.0f%%, PPV(H)=%.0f%%, PPV(Hc)=%.0f%%",
      100 * find["Any"],
      100 * find["Exact"],
      100 * sens["sens_H"],
      100 * sens["sens_Hc"],
      100 * sens["ppv_H"],
      100 * sens["ppv_Hc"]
    )

    estimates_gt <- estimates_gt |>
      gt::tab_source_note(
        source_note = gt::md(cv_text)
      )
  }

  # Add subgroup definition
  if (!is.null(sg_definition)) {
    estimates_gt <- estimates_gt |>
      gt::tab_source_note(
        source_note = gt::md(paste0("**Identified subgroup (H):** ", sg_definition))
      )
  }

  return(estimates_gt)
}
