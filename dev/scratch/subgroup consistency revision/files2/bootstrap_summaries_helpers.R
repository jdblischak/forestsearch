#' Format Bootstrap Results Table with gt
#'
#' Creates a publication-ready table from ForestSearch bootstrap results,
#' with bias-corrected confidence intervals, informative formatting, and
#' optional subgroup definition footnote.
#'
#' @param FSsg_tab Data frame or matrix from forestsearch_bootstrap_dofuture()$FSsg_tab
#' @param nb_boots Integer. Number of bootstrap iterations performed
#' @param est.scale Character. "hr" or "1/hr" for effect scale
#' @param boot_success_rate Numeric. Proportion of bootstraps that found subgroups
#' @param sg_definition Character. Subgroup definition string to display as footnote
#'   (e.g., "{age>=50} & {nodes>=3}"). If NULL, no subgroup footnote is added.
#' @param title Character. Custom title (optional)
#' @param subtitle Character. Custom subtitle (optional)
#'
#' @return A gt table object
#'
#' @details
#' The table includes:
#' \itemize{
#'   \item Sample size columns (N, N_T, Events)
#'   \item Survival columns (Median T, Median C, RMST difference)
#'   \item Treatment effect columns (unadjusted HR, bias-corrected HR)
#'   \item Footnotes explaining methodology
#'   \item Optional subgroup definition footnote attached to the "H" row
#' }
#'
#' @seealso
#' \code{\link{summarize_bootstrap_results}} which calls this function
#' \code{\link{sg_tables}} for analogous table formatting in main analysis
#'
#' @importFrom gt gt tab_header tab_spanner tab_footnote tab_source_note md
#'   cols_label tab_style cell_fill cell_text cells_body cells_column_labels
#' @importFrom dplyr all_of
#' @export
format_bootstrap_table <- function(FSsg_tab,
                                   nb_boots,
                                   est.scale = "hr",
                                   boot_success_rate = NULL,
                                   sg_definition = NULL,
                                   title = NULL,
                                   subtitle = NULL) {

 if (!requireNamespace("gt", quietly = TRUE)) {
    stop(
      "Package 'gt' is required for table formatting. ",
      "Install with: install.packages('gt')"
    )
  }

  # ---------------------------------------------------------------------------
  # INPUT VALIDATION AND CONVERSION
  # ---------------------------------------------------------------------------

  # Convert matrix to data frame if needed
  if (is.matrix(FSsg_tab)) {
    FSsg_tab <- as.data.frame(FSsg_tab, stringsAsFactors = FALSE)
  }

  if (!is.data.frame(FSsg_tab)) {
    FSsg_tab <- as.data.frame(FSsg_tab, stringsAsFactors = FALSE)
  }

  # ---------------------------------------------------------------------------
  # DEFAULT TITLE AND SUBTITLE
  # ---------------------------------------------------------------------------

  if (is.null(title)) {
    title <- "Treatment Effect by Subgroup"
  }

  if (is.null(subtitle)) {
    effect_label <- ifelse(
      est.scale == "hr",
      "Hazard Ratio",
      "Inverse Hazard Ratio (1/HR)"
    )
    subtitle <- sprintf(
      "Bootstrap bias-corrected estimates (%d iterations)",
      nb_boots
    )
  }

  # ---------------------------------------------------------------------------
  # COLUMN LABELS CONFIGURATION
  # ---------------------------------------------------------------------------

  col_names <- colnames(FSsg_tab)

  # Create base labels list
  labels_list <- list(Subgroup = "Subgroup")

  # Add labels for columns that exist
  if ("n" %in% col_names) labels_list$n <- "N"
  if ("n1" %in% col_names) labels_list$n1 <- gt::md("N<sub>T</sub>")
  if ("events" %in% col_names) labels_list$events <- "Events"
  if ("m1" %in% col_names) labels_list$m1 <- gt::md("Med<sub>T</sub>")
  if ("m0" %in% col_names) labels_list$m0 <- gt::md("Med<sub>C</sub>")
  if ("RMST" %in% col_names) labels_list$RMST <- gt::md("RMST<sub>d</sub>")

  # Handle HR column (might be "HR (95% CI)" or similar)
  hr_col <- grep("HR.*CI", col_names, value = TRUE)[1]
  if (!is.na(hr_col) && length(hr_col) > 0) {
    labels_list[[hr_col]] <- gt::md("HR<br/>(95% CI)<sup>\u2020</sup>")
  }

  # Handle adjusted HR column (might be "HR*" or similar)
  hr_adj_col <- grep("HR\\*", col_names, value = TRUE)[1]
  if (!is.na(hr_adj_col) && length(hr_adj_col) > 0) {
    labels_list[[hr_adj_col]] <- gt::md("HR<sup>\u2021</sup><br/>(95% CI)")
  }

  # ---------------------------------------------------------------------------
  # CREATE BASE GT TABLE
  # ---------------------------------------------------------------------------

  tbl <- FSsg_tab |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md(paste0("**", title, "**")),
      subtitle = subtitle
    ) |>
    gt::cols_label(.list = labels_list)

  # ---------------------------------------------------------------------------
  # ADD COLUMN SPANNERS
  # ---------------------------------------------------------------------------

  sample_size_cols <- intersect(c("n", "n1", "events"), col_names)
  if (length(sample_size_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(
        label = "Sample Size",
        columns = sample_size_cols
      )
  }

  survival_cols <- intersect(c("m1", "m0", "RMST"), col_names)
  if (length(survival_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(
        label = "Survival",
        columns = survival_cols
      )
  }

  hr_cols <- grep("HR", col_names, value = TRUE)
  if (length(hr_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(
        label = "Treatment Effect",
        columns = gt::starts_with("HR")
      )
  }

  # ---------------------------------------------------------------------------
  # ADD METHODOLOGY FOOTNOTES
  # ---------------------------------------------------------------------------

  if (!is.na(hr_col) && length(hr_col) > 0) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md(
          "**Unadjusted HR**: Standard Cox regression hazard ratio with ",
          "robust standard errors"
        ),
        locations = gt::cells_column_labels(columns = dplyr::all_of(hr_col))
      )
  }

  if (!is.na(hr_adj_col) && length(hr_adj_col) > 0) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md(sprintf(
          "**Bias-corrected HR**: Bootstrap-adjusted estimate using ",
          "infinitesimal jackknife method (%d iterations). Corrects for ",
          "optimism in subgroup selection.",
          nb_boots
        )),
        locations = gt::cells_column_labels(columns = dplyr::all_of(hr_adj_col))
      )
  }

  # ---------------------------------------------------------------------------
  # ADD SUBGROUP DEFINITION FOOTNOTE (analogous to sg_tables tab_estimates)
  # ---------------------------------------------------------------------------

  if (!is.null(sg_definition) && nzchar(sg_definition)) {
    # Find the row containing "H" or "Questionable" subgroup
    # This is analogous to sg_tables() which targets the questionable subgroup row
    subgroup_col <- "Subgroup"

    if (subgroup_col %in% col_names) {
      # Identify rows to attach footnote to (H subgroup row)
      h_row_pattern <- "^H$|Questionable|H \\(|Harm"

      tbl <- tbl |>
        gt::tab_footnote(
          footnote = gt::md(paste0("**Identified subgroup (H)**: ", sg_definition)),
          locations = gt::cells_body(
            columns = subgroup_col,
            rows = grepl(h_row_pattern, FSsg_tab[[subgroup_col]])
          )
        )
    }
  }

  # ---------------------------------------------------------------------------
  # ADD SOURCE NOTE
  # ---------------------------------------------------------------------------

  source_note_text <- paste0(
    "*Note*: Med = Median survival time (months). ",
    "RMST<sub>d</sub> = Restricted mean survival time difference."
  )

  if (!is.null(boot_success_rate)) {
    source_note_text <- paste0(
      source_note_text,
      sprintf(
        " Subgroup identified in %.1f%% of bootstrap samples.",
        boot_success_rate * 100
      )
    )
  }

  tbl <- tbl |>
    gt::tab_source_note(source_note = gt::md(source_note_text))

  # ---------------------------------------------------------------------------
  # ADD INTERPRETIVE CAPTION
  # ---------------------------------------------------------------------------

  effect_name <- ifelse(est.scale == "hr", "Hazard Ratios", "Inverse Hazard Ratios")

  caption <- sprintf(
    paste0(
      "%s less than 1.0 indicate benefit from treatment. ",
      "Bias correction typically shifts estimates toward the null, ",
      "reflecting the optimism inherent in data-driven subgroup identification. ",
      "Bootstrap analysis based on %d iterations (%.1f%% successful)."
    ),
    effect_name,
    nb_boots,
    if (!is.null(boot_success_rate)) boot_success_rate * 100 else NA
  )

  tbl <- tbl |>
    gt::tab_source_note(source_note = gt::md(paste0("*", caption, "*")))

  # ---------------------------------------------------------------------------
  # FINAL STYLING
  # ---------------------------------------------------------------------------

  tbl <- tbl |>
    gt::tab_options(
      table.border.top.style = "solid",
      table.border.top.width = gt::px(3),
      table.border.top.color = "#333333",
      heading.border.bottom.style = "solid",
      heading.border.bottom.width = gt::px(2),
      heading.border.bottom.color = "#333333",
      column_labels.border.bottom.style = "solid",
      column_labels.border.bottom.width = gt::px(2),
      column_labels.border.bottom.color = "#333333",
      table.font.size = gt::px(13),
      heading.align = "left"
    )

  return(tbl)
}
