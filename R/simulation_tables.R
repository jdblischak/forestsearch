# =============================================================================
# R/simulation_tables.R
# Publication-quality tables for simulation operating characteristics
# =============================================================================

# NOTE: Global variable declarations for data.table columns (Group, Scenario,
# Subgroup) are consolidated in R/globals.R per package convention.

# =============================================================================
# Classification / Identification Rate Table
# =============================================================================

#' Build Classification Rate Table from Simulation Results
#'
#' Constructs a publication-quality \code{gt} table summarizing subgroup
#' identification and classification rates across one or more data generation
#' scenarios and analysis methods. The layout mirrors Table 4 of
#' Leon et al. (2024) with metrics grouped by model scenario (null / alt)
#' and columns for each analysis method.
#'
#' @param scenario_results Named list. Each element is itself a list with:
#'   \describe{
#'     \item{results}{\code{data.table} from
#'       \code{\link{run_simulation_analysis}}.}
#'     \item{label}{Character scenario label, e.g., \code{"M1"}.}
#'     \item{n_sample}{Integer sample size.}
#'     \item{dgm}{DGM object (for true HRs and subgroup prevalence).}
#'     \item{hypothesis}{Character: \code{"null"} or \code{"alt"}.}
#'   }
#' @param analyses Character vector of analysis labels to include
#'   (e.g., \code{c("FS", "FSlg", "GRF")}). When \code{NULL}, all unique
#'   values of \code{results$analysis} across scenarios are used.
#' @param digits Integer. Decimal places for proportions. Default: 2.
#' @param title Character. Table title. Default:
#'   \code{"Subgroup Identification and Classification Rates"}.
#' @param n_sims Integer. Number of simulations (for subtitle). Default:
#'   \code{NULL}.
#' @param bold_threshold Numeric. Type I error threshold above which the
#'   \code{any(H)} value is shown in bold. Set \code{NULL} to disable.
#'   Default: 0.05.
#'
#' @return A \code{gt} table object.
#'
#' @details
#' For each scenario the function computes:
#' \itemize{
#'   \item \code{any(H)}: Proportion of simulations identifying any subgroup.
#'   \item \code{sens(H)}: Mean sensitivity (only under alternative).
#'   \item \code{sens(Hc)}: Mean specificity.
#'   \item \code{ppv(H)}: Mean positive predictive value (only under
#'     alternative).
#'   \item \code{ppv(Hc)}: Mean negative predictive value.
#'   \item \code{avg|H|}: Mean size of identified subgroup (when found).
#' }
#'
#' Under the null hypothesis the rows are reduced to \code{any(H)},
#' \code{sens(Hc)}, \code{ppv(Hc)}, and \code{avg|H|}.
#'
#' @examples
#' \dontrun{
#' # Assemble results from H0 and H1 simulations
#' scenarios <- list(
#'   null = list(
#'     results = results_null, label = "M1",
#'     n_sample = 700, dgm = dgm_null, hypothesis = "null"
#'   ),
#'   alt = list(
#'     results = results_alt, label = "M1",
#'     n_sample = 700, dgm = dgm_calibrated, hypothesis = "alt"
#'   )
#' )
#'
#' build_classification_table(scenarios, n_sims = 100)
#' }
#'
#' @seealso \code{\link{format_oc_results}},
#'   \code{\link{summarize_simulation_results}}
#'
#' @importFrom data.table as.data.table rbindlist
#' @importFrom gt gt tab_header cols_hide cols_label tab_style cell_text
#'   tab_options px cells_body cells_row_groups
#' @export
build_classification_table <- function(
    scenario_results,
    analyses = NULL,
    digits = 2,
    title = "Subgroup Identification and Classification Rates",
    n_sims = NULL,
    bold_threshold = 0.05
) {

  rows_list <- list()

 for (sc_name in names(scenario_results)) {
    sc  <- scenario_results[[sc_name]]
    res <- data.table::as.data.table(sc$results)
    dgm <- sc$dgm
    hyp <- sc$hypothesis
    n   <- sc$n_sample
    lab <- sc$label

    # Determine analyses if not specified
    if (is.null(analyses)) {
      analyses_use <- sort(unique(res$analysis))
    } else {
      analyses_use <- analyses
    }

    # ── Section header label ────────────────────────────────────────────────
    if (hyp == "null") {
      hr_itt <- round(dgm$hr_causal, 2)
      section_label <- sprintf(
        "%s Null: N=%d, theta(ITT) = %s",
        lab, n, hr_itt
      )
    } else {
      prop_H <- round(
        100 * mean(dgm$df_super_rand$flag.harm, na.rm = TRUE), 0
      )
      hr_H   <- round(dgm$hr_H_true, 2)
      hr_Hc  <- round(dgm$hr_Hc_true, 2)
      hr_itt <- round(dgm$hr_causal, 2)
      section_label <- sprintf(
        "%s Alt: N=%d, p_H=%d%%, theta(H)=%s, theta(Hc)=%s, theta(ITT)=%s",
        lab, n, prop_H, hr_H, hr_Hc, hr_itt
      )
    }

    # ── Metric names depend on hypothesis ───────────────────────────────────
    metric_names <- if (hyp == "null") {
      c("any(H)", "sens(Hc)", "ppv(Hc)", "avg|H|")
    } else {
      c("any(H)", "sens(H)", "sens(Hc)", "ppv(H)", "ppv(Hc)", "avg|H|")
    }

    for (metric in metric_names) {
      row_vals <- list(Scenario = section_label, Metric = metric)

      for (a in analyses_use) {
        r <- res[res$analysis == a, ]
        r_found <- r[r$any.H == 1, ]

        val <- switch(
          metric,
          "any(H)"  = mean(r$any.H, na.rm = TRUE),
          "sens(H)"  = mean(r$sens, na.rm = TRUE),
          "sens(Hc)" = mean(r$spec, na.rm = TRUE),
          "ppv(H)"   = mean(r$ppv, na.rm = TRUE),
          "ppv(Hc)"  = mean(r$npv, na.rm = TRUE),
          "avg|H|"   = if (nrow(r_found) > 0) {
            round(mean(r_found$size.H, na.rm = TRUE), 0)
          } else {
            NA_real_
          }
        )

        row_vals[[a]] <- if (metric == "avg|H|") val else round(val, digits)
      }

      rows_list[[length(rows_list) + 1]] <- as.data.frame(
        row_vals, stringsAsFactors = FALSE
      )
    }
  }

  table_df <- data.table::rbindlist(rows_list, fill = TRUE)

  # ── Grouping variable ───────────────────────────────────────────────────
  table_df[, Group := Scenario]

  gt_tbl <- gt::gt(table_df, groupname_col = "Group") |>
    gt::cols_hide(columns = "Scenario") |>
    gt::tab_header(
      title = title,
      subtitle = if (!is.null(n_sims)) {
        sprintf(
          "Across %s simulations per scenario",
          format(n_sims, big.mark = ",")
        )
      }
    ) |>
    gt::cols_label(Metric = "") |>
    gt::tab_style(
      style = gt::cell_text(size = "small"),
      locations = gt::cells_body()
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold", size = "small"),
      locations = gt::cells_row_groups()
    ) |>
    gt::tab_options(
      table.font.size = gt::px(12),
      row_group.padding = gt::px(4)
    )

  # ── Bold inflated Type I error values ───────────────────────────────────
  if (!is.null(bold_threshold)) {
    analyses_in_tbl <- setdiff(
      names(table_df),
      c("Scenario", "Metric", "Group")
    )

    for (col_name in analyses_in_tbl) {
      bold_rows <- which(
        table_df$Metric == "any(H)" &
          !is.na(table_df[[col_name]]) &
          is.numeric(table_df[[col_name]]) &
          table_df[[col_name]] > bold_threshold
      )
      if (length(bold_rows) > 0) {
        gt_tbl <- gt::tab_style(
          gt_tbl,
          style = gt::cell_text(weight = "bold"),
          locations = gt::cells_body(columns = col_name, rows = bold_rows)
        )
      }
    }
  }

  gt_tbl
}


# =============================================================================
# Estimation Properties Table
# =============================================================================

#' Build Estimation Properties Table from Simulation Results
#'
#' Constructs a publication-quality \code{gt} table summarizing estimation
#' properties for hazard ratios in the identified subgroup and its complement.
#' The layout mirrors the estimation table in Leon et al. (2024), showing
#' average estimate, empirical SD, min, max, and relative bias for each
#' estimator.
#'
#' @param results \code{data.table} or \code{data.frame}. Simulation results
#'   from \code{\link{run_simulation_analysis}}, optionally enriched with
#'   bootstrap bias-corrected columns (\code{hr.H.bc}, \code{hr.Hc.bc}).
#' @param dgm DGM object. Used for true parameter values (\code{hr_H_true},
#'   \code{hr_Hc_true}).
#' @param analysis_method Character. Which analysis method to tabulate
#'   (e.g., \code{"FSlg"}). Default: \code{"FSlg"}.
#' @param n_boots Integer. Number of bootstraps (for subtitle). Default: 300.
#' @param digits Integer. Decimal places. Default: 2.
#' @param title Character. Table title.
#'
#' @return A \code{gt} table object, or \code{NULL} if no estimable
#'   realizations exist.
#'
#' @details
#' For each subgroup (H and Hc) the function reports:
#' \itemize{
#'   \item \strong{Avg}: Mean of the estimates across estimable simulations.
#'   \item \strong{SD}: Empirical standard deviation.
#'   \item \strong{Min / Max}: Range.
#'   \item \strong{Rel Bias}: Relative bias in percent, \code{100 * (Avg - true) / true}.
#' }
#'
#' When bootstrap-corrected columns (\code{hr.H.bc}, \code{hr.Hc.bc}) are
#' present in \code{results}, an additional bias-corrected row
#' (\code{"theta-hat*(H)"} or \code{"theta-hat*(Hc)"}) is added per subgroup.
#'
#' @examples
#' \dontrun{
#' build_estimation_table(
#'   results = results_alt,
#'   dgm = dgm_calibrated,
#'   analysis_method = "FSlg"
#' )
#' }
#'
#' @seealso \code{\link{build_classification_table}},
#'   \code{\link{format_oc_results}}
#'
#' @importFrom data.table as.data.table rbindlist
#' @importFrom gt gt tab_header cols_label tab_style cell_text tab_options px
#'   cells_body cells_row_groups
#' @export
build_estimation_table <- function(
    results,
    dgm,
    analysis_method = "FSlg",
    n_boots = 300,
    digits = 2,
    title = "Estimation Properties"
) {

  res <- data.table::as.data.table(results)

  if ("analysis" %in% names(res)) {
    res <- res[res$analysis == analysis_method, ]
  }

  res_found <- res[res$any.H == 1, ]
  n_estimable <- nrow(res_found)

  if (n_estimable == 0) {
    message("No estimable realizations for analysis = ", analysis_method)
    return(NULL)
  }

  # ── True values from DGM ──────────────────────────────────────────────────
  theta_H_true  <- dgm$hr_H_true
  theta_Hc_true <- dgm$hr_Hc_true
  avg_size_H    <- round(mean(res_found$size.H, na.rm = TRUE), 0)

  # ── Helper: compute one estimation row ────────────────────────────────────
  make_est_row <- function(estimates, true_val, label) {
    est <- estimates[!is.na(estimates)]
    if (length(est) == 0) return(NULL)

    avg_est <- mean(est)
    sd_est  <- sd(est)
    min_est <- min(est)
    max_est <- max(est)
    rel_bias <- 100 * (avg_est - true_val) / true_val

    data.frame(
      Estimator      = label,
      Avg            = round(avg_est, digits),
      SD             = round(sd_est, digits),
      Min            = round(min_est, digits),
      Max            = round(max_est, digits),
      `Rel Bias (%)` = round(rel_bias, digits),
      check.names    = FALSE,
      stringsAsFactors = FALSE
    )
  }

  # ── Rows for H (harm subgroup) ────────────────────────────────────────────
  rows_H <- list()

  if ("hr.H.hat" %in% names(res_found)) {
    rows_H[[1]] <- make_est_row(
      res_found$hr.H.hat, theta_H_true, "theta-hat(H-hat)"
    )
  }

  if ("hr.H.bc" %in% names(res_found)) {
    rows_H[[length(rows_H) + 1]] <- make_est_row(
      res_found$hr.H.bc, theta_H_true, "theta-hat*(H-hat)"
    )
  }

  # ── Rows for Hc (complement) ──────────────────────────────────────────────
  rows_Hc <- list()

  if ("hr.Hc.hat" %in% names(res_found)) {
    rows_Hc[[1]] <- make_est_row(
      res_found$hr.Hc.hat, theta_Hc_true, "theta-hat(Hc-hat)"
    )
  }

  if ("hr.Hc.bc" %in% names(res_found)) {
    rows_Hc[[length(rows_Hc) + 1]] <- make_est_row(
      res_found$hr.Hc.bc, theta_Hc_true, "theta-hat*(Hc-hat)"
    )
  }

  # ── Combine ───────────────────────────────────────────────────────────────
  df_H  <- if (length(rows_H) > 0) {
    data.table::rbindlist(rows_H, fill = TRUE)
  }
  df_Hc <- if (length(rows_Hc) > 0) {
    data.table::rbindlist(rows_Hc, fill = TRUE)
  }

  if (!is.null(df_H)) {
    df_H[, Subgroup := sprintf(
      "H: %d estimable, avg |H| = %d, theta(H) = %s",
      n_estimable, avg_size_H, round(theta_H_true, 2)
    )]
  }
  if (!is.null(df_Hc)) {
    df_Hc[, Subgroup := sprintf(
      "Hc: avg |Hc| = %d, theta(Hc) = %s",
      round(mean(res_found$size.Hc, na.rm = TRUE), 0),
      round(theta_Hc_true, 2)
    )]
  }

  table_df <- data.table::rbindlist(list(df_H, df_Hc), fill = TRUE)

  if (nrow(table_df) == 0) {
    message("No estimation columns found in results.")
    return(NULL)
  }

  # ── gt table ──────────────────────────────────────────────────────────────
  gt_tbl <- gt::gt(table_df, groupname_col = "Subgroup") |>
    gt::tab_header(
      title = title,
      subtitle = sprintf(
        "%s: %d estimable realizations (B = %d bootstraps)",
        analysis_method, n_estimable, n_boots
      )
    ) |>
    gt::cols_label(Estimator = "") |>
    gt::tab_style(
      style = gt::cell_text(size = "small"),
      locations = gt::cells_body()
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold", size = "small"),
      locations = gt::cells_row_groups()
    ) |>
    gt::tab_options(
      table.font.size = gt::px(12),
      row_group.padding = gt::px(4)
    )

  gt_tbl
}


# =============================================================================
# Render LaTeX Table Data as gt
# =============================================================================

#' Render Reference Simulation Table as gt
#'
#' Converts a data frame of pre-computed reference simulation results (e.g.,
#' digitized from a published LaTeX table) into a styled \code{gt} table. This
#' is useful for displaying published benchmark results alongside new
#' simulation output within vignettes or reports.
#'
#' @param ref_df \code{data.frame}. Must contain a \code{Metric} column and
#'   one column per analysis method, with a \code{Scenario} column for row
#'   grouping.
#' @param title Character. Table title.
#' @param subtitle Character. Table subtitle. Default: \code{NULL}.
#' @param bold_threshold Numeric. Values in \code{any(H)} rows exceeding this
#'   threshold are shown in bold. Set \code{NULL} to disable. Default: 0.05.
#'
#' @return A \code{gt} table object.
#'
#' @examples
#' \dontrun{
#' ref <- data.frame(
#'   Scenario = "M1 Null: N=700",
#'   Metric   = "any(H)",
#'   FS       = 0.02,
#'   FSlg     = 0.03,
#'   GRF      = 0.25
#' )
#' render_reference_table(ref, title = "Reference Results")
#' }
#'
#' @importFrom gt gt tab_header cols_label tab_style cell_text tab_options px
#'   cells_body cells_row_groups
#' @export
render_reference_table <- function(
    ref_df,
    title = "Reference Simulation Results",
    subtitle = NULL,
    bold_threshold = 0.05
) {

  ref_df <- as.data.frame(ref_df, stringsAsFactors = FALSE)

  gt_tbl <- gt::gt(ref_df, groupname_col = "Scenario") |>
    gt::tab_header(title = title, subtitle = subtitle) |>
    gt::cols_label(Metric = "") |>
    gt::tab_style(
      style = gt::cell_text(size = "small"),
      locations = gt::cells_body()
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold", size = "small"),
      locations = gt::cells_row_groups()
    ) |>
    gt::tab_options(
      table.font.size = gt::px(12),
      row_group.padding = gt::px(4)
    )

  # Bold inflated Type I error
  if (!is.null(bold_threshold)) {
    analysis_cols <- setdiff(names(ref_df), c("Scenario", "Metric"))
    for (col_name in analysis_cols) {
      bold_rows <- which(
        ref_df$Metric == "any(H)" &
          !is.na(ref_df[[col_name]]) &
          is.numeric(ref_df[[col_name]]) &
          ref_df[[col_name]] > bold_threshold
      )
      if (length(bold_rows) > 0) {
        gt_tbl <- gt::tab_style(
          gt_tbl,
          style = gt::cell_text(weight = "bold"),
          locations = gt::cells_body(columns = col_name, rows = bold_rows)
        )
      }
    }
  }

  gt_tbl
}
