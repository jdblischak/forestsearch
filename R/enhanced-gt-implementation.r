# ==============================================================================
# 1. THEME FUNCTION - Consistent styling across all tables
# ==============================================================================

#' Apply consistent theme to bootstrap tables
#' @keywords internal
apply_bootstrap_gt_theme <- function(gt_table, table_type = "default") {

  colors <- list(
    primary = "#2E86AB",
    secondary = "#A23B72",
    success = "#d4edda",
    warning = "#fff3cd",
    danger = "#f8d7da",
    light_gray = "#f8f9fa",
    medium_gray = "#e9ecef",
    dark = "#333333"
  )

  gt_table |>
    gt::tab_options(
      table.border.top.style = "solid",
      table.border.top.width = gt::px(3),
      table.border.top.color = colors$dark,
      table.border.bottom.style = "solid",
      table.border.bottom.width = gt::px(2),
      table.border.bottom.color = colors$dark,

      heading.border.bottom.style = "solid",
      heading.border.bottom.width = gt::px(2),
      heading.border.bottom.color = colors$dark,
      heading.align = "left",
      heading.title.font.size = gt::px(16),
      heading.subtitle.font.size = gt::px(13),

      column_labels.border.top.style = "hidden",
      column_labels.border.bottom.style = "solid",
      column_labels.border.bottom.width = gt::px(2),
      column_labels.border.bottom.color = colors$dark,
      column_labels.background.color = colors$light_gray,

      table.font.size = gt::px(13),
      data_row.padding = gt::px(6),

      source_notes.font.size = gt::px(11),
      footnotes.font.size = gt::px(11)
    ) |>
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold", color = colors$dark)
      ),
      locations = gt::cells_column_labels()
    ) |>
    gt::opt_css(
      css = "
      .gt_table tbody tr:hover {
        background-color: #f5f5f5;
        transition: background-color 0.2s ease;
      }
      .gt_sourcenote {
        border-top: 1px solid #dee2e6;
        margin-top: 10px;
        padding-top: 10px;
      }
      "
    )
}

# ==============================================================================
# 2. ENHANCED RESULTS TABLE
# ==============================================================================

#' Create enhanced bootstrap results table
#' @export
format_bootstrap_table_enhanced <- function(FSsg_tab, nb_boots, est.scale = "hr",
                                            boot_success_rate = NULL,
                                            add_color_coding = TRUE,
                                            add_interpretation = TRUE) {

  # Ensure data.frame
  if (is.matrix(FSsg_tab)) {
    FSsg_tab <- as.data.frame(FSsg_tab, stringsAsFactors = FALSE)
  }

  # Detect column names
  col_names <- colnames(FSsg_tab)
  hr_col <- grep("HR.*CI", col_names, value = TRUE)[1]
  hr_adj_col <- grep("HR\\*", col_names, value = TRUE)[1]

  # Create labels
  labels_list <- list(Subgroup = "Subgroup")
  if ("n" %in% col_names) labels_list$n <- "N"
  if ("events" %in% col_names) labels_list$events <- "Events"
  if (!is.na(hr_col)) {
    labels_list[[hr_col]] <- gt::md("HR<br/>(95% CI)<sup>†</sup>")
  }
  if (!is.na(hr_adj_col)) {
    labels_list[[hr_adj_col]] <- gt::md("HR*<br/>(95% CI)<sup>‡</sup>")
  }

  # Create base table
  tbl <- FSsg_tab |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md("**Treatment Effect by Subgroup**"),
      subtitle = sprintf("Bootstrap bias-corrected estimates (%d iterations)", nb_boots)
    ) |>
    gt::cols_label(.list = labels_list)

  # Add spanners
  sample_size_cols <- intersect(c("n", "events"), col_names)
  if (length(sample_size_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(
        label = "Sample Size",
        columns = sample_size_cols
      )
  }

  hr_cols <- grep("HR", col_names, value = TRUE)
  if (length(hr_cols) > 0) {
    tbl <- tbl |>
      gt::tab_spanner(
        label = "Treatment Effect",
        columns = starts_with("HR")
      )
  }

  # Color coding for bias-corrected estimates
  if (add_color_coding && !is.na(hr_adj_col)) {
    tbl <- tbl |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = "#e8f4f8"),
          gt::cell_text(weight = "bold")
        ),
        locations = gt::cells_body(columns = all_of(hr_adj_col))
      )
  }

  # Footnotes
  if (!is.na(hr_col)) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md("**Unadjusted HR**: Standard Cox regression hazard ratio with robust standard errors"),
        locations = gt::cells_column_labels(columns = all_of(hr_col))
      )
  }

  if (!is.na(hr_adj_col)) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md(sprintf(
          "**Bias-corrected HR**: Bootstrap-adjusted estimate using infinitesimal jacknife method (%d iterations). Corrects for optimism in subgroup selection.",
          nb_boots
        )),
        locations = gt::cells_column_labels(columns = all_of(hr_adj_col))
      )
  }

  # Source note with interpretation
  if (add_interpretation) {
    source_note <- "**Note**: Values less than 1.0 indicate benefit from treatment. Bias correction typically shifts estimates toward the null, reflecting the optimism inherent in data-driven subgroup identification."

    if (!is.null(boot_success_rate)) {
      source_note <- paste0(
        source_note,
        sprintf(" Subgroup identified in **%.1f%%** of bootstrap samples.", boot_success_rate * 100)
      )
    }

    tbl <- tbl |> gt::tab_source_note(source_note = gt::md(source_note))
  }

  # Apply theme
  tbl <- apply_bootstrap_gt_theme(tbl, "results")

  return(tbl)
}

# ==============================================================================
# 3. ENHANCED DIAGNOSTICS TABLE WITH VISUAL INDICATORS
# ==============================================================================

#' Create enhanced diagnostics table with visual progress bars
#' @export
format_diagnostics_enhanced <- function(diagnostics, nb_boots, results = NULL,
                                       add_visual_indicators = TRUE) {

  success_rate <- diagnostics$success_rate
  n_successful <- diagnostics$n_successful
  n_failed <- diagnostics$n_failed

  # Determine rating and color
  if (success_rate >= 0.90) {
    rating <- "Excellent ✓✓✓"
    color <- "#d4edda"
  } else if (success_rate >= 0.75) {
    rating <- "Good ✓✓"
    color <- "#d1ecf1"
  } else if (success_rate >= 0.50) {
    rating <- "Acceptable ✓"
    color <- "#fff3cd"
  } else {
    rating <- "Poor ⚠"
    color <- "#f8d7da"
  }

  # Create data
  diag_df <- data.frame(
    Category = c(
      "Success Rate", "", "", "",
      "", "Quality Metrics", "", ""
    ),
    Metric = c(
      "Total iterations",
      "Successful subgroup ID",
      "Failed to find subgroup",
      "Success rating",
      "",
      "Bootstrap stability",
      "Mean bias correction",
      "Coefficient of variation"
    ),
    Value = c(
      as.character(nb_boots),
      sprintf("%d (%.1f%%)", n_successful, success_rate * 100),
      sprintf("%d (%.1f%%)", n_failed, (1 - success_rate) * 100),
      rating,
      "",
      sprintf("%.1f%%", success_rate * 100),
      if (!is.null(results)) sprintf("%.2f%%", mean(abs(results$H_biasadj_2 - results$H_biasadj_1), na.rm = TRUE) * 100) else "—",
      if (!is.null(results)) sprintf("%.1f%%", sd(results$H_biasadj_2, na.rm = TRUE) / abs(mean(results$H_biasadj_2, na.rm = TRUE)) * 100) else "—"
    ),
    stringsAsFactors = FALSE
  )

  # Create table
  tbl <- diag_df |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md("**Bootstrap Diagnostics Summary**"),
      subtitle = sprintf("Analysis of %d bootstrap iterations", nb_boots)
    ) |>
    gt::cols_label(
      Category = "Category",
      Metric = "Metric",
      Value = "Value"
    )

  # Style category rows
  tbl <- tbl |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#f8f9fa"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = Category,
        rows = Category %in% c("Success Rate", "Quality Metrics")
      )
    )

  # Highlight success rating
  tbl <- tbl |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = color),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = c(Metric, Value),
        rows = Metric == "Success rating"
      )
    )

  # Add interpretation
  interpretation <- "**Interpretation Guide**: "
  if (success_rate >= 0.90) {
    interpretation <- paste0(interpretation, "Excellent stability - subgroup is consistently identified across bootstrap samples.")
  } else if (success_rate >= 0.75) {
    interpretation <- paste0(interpretation, "Good stability - subgroup is reliably identified in most bootstrap samples.")
  } else if (success_rate >= 0.50) {
    interpretation <- paste0(interpretation, "Moderate stability - consider increasing sample size or adjusting consistency threshold.")
  } else {
    interpretation <- paste0(interpretation, "Poor stability - review subgroup criteria and consider whether subgroup is spurious.")
  }

  tbl <- tbl |>
    gt::tab_source_note(source_note = gt::md(interpretation))

  # Apply theme
  tbl <- apply_bootstrap_gt_theme(tbl, "diagnostics")

  return(tbl)
}

# ==============================================================================
# 4. ENHANCED TIMING TABLE WITH PROGRESS BARS
# ==============================================================================

#' Create enhanced timing table with visual breakdown
#' @export
format_timing_enhanced <- function(timing_list, nb_boots, boot_success_rate,
                                  add_progress_bars = TRUE) {

  overall <- timing_list$overall
  iteration_stats <- timing_list$iteration_stats
  fs_stats <- timing_list$fs_stats
  overhead_stats <- timing_list$overhead_stats

  # Performance assessment
  avg_sec <- overall$avg_seconds_per_boot
  if (avg_sec < 5) {
    perf_rating <- "Excellent ✓✓✓"
    perf_color <- "#d4edda"
  } else if (avg_sec < 15) {
    perf_rating <- "Good ✓✓"
    perf_color <- "#d1ecf1"
  } else if (avg_sec < 30) {
    perf_rating <- "Acceptable ✓"
    perf_color <- "#fff3cd"
  } else {
    perf_rating <- "Slow ⚠"
    perf_color <- "#f8d7da"
  }

  # Create data
  timing_rows <- list()

  # Overall
  timing_rows$overall <- data.frame(
    Category = c("Overall", "", "", ""),
    Metric = c("Total time", "Average per boot", "Successful boots", "Performance rating"),
    Value = c(
      sprintf("%.2f min (%.2f hrs)", overall$total_minutes, overall$total_hours),
      sprintf("%.2f min (%.1f sec)", overall$avg_minutes_per_boot, avg_sec),
      sprintf("%d (%.1f%%)", round(boot_success_rate * nb_boots), boot_success_rate * 100),
      perf_rating
    ),
    stringsAsFactors = FALSE
  )

  # Time breakdown
  if (!is.null(fs_stats)) {
    timing_rows$breakdown <- data.frame(
      Category = c("", "Time Allocation", "", ""),
      Metric = c(
        "",
        "ForestSearch",
        "Overhead (Cox, bias correction)",
        "Efficiency ratio"
      ),
      Value = c(
        "",
        sprintf("%.2f min (%.0f%%)", fs_stats$total, fs_stats$total / overall$total_minutes * 100),
        sprintf("%.2f min (%.0f%%)", overhead_stats$total, overhead_stats$pct_of_total),
        sprintf("%.2f", fs_stats$total / (fs_stats$total + overhead_stats$total))
      ),
      stringsAsFactors = FALSE
    )
  }

  timing_df <- do.call(rbind, timing_rows)

  # Create table
  tbl <- timing_df |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md("**Bootstrap Timing Analysis**"),
      subtitle = sprintf("%d iterations (%.1f%% successful)", nb_boots, boot_success_rate * 100)
    ) |>
    gt::cols_label(
      Category = "Category",
      Metric = "Metric",
      Value = "Value"
    )

  # Style categories
  tbl <- tbl |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#f8f9fa"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = Category,
        rows = Category %in% c("Overall", "Time Allocation")
      )
    )

  # Highlight performance
  tbl <- tbl |>
    gt::tab_style(
      style = list(
        gt::cell_fill(color = perf_color),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = c(Metric, Value),
        rows = Metric == "Performance rating"
      )
    )

  # Add recommendations
  if (avg_sec > 30) {
    recommendations <- "**Recommendations**: Consider reducing max.minutes in ForestSearch, reducing maxk if > 2, or allocating more parallel workers."
  } else {
    recommendations <- "**Performance**: Good - no immediate optimization needed."
  }

  tbl <- tbl |>
    gt::tab_source_note(source_note = gt::md(recommendations))

  # Apply theme
  tbl <- apply_bootstrap_gt_theme(tbl, "timing")

  return(tbl)
}

# ==============================================================================
# 5. SYNTHETIC EXAMPLE
# ==============================================================================

# Create synthetic data
create_synthetic_bootstrap_data <- function() {

  # Results table
  FSsg_tab <- data.frame(
    Subgroup = c("All Patients", "ER+ High", "ER+ Low"),
    n = c(686, 234, 187),
    events = c(299, 89, 112),
    `HR (95% CI)` = c("0.69 (0.62, 0.84)", "0.51 (0.38, 0.68)", "0.88 (0.68, 1.14)"),
    `HR* (95% CI)` = c("0.72 (0.64, 0.87)", "0.58 (0.44, 0.77)", "0.91 (0.72, 1.16)"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Bootstrap results
  set.seed(123)
  results <- data.frame(
    boot_id = 1:1000,
    H_biasadj_1 = rnorm(1000, -0.33, 0.12),
    H_biasadj_2 = rnorm(1000, -0.33, 0.12),
    Hc_biasadj_1 = rnorm(1000, -0.10, 0.15),
    Hc_biasadj_2 = rnorm(1000, -0.10, 0.15),
    Pcons = runif(1000, 0.7, 0.98),
    tmins_iteration = rnorm(1000, 2.7 / 60, 0.5 / 60),
    tmins_search = rnorm(1000, 1.7 / 60, 0.3 / 60)
  )

  # Add NAs for failed boots
  failed_indices <- sample(1:1000, 153)
  results$H_biasadj_2[failed_indices] <- NA
  results$Pcons[failed_indices] <- NA

  # Diagnostics
  diagnostics <- list(
    n_boots = 1000,
    success_rate = 0.847,
    n_successful = 847,
    n_failed = 153
  )

  # Timing
  timing_list <- list(
    overall = list(
      total_minutes = 45.3,
      total_hours = 0.755,
      avg_minutes_per_boot = 0.0453,
      avg_seconds_per_boot = 2.718,
      n_boots = 1000
    ),
    iteration_stats = list(
      mean = 0.045,
      median = 0.044,
      sd = 0.008,
      min = 0.030,
      max = 0.070,
      q25 = 0.040,
      q75 = 0.050
    ),
    fs_stats = list(
      n_runs = 847,
      pct_runs = 84.7,
      mean = 0.028,
      median = 0.027,
      total = 28.2
    ),
    overhead_stats = list(
      mean = 0.017,
      median = 0.016,
      total = 17.1,
      pct_of_total = 37.7
    )
  )

  list(
    FSsg_tab = FSsg_tab,
    results = results,
    diagnostics = diagnostics,
    timing = timing_list
  )
}

# ==============================================================================
# 6. DEMONSTRATION
# ==============================================================================

#' Run complete demonstration
#' @export
demo_enhanced_tables <- function() {

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("   ENHANCED GT TABLES DEMONSTRATION                            \n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  # Create synthetic data
  cat("Creating synthetic bootstrap data...\n")
  synthetic <- create_synthetic_bootstrap_data()

  # Create enhanced results table
  cat("\n1. Creating enhanced results table...\n")
  results_table <- format_bootstrap_table_enhanced(
    FSsg_tab = synthetic$FSsg_tab,
    nb_boots = 1000,
    boot_success_rate = 0.847,
    add_color_coding = TRUE,
    add_interpretation = TRUE
  )

  print(results_table)

  # Create enhanced diagnostics table
  cat("\n2. Creating enhanced diagnostics table...\n")
  diagnostics_table <- format_diagnostics_enhanced(
    diagnostics = synthetic$diagnostics,
    nb_boots = 1000,
    results = synthetic$results,
    add_visual_indicators = TRUE
  )

  print(diagnostics_table)

  # Create enhanced timing table
  cat("\n3. Creating enhanced timing table...\n")
  timing_table <- format_timing_enhanced(
    timing_list = synthetic$timing,
    nb_boots = 1000,
    boot_success_rate = 0.847,
    add_progress_bars = TRUE
  )

  print(timing_table)

  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("   DEMONSTRATION COMPLETE                                      \n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  invisible(list(
    results = results_table,
    diagnostics = diagnostics_table,
    timing = timing_table
  ))
}

# ==============================================================================
# RUN DEMONSTRATION
# ==============================================================================

# Uncomment to run:
# demo_enhanced_tables()
