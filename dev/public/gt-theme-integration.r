      overhead_stats = if (all(c("tmins_iteration", "tmins_search") %in% names(results))) {
        overhead_times <- results$tmins_iteration - results$tmins_search
        overhead_times <- overhead_times[!is.na(overhead_times)]
        if (length(overhead_times) > 0) {
          list(
            mean = mean(overhead_times),
            median = median(overhead_times),
            total = sum(overhead_times),
            pct_of_total = sum(overhead_times) / attr(results, "timing")$total_minutes * 100
          )
        } else NULL
      } else NULL
    )
    
    timing_table <- format_bootstrap_timing_table_themed(
      timing_list = timing_list,
      nb_boots = nb_boots,
      boot_success_rate = boot_success_rate
    )
  }
  
  # Compile output
  output <- list(
    table = formatted_table,
    diagnostics = diagnostics,
    diagnostics_table = diagnostics_table,
    timing_table = timing_table,
    plots = if (create_plots) create_bootstrap_diagnostic_plots(results, H_estimates, Hc_estimates) else NULL
  )
  
  # Print summary
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("           BOOTSTRAP ANALYSIS SUMMARY (THEMED)                 \n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  cat("SUCCESS METRICS:\n")
  cat(sprintf("  Total iterations:        %d\n", nb_boots))
  cat(sprintf("  Successful:              %d (%.1f%%)\n", 
              diagnostics$n_successful, boot_success_rate * 100))
  cat(sprintf("  Failed:                  %d (%.1f%%)\n",
              diagnostics$n_failed, (1 - boot_success_rate) * 100))
  cat("\n")
  
  cat("Themed tables created with consistent styling:\n")
  cat("  ✓ Results table (apply_results_theme)\n")
  cat("  ✓ Diagnostics table (apply_diagnostics_theme)\n")
  if (!is.null(timing_table)) {
    cat("  ✓ Timing table (apply_timing_theme)\n")
  }
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  invisible(output)
}

# ==============================================================================
# STEP 4: Create Convenience Function for Theme Switching
# ==============================================================================

#' Apply Custom Theme to Bootstrap Table
#' 
#' Allows users to apply custom themes or modify existing themes.
#' Useful for matching institutional branding or publication requirements.
#' 
#' @param gt_table A gt table object
#' @param theme Character. One of "results", "diagnostics", "timing", "subgroup", "custom"
#' @param custom_colors Named list of colors (optional, for theme="custom")
#' @return Modified gt table object
#' @export
#' 
#' @examples
#' \dontrun{
#' # Use built-in theme
#' my_table |> apply_bootstrap_theme("results")
#' 
#' # Use custom colors
#' custom_colors <- list(
#'   primary = "#003366",
#'   success = "#00AA00"
#' )
#' my_table |> apply_bootstrap_theme("custom", custom_colors = custom_colors)
#' }
apply_bootstrap_theme <- function(gt_table, theme = "results", custom_colors = NULL) {
  
  theme <- match.arg(theme, c("results", "diagnostics", "timing", "subgroup", "custom", "base"))
  
  if (theme == "custom") {
    if (is.null(custom_colors)) {
      warning("No custom_colors provided. Using base theme instead.")
      return(apply_bootstrap_base_theme(gt_table))
    }
    # Apply base theme then override with custom colors
    return(apply_custom_theme(gt_table, custom_colors))
  }
  
  switch(theme,
    base = apply_bootstrap_base_theme(gt_table),
    results = apply_results_theme(gt_table),
    diagnostics = apply_diagnostics_theme(gt_table),
    timing = apply_timing_theme(gt_table),
    subgroup = apply_subgroup_theme(gt_table)
  )
}

#' Apply Custom Theme with User-Specified Colors
#' @keywords internal
apply_custom_theme <- function(gt_table, custom_colors) {
  
  # Start with base theme
  tbl <- apply_bootstrap_base_theme(gt_table)
  
  # Override colors if specified
  if (!is.null(custom_colors$primary)) {
    tbl <- tbl |>
      gt::tab_options(
        table.border.top.color = custom_colors$primary
      )
  }
  
  if (!is.null(custom_colors$success)) {
    # Would need to update specific cell styles
    # This is a simplified example
  }
  
  return(tbl)
}

# ==============================================================================
# STEP 5: Create Theme Preview Function
# ==============================================================================

#' Preview All Bootstrap Themes
#' 
#' Creates a simple table with all theme variations for comparison.
#' Useful for deciding which theme to use or for documentation.
#' 
#' @return List of gt table objects, one for each theme
#' @export
#' 
#' @examples
#' \dontrun{
#' theme_examples <- preview_bootstrap_themes()
#' theme_examples$results
#' theme_examples$diagnostics
#' }
preview_bootstrap_themes <- function() {
  
  # Create sample data
  sample_data <- data.frame(
    Category = c("Metric 1", "Metric 2", "Metric 3"),
    Value1 = c("100", "85.5%", "Good ✓"),
    Value2 = c("200", "92.3%", "Excellent ✓✓")
  )
  
  # Create base table
  base_table <- sample_data |>
    gt::gt() |>
    gt::tab_header(
      title = "Theme Preview",
      subtitle = "Sample bootstrap analysis table"
    )
  
  # Apply each theme
  list(
    base = base_table |> apply_bootstrap_base_theme(),
    results = base_table |> apply_results_theme(),
    diagnostics = base_table |> apply_diagnostics_theme(),
    timing = base_table |> apply_timing_theme(),
    subgroup = base_table |> apply_subgroup_theme()
  )
}

# ==============================================================================
# STEP 6: Integration with Existing Code - Backward Compatibility
# ==============================================================================

#' Wrapper to Maintain Backward Compatibility
#' 
#' Allows existing code to continue working while optionally using new themes.
#' 
#' @param FSsg_tab Data frame or matrix
#' @param nb_boots Integer
#' @param est.scale Character
#' @param boot_success_rate Numeric
#' @param use_theme Logical. Use new theme system? (default: TRUE)
#' @param ... Additional arguments passed to formatting function
#' @return gt table object
#' @export
format_bootstrap_table_compat <- function(FSsg_tab, nb_boots, est.scale = "hr",
                                          boot_success_rate = NULL,
                                          use_theme = TRUE, ...) {
  
  if (use_theme) {
    format_bootstrap_table_themed(
      FSsg_tab = FSsg_tab,
      nb_boots = nb_boots,
      est.scale = est.scale,
      boot_success_rate = boot_success_rate,
      ...
    )
  } else {
    # Call original function (if you want to keep old version)
    format_bootstrap_table(
      FSsg_tab = FSsg_tab,
      nb_boots = nb_boots,
      est.scale = est.scale,
      boot_success_rate = boot_success_rate,
      ...
    )
  }
}

# ==============================================================================
# STEP 7: Complete Working Example with Current Codebase
# ==============================================================================

#' Complete Example Using Themed Tables
#' 
#' Demonstrates how to use the theme system with real bootstrap results
#' 
#' @export
example_themed_bootstrap_workflow <- function() {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("   THEMED BOOTSTRAP TABLES - COMPLETE WORKFLOW                 \n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  # Simulate bootstrap results (in practice, these come from your analysis)
  cat("1. Creating synthetic bootstrap results...\n")
  
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
  
  # Bootstrap iterations
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
  
  # Simulate some failures
  failed_idx <- sample(1:1000, 153)
  results$H_biasadj_2[failed_idx] <- NA
  results$Pcons[failed_idx] <- NA
  
  # Add timing metadata
  attr(results, "timing") <- list(
    total_minutes = 45.3,
    total_hours = 0.755,
    avg_minutes_per_boot = 0.0453,
    avg_seconds_per_boot = 2.718,
    n_boots = 1000
  )
  
  # Create bootstrap results object
  boot_results <- list(
    FSsg_tab = FSsg_tab,
    results = results,
    H_estimates = list(H0 = 0.69, H2 = 0.72),
    Hc_estimates = list(H0 = 0.88, H2 = 0.91)
  )
  
  cat("   ✓ Synthetic data created\n\n")
  
  # ===========================================================================
  # METHOD 1: Using summarize_bootstrap_results_themed (recommended)
  # ===========================================================================
  
  cat("2. Creating themed summary (Method 1 - Recommended)...\n")
  
  summary <- summarize_bootstrap_results_themed(
    sgharm = NULL,  # Original subgroup definition (optional)
    boot_results = boot_results,
    create_plots = FALSE,
    est.scale = "hr"
  )
  
  cat("   ✓ All themed tables created\n\n")
  
  # Display tables
  cat("3. Displaying themed tables:\n\n")
  
  cat("   Results Table:\n")
  print(summary$table)
  
  cat("\n   Diagnostics Table:\n")
  print(summary$diagnostics_table)
  
  if (!is.null(summary$timing_table)) {
    cat("\n   Timing Table:\n")
    print(summary$timing_table)
  }
  
  # ===========================================================================
  # METHOD 2: Using individual themed functions
  # ===========================================================================
  
  cat("\n4. Alternative approach - individual themed functions:\n\n")
  
  # Create each table individually
  table_results <- format_bootstrap_table_themed(
    FSsg_tab = FSsg_tab,
    nb_boots = 1000,
    boot_success_rate = 0.847
  )
  
  diagnostics <- list(
    n_boots = 1000,
    success_rate = 0.847,
    n_successful = 847,
    n_failed = 153
  )
  
  table_diagnostics <- format_bootstrap_diagnostics_table_themed(
    diagnostics = diagnostics,
    nb_boots = 1000,
    results = results
  )
  
  cat("   ✓ Individual themed tables created\n\n")
  
  # ===========================================================================
  # METHOD 3: Applying themes to existing tables
  # ===========================================================================
  
  cat("5. Applying themes to existing tables:\n\n")
  
  # Create a basic table
  basic_table <- FSsg_tab |>
    gt::gt() |>
    gt::tab_header(title = "Basic Table")
  
  # Apply different themes
  cat("   Applying results theme...\n")
  themed_table <- basic_table |> apply_bootstrap_theme("results")
  
  cat("   ✓ Theme applied\n\n")
  
  # ===========================================================================
  # METHOD 4: Preview all themes
  # ===========================================================================
  
  cat("6. Previewing all available themes:\n\n")
  
  themes <- preview_bootstrap_themes()
  
  cat("   Available themes:\n")
  cat("   • base        - Foundational styling\n")
  cat("   • results     - For treatment effect tables\n")
  cat("   • diagnostics - For quality metrics\n")
  cat("   • timing      - For performance analysis\n")
  cat("   • subgroup    - For subgroup summaries\n")
  cat("   • custom      - User-defined colors\n\n")
  
  # ===========================================================================
  # METHOD 5: Exporting themed tables
  # ===========================================================================
  
  cat("7. Exporting themed tables:\n\n")
  
  # Save to HTML (recommended for sharing)
  cat("   Saving tables to HTML...\n")
  # gt::gtsave(summary$table, "bootstrap_results_themed.html")
  cat("   • bootstrap_results_themed.html (commented out for demo)\n")
  
  # Save to PNG (for manuscripts)
  # gt::gtsave(summary$table, "bootstrap_results_themed.png")
  cat("   • bootstrap_results_themed.png (commented out for demo)\n\n")
  
  # ===========================================================================
  # Summary
  # ===========================================================================
  
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("   WORKFLOW COMPLETE                                           \n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  cat("Key Features Demonstrated:\n")
  cat("  ✓ Consistent color palette across all tables\n")
  cat("  ✓ Professional borders and spacing\n")
  cat("  ✓ Hover effects for interactivity\n")
  cat("  ✓ Success rate color coding\n")
  cat("  ✓ Category row highlighting\n")
  cat("  ✓ Footnote and source note styling\n")
  cat("  ✓ Accessibility features (focus indicators)\n")
  cat("  ✓ Print-friendly CSS\n\n")
  
  cat("Integration Notes:\n")
  cat("  • Backward compatible with existing code\n")
  cat("  • Use summarize_bootstrap_results_themed() for simplest integration\n")
  cat("  • Individual theme functions available for customization\n")
  cat("  • Theme system easily extensible for new table types\n\n")
  
  invisible(list(
    summary = summary,
    themes = themes,
    individual_tables = list(
      results = table_results,
      diagnostics = table_diagnostics
    )
  ))
}

# ==============================================================================
# STEP 8: Quick Start Guide
# ==============================================================================

#' Quick Start: Apply Themes to Your Bootstrap Results
#' 
#' @examples
#' \dontrun{
#' # OPTION 1: Simplest - use in existing workflow
#' # Replace your existing summarize_bootstrap_results() call with:
#' 
#' summary <- summarize_bootstrap_results_themed(
#'   sgharm = original_subgroup,
#'   boot_results = boot_results,
#'   create_plots = TRUE
#' )
#' 
#' # View themed tables
#' summary$table              # Results with results theme
#' summary$diagnostics_table  # Diagnostics with diagnostics theme
#' summary$timing_table       # Timing with timing theme
#' 
#' 
#' # OPTION 2: Apply themes to existing tables
#' 
#' # Create table as usual
#' my_table <- format_bootstrap_table(FSsg_tab, nb_boots = 1000)
#' 
#' # Then apply theme
#' themed_table <- my_table |> apply_bootstrap_theme("results")
#' 
#' 
#' # OPTION 3: Use themed functions directly
#' 
#' results_table <- format_bootstrap_table_themed(
#'   FSsg_tab = FSsg_tab,
#'   nb_boots = 1000,
#'   boot_success_rate = 0.85
#' )
#' 
#' diagnostics_table <- format_bootstrap_diagnostics_table_themed(
#'   diagnostics = diagnostics,
#'   nb_boots = 1000,
#'   results = results
#' )
#' 
#' 
#' # OPTION 4: Custom colors for institutional branding
#' 
#' custom_colors <- list(
#'   primary = "#003366",    # University blue
#'   success = "#00AA00",    # Custom green
#'   warning = "#FF9900"     # Custom orange
#' )
#' 
#' my_table |> apply_bootstrap_theme("custom", custom_colors = custom_colors)
#' 
#' 
#' # Export themed tables
#' gt::gtsave(themed_table, "results.html")  # For sharing
#' gt::gtsave(themed_table, "results.png")   # For manuscripts
#' }
#' 
#' @name quick_start_themes
NULL

# ==============================================================================
# STEP 9: Migration Guide
# ==============================================================================

#' Migration Guide: Updating Existing Code to Use Themes
#' 
#' @section Before:
#' \preformatted{
#' # Old code (still works!)
#' summary <- summarize_bootstrap_results(
#'   sgharm = original_sg,
#'   boot_results = boot_results
#' )
#' 
#' results_table <- format_bootstrap_table(
#'   FSsg_tab = FSsg_tab,
#'   nb_boots = 1000
#' )
#' }
#' 
#' @section After:
#' \preformatted{
#' # New code with themes
#' summary <- summarize_bootstrap_results_themed(
#'   sgharm = original_sg,
#'   boot_results = boot_results
#' )
#' 
#' results_table <- format_bootstrap_table_themed(
#'   FSsg_tab = FSsg_tab,
#'   nb_boots = 1000
#' )
#' }
#' 
#' @section Backward Compatibility:
#' \preformatted{
#' # Use compatibility wrapper to try themes without breaking existing code
#' results_table <- format_bootstrap_table_compat(
#'   FSsg_tab = FSsg_tab,
#'   nb_boots = 1000,
#'   use_theme = TRUE  # Set to FALSE to use old version
#' )
#' }
#' 
#' @section Testing Themes:
#' \preformatted{
#' # Preview all themes before deciding
#' theme_examples <- preview_bootstrap_themes()
#' 
#' # View each theme
#' theme_examples$base
#' theme_examples$results
#' theme_examples$diagnostics
#' }
#' 
#' @name migration_guide
NULL

# ==============================================================================
# STEP 10: Package-Level Documentation
# ==============================================================================

#' GT Theme System for Bootstrap Analysis
#' 
#' @description
#' This theme system provides consistent, publication-ready styling for all
#' bootstrap analysis tables. It includes:
#' 
#' - Consistent color palette
#' - Professional borders and spacing
#' - Interactive hover effects
#' - Accessibility features
#' - Multiple specialized themes
#' 
#' @section Available Themes:
#' 
#' \describe{
#'   \item{base}{Foundational theme inherited by all others}
#'   \item{results}{For treatment effect and subgroup tables}
#'   \item{diagnostics}{For bootstrap quality metrics}
#'   \item{timing}{For performance and timing analysis}
#'   \item{subgroup}{For subgroup characteristic summaries}
#'   \item{custom}{User-defined color scheme}
#' }
#' 
#' @section Key Functions:
#' 
#' \describe{
#'   \item{apply_bootstrap_theme()}{Apply theme to any gt table}
#'   \item{format_bootstrap_table_themed()}{Create themed results table}
#'   \item{format_bootstrap_diagnostics_table_themed()}{Create themed diagnostics}
#'   \item{format_bootstrap_timing_table_themed()}{Create themed timing table}
#'   \item{summarize_bootstrap_results_themed()}{Complete themed workflow}
#'   \item{preview_bootstrap_themes()}{Preview all available themes}
#' }
#' 
#' @section Color Palette:
#' 
#' The theme system uses a carefully chosen color palette:
#' 
#' \describe{
#'   \item{Primary}{#2E86AB (Blue) - Main brand color}
#'   \item{Secondary}{#A23B72 (Purple) - Accent color}
#'   \item{Tertiary}{#F18F01 (Orange) - Highlight color}
#'   \item{Success}{#28a745 (Green) - Positive indicators}
#'   \item{Warning}{#ffc107 (Yellow) - Caution indicators}
#'   \item{Danger}{#dc3545 (Red) - Problem indicators}
#' }
#' 
#' @section Usage:
#' 
#' Basic usage:
#' \preformatted{
#' # Complete workflow with themes
#' summary <- summarize_bootstrap_results_themed(
#'   sgharm = original_sg,
#'   boot_results = boot_results
#' )
#' 
#' # Individual themed table
#' table <- format_bootstrap_table_themed(
#'   FSsg_tab = results,
#'   nb_boots = 1000
#' )
#' 
#' # Apply theme to existing table
#' themed_table <- my_table |> apply_bootstrap_theme("results")
#' }
#' 
#' @section Integration:
#' 
#' The theme system integrates seamlessly with existing code:
#' 
#' 1. **No Breaking Changes**: Original functions still work
#' 2. **Easy Migration**: Add "_themed" suffix to function names
#' 3. **Backward Compatible**: Use compatibility wrappers if needed
#' 4. **Extensible**: Easy to add new themes or customize colors
#' 
#' @examples
#' \dontrun{
#' # Run complete example
#' example_themed_bootstrap_workflow()
#' 
#' # Preview themes
#' preview_bootstrap_themes()
#' 
#' # Apply theme
#' my_table |> apply_bootstrap_theme("results")
#' }
#' 
#' @name bootstrap_gt_themes
#' @docType package
NULL# ==============================================================================
# GT THEME SYSTEM FOR BOOTSTRAP ANALYSIS
# Integration with existing R/bootstrap_summaries_helpers.R
# ==============================================================================

# ==============================================================================
# STEP 1: Create R/gt_themes.R (NEW FILE)
# ==============================================================================

#' Bootstrap Analysis GT Theme System
#' 
#' Provides consistent styling across all bootstrap analysis tables.
#' Includes themes for results, diagnostics, timing, and subgroup tables.
#'
#' @name gt_themes
NULL

# ------------------------------------------------------------------------------
# Color Palette Definition
# ------------------------------------------------------------------------------

#' Get bootstrap analysis color palette
#' @keywords internal
get_bootstrap_colors <- function() {
  list(
    # Primary colors
    primary = "#2E86AB",
    secondary = "#A23B72",
    tertiary = "#F18F01",
    
    # Status colors
    success = "#28a745",
    success_light = "#d4edda",
    warning = "#ffc107",
    warning_light = "#fff3cd",
    danger = "#dc3545",
    danger_light = "#f8d7da",
    info = "#17a2b8",
    info_light = "#d1ecf1",
    
    # Neutral colors
    dark = "#333333",
    gray_dark = "#495057",
    gray = "#6c757d",
    gray_light = "#dee2e6",
    gray_lighter = "#e9ecef",
    gray_lightest = "#f8f9fa",
    white = "#ffffff",
    
    # Accent colors (for data visualization)
    accent1 = "#667eea",
    accent2 = "#764ba2",
    accent3 = "#f093fb",
    accent4 = "#4facfe"
  )
}

# ------------------------------------------------------------------------------
# Base Theme Function
# ------------------------------------------------------------------------------

#' Apply base bootstrap GT theme
#' 
#' Provides foundational styling that all bootstrap tables inherit.
#' Includes borders, fonts, padding, and hover effects.
#' 
#' @param gt_table A gt table object
#' @return Modified gt table object
#' @keywords internal
#' @importFrom gt tab_options tab_style cell_text cells_column_labels opt_css px
apply_bootstrap_base_theme <- function(gt_table) {
  
  colors <- get_bootstrap_colors()
  
  gt_table |>
    # Overall table styling
    gt::tab_options(
      # Top border
      table.border.top.style = "solid",
      table.border.top.width = gt::px(3),
      table.border.top.color = colors$primary,
      
      # Bottom border
      table.border.bottom.style = "solid",
      table.border.bottom.width = gt::px(2),
      table.border.bottom.color = colors$dark,
      
      # Heading borders
      heading.border.bottom.style = "solid",
      heading.border.bottom.width = gt::px(2),
      heading.border.bottom.color = colors$dark,
      heading.align = "left",
      heading.title.font.size = gt::px(16),
      heading.title.font.weight = "bold",
      heading.subtitle.font.size = gt::px(13),
      heading.subtitle.font.weight = "normal",
      
      # Column labels
      column_labels.border.top.style = "hidden",
      column_labels.border.bottom.style = "solid",
      column_labels.border.bottom.width = gt::px(2),
      column_labels.border.bottom.color = colors$dark,
      column_labels.background.color = colors$gray_lightest,
      column_labels.font.size = gt::px(12),
      column_labels.font.weight = "bold",
      
      # Data rows
      table.font.size = gt::px(13),
      data_row.padding = gt::px(6),
      
      # Footnotes and source notes
      source_notes.font.size = gt::px(11),
      source_notes.padding = gt::px(8),
      source_notes.border.top.style = "solid",
      source_notes.border.top.width = gt::px(1),
      source_notes.border.top.color = colors$gray_light,
      
      footnotes.font.size = gt::px(11),
      footnotes.padding = gt::px(8),
      footnotes.border.bottom.style = "solid",
      footnotes.border.bottom.width = gt::px(1),
      footnotes.border.bottom.color = colors$gray_light
    ) |>
    
    # Column label styling
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold", color = colors$dark)
      ),
      locations = gt::cells_column_labels()
    ) |>
    
    # Add hover effects and accessibility CSS
    gt::opt_css(
      css = sprintf("
      /* Hover effect for data rows */
      .gt_table tbody tr:hover {
        background-color: %s;
        transition: background-color 0.2s ease;
      }
      
      /* Source note styling */
      .gt_sourcenote {
        border-top: 1px solid %s;
        margin-top: 10px;
        padding-top: 10px;
        font-style: italic;
      }
      
      /* Footnote styling */
      .gt_footnote {
        border-bottom: 1px solid %s;
        margin-bottom: 8px;
        padding-bottom: 8px;
      }
      
      /* Keyboard focus for accessibility */
      .gt_table *:focus {
        outline: 2px solid %s;
        outline-offset: 2px;
      }
      
      /* Improve print appearance */
      @media print {
        .gt_table tbody tr:hover {
          background-color: transparent;
        }
      }
      ",
      colors$gray_lightest,
      colors$gray_light,
      colors$gray_light,
      colors$primary
      )
    )
}

# ------------------------------------------------------------------------------
# Specialized Theme Functions
# ------------------------------------------------------------------------------

#' Apply results table theme
#' 
#' Specialized theme for bootstrap results tables showing treatment effects.
#' Highlights bias-corrected estimates and adds interpretation aids.
#' 
#' @param gt_table A gt table object
#' @return Modified gt table object
#' @export
apply_results_theme <- function(gt_table) {
  
  colors <- get_bootstrap_colors()
  
  gt_table |>
    apply_bootstrap_base_theme() |>
    
    # Highlight first row (usually "All Patients" / ITT)
    gt::tab_style(
      style = list(
        gt::cell_fill(color = colors$info_light),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(rows = 1)
    ) |>
    
    # Additional CSS for results
    gt::opt_css(
      css = sprintf("
      /* Spanner styling */
      .gt_column_spanner {
        background-color: %s;
        font-weight: bold;
        border-bottom: 2px solid %s;
      }
      
      /* Make HR confidence intervals slightly smaller */
      .gt_table td:nth-child(n+4) {
        font-size: 12px;
      }
      ",
      colors$gray_lighter,
      colors$dark
      )
    )
}

#' Apply diagnostics table theme
#' 
#' Specialized theme for diagnostic tables showing bootstrap quality metrics.
#' Emphasizes success rates and quality indicators.
#' 
#' @param gt_table A gt table object
#' @return Modified gt table object
#' @export
apply_diagnostics_theme <- function(gt_table) {
  
  colors <- get_bootstrap_colors()
  
  gt_table |>
    apply_bootstrap_base_theme() |>
    
    # Additional options for diagnostics
    gt::tab_options(
      table.border.top.color = colors$secondary,  # Different top border
      row_group.background.color = colors$gray_lightest,
      row_group.font.weight = "bold"
    ) |>
    
    # Additional CSS for diagnostics
    gt::opt_css(
      css = sprintf("
      /* Category row styling */
      .gt_table tbody td:first-child {
        font-weight: 500;
      }
      
      /* Empty rows for spacing */
      .gt_table tbody tr:has(td[style*='']) {
        height: 8px;
      }
      ",
      colors$gray_lighter
      )
    )
}

#' Apply timing table theme
#' 
#' Specialized theme for timing analysis tables.
#' Emphasizes performance metrics and time allocations.
#' 
#' @param gt_table A gt table object
#' @return Modified gt table object
#' @export
apply_timing_theme <- function(gt_table) {
  
  colors <- get_bootstrap_colors()
  
  gt_table |>
    apply_bootstrap_base_theme() |>
    
    # Additional options for timing
    gt::tab_options(
      table.border.top.color = colors$tertiary,  # Different top border
      row_group.background.color = colors$warning_light
    ) |>
    
    # Additional CSS for timing
    gt::opt_css(
      css = sprintf("
      /* Highlight performance rating rows */
      .gt_table tbody tr:has(td:contains('Performance')) {
        font-weight: bold;
      }
      
      /* Time value formatting */
      .gt_table tbody td:nth-child(3) {
        font-family: 'Courier New', monospace;
      }
      "
      )
    )
}

#' Apply subgroup summary theme
#' 
#' Specialized theme for subgroup summary tables.
#' Highlights factor frequencies and agreements.
#' 
#' @param gt_table A gt table object
#' @return Modified gt table object
#' @export
apply_subgroup_theme <- function(gt_table) {
  
  colors <- get_bootstrap_colors()
  
  gt_table |>
    apply_bootstrap_base_theme() |>
    
    # Additional options for subgroup summaries
    gt::tab_options(
      row_group.background.color = colors$accent1,
      row_group.font.weight = "bold",
      row_group.text_transform = "uppercase"
    )
}

# ==============================================================================
# STEP 2: Modify R/bootstrap_summaries_helpers.R
# ==============================================================================

# Update format_bootstrap_table() to use theme system
#' Format Bootstrap Results Table with Consistent Theme
#' @export
format_bootstrap_table_themed <- function(FSsg_tab, nb_boots, est.scale = "hr",
                                          boot_success_rate = NULL,
                                          title = NULL, subtitle = NULL) {
  
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required")
  }
  
  # Convert matrix to data.frame if needed
  if (is.matrix(FSsg_tab)) {
    FSsg_tab <- as.data.frame(FSsg_tab, stringsAsFactors = FALSE)
  }
  
  # Default title and subtitle
  if (is.null(title)) {
    title <- "Treatment Effect by Subgroup"
  }
  
  if (is.null(subtitle)) {
    effect_label <- ifelse(est.scale == "hr", "Hazard Ratio", "Inverse Hazard Ratio (1/HR)")
    subtitle <- sprintf("Bootstrap bias-corrected estimates (%d iterations)", nb_boots)
  }
  
  # Get column names
  col_names <- colnames(FSsg_tab)
  
  # Create labels
  labels_list <- list(Subgroup = "Subgroup")
  
  if ("n" %in% col_names) labels_list$n <- "N"
  if ("n1" %in% col_names) labels_list$n1 <- gt::md("N<sub>T</sub>")
  if ("events" %in% col_names) labels_list$events <- "Events"
  if ("m1" %in% col_names) labels_list$m1 <- gt::md("Med<sub>T</sub>")
  if ("m0" %in% col_names) labels_list$m0 <- gt::md("Med<sub>C</sub>")
  if ("RMST" %in% col_names) labels_list$RMST <- gt::md("RMST<sub>d</sub>")
  
  # Handle HR columns
  hr_col <- grep("HR.*CI", col_names, value = TRUE)[1]
  if (!is.na(hr_col) && length(hr_col) > 0) {
    labels_list[[hr_col]] <- gt::md("HR<br/>(95% CI)<sup>†</sup>")
  }
  
  hr_adj_col <- grep("HR\\*", col_names, value = TRUE)[1]
  if (!is.na(hr_adj_col) && length(hr_adj_col) > 0) {
    labels_list[[hr_adj_col]] <- gt::md("HR*<br/>(95% CI)<sup>‡</sup>")
  }
  
  # Create the gt table
  tbl <- FSsg_tab |>
    gt::gt() |>
    
    # Title and subtitle
    gt::tab_header(
      title = gt::md(paste0("**", title, "**")),
      subtitle = subtitle
    ) |>
    
    # Column labels
    gt::cols_label(.list = labels_list)
  
  # Add spanners
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
  
  # Add footnotes
  if (!is.na(hr_col) && length(hr_col) > 0) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md("**Unadjusted HR**: Standard Cox regression hazard ratio with robust standard errors"),
        locations = gt::cells_column_labels(columns = dplyr::all_of(hr_col))
      )
  }
  
  if (!is.na(hr_adj_col) && length(hr_adj_col) > 0) {
    tbl <- tbl |>
      gt::tab_footnote(
        footnote = gt::md(sprintf(
          "**Bias-corrected HR**: Bootstrap-adjusted estimate using infinitesimal jacknife method (%d iterations). Corrects for optimism in subgroup selection.",
          nb_boots
        )),
        locations = gt::cells_column_labels(columns = dplyr::all_of(hr_adj_col))
      )
  }
  
  # Source note
  source_note_text <- "*Note*: Med = Median survival time (months). RMST<sub>d</sub> = Restricted mean survival time difference."
  if (!is.null(boot_success_rate)) {
    source_note_text <- paste0(
      source_note_text,
      sprintf(" Subgroup identified in %.1f%% of bootstrap samples.", boot_success_rate * 100)
    )
  }
  
  tbl <- tbl |>
    gt::tab_source_note(source_note = gt::md(source_note_text))
  
  # Highlight bias-corrected column if it exists
  if (!is.na(hr_adj_col) && length(hr_adj_col) > 0) {
    colors <- get_bootstrap_colors()
    tbl <- tbl |>
      gt::tab_style(
        style = list(
          gt::cell_fill(color = colors$info_light),
          gt::cell_text(weight = "bold")
        ),
        locations = gt::cells_body(columns = dplyr::all_of(hr_adj_col))
      )
  }
  
  # Apply consistent theme
  tbl <- apply_results_theme(tbl)
  
  return(tbl)
}

# Update format_bootstrap_diagnostics_table() to use theme
#' Format Bootstrap Diagnostics Table with Consistent Theme
#' @export
format_bootstrap_diagnostics_table_themed <- function(diagnostics, nb_boots, results,
                                                      H_estimates = NULL, Hc_estimates = NULL) {
  
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required")
  }
  
  if (is.null(diagnostics)) {
    stop("diagnostics cannot be NULL")
  }
  
  colors <- get_bootstrap_colors()
  
  # Extract metrics
  success_rate <- diagnostics$success_rate
  n_successful <- diagnostics$n_successful
  n_failed <- diagnostics$n_failed
  
  # Determine success rating
  if (success_rate >= 0.90) {
    success_color <- colors$success_light
    success_rating <- "Excellent ✓✓✓"
  } else if (success_rate >= 0.75) {
    success_color <- colors$info_light
    success_rating <- "Good ✓✓"
  } else if (success_rate >= 0.50) {
    success_color <- colors$warning_light
    success_rating <- "Acceptable ✓"
  } else {
    success_color <- colors$danger_light
    success_rating <- "Poor ⚠"
  }
  
  # Build data frame
  diagnostics_rows <- list()
  
  # Success metrics
  diagnostics_rows$success <- data.frame(
    Category = c("Success Rate", "", "", ""),
    Metric = c(
      "Total iterations",
      "Successful subgroup ID",
      "Failed to find subgroup",
      "Success rating"
    ),
    Value = c(
      sprintf("%d", nb_boots),
      sprintf("%d (%.1f%%)", n_successful, success_rate * 100),
      sprintf("%d (%.1f%%)", n_failed, (1 - success_rate) * 100),
      success_rating
    ),
    stringsAsFactors = FALSE
  )
  
  # Bootstrap quality metrics
  if (!is.null(results)) {
    H_valid <- results$H_biasadj_2[!is.na(results$H_biasadj_2)]
    if (length(H_valid) > 0) {
      H_mean <- mean(H_valid)
      H_sd <- sd(H_valid)
      H_cv <- (H_sd / abs(H_mean)) * 100
      
      diagnostics_rows$quality <- data.frame(
        Category = c("", "Bootstrap Quality", "", ""),
        Metric = c(
          "",
          "Valid iterations",
          "Mean (SD)",
          "Coefficient of variation"
        ),
        Value = c(
          "",
          sprintf("%d", length(H_valid)),
          sprintf("%.2f (%.2f)", H_mean, H_sd),
          sprintf("%.1f%%", H_cv)
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine all rows
  diagnostics_df <- do.call(rbind, diagnostics_rows)
  rownames(diagnostics_df) <- NULL
  
  # Create gt table
  tbl <- diagnostics_df |>
    gt::gt() |>
    
    gt::tab_header(
      title = gt::md("**Bootstrap Diagnostics Summary**"),
      subtitle = sprintf("Analysis of %d bootstrap iterations", nb_boots)
    ) |>
    
    gt::cols_label(
      Category = "Category",
      Metric = "Metric",
      Value = "Value"
    ) |>
    
    # Style category rows
    gt::tab_style(
      style = list(
        gt::cell_fill(color = colors$gray_lightest),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = Category,
        rows = Category %in% c("Success Rate", "Bootstrap Quality")
      )
    ) |>
    
    # Highlight success rating
    gt::tab_style(
      style = list(
        gt::cell_fill(color = success_color),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = c(Metric, Value),
        rows = Metric == "Success rating"
      )
    )
  
  # Add interpretation
  interpretation <- "**Interpretation**: "
  if (success_rate >= 0.90) {
    interpretation <- paste0(interpretation, "Excellent stability - subgroup is consistently identified across bootstrap samples.")
  } else if (success_rate >= 0.75) {
    interpretation <- paste0(interpretation, "Good stability - subgroup is reliably identified in most bootstrap samples.")
  } else if (success_rate >= 0.50) {
    interpretation <- paste0(interpretation, "Moderate stability - consider increasing sample size or adjusting threshold.")
  } else {
    interpretation <- paste0(interpretation, "Poor stability - review subgroup criteria and consider if subgroup is spurious.")
  }
  
  tbl <- tbl |>
    gt::tab_source_note(source_note = gt::md(interpretation))
  
  # Apply consistent theme
  tbl <- apply_diagnostics_theme(tbl)
  
  return(tbl)
}

# Update format_bootstrap_timing_table() to use theme
#' Format Bootstrap Timing Table with Consistent Theme
#' @export
format_bootstrap_timing_table_themed <- function(timing_list, nb_boots, boot_success_rate) {
  
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required")
  }
  
  if (is.null(timing_list)) {
    stop("timing_list cannot be NULL")
  }
  
  colors <- get_bootstrap_colors()
  
  # Extract components
  overall <- timing_list$overall
  iteration_stats <- timing_list$iteration_stats
  fs_stats <- timing_list$fs_stats
  overhead_stats <- timing_list$overhead_stats
  
  # Build timing rows
  timing_rows <- list()
  
  # Overall timing
  if (!is.null(overall)) {
    timing_rows$overall <- data.frame(
      Category = c("Overall", "", "", ""),
      Metric = c(
        "Total time",
        "Average per boot",
        "Successful boots",
        "Projected for 1000"
      ),
      Value = c(
        sprintf("%.2f min (%.2f hrs)", overall$total_minutes, overall$total_hours),
        sprintf("%.2f min (%.1f sec)", overall$avg_minutes_per_boot, overall$avg_seconds_per_boot),
        sprintf("%d (%.1f%%)", round(boot_success_rate * nb_boots), boot_success_rate * 100),
        if (nb_boots < 1000) {
          projected <- overall$avg_minutes_per_boot * 1000
          sprintf("%.2f min (%.2f hrs)", projected, projected / 60)
        } else {
          "—"
        }
      ),
      stringsAsFactors = FALSE
    )
  }
  
  # Per-iteration stats
  if (!is.null(iteration_stats)) {
    timing_rows$iteration <- data.frame(
      Category = c("", "Per-Iteration", "", ""),
      Metric = c(
        "",
        "Mean",
        "Median",
        "Range"
      ),
      Value = c(
        "",
        sprintf("%.2f min (%.1f sec)", iteration_stats$mean, iteration_stats$mean * 60),
        sprintf("%.2f min (%.1f sec)", iteration_stats$median, iteration_stats$median * 60),
        sprintf("[%.2f, %.2f] min", iteration_stats$min, iteration_stats$max)
      ),
      stringsAsFactors = FALSE
    )
  }
  
  # ForestSearch timing
  if (!is.null(fs_stats)) {
    timing_rows$fs <- data.frame(
      Category = c("", "ForestSearch", "", ""),
      Metric = c(
        "",
        "Iterations with FS",
        "Mean FS time",
        "Total FS time"
      ),
      Value = c(
        "",
        sprintf("%d (%.1f%%)", fs_stats$n_runs, fs_stats$pct_runs),
        sprintf("%.2f min", fs_stats$mean),
        sprintf("%.2f min (%.1f%% of total)", fs_stats$total, 
                fs_stats$total / overall$total_minutes * 100)
      ),
      stringsAsFactors = FALSE
    )
  }
  
  # Performance rating
  if (!is.null(iteration_stats)) {
    avg_sec <- overall$avg_seconds_per_boot
    
    if (avg_sec < 5) {
      perf_rating <- "Excellent ✓✓✓"
      perf_color <- colors$success_light
    } else if (avg_sec < 15) {
      perf_rating <- "Good ✓✓"
      perf_color <- colors$info_light
    } else if (avg_sec < 30) {
      perf_rating <- "Acceptable ✓"
      perf_color <- colors$warning_light
    } else {
      perf_rating <- "Slow ⚠"
      perf_color <- colors$danger_light
    }
    
    timing_rows$performance <- data.frame(
      Category = c("", "Performance", ""),
      Metric = c(
        "",
        "Rating",
        "Speed"
      ),
      Value = c(
        "",
        perf_rating,
        sprintf("%.1f sec/iteration", avg_sec)
      ),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all rows
  timing_df <- do.call(rbind, timing_rows)
  rownames(timing_df) <- NULL
  
  # Create gt table
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
    ) |>
    
    # Style category rows
    gt::tab_style(
      style = list(
        gt::cell_fill(color = colors$gray_lightest),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = Category,
        rows = Category %in% c("Overall", "Per-Iteration", "ForestSearch", "Performance")
      )
    ) |>
    
    # Highlight performance rating
    gt::tab_style(
      style = list(
        gt::cell_fill(color = perf_color),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = c(Metric, Value),
        rows = Metric == "Rating"
      )
    )
  
  # Add recommendations if slow
  if (!is.null(iteration_stats) && avg_sec > 30) {
    recommendations <- "**Recommendations**: Consider reducing max.minutes, reducing maxk if > 2, or allocating more parallel workers."
  } else {
    recommendations <- "**Performance**: Good - no immediate optimization needed."
  }
  
  tbl <- tbl |>
    gt::tab_source_note(source_note = gt::md(recommendations))
  
  # Apply consistent theme
  tbl <- apply_timing_theme(tbl)
  
  return(tbl)
}

# ==============================================================================
# STEP 3: Update R/summarize_bootstrap_results.R
# ==============================================================================

#' Summarize Bootstrap Results with Themed Tables
#' 
#' Updated version that uses consistent theme system
#' 
#' @export
summarize_bootstrap_results_themed <- function(sgharm, boot_results, create_plots = FALSE,
                                               est.scale = "hr") {
  
  # Extract components (existing logic)
  FSsg_tab <- boot_results$FSsg_tab
  results <- boot_results$results
  H_estimates <- boot_results$H_estimates
  Hc_estimates <- boot_results$Hc_estimates
  
  # Calculate success rate
  boot_success_rate <- mean(!is.na(results$H_biasadj_2))
  nb_boots <- nrow(results)
  
  # Create themed tables
  formatted_table <- format_bootstrap_table_themed(
    FSsg_tab = FSsg_tab,
    nb_boots = nb_boots,
    est.scale = est.scale,
    boot_success_rate = boot_success_rate
  )
  
  # Diagnostics
  diagnostics <- list(
    n_boots = nb_boots,
    success_rate = boot_success_rate,
    n_successful = sum(!is.na(results$H_biasadj_2)),
    n_failed = sum(is.na(results$H_biasadj_2))
  )
  
  diagnostics_table <- format_bootstrap_diagnostics_table_themed(
    diagnostics = diagnostics,
    nb_boots = nb_boots,
    results = results,
    H_estimates = H_estimates,
    Hc_estimates = Hc_estimates
  )
  
  # Timing (if available)
  timing_table <- NULL
  if (!is.null(attr(results, "timing"))) {
    timing_list <- list(
      overall = attr(results, "timing"),
      iteration_stats = if ("tmins_iteration" %in% names(results)) {
        list(
          mean = mean(results$tmins_iteration, na.rm = TRUE),
          median = median(results$tmins_iteration, na.rm = TRUE),
          sd = sd(results$tmins_iteration, na.rm = TRUE),
          min = min(results$tmins_iteration, na.rm = TRUE),
          max = max(results$tmins_iteration, na.rm = TRUE)
        )
      } else NULL,
      fs_stats = if ("tmins_search" %in% names(results)) {
        fs_times <- results$tmins_search[!is.na(results$tmins_search)]
        list(
          n_runs = length(fs_times),
          pct_runs = length(fs_times) / nb_boots * 100,
          mean = mean(fs_times),
          median = median(fs_times),
          total = sum(fs_times)
        )
      } else NULL,