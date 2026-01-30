# ==============================================================================
# GT THEME SYSTEM FOR BOOTSTRAP ANALYSIS
# File: R/gt_themes.R
# ==============================================================================
#
# This file provides a consistent theming system for all bootstrap analysis
# tables. It includes:
#   - Color palette management
#   - Base theme with common styling
#   - Specialized themes for different table types
#   - Theme application functions
#   - Custom theme support
#
# ==============================================================================

#' Bootstrap Analysis GT Theme System
#' 
#' @description
#' Provides consistent styling across all bootstrap analysis tables.
#' Includes themes for results, diagnostics, timing, and subgroup tables.
#'
#' @name gt_themes
#' @keywords internal
NULL

# ------------------------------------------------------------------------------
# COLOR PALETTE
# ------------------------------------------------------------------------------

#' Get Bootstrap Analysis Color Palette
#' 
#' Returns a named list of colors used throughout the theme system.
#' Modify this function to change colors globally across all themes.
#' 
#' @return Named list of color codes
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' colors <- get_bootstrap_colors()
#' colors$primary  # "#2E86AB"
#' }
get_bootstrap_colors <- function() {
  list(
    # Primary brand colors
    primary = "#2E86AB",      # Blue - main brand color
    secondary = "#A23B72",    # Purple - accent color
    tertiary = "#F18F01",     # Orange - highlight color
    
    # Status colors (semantic)
    success = "#28a745",           # Green
    success_light = "#d4edda",     # Light green background
    success_border = "#c3e6cb",    # Green border
    
    warning = "#ffc107",           # Yellow/Orange
    warning_light = "#fff3cd",     # Light yellow background
    warning_border = "#ffeaa7",    # Yellow border
    
    danger = "#dc3545",            # Red
    danger_light = "#f8d7da",      # Light red background
    danger_border = "#f5c6cb",     # Red border
    
    info = "#17a2b8",              # Cyan
    info_light = "#d1ecf1",        # Light cyan background
    info_border = "#bee5eb",       # Cyan border
    
    # Neutral colors
    dark = "#333333",              # Very dark gray (text)
    gray_dark = "#495057",         # Dark gray
    gray = "#6c757d",              # Medium gray
    gray_light = "#dee2e6",        # Light gray (borders)
    gray_lighter = "#e9ecef",      # Lighter gray (backgrounds)
    gray_lightest = "#f8f9fa",     # Lightest gray (subtle backgrounds)
    white = "#ffffff",             # White
    
    # Accent colors for data visualization
    accent1 = "#667eea",           # Purple-blue
    accent2 = "#764ba2",           # Deep purple
    accent3 = "#f093fb",           # Pink
    accent4 = "#4facfe",           # Sky blue
    accent5 = "#00f2fe",           # Cyan
    accent6 = "#43e97b"            # Green
  )
}

# ------------------------------------------------------------------------------
# BASE THEME
# ------------------------------------------------------------------------------

#' Apply Base Bootstrap GT Theme
#' 
#' Provides foundational styling that all bootstrap tables inherit.
#' Includes borders, fonts, padding, hover effects, and accessibility features.
#' 
#' @param gt_table A gt table object
#' @return Modified gt table object with base theme applied
#' @keywords internal
#' @importFrom gt tab_options tab_style cell_text cells_column_labels opt_css px
#' 
#' @examples
#' \dontrun{
#' library(gt)
#' my_table <- data.frame(x = 1:3, y = 4:6) |>
#'   gt() |>
#'   apply_bootstrap_base_theme()
#' }
apply_bootstrap_base_theme <- function(gt_table) {
  
  colors <- get_bootstrap_colors()
  
  gt_table |>
    # Overall table styling
    gt::tab_options(
      # Top border (distinctive)
      table.border.top.style = "solid",
      table.border.top.width = gt::px(3),
      table.border.top.color = colors$primary,
      
      # Bottom border
      table.border.bottom.style = "solid",
      table.border.bottom.width = gt::px(2),
      table.border.bottom.color = colors$dark,
      
      # Heading section
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
      table.font.color = colors$dark,
      data_row.padding = gt::px(6),
      
      # Footnotes
      footnotes.font.size = gt::px(11),
      footnotes.padding = gt::px(8),
      footnotes.border.bottom.style = "solid",
      footnotes.border.bottom.width = gt::px(1),
      footnotes.border.bottom.color = colors$gray_light,
      
      # Source notes
      source_notes.font.size = gt::px(11),
      source_notes.padding = gt::px(8),
      source_notes.border.top.style = "solid",
      source_notes.border.top.width = gt::px(1),
      source_notes.border.top.color = colors$gray_light
    ) |>
    
    # Column label styling
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold", color = colors$dark)
      ),
      locations = gt::cells_column_labels()
    ) |>
    
    # Add CSS for hover effects and accessibility
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
      
      /* Keyboard focus for accessibility (WCAG) */
      .gt_table *:focus {
        outline: 2px solid %s;
        outline-offset: 2px;
      }
      
      /* Improve print appearance */
      @media print {
        .gt_table tbody tr:hover {
          background-color: transparent;
        }
        .gt_table {
          page-break-inside: avoid;
        }
      }
      
      /* Spanner styling */
      .gt_column_spanner {
        background-color: %s;
        font-weight: bold;
        border-bottom: 2px solid %s;
        padding: 8px;
      }
      ",
      colors$gray_lightest,     # Hover color
      colors$gray_light,        # Source note border
      colors$gray_light,        # Footnote border
      colors$primary,           # Focus outline
      colors$gray_lighter,      # Spanner background
      colors$dark               # Spanner border
      )
    )
}

# ------------------------------------------------------------------------------
# SPECIALIZED THEMES
# ------------------------------------------------------------------------------

#' Apply Results Table Theme
#' 
#' Specialized theme for bootstrap results tables showing treatment effects.
#' Highlights bias-corrected estimates and adds interpretation aids.
#' 
#' @param gt_table A gt table object
#' @return Modified gt table object with results theme applied
#' @export
#' 
#' @examples
#' \dontrun{
#' results_table <- format_bootstrap_table(...) |>
#'   apply_results_theme()
#' }
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
    
    # Additional CSS for results-specific styling
    gt::opt_css(
      css = sprintf("
      /* Spanner styling for results tables */
      .gt_column_spanner {
        background-color: %s;
        font-weight: bold;
        border-bottom: 2px solid %s;
      }
      
      /* Make confidence intervals slightly smaller font */
      .gt_table td:nth-child(n+4) {
        font-size: 12px;
      }
      
      /* Highlight bias-corrected columns */
      .gt_table td[style*='background-color: %s'] {
        border-left: 3px solid %s;
      }
      ",
      colors$gray_lighter,   # Spanner background
      colors$dark,           # Spanner border
      colors$info_light,     # Bias-corrected background
      colors$info            # Bias-corrected border
      )
    )
}

#' Apply Diagnostics Table Theme
#' 
#' Specialized theme for diagnostic tables showing bootstrap quality metrics.
#' Emphasizes success rates and quality indicators.
#' 
#' @param gt_table A gt table object
#' @return Modified gt table object with diagnostics theme applied
#' @export
#' 
#' @examples
#' \dontrun{
#' diagnostics_table <- format_bootstrap_diagnostics_table(...) |>
#'   apply_diagnostics_theme()
#' }
apply_diagnostics_theme <- function(gt_table) {
  
  colors <- get_bootstrap_colors()
  
  gt_table |>
    apply_bootstrap_base_theme() |>
    
    # Different top border color for diagnostics
    gt::tab_options(
      table.border.top.color = colors$secondary,
      row_group.background.color = colors$gray_lightest,
      row_group.font.weight = "bold",
      row_group.border.top.style = "solid",
      row_group.border.top.width = gt::px(2),
      row_group.border.top.color = colors$gray_light,
      row_group.border.bottom.style = "solid",
      row_group.border.bottom.width = gt::px(1),
      row_group.border.bottom.color = colors$gray_light
    ) |>
    
    # Additional CSS for diagnostics-specific styling
    gt::opt_css(
      css = sprintf("
      /* Category row styling */
      .gt_table tbody td:first-child {
        font-weight: 500;
      }
      
      /* Empty rows for spacing */
      .gt_table tbody tr:has(td:empty) {
        height: 8px;
        background-color: transparent;
      }
      
      /* Highlight metric rows */
      .gt_table tbody tr td:nth-child(2) {
        padding-left: 20px;
      }
      "
      )
    )
}

#' Apply Timing Table Theme
#' 
#' Specialized theme for timing analysis tables.
#' Emphasizes performance metrics and time allocations.
#' 
#' @param gt_table A gt table object
#' @return Modified gt table object with timing theme applied
#' @export
#' 
#' @examples
#' \dontrun{
#' timing_table <- format_bootstrap_timing_table(...) |>
#'   apply_timing_theme()
#' }
apply_timing_theme <- function(gt_table) {
  
  colors <- get_bootstrap_colors()
  
  gt_table |>
    apply_bootstrap_base_theme() |>
    
    # Different top border color for timing
    gt::tab_options(
      table.border.top.color = colors$tertiary,
      row_group.background.color = colors$warning_light,
      row_group.font.weight = "bold"
    ) |>
    
    # Additional CSS for timing-specific styling
    gt::opt_css(
      css = sprintf("
      /* Monospace font for time values */
      .gt_table tbody td:nth-child(3) {
        font-family: 'Courier New', monospace;
        font-size: 12px;
      }
      
      /* Highlight performance rows */
      .gt_table tbody tr:has(td:contains('Performance')) {
        font-weight: bold;
      }
      "
      )
    )
}

#' Apply Subgroup Summary Theme
#' 
#' Specialized theme for subgroup summary tables.
#' Highlights factor frequencies and agreements.
#' 
#' @param gt_table A gt table object
#' @return Modified gt table object with subgroup theme applied
#' @export
#' 
#' @examples
#' \dontrun{
#' subgroup_table <- format_subgroup_summary_tables(...) |>
#'   apply_subgroup_theme()
#' }
apply_subgroup_theme <- function(gt_table) {
  
  colors <- get_bootstrap_colors()
  
  gt_table |>
    apply_bootstrap_base_theme() |>
    
    # Styling for subgroup summaries
    gt::tab_options(
      row_group.background.color = colors$accent1,
      row_group.font.weight = "bold",
      row_group.text_transform = "uppercase",
      row_group.font.size = gt::px(11)
    )
}

# ------------------------------------------------------------------------------
# THEME APPLICATION FUNCTION
# ------------------------------------------------------------------------------

#' Apply Bootstrap Theme to GT Table
#' 
#' Main function for applying themes. Allows users to choose from predefined
#' themes or use custom colors.
#' 
#' @param gt_table A gt table object
#' @param theme Character. One of "base", "results", "diagnostics", "timing", 
#'   "subgroup", or "custom"
#' @param custom_colors Named list of colors (optional, for theme="custom")
#' 
#' @return Modified gt table object with theme applied
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
apply_bootstrap_theme <- function(gt_table, 
                                  theme = c("results", "diagnostics", "timing", 
                                            "subgroup", "base", "custom"),
                                  custom_colors = NULL) {
  
  theme <- match.arg(theme)
  
  if (theme == "custom") {
    if (is.null(custom_colors)) {
      warning("No custom_colors provided. Using base theme instead.")
      return(apply_bootstrap_base_theme(gt_table))
    }
    return(apply_custom_theme(gt_table, custom_colors))
  }
  
  switch(theme,
    base = apply_bootstrap_base_theme(gt_table),
    results = apply_results_theme(gt_table),
    diagnostics = apply_diagnostics_theme(gt_table),
    timing = apply_timing_theme(gt_table),
    subgroup = apply_subgroup_theme(gt_table),
    apply_bootstrap_base_theme(gt_table)  # Default fallback
  )
}

# ------------------------------------------------------------------------------
# CUSTOM THEME SUPPORT
# ------------------------------------------------------------------------------

#' Apply Custom Theme with User-Specified Colors
#' 
#' Allows users to override default colors while maintaining base styling.
#' 
#' @param gt_table A gt table object
#' @param custom_colors Named list of color overrides
#' @return Modified gt table object
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' custom_colors <- list(
#'   primary = "#003366",
#'   success = "#00AA00",
#'   warning = "#FFB000"
#' )
#' my_table |> apply_custom_theme(custom_colors)
#' }
apply_custom_theme <- function(gt_table, custom_colors) {
  
  # Start with base theme
  tbl <- apply_bootstrap_base_theme(gt_table)
  
  # Override specific colors if provided
  if (!is.null(custom_colors$primary)) {
    tbl <- tbl |>
      gt::tab_options(
        table.border.top.color = custom_colors$primary
      )
  }
  
  if (!is.null(custom_colors$success)) {
    # Would apply success color to relevant elements
    # Implementation depends on specific table structure
  }
  
  # Add custom CSS if needed
  if (!is.null(custom_colors)) {
    css_overrides <- generate_custom_css(custom_colors)
    tbl <- tbl |> gt::opt_css(css = css_overrides)
  }
  
  return(tbl)
}

#' Generate Custom CSS from Color Overrides
#' @keywords internal
generate_custom_css <- function(custom_colors) {
  
  css <- ""
  
  if (!is.null(custom_colors$primary)) {
    css <- paste0(css, sprintf("
      .gt_table {
        border-top-color: %s !important;
      }
      .gt_table *:focus {
        outline-color: %s !important;
      }
    ", custom_colors$primary, custom_colors$primary))
  }
  
  if (!is.null(custom_colors$hover)) {
    css <- paste0(css, sprintf("
      .gt_table tbody tr:hover {
        background-color: %s !important;
      }
    ", custom_colors$hover))
  }
  
  return(css)
}

# ------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ------------------------------------------------------------------------------

#' Preview All Bootstrap Themes
#' 
#' Creates sample tables with all theme variations for comparison.
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
#' theme_examples$timing
#' }
preview_bootstrap_themes <- function() {
  
  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required")
  }
  
  # Create sample data
  sample_data <- data.frame(
    Category = c("Metric 1", "Metric 2", "Metric 3"),
    Value1 = c("100", "85.5%", "Good ✓"),
    Value2 = c("200", "92.3%", "Excellent ✓✓"),
    stringsAsFactors = FALSE
  )
  
  # Create base table
  base_table <- sample_data |>
    gt::gt() |>
    gt::tab_header(
      title = gt::md("**Theme Preview**"),
      subtitle = "Sample bootstrap analysis table"
    ) |>
    gt::cols_label(
      Category = "Category",
      Value1 = "Value 1",
      Value2 = "Value 2"
    )
  
  # Apply each theme and return as list
  list(
    base = base_table |> apply_bootstrap_base_theme(),
    results = base_table |> apply_results_theme(),
    diagnostics = base_table |> apply_diagnostics_theme(),
    timing = base_table |> apply_timing_theme(),
    subgroup = base_table |> apply_subgroup_theme()
  )
}

#' Get Current Theme Colors
#' 
#' Returns the color palette being used by the theme system.
#' Useful for creating consistent custom visualizations.
#' 
#' @return Named list of colors
#' @export
#' 
#' @examples
#' \dontrun{
#' colors <- get_theme_colors()
#' ggplot(...) + scale_color_manual(values = c(colors$primary, colors$secondary))
#' }
get_theme_colors <- function() {
  get_bootstrap_colors()
}

# ------------------------------------------------------------------------------
# PACKAGE DOCUMENTATION
# ------------------------------------------------------------------------------

#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
