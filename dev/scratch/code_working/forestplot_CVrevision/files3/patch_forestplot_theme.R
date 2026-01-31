# =============================================================================
# PATCH for plot_subgroup_results_forestplot.R
# =============================================================================
# Apply these changes to add theme support with proper coloring and CV font size
# =============================================================================

# =============================================================================
# CHANGE 1: Add roxygen documentation for theme parameter
# =============================================================================
# FIND the @param conf.level line and ADD AFTER IT:

#' @param theme An fs_forest_theme object from \code{create_forest_theme()}.
#'   Use this to control plot sizing (fonts, row height, CI appearance,
#'   CV annotation font size). Default: NULL (uses default theme).


# =============================================================================
# CHANGE 2: Add theme = NULL to function signature
# =============================================================================
# FIND:
#     ci_column_spaces = 20,
#     conf.level = 0.95
# ) {

# REPLACE WITH:
#     ci_column_spaces = 20,
#     conf.level = 0.95,
#     theme = NULL
# ) {


# =============================================================================
# CHANGE 3: Replace the ENTIRE theme and forest plot creation section
# =============================================================================
# FIND the section that starts with:
#   # ==========================================================================
#   # Create Forest Plot
#   # ==========================================================================
#
#   # Apply theme
#   tm <- forestploter::forest_theme(
#     ...
#   )
#
# And ENDS just before:
#   # ==========================================================================
#   # Return Results
#   # ==========================================================================
#
# REPLACE THE ENTIRE SECTION WITH THE CODE BELOW:

  # ==========================================================================
  # Create Forest Plot
  # ==========================================================================

  # Extract theme parameters or use defaults
  if (!is.null(theme) && inherits(theme, "fs_forest_theme")) {
    # Use custom theme parameters
    base_size <- theme$base_size
    row_padding <- theme$row_padding
    body_fontsize <- theme$body_fontsize
    header_fontsize <- theme$header_fontsize
    footnote_fontsize <- theme$footnote_fontsize
    footnote_col <- theme$footnote_col
    cv_fontsize <- theme$cv_fontsize
    cv_col <- theme$cv_col
    ci_pch <- theme$ci_pch
    ci_lwd <- theme$ci_lwd
    ci_Theight <- theme$ci_Theight
    ci_col <- theme$ci_col
    refline_lwd <- theme$refline_lwd
    refline_lty <- theme$refline_lty
    refline_col <- theme$refline_col
  } else {
    # Default values
    base_size <- 10
    row_padding <- c(4, 4)
    body_fontsize <- 10
    header_fontsize <- 11
    footnote_fontsize <- 9
    footnote_col <- "darkcyan"
    cv_fontsize <- 7
    cv_col <- "gray30"
    ci_pch <- 15
    ci_lwd <- 1.5
    ci_Theight <- 0.2
    ci_col <- "black"
    refline_lwd <- 1
    refline_lty <- "dashed"
    refline_col <- "gray30"
  }

  # Create forestploter theme with row-specific colors
  tm <- forestploter::forest_theme(
    base_size = base_size,
    core = list(
      fg_params = list(
        hjust = 1,
        x = 0.9,
        fontsize = body_fontsize
      ),
      bg_params = list(
        fill = sg_colors
      ),
      padding = grid::unit(row_padding, "mm")
    ),
    colhead = list(
      fg_params = list(
        hjust = 0.5,
        x = 0.5,
        fontsize = header_fontsize,
        fontface = "bold"
      )
    ),
    ci_pch = ci_pch,
    ci_col = ci_col,
    ci_lwd = ci_lwd,
    ci_Theight = ci_Theight,
    refline_lwd = refline_lwd,
    refline_lty = refline_lty,
    refline_col = refline_col,
    footnote_gp = grid::gpar(
      fontsize = footnote_fontsize,
      fontface = "italic",
      col = footnote_col
    )
  )

  # Add spacing column for CI plot (width controlled by ci_column_spaces)
  dt$` ` <- paste(rep(" ", ci_column_spaces), collapse = " ")

  # Create HR (xx% CI) text column with dynamic label
  dt[[ci_label]] <- ifelse(is.na(dt$se), "",
                           sprintf("%.2f (%.2f to %.2f)", dt$est, dt$low, dt$hi))

  # Generate the forest plot
  p <- forestploter::forest(
    dt[, c("Subgroup", E.name, C.name, " ", ci_label)],
    title = title_text,
    est = dt$est,
    lower = dt$low,
    upper = dt$hi,
    sizes = 0.4,
    ci_column = 4,
    ref_line = 1,
    arrow_lab = arrow_text,
    xlim = xlim,
    ticks_at = ticks_at,
    xlog = xlog,
    footnote = footnote_text,
    theme = tm
  )

  # Add CV annotation text using insert_text (spans across first 3 columns)
  if (length(cv_row_positions) > 0) {
    for (i in seq_along(cv_texts)) {
      p <- forestploter::insert_text(
        p,
        text = paste0("  ", cv_texts[[i]]),
        row = cv_row_positions[[i]],
        col = 1:3,
        part = "body",
        just = "center",
        gp = grid::gpar(fontsize = cv_fontsize, fontface = "italic", col = cv_col)
      )
    }
  }


# =============================================================================
# END OF CHANGES
# =============================================================================
# After making these changes:
# 1. Run devtools::document()
# 2. Run devtools::load_all()
#
# Then use like this:
#
# large_theme <- create_forest_theme(
#   base_size = 14,
#   row_padding = c(6, 4),
#   cv_fontsize = 11,  # Larger CV text
#   ci_lwd = 2
# )
#
# result <- plot_subgroup_results_forestplot(
#   fs_results = list(fs.est = fs, fs_bc = fs_bc, fs_kfold = fs_kfold),
#   df_analysis = df.analysis,
#   subgroup_list = subgroups,
#   outcome.name = outcome.name,
#   event.name = event.name,
#   treat.name = treat.name,
#   theme = large_theme
# )
#
# render_forestplot(result)
# =============================================================================
