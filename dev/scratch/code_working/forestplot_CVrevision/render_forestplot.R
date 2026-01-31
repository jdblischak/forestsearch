# =============================================================================
# render_forestplot.R - Helper for Rendering ForestSearch Forest Plots
# =============================================================================
#
# Provides functions to properly render forestploter plots with theme control
# for sizing in Quarto/RMarkdown documents.
#
# =============================================================================

#' Create Forest Plot Theme with Size Controls
#'
#' Creates a forestploter theme with parameters that control overall plot
#' sizing and appearance. This is the primary way to control how large
#' the forest plot renders.
#'
#' @param base_size Numeric. Base font size in points. Increase for larger
#'   overall plot. Default: 10.
#' @param row_padding Numeric vector of length 2. Padding around row content
#'   in mm as c(vertical, horizontal). Increase vertical for taller rows.
#'   Default: c(4, 4).
#' @param ci_size Numeric. Size multiplier for CI points. Default: 0.4.
#' @param ci_lwd Numeric. Line width for CI lines. Default: 1.5.
#' @param ci_pch Integer. Point character for CI. 15=square, 16=circle,
#'   18=diamond. Default: 15.
#' @param ci_Theight Numeric. Height of T-bar ends on CI. Default: 0.2.
#' @param header_fontsize Numeric. Font size for column headers.
#'
#'   Default: base_size + 1.
#' @param body_fontsize Numeric. Font size for body text.
#'   Default: base_size.
#' @param footnote_fontsize Numeric. Font size for footnotes.
#'   Default: base_size - 1.
#' @param title_fontsize Numeric. Font size for title.
#'   Default: base_size + 4.
#' @param row_colors Character vector. Colors for row backgrounds.
#'   If NULL, uses default from sg_colors parameter in main function.
#'   Default: NULL.
#' @param refline_lwd Numeric. Reference line width. Default: 1.
#' @param refline_col Character. Reference line color. Default: "gray30".
#'
#' @return A forestploter theme object.
#'
#' @details
#' The most important parameters for controlling plot size are:
#' \itemize{
#'   \item \code{base_size}: Scales all fonts proportionally
#'   \item \code{row_padding}: Controls row height (first element is vertical)
#'   \item \code{ci_size}: Controls point size in CI column
#' }
#'
#' @examples
#' \dontrun{
#' # Create a larger theme
#' large_theme <- create_forest_theme(
#'   base_size = 14,
#'   row_padding = c(6, 4),
#'   ci_size = 0.5,
#'   ci_lwd = 2
#' )
#'
#' # Use with plot_subgroup_results_forestplot
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(fs.est = fs, fs_bc = fs_bc),
#'   df_analysis = df.analysis,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment",
#'   theme = large_theme
#' )
#' }
#'
#' @seealso \code{\link{plot_subgroup_results_forestplot}}
#'
#' @importFrom grid unit gpar
#' @export
create_forest_theme <- function(
    base_size = 10,
    row_padding = c(4, 4),
    ci_size = 0.4,
    ci_lwd = 1.5,
    ci_pch = 15,
    ci_Theight = 0.2,
    header_fontsize = NULL,
    body_fontsize = NULL,
    footnote_fontsize = NULL,
    title_fontsize = NULL,
    row_colors = NULL,
    refline_lwd = 1,
    refline_col = "gray30"
) {

  if (!requireNamespace("forestploter", quietly = TRUE)) {
    stop("Package 'forestploter' is required. Install with: install.packages('forestploter')")
  }

  # Calculate derived font sizes
  if (is.null(header_fontsize)) header_fontsize <- base_size + 1
  if (is.null(body_fontsize)) body_fontsize <- base_size
  if (is.null(footnote_fontsize)) footnote_fontsize <- base_size - 1
  if (is.null(title_fontsize)) title_fontsize <- base_size + 4

  # Default row colors if not provided
  if (is.null(row_colors)) {
    row_colors <- c("white", "gray95")
  }

  # Create theme
  tm <- forestploter::forest_theme(
    base_size = base_size,

    # Core/body cells
    core = list(
      fg_params = list(
        hjust = 0,
        x = 0.05,
        fontsize = body_fontsize
      ),
      bg_params = list(
        fill = row_colors
      ),
      padding = grid::unit(row_padding, "mm")
    ),

    # Column headers
    colhead = list(
      fg_params = list(
        hjust = 0.5,
        x = 0.5,
        fontsize = header_fontsize,
        fontface = "bold"
      ),
      bg_params = list(
        fill = "gray90"
      )
    ),

    # CI appearance
    ci_pch = ci_pch,
    ci_lwd = ci_lwd,
    ci_Theight = ci_Theight,

    # Reference line
    refline_lwd = refline_lwd,
    refline_lty = "dashed",
    refline_col = refline_col,

    # Footnote
    footnote_gp = grid::gpar(
      fontsize = footnote_fontsize,
      fontface = "italic",
      col = "gray40"
    ),

    # Title
    title_gp = grid::gpar(
      fontsize = title_fontsize,
      fontface = "bold"
    )
  )

  # Store parameters for reference
  attr(tm, "fs_params") <- list(
    base_size = base_size,
    row_padding = row_padding,
    ci_size = ci_size,
    ci_lwd = ci_lwd
  )

  return(tm)
}


#' Render ForestSearch Forest Plot
#'
#' Renders a forest plot. For size control, use \code{create_forest_theme()}
#' and pass to \code{plot_subgroup_results_forestplot(theme = ...)}.
#'
#' @param x An fs_forestplot object from \code{plot_subgroup_results_forestplot()}.
#' @param newpage Logical. Call grid.newpage() before drawing. Default: TRUE.
#'
#' @return Invisibly returns the grob object.
#'
#' @details
#' To control plot sizing, create a custom theme and pass it to
#' \code{plot_subgroup_results_forestplot()}:
#'
#' \code{my_theme <- create_forest_theme(base_size = 14, row_padding = c(6, 4))}
#'
#' \code{result <- plot_subgroup_results_forestplot(..., theme = my_theme)}
#'
#' @examples
#' \dontrun{
#' # For larger plot, use theme parameter in plot_subgroup_results_forestplot:
#' large_theme <- create_forest_theme(base_size = 14, row_padding = c(6, 4))
#'
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(fs.est = fs, fs_bc = fs_bc),
#'   df_analysis = df.analysis,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment",
#'   theme = large_theme
#' )
#'
#' render_forestplot(result)
#' }
#'
#' @importFrom grid grid.newpage grid.draw
#' @export
render_forestplot <- function(x, newpage = TRUE) {

 if (!inherits(x, "fs_forestplot")) {
    stop("x must be an fs_forestplot object from plot_subgroup_results_forestplot()")
  }

  if (newpage) grid::grid.newpage()
  grid::grid.draw(x$plot)

  invisible(x$plot)
}


#' Render Forest Plot with Preset Size
#'
#' Note: This function requires creating the plot with a custom theme.
#' Use \code{create_forest_theme()} with \code{plot_subgroup_results_forestplot()}.
#'
#' @param x An fs_forestplot object.
#' @param size Character. Ignored - use create_forest_theme() instead.
#'
#' @return Invisibly returns the grob.
#'
#' @examples
#' \dontrun{
#' # To get a larger plot, use theme when creating:
#' large_theme <- create_forest_theme(base_size = 14, row_padding = c(6, 4))
#' result <- plot_subgroup_results_forestplot(..., theme = large_theme)
#' render_forestplot(result)
#' }
#'
#' @export
render_forestplot_large <- function(x, size = c("large", "medium", "xlarge")) {
  message("Note: For size control, use create_forest_theme() with plot_subgroup_results_forestplot()")
  message("Example: theme <- create_forest_theme(base_size = 14, row_padding = c(6, 4))")
  message("         result <- plot_subgroup_results_forestplot(..., theme = theme)")
  render_forestplot(x)
}


#' Save ForestSearch Forest Plot to File
#'
#' Saves a forest plot to a file (PDF, PNG, etc.) with explicit dimensions.
#' For size control, use \code{create_forest_theme()} when creating the plot.
#'
#' @param x An fs_forestplot object.
#' @param filename Character. Output filename. Extension determines format.
#' @param width Numeric. Plot width in inches. Default: 12.
#' @param height Numeric. Plot height in inches. Default: 10.
#' @param dpi Numeric. Resolution for raster formats. Default: 300.
#' @param bg Character. Background color. Default: "white".
#'
#' @return Invisibly returns the filename.
#'
#' @examples
#' \dontrun{
#' # For larger plot, use theme when creating:
#' large_theme <- create_forest_theme(base_size = 14, row_padding = c(6, 4))
#'
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(fs.est = fs, fs_bc = fs_bc),
#'   df_analysis = df.analysis,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment",
#'   theme = large_theme
#' )
#'
#' # Save to file
#' save_forestplot(result, "forest_plot.pdf", width = 14, height = 12)
#' }
#'
#' @importFrom grDevices pdf png svg jpeg dev.off
#' @export
save_forestplot <- function(
    x,
    filename,
    width = 12,
    height = 10,
    dpi = 300,
    bg = "white"
) {

  if (!inherits(x, "fs_forestplot")) {
    stop("x must be an fs_forestplot object")
  }

  # Determine format from extension
  ext <- tolower(tools::file_ext(filename))

  # Open appropriate device
  if (ext == "pdf") {
    grDevices::pdf(filename, width = width, height = height, bg = bg)
  } else if (ext == "png") {
    grDevices::png(filename, width = width, height = height, units = "in",
                   res = dpi, bg = bg)
  } else if (ext == "svg") {
    if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(filename, width = width, height = height, bg = bg)
    } else {
      grDevices::svg(filename, width = width, height = height, bg = bg)
    }
  } else if (ext %in% c("jpg", "jpeg")) {
    grDevices::jpeg(filename, width = width, height = height, units = "in",
                    res = dpi, bg = bg, quality = 95)
  } else {
    stop("Unsupported file format: ", ext, ". Use .pdf, .png, .svg, or .jpg")
  }

  # Draw the plot
  grid::grid.newpage()
  grid::grid.draw(x$plot)

  # Close device
  grDevices::dev.off()

  message("Saved forest plot to: ", filename)
  invisible(filename)
}


#' Print Method for ForestSearch Forest Plot
#'
#' @param x An fs_forestplot object
#' @param ... Additional arguments (ignored)
#' @export
print.fs_forestplot <- function(x, ...) {
  cat("ForestSearch Subgroup Results Forest Plot\n")
  cat("=========================================\n")
  cat("Number of rows:", nrow(x$data), "\n")
  cat("Row types:", paste(unique(x$row_types), collapse = ", "), "\n")

  if (length(x$cv_metrics) > 0) {
    cat("\nCross-validation metrics:\n")
    for (i in seq_along(x$cv_metrics)) {
      cat("  -", x$cv_metrics[[i]], "\n")
    }
  }

  cat("\nTo display: render_forestplot(x) or plot(x)\n")
  cat("To save:    save_forestplot(x, 'plot.pdf', width = 14, height = 12)\n")
  cat("\nFor larger plot, recreate with custom theme:\n")
  cat("  theme <- create_forest_theme(base_size = 14, row_padding = c(6, 4))\n")
  cat("  result <- plot_subgroup_results_forestplot(..., theme = theme)\n")

  invisible(x)
}


#' Plot Method for ForestSearch Forest Plot
#'
#' @param x An fs_forestplot object
#' @param ... Additional arguments (ignored)
#' @export
plot.fs_forestplot <- function(x, ...) {
  render_forestplot(x)
  invisible(x)
}

