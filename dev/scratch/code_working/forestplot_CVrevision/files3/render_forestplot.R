# =============================================================================
# render_forestplot.R - Helper for Rendering ForestSearch Forest Plots
# =============================================================================
#
# Provides functions to create custom themes and render forestploter plots
# with size control in Quarto/RMarkdown documents.
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
#' @param ci_pch Integer. Point character for CI. 15=square, 16=circle,
#'   18=diamond. Default: 15.
#' @param ci_lwd Numeric. Line width for CI lines. Default: 1.5.
#' @param ci_Theight Numeric. Height of T-bar ends on CI. Default: 0.2.
#' @param ci_col Character. Color for CI lines and points. Default: "black".
#' @param header_fontsize Numeric. Font size for column headers.
#'   Default: base_size + 1.
#' @param body_fontsize Numeric. Font size for body text.
#'   Default: base_size.
#' @param footnote_fontsize Numeric. Font size for footnotes.
#'   Default: base_size - 1.
#' @param footnote_col Character. Color for footnote text. Default: "darkcyan".
#' @param title_fontsize Numeric. Font size for title.
#'   Default: base_size + 4.
#' @param cv_fontsize Numeric. Font size for CV annotation text.
#'   Default: 9. Set higher (e.g., 10-12) for larger CV text.
#' @param cv_col Character. Color for CV annotation text. Default: "gray30".
#' @param refline_lwd Numeric. Reference line width. Default: 1.
#' @param refline_lty Character. Reference line type. Default: "dashed".
#' @param refline_col Character. Reference line color. Default: "gray30".
#' @param vertline_lwd Numeric. Vertical line width. Default: 1.
#' @param vertline_lty Character. Vertical line type. Default: "dashed".
#' @param vertline_col Character. Vertical line color. Default: "gray20".
#' @param arrow_type Character. Arrow type: "open" or "closed". Default: "closed".
#' @param arrow_col Character. Arrow color. Default: "black".
#' @param summary_fill Character. Fill color for summary diamonds. Default: "black".
#' @param summary_col Character. Border color for summary diamonds. Default: "black".
#'
#' @return A list containing a forestploter theme object and additional parameters.
#'   The list has class "fs_forest_theme" and includes:
#'   \describe{
#'     \item{theme}{The forestploter theme object}
#'     \item{cv_fontsize}{Font size for CV annotations}
#'     \item{cv_col}{Color for CV annotations}
#'     \item{params}{All input parameters for reference}
#'   }
#'
#' @details
#' The most important parameters for controlling plot size are:
#' \itemize{
#'   \item \code{base_size}: Scales all fonts proportionally
#'   \item \code{row_padding}: Controls row height (first element is vertical)
#'   \item \code{cv_fontsize}: Controls the CV annotation text size
#' }
#'
#' The theme does NOT set row background colors - those are determined
#' automatically by \code{plot_subgroup_results_forestplot()} based on
#' row types (ITT, reference, posthoc, etc.).
#'
#' @examples
#' \dontrun{
#' # Create a larger theme with bigger CV text
#' large_theme <- create_forest_theme(
#'   base_size = 14,
#'   row_padding = c(6, 4),
#'   cv_fontsize = 11,
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
#'
#' render_forestplot(result)
#' }
#'
#' @seealso \code{\link{plot_subgroup_results_forestplot}}, \code{\link{render_forestplot}}
#'
#' @importFrom grid unit gpar
#' @export
create_forest_theme <- function(
    base_size = 10,
    row_padding = c(4, 4),
    ci_pch = 15,
    ci_lwd = 1.5,
    ci_Theight = 0.2,
    ci_col = "black",
    header_fontsize = NULL,
    body_fontsize = NULL,
    footnote_fontsize = NULL,
    footnote_col = "darkcyan",
    title_fontsize = NULL,
    cv_fontsize = 9,
    cv_col = "gray30",
    refline_lwd = 1,
    refline_lty = "dashed",
    refline_col = "gray30",
    vertline_lwd = 1,
    vertline_lty = "dashed",
    vertline_col = "gray20",
    arrow_type = "closed",
    arrow_col = "black",
    summary_fill = "black",
    summary_col = "black"
) {

  if (!requireNamespace("forestploter", quietly = TRUE)) {
    stop("Package 'forestploter' is required. Install with: install.packages('forestploter')")
  }

  # Calculate derived font sizes if not provided
  if (is.null(header_fontsize)) header_fontsize <- base_size + 1
  if (is.null(body_fontsize)) body_fontsize <- base_size
  if (is.null(footnote_fontsize)) footnote_fontsize <- max(base_size - 1, 8)
  if (is.null(title_fontsize)) title_fontsize <- base_size + 4

  # Create result object with parameters
  # The actual forestploter theme will be created in plot_subgroup_results_forestplot()
 # where we have access to sg_colors
  result <- list(
    base_size = base_size,
    row_padding = row_padding,
    body_fontsize = body_fontsize,
    header_fontsize = header_fontsize,
    footnote_fontsize = footnote_fontsize,
    footnote_col = footnote_col,
    title_fontsize = title_fontsize,
    cv_fontsize = cv_fontsize,
    cv_col = cv_col,
    ci_pch = ci_pch,
    ci_lwd = ci_lwd,
    ci_Theight = ci_Theight,
    ci_col = ci_col,
    refline_lwd = refline_lwd,
    refline_lty = refline_lty,
    refline_col = refline_col,
    vertline_lwd = vertline_lwd,
    vertline_lty = vertline_lty,
    vertline_col = vertline_col,
    arrow_type = arrow_type,
    arrow_col = arrow_col,
    summary_fill = summary_fill,
    summary_col = summary_col
  )

  class(result) <- c("fs_forest_theme", "list")

  return(result)
}


#' Print Method for ForestSearch Forest Theme
#'
#' @param x An fs_forest_theme object
#' @param ... Additional arguments (ignored)
#' @export
print.fs_forest_theme <- function(x, ...) {
  cat("ForestSearch Forest Plot Theme\n")
  cat("==============================\n")
  cat("Base size:       ", x$base_size, "\n")
  cat("Row padding:     ", paste(x$row_padding, collapse = ", "), " mm\n")
  cat("Body font size:  ", x$body_fontsize, "\n")
  cat("Header font size:", x$header_fontsize, "\n")
  cat("CV font size:    ", x$cv_fontsize, "\n")
  cat("CI line width:   ", x$ci_lwd, "\n")
  cat("\nUse with: plot_subgroup_results_forestplot(..., theme = x)\n")
  invisible(x)
}


#' Render ForestSearch Forest Plot
#'
#' Renders a forest plot from \code{plot_subgroup_results_forestplot()}.
#'
#' @param x An fs_forestplot object from \code{plot_subgroup_results_forestplot()}.
#' @param newpage Logical. Call grid.newpage() before drawing. Default: TRUE.
#'
#' @return Invisibly returns the grob object.
#'
#' @details
#' To control plot sizing, create a custom theme using \code{create_forest_theme()}
#' and pass it to \code{plot_subgroup_results_forestplot()}:
#'
#' \code{my_theme <- create_forest_theme(base_size = 14, row_padding = c(6, 4))}
#'
#' \code{result <- plot_subgroup_results_forestplot(..., theme = my_theme)}
#'
#' \code{render_forestplot(result)}
#'
#' @examples
#' \dontrun{
#' # For larger plot, use theme parameter in plot_subgroup_results_forestplot:
#' large_theme <- create_forest_theme(
#'   base_size = 14,
#'   row_padding = c(6, 4),
#'   cv_fontsize = 11
#' )
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


#' Save ForestSearch Forest Plot to File
#'
#' Saves a forest plot to a file (PDF, PNG, etc.) with explicit dimensions.
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
#' # Create plot with custom theme
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
  cat("\nFor larger plot or bigger CV text, recreate with custom theme:\n")
  cat("  theme <- create_forest_theme(base_size = 14, cv_fontsize = 11)\n")
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
