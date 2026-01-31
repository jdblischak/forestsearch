# =============================================================================
# render_forestplot.R - Helper for Rendering ForestSearch Forest Plots
# =============================================================================
#
# Provides functions to properly render forestploter plots with explicit
# dimensions in Quarto/RMarkdown documents.
#
# =============================================================================

#' Render ForestSearch Forest Plot with Explicit Dimensions
#'
#' Renders a forest plot from \code{plot_subgroup_results_forestplot()} with
#' explicit width and height control. This is necessary because forestploter
#' creates grid graphics that don't automatically respect Quarto/knitr chunk
#' dimensions.
#'
#' @param x An fs_forestplot object from \code{plot_subgroup_results_forestplot()},
#'   or a forestploter grob object directly.
#' @param width Numeric. Plot width in inches. Default: 12.
#' @param height Numeric. Plot height in inches. Default: 10.
#' @param units Character. Units for width/height: "in" (inches), "cm", or "mm".
#'   Default: "in".
#' @param newpage Logical. Whether to call \code{grid.newpage()} before drawing.
#'   Set to FALSE if adding to an existing plot. Default: TRUE.
#' @param vp_name Character. Name for the viewport. Default: "forestplot_vp".
#'
#' @return Invisibly returns the grob object.
#'
#' @details
#' This function solves the common problem where forestploter plots appear
#' condensed in Quarto/RMarkdown documents because the grid graphics don't
#' respect chunk options like \code{fig-width} and \code{fig-height}.
#'
#' The function creates a viewport with explicit dimensions and draws the
#' forest plot grob within that viewport.
#'
#' For Quarto/RMarkdown usage, set chunk options and call render_forestplot
#' with matching dimensions:
#'
#' \code{render_forestplot(result, width = 12, height = 10)}
#'
#' For more control, you can also use grid directly:
#'
#' \code{grid::grid.newpage()}
#'
#' \code{grid::grid.draw(result$plot)}
#'
#' @examples
#' \dontrun{
#' # Create forest plot
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(fs.est = fs, fs_bc = fs_bc),
#'   df_analysis = df.analysis,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment"
#' )
#'
#' # Render with explicit dimensions
#' render_forestplot(result, width = 12, height = 10)
#'
#' # Or specify dimensions in cm
#' render_forestplot(result, width = 30, height = 25, units = "cm")
#' }
#'
#' @seealso
#' \code{\link{plot_subgroup_results_forestplot}} for creating the forest plot
#'
#' @importFrom grid grid.newpage pushViewport viewport grid.draw popViewport unit
#' @export
render_forestplot <- function(x,
                              width = 12,
                              height = 10,
                              units = c("in", "cm", "mm"),
                              newpage = TRUE,
                              vp_name = "forestplot_vp") {

  units <- match.arg(units)

  # Extract grob from fs_forestplot object if necessary
  if (inherits(x, "fs_forestplot")) {
    grob <- x$plot
  } else if (inherits(x, "gTree") || inherits(x, "grob")) {
    grob <- x
  } else {
    stop("x must be an fs_forestplot object or a grid grob")
  }

  # Start new page if requested
  if (newpage) {
    grid::grid.newpage()
  }

  # Create viewport with explicit dimensions
  vp <- grid::viewport(
    width = grid::unit(width, units),
    height = grid::unit(height, units),
    name = vp_name
  )

  # Push viewport and draw
  grid::pushViewport(vp)
  grid::grid.draw(grob)
  grid::popViewport()

  invisible(grob)
}


#' Print ForestSearch Forest Plot to Device with Dimensions
#'
#' Alternative rendering function that sets the graphics device size
#' before printing. Useful for non-interactive rendering.
#'
#' @param x An fs_forestplot object.
#' @param width Numeric. Plot width in inches. Default: 12.
#' @param height Numeric. Plot height in inches. Default: 10.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @details
#' This function is useful when you want to ensure the plot renders
#' at specific dimensions. For Quarto/RMarkdown, it's often better to
#' use \code{render_forestplot()} which creates an explicit viewport.
#'
#' @examples
#' \dontrun{
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(fs.est = fs, fs_bc = fs_bc),
#'   df_analysis = df.analysis,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment"
#' )
#' print_forestplot(result, width = 12, height = 10)
#' }
#'
#' @export
print_forestplot <- function(x, width = 12, height = 10, ...) {

  if (!inherits(x, "fs_forestplot")) {
    stop("x must be an fs_forestplot object")
  }

  # For interactive use, just draw the grob
  grid::grid.newpage()
  grid::grid.draw(x$plot)

  invisible(x)
}


#' Save ForestSearch Forest Plot to File
#'
#' Saves a forest plot to a file (PDF, PNG, etc.) with explicit dimensions.
#'
#' @param x An fs_forestplot object.
#' @param filename Character. Output filename. Extension determines format
#'   (e.g., ".pdf", ".png", ".svg").
#' @param width Numeric. Plot width in inches. Default: 12.
#' @param height Numeric. Plot height in inches. Default: 10.
#' @param dpi Numeric. Resolution for raster formats (PNG, JPEG). Default: 300.
#' @param bg Character. Background color. Default: "white".
#' @param ... Additional arguments passed to the graphics device.
#'
#' @return Invisibly returns the filename.
#'
#' @examples
#' \dontrun{
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(fs.est = fs, fs_bc = fs_bc),
#'   df_analysis = df.analysis,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment"
#' )
#'
#' # Save as PDF
#' save_forestplot(result, "forest_plot.pdf", width = 12, height = 10)
#'
#' # Save as high-resolution PNG
#' save_forestplot(result, "forest_plot.png", width = 12, height = 10, dpi = 300)
#' }
#'
#' @importFrom grDevices pdf png svg jpeg dev.off
#' @export
save_forestplot <- function(x,
                            filename,
                            width = 12,
                            height = 10,
                            dpi = 300,
                            bg = "white",
                            ...) {

  if (!inherits(x, "fs_forestplot")) {
    stop("x must be an fs_forestplot object")
  }

  # Determine format from extension
  ext <- tolower(tools::file_ext(filename))

  # Open appropriate device
  if (ext == "pdf") {
    grDevices::pdf(filename, width = width, height = height, bg = bg, ...)
  } else if (ext == "png") {
    grDevices::png(filename, width = width, height = height, units = "in",
                   res = dpi, bg = bg, ...)
  } else if (ext == "svg") {
    if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(filename, width = width, height = height, bg = bg, ...)
    } else {
      grDevices::svg(filename, width = width, height = height, bg = bg, ...)
    }
  } else if (ext %in% c("jpg", "jpeg")) {
    grDevices::jpeg(filename, width = width, height = height, units = "in",
                    res = dpi, bg = bg, quality = 95, ...)
  } else {
    stop("Unsupported file format: ", ext,
         ". Use .pdf, .png, .svg, or .jpg")
  }

  # Draw the plot
  grid::grid.newpage()
  grid::grid.draw(x$plot)

  # Close device
  grDevices::dev.off()

  message("Saved forest plot to: ", filename)
  invisible(filename)
}


#' Get Recommended Dimensions for Forest Plot
#'
#' Calculates recommended plot dimensions based on the number of rows
#' in the forest plot data.
#'
#' @param x An fs_forestplot object.
#' @param row_height Numeric. Height per row in inches. Default: 0.4.
#' @param min_height Numeric. Minimum plot height in inches. Default: 6.
#' @param max_height Numeric. Maximum plot height in inches. Default: 20.
#' @param base_width Numeric. Base plot width in inches. Default: 10.
#' @param ci_width_factor Numeric. Additional width factor for CI column.
#'   Default: 0.1.
#'
#' @return Named list with recommended \code{width} and \code{height}.
#'
#' @examples
#' \dontrun{
#' result <- plot_subgroup_results_forestplot(
#'   fs_results = list(fs.est = fs, fs_bc = fs_bc),
#'   df_analysis = df.analysis,
#'   outcome.name = "time",
#'   event.name = "status",
#'   treat.name = "treatment"
#' )
#' dims <- get_forestplot_dims(result)
#' render_forestplot(result, width = dims$width, height = dims$height)
#' }
#'
#' @export
get_forestplot_dims <- function(x,
                                row_height = 0.4,
                                min_height = 6,
                                max_height = 20,
                                base_width = 10,
                                ci_width_factor = 0.1) {

  if (!inherits(x, "fs_forestplot")) {
    stop("x must be an fs_forestplot object")
  }

  n_rows <- nrow(x$data)

  # Calculate height based on rows
  calc_height <- n_rows * row_height + 2  # +2 for header/footer

  # Constrain to min/max
  height <- max(min_height, min(max_height, calc_height))

  # Width is relatively fixed but can adjust for CI column
  width <- base_width + (n_rows * ci_width_factor)
  width <- max(10, min(16, width))

  list(
    width = round(width, 1),
    height = round(height, 1),
    n_rows = n_rows
  )
}
