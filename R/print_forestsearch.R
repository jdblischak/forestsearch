#' Print Method for forestsearch Objects
#'
#' @param x A forestsearch object to print
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the input object
#' @export
print.forestsearch <- function(x, ...) {
  cat("Forest Search Results\n")
  cat("====================\n")
  if (!is.null(x$sg.harm)) {
    cat("Selected subgroup:", x$sg.harm$definition, "\n")
    cat("Sample size:", x$sg.harm$n, "\n")
    cat("Hazard ratio:", round(x$sg.harm$hr, 3), "\n")
  }
  cat("\nUse str() to see full object structure\n")
  invisible(x)
}
