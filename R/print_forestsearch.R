#' Print Method for forestsearch Objects
#'
#' @param x A forestsearch object to print
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the input object
#' @export
print.forestsearch <- function(x, ...) {
  cat("ForestSearch Results\n")
  cat("====================\n")
  cat("Selection criterion (sg_focus):", x$parameters$sg_focus, "\n")
  cat("Subgroup identified:", if(!is.null(x$sg.harm)) "Yes" else "No", "\n")

  if (!is.null(x$sg.harm)) {
    cat("\nSelected Subgroup Characteristics:\n")
    cat("  Definition:", x$subgroup_definition, "\n")
    cat("  Sample size:", x$sg.harm$N, "\n")
    cat("  Hazard ratio:", sprintf("%.2f", x$sg.harm$hr), "\n")
    cat("  Consistency:", sprintf("%.1f%%", x$sg.harm$Pcons * 100), "\n")
    cat("  Control events:", x$sg.harm$d0, "\n")
    cat("  Treatment events:", x$sg.harm$d1, "\n")
  }

  cat("\nComputation time:", sprintf("%.1f", x$computation_time), "seconds\n")
  invisible(x)
}
