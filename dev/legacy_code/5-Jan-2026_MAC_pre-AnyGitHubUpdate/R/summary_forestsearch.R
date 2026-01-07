#' Summary Method for forestsearch Objects
#'
#' @param object A forestsearch object to summarize
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the input object
#' @export
summary.forestsearch <- function(object, ...) {
  cat("ForestSearch Summary\n")
  cat("====================\n\n")

  # Parameters section
  cat("Analysis Parameters:\n")
  cat("  sg_focus:", object$parameters$sg_focus, "\n")
  cat("  hr.threshold:", object$parameters$hr.threshold, "\n")
  cat("  hr.consistency:", object$parameters$hr.consistency, "\n")
  cat("  pconsistency.threshold:", object$parameters$pconsistency.threshold, "\n")
  cat("  n.min:", object$parameters$n.min, "\n")
  cat("  n.splits:", object$parameters$n.splits, "\n")
  cat("  maxk:", object$parameters$maxk, "\n\n")

  # Variable selection
  cat("Variable Selection:\n")
  cat("  GRF variables:", length(object$grf_importance), "\n")
  cat("  LASSO selected:", length(object$lasso_selected), "\n\n")

  # Consistency results
  cat("Consistency Evaluation:\n")
  cat("  Candidates evaluated:", nrow(object$grp.consistency$all_candidates), "\n")
  cat("  Meeting hr.threshold:",
      sum(object$grp.consistency$all_candidates$hr > object$parameters$hr.threshold), "\n")
  cat("  Meeting consistency:",
      sum(object$grp.consistency$all_candidates$Pcons >= object$parameters$pconsistency.threshold), "\n\n")

  # Results by sg_focus
  cat("Results by sg_focus Option:\n")
  for (focus in c("hr", "maxSG", "minSG")) {
    result_name <- paste0("out_", focus)
    if (!is.null(object$grp.consistency[[result_name]])) {
      sg <- object$grp.consistency[[result_name]]$sg.harm
      if (!is.null(sg)) {
        cat(sprintf("  %-8s: N=%d, HR=%.2f, Pcons=%.1f%%\n",
                    focus, sg$N, sg$hr, sg$Pcons * 100))
      }
    }
  }

  invisible(object)
}
