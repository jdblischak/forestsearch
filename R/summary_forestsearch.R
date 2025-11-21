#' Summary Method for forestsearch Objects
#'
#' @param object A forestsearch object to summarize
#' @param ... Additional arguments (currently unused)
#'
#' @return A list with summary information
#' @export
summary.forestsearch <- function(object, ...) {
  summary_list <- list(
    n_candidates = length(object$grp.consistency),
    selected_subgroup = object$sg.harm$definition,
    focus_type = object$parameters$sg_focus
  )
  class(summary_list) <- "summary.forestsearch"
  return(summary_list)
}
