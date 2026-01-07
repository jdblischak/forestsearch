# List of required packages for ForestSearch analysis

required_packages <- c("grf","policytree","data.table","randomForest","survival","weightedsurv","future.apply")
missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if(length(missing) > 0) stop("Missing required packages: ", paste(missing, collapse = ", "))


#' Add ID Column to Data Frame
#'
#' Ensures that a data frame has a unique ID column. If \code{id.name} is not provided,
#' a column named "id" is added. If \code{id.name} is provided but does not exist in the data frame,
#' it is created with unique integer values.
#'
#' @param df.analysis Data frame to which the ID column will be added.
#' @param id.name Character. Name of the ID column to add (default is \code{NULL}, which uses "id").
#'
#' @return Data frame with the ID column added if necessary.
#' @export

add_id_column <- function(df.analysis, id.name = NULL) {
  if (is.null(id.name)) {
    df.analysis$id <- seq_len(nrow(df.analysis))
    id.name <- "id"
  } else if (!(id.name %in% names(df.analysis))) {
    df.analysis[[id.name]] <- seq_len(nrow(df.analysis))
  }
  return(df.analysis)
}


#' Generate Prediction Dataset with Subgroup Treatment Recommendation
#'
#' Creates a prediction dataset with a treatment recommendation flag based on subgroup definition.
#'
#' @param df.predict Data frame for prediction (test or validation set).
#' @param sg.harm Character vector of subgroup-defining covariate names.
#' @param version Integer; 1 uses \code{dummy()}, 2 uses \code{dummy2()} for factor encoding.
#'
#' @return Data frame with treatment recommendation flag (\code{treat.recommend}).
#' @export

get_dfpred <- function(df.predict, sg.harm, version = 1) {
  if (version == 1) df.pred <- dummy(df.predict)
  if (version == 2) df.pred <- dummy2(df.predict)
  df.pred$treat.recommend <- NA
  id.harm <- paste(sg.harm, collapse = "==1 & ")
  id.harm <- paste(id.harm, "==1")
  df.pred.0 <- subset(df.pred, eval(parse(text = id.harm)))
  if (nrow(df.pred.0) > 0) {
    df.pred.0$treat.recommend <- 0
  }
  id.noharm <- paste(sg.harm, collapse = "!=1 | ")
  id.noharm <- paste(id.noharm, "!=1")
  df.pred.1 <- subset(df.pred, eval(parse(text = id.noharm)))
  if (nrow(df.pred.1) > 0) {
    df.pred.1$treat.recommend <- 1
  }
  if (nrow(df.pred.0) > 0 && nrow(df.pred.1) > 0) df.pred.out <- data.frame(rbind(df.pred.0, df.pred.1))
  if (nrow(df.pred.0) == 0 && nrow(df.pred.1) > 0) df.pred.out <- data.frame(df.pred.1)
  if (nrow(df.pred.0) > 0 && nrow(df.pred.1) == 0) df.pred.out <- data.frame(df.pred.0)
  return(df.pred.out)
}

#' Get parameter with default fallback
#' @keywords internal
get_param <- function(args_list, param_name, default_value) {
  if (hasName(args_list, param_name) && !is.null(args_list[[param_name]])) {
    return(args_list[[param_name]])
  }
  return(default_value)
}


#' @rdname forestsearch
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

#' @rdname forestsearch
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

#' Plot ForestSearch Results
#'
#' @description
#' Creates visualization of forestsearch results including subgroup characteristics,
#' consistency metrics, and treatment effect estimates.
#'
#' @param x A forestsearch object
#' @param type Character. Type of plot: "consistency", "effects", "selection", or "all"
#' @param ... Additional plotting parameters
#'
#' @export
plot.forestsearch <- function(x, type = c("consistency", "effects", "selection", "all"), ...) {
  type <- match.arg(type)

  # Implementation would create appropriate visualizations
  # This is a placeholder for the documentation

  invisible(x)
}


