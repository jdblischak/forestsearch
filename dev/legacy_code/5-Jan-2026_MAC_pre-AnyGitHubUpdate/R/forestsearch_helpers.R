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


