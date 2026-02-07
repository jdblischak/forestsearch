# # List of required packages for ForestSearch analysis
# required_packages <- c("grf","policytree","data.table","randomForest","survival","weightedsurv","future.apply")
# missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
# if(length(missing) > 0) stop("Missing required packages: ", paste(missing, collapse = ", "))


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
#' Creates a prediction dataset with a treatment recommendation flag based
#' on the subgroup definition. Supports both q-coded column names (for
#' internally encoded data) and label expressions (for raw prediction data).
#'
#' @param df.predict Data frame for prediction (test or validation set).
#' @param sg.harm Character vector of subgroup-defining labels. If named,
#'   values are used as label expressions (e.g., "er <= 0", "size <= 35").
#' @param version Integer; encoding version (maintained for backward
#'   compatibility). Default: 1.
#'
#' @return Data frame with treatment recommendation flag
#'   (\code{treat.recommend}): 0 for harm subgroup, 1 for complement.
#' @export
get_dfpred <- function(df.predict, sg.harm, version = 1) {

  df.pred <- df.predict
  df.pred$treat.recommend <- NA

  # Extract label expressions, stripping braces and negation
  # sg.harm values look like "{er <= 0}" or "!{er <= 0}"
  labels <- if (!is.null(names(sg.harm))) unname(sg.harm) else sg.harm

  # Build membership indicator for each factor
  in_harm <- rep(TRUE, nrow(df.pred))

  for (lab in labels) {
    # Strip braces: "{er <= 0}" -> "er <= 0"
    clean <- gsub("^!?\\{(.*)\\}$", "\\1", lab)
    is_negated <- grepl("^!", lab)

    # Evaluate the expression on the data
    member <- tryCatch(
      eval(parse(text = clean), envir = df.pred),
      error = function(e) {
        warning("Could not evaluate: ", clean, " - ", e$message)
        rep(NA, nrow(df.pred))
      }
    )

    # Apply negation if needed
    if (is_negated) member <- !member

    in_harm <- in_harm & member
  }

  df.pred$treat.recommend <- ifelse(in_harm, 0L, 1L)
  df.pred
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


