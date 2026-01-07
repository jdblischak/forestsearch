#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' Package-level imports for ForestSearch
#'
#' @importFrom graphics par
#' @importFrom stats as.formula coef na.exclude na.omit predict setNames
#' @importFrom utils hasName head tail
#' @importFrom data.table := .SD .GRP .I .N rbindlist setnames setcolorder is.data.table
#' @importFrom future multisession multicore sequential
NULL

# Note: All globalVariables() declarations are consolidated in globals.R
# Note: install.packages is intentionally NOT imported as it's generally not
# appropriate for package code. Use requireNamespace() instead.
