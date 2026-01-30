# =============================================================================
# Test Script for plot_sg_weighted_km()
# =============================================================================
# Run this script to debug issues with plot_sg_weighted_km()
#
# Usage: source("test_plot_sg_weighted_km.R")
# =============================================================================

library(survival)

# Source the function (adjust path as needed)
# source("plot_sg_weighted_km.R")

# -----------------------------------------------------------------------------
# Create synthetic dataset similar to forestsearch output
# -----------------------------------------------------------------------------
set.seed(123)
n <- 200

df_test <- data.frame(
  id = 1:n,
  time_months = rexp(n, rate = 0.05),
  status = rbinom(n, 1, 0.7),
  hormon = rbinom(n, 1, 0.5),
  treat.recommend = c(rep(0, 80), rep(1, 120))
)

# Create a mock forestsearch object
fs_test <- list(
  df.est = df_test,
  sg.harm = c("age_high", "grade_3"),
  grp.consistency = NULL
)
class(fs_test) <- c("forestsearch", "list")

# -----------------------------------------------------------------------------
# Print diagnostic info
# -----------------------------------------------------------------------------
cat("=== Test Data Summary ===\n")
cat("  N total:", nrow(df_test), "\n")
cat("  H (treat.recommend=0):", sum(df_test$treat.recommend == 0), "\n")
cat("  Hc (treat.recommend=1):", sum(df_test$treat.recommend == 1), "\n")
cat("  Events:", sum(df_test$status), "\n")
cat("  Treatment=1:", sum(df_test$hormon), "\n\n")

# Check if weightedsurv is available
cat("=== Package Check ===\n")
ws_available <- requireNamespace("weightedsurv", quietly = TRUE)
cat("  weightedsurv available:", ws_available, "\n")
if (ws_available) {
  cat("  weightedsurv version:", as.character(packageVersion("weightedsurv")), "\n")
}
cat("\n")

# -----------------------------------------------------------------------------
# Test 1: Basic call with verbose output
# -----------------------------------------------------------------------------
cat("=== Test 1: Basic call with verbose=TRUE ===\n")
result1 <- tryCatch({
  plot_sg_weighted_km(
    fs.est = fs_test,
    outcome.name = "time_months",
    event.name = "status",
    treat.name = "hormon",
    verbose = TRUE
  )
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})
cat(if (!is.null(result1)) "SUCCESS!\n\n" else "FAILED\n\n")

# -----------------------------------------------------------------------------
# Test 2: Force base survival plotting
# -----------------------------------------------------------------------------
cat("=== Test 2: Force base survival (use_weightedsurv=FALSE) ===\n")
result2 <- tryCatch({
  plot_sg_weighted_km(
    fs.est = fs_test,
    outcome.name = "time_months",
    event.name = "status",
    treat.name = "hormon",
    use_weightedsurv = TRUE,
    verbose = TRUE
  )
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})
cat(if (!is.null(result2)) "SUCCESS!\n\n" else "FAILED\n\n")

# -----------------------------------------------------------------------------
# Test 3: Disable HR annotation
# -----------------------------------------------------------------------------
cat("=== Test 3: Disable HR annotation (show_hr=FALSE) ===\n")
result3 <- tryCatch({
  plot_sg_weighted_km(
    fs.est = fs_test,
    outcome.name = "time_months",
    event.name = "status",
    treat.name = "hormon",
    show_hr = FALSE,
    verbose = TRUE
  )
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  NULL
})
cat(if (!is.null(result3)) "SUCCESS!\n\n" else "FAILED\n\n")

cat("=== Tests Complete ===\n")
