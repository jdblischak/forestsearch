# Demonstration of Enhanced AFT DGM with Spline Model Option
# This script shows how to use the generate_aft_dgm_flex function with:
# 1. Standard subgroup models (alt/null)
# 2. New spline-based hazard ratio models

# Load required libraries
library(survival)
library(data.table)

# Source the enhanced functions
source("generate_aft_dgm_flex_with_spline.R")
source("simulate_from_dgm_flex.R")

# For demonstration, we'll create a simple dataset
# In practice, you would use your actual data (e.g., GBSG data)
set.seed(42)
n <- 500

# Create example data
example_data <- data.frame(
  # Continuous variables
  age = rnorm(n, 50, 10),
  size = rlnorm(n, log(30), 0.5),
  nodes = rpois(n, 3),
  pgr = rlnorm(n, log(50), 1),
  er = rlnorm(n, log(40), 1),

  # Factor variables
  meno = sample(0:1, n, replace = TRUE),
  grade = sample(1:3, n, replace = TRUE),

  # Treatment
  hormon = rbinom(n, 1, 0.5)
)

# Generate survival times
# For demonstration, we'll create some artificial survival data
example_data$rfstime <- rweibull(n, shape = 1.5, scale = 1000)
example_data$status <- rbinom(n, 1, 0.7)  # 70% event rate

# ============================================================================
# Example 1: Standard Alt Model with Flexible Subgroups
# ============================================================================

cat("\n========================================\n")
cat("Example 1: Standard Alt Model\n")
cat("========================================\n")

dgm_standard <- generate_aft_dgm_flex(
  data = example_data,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.25),  # Lower quartile
    meno = 0  # Pre-menopausal
  ),
  model = "alt",
  k_treat = 1.2,  # Enhance treatment effect
  k_inter = 1.5,  # Enhance interaction
  n_super = 2000,
  verbose = TRUE
)

# Print summary
print(dgm_standard)

# Simulate data from standard model
sim_standard <- simulate_from_dgm(
  dgm = dgm_standard,
  n_sim = 1000,
  draw_treatment = TRUE,
  seed = 123,
  verbose = TRUE
)

# Analyze results
cat("\nAnalysis of Standard Model Simulation:\n")
analyze_simulation(sim_standard, dgm = dgm_standard, verbose = TRUE)

# ============================================================================
# Example 2: Spline Model with Smooth Hazard Ratio Function
# ============================================================================

cat("\n\n========================================\n")
cat("Example 2: Spline Model\n")
cat("========================================\n")

dgm_spline <- generate_aft_dgm_flex(
  data = example_data,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  model = "spline",
  spline_params = list(
    knot = 5,                           # Knot at nodes = 5
    zeta = 10,                          # Evaluation point at nodes = 10
    log_hrs = log(c(0.5, 0.75, 1.2)),  # HR increases with biomarker
    biomarker_var = "nodes"            # Use nodes as biomarker
  ),
  n_super = 2000,
  verbose = TRUE
)

# Print summary
print(dgm_spline)

# Simulate data from spline model
sim_spline <- simulate_from_dgm(
  dgm = dgm_spline,
  n_sim = 1000,
  draw_treatment = TRUE,
  seed = 456,
  verbose = TRUE
)

# Analyze results
cat("\nAnalysis of Spline Model Simulation:\n")
analyze_simulation(sim_spline, dgm = dgm_spline, verbose = TRUE)

# ============================================================================
# Example 3: Visualizing Spline Model Hazard Ratios
# ============================================================================

cat("\n\n========================================\n")
cat("Example 3: Visualizing Spline Effects\n")
cat("========================================\n")

# For spline model, we can visualize how HR changes with biomarker
if (dgm_spline$model_type == "spline") {

  # Extract parameters
  knot <- dgm_spline$spline_info$knot
  b0 <- dgm_spline$model_params$b0_weibull

  # Create biomarker range
  z_range <- seq(0, 20, by = 0.5)

  # Calculate log HR at each point
  log_hr <- b0[1] +
            b0[3] * z_range +
            b0[5] * (z_range - knot) * ifelse(z_range > knot, 1, 0)

  # Convert to HR
  hr <- exp(log_hr)

  # Create plot
  par(mfrow = c(1, 2))

  # Plot 1: Log Hazard Ratio
  plot(z_range, log_hr, type = "l", lwd = 2,
       xlab = "Biomarker (nodes)",
       ylab = "Log Hazard Ratio",
       main = "Log HR vs Biomarker",
       col = "blue")
  abline(v = knot, lty = 2, col = "red")
  abline(h = 0, lty = 3, col = "gray")
  text(knot, min(log_hr), paste("Knot =", knot), pos = 4)

  # Plot 2: Hazard Ratio
  plot(z_range, hr, type = "l", lwd = 2,
       xlab = "Biomarker (nodes)",
       ylab = "Hazard Ratio",
       main = "HR vs Biomarker",
       col = "darkgreen")
  abline(v = knot, lty = 2, col = "red")
  abline(h = 1, lty = 3, col = "gray")
  text(knot, min(hr), paste("Knot =", knot), pos = 4)

  # Add data distribution
  rug(sim_spline$z, col = rgb(0, 0, 0, 0.3))

  par(mfrow = c(1, 1))
}

# ============================================================================
# Example 4: Comparing Multiple Spline Configurations
# ============================================================================

cat("\n\n========================================\n")
cat("Example 4: Comparing Spline Configurations\n")
cat("========================================\n")

# Define different spline configurations
spline_configs <- list(
  protective = list(
    knot = 5,
    zeta = 10,
    log_hrs = log(c(0.5, 0.6, 0.7)),  # Protective effect
    biomarker_var = "nodes"
  ),
  harmful = list(
    knot = 5,
    zeta = 10,
    log_hrs = log(c(1.2, 1.5, 2.0)),  # Harmful effect
    biomarker_var = "nodes"
  ),
  u_shaped = list(
    knot = 5,
    zeta = 10,
    log_hrs = log(c(1.5, 0.7, 1.8)),  # U-shaped relationship
    biomarker_var = "nodes"
  )
)

# Generate DGMs for each configuration
results_comparison <- list()

for (config_name in names(spline_configs)) {
  cat("\n--- Configuration:", config_name, "---\n")

  dgm_temp <- generate_aft_dgm_flex(
    data = example_data,
    continuous_vars = c("age", "size", "nodes", "pgr", "er"),
    factor_vars = c("meno", "grade"),
    outcome_var = "rfstime",
    event_var = "status",
    treatment_var = "hormon",
    model = "spline",
    spline_params = spline_configs[[config_name]],
    n_super = 1000,
    verbose = FALSE
  )

  # Simulate and analyze
  sim_temp <- simulate_from_dgm(
    dgm = dgm_temp,
    n_sim = 500,
    seed = which(names(spline_configs) == config_name) * 100,
    verbose = FALSE
  )

  analysis_temp <- analyze_simulation(sim_temp, dgm = dgm_temp, verbose = FALSE)

  results_comparison[[config_name]] <- list(
    dgm = dgm_temp,
    sim = sim_temp,
    analysis = analysis_temp,
    log_hrs = spline_configs[[config_name]]$log_hrs
  )

  cat("  Log HRs specified:", paste(round(spline_configs[[config_name]]$log_hrs, 2),
                                    collapse = ", "), "\n")
  cat("  Estimated HR:", round(analysis_temp$hr_estimate, 3), "\n")
  cat("  95% CI: [", round(analysis_temp$hr_ci_lower, 3), ",",
      round(analysis_temp$hr_ci_upper, 3), "]\n")
}

# ============================================================================
# Example 5: Multiple Simulations for Operating Characteristics
# ============================================================================

cat("\n\n========================================\n")
cat("Example 5: Operating Characteristics\n")
cat("========================================\n")

# Function to run multiple simulations
run_simulation_study <- function(dgm, n_sims = 100, n_per_sim = 500, verbose = FALSE) {

  results <- data.frame(
    sim = 1:n_sims,
    hr_est = numeric(n_sims),
    hr_lower = numeric(n_sims),
    hr_upper = numeric(n_sims),
    p_value = numeric(n_sims),
    covers_truth = logical(n_sims)
  )

  true_hr <- dgm$hazard_ratios$overall

  for (i in 1:n_sims) {
    if (verbose && i %% 10 == 0) cat("Simulation", i, "of", n_sims, "\n")

    # Simulate data
    sim_data <- simulate_from_dgm(
      dgm = dgm,
      n_sim = n_per_sim,
      seed = i * 1000,
      verbose = FALSE
    )

    # Analyze
    analysis <- analyze_simulation(sim_data, dgm = NULL, verbose = FALSE)

    # Store results
    results$hr_est[i] <- analysis$hr_estimate
    results$hr_lower[i] <- analysis$hr_ci_lower
    results$hr_upper[i] <- analysis$hr_ci_upper
    results$p_value[i] <- analysis$p_value

    if (!is.null(true_hr)) {
      results$covers_truth[i] <- (true_hr >= analysis$hr_ci_lower) &
                                 (true_hr <= analysis$hr_ci_upper)
    }
  }

  return(results)
}

# Run simulation study for spline model
cat("Running 50 simulations for spline model...\n")
sim_study_results <- run_simulation_study(
  dgm = dgm_spline,
  n_sims = 50,
  n_per_sim = 500,
  verbose = TRUE
)

# Summarize operating characteristics
cat("\n--- Operating Characteristics ---\n")
cat("Mean HR estimate:", round(mean(sim_study_results$hr_est), 3), "\n")
cat("SD of HR estimates:", round(sd(sim_study_results$hr_est), 3), "\n")
cat("Mean CI width:", round(mean(sim_study_results$hr_upper -
                                 sim_study_results$hr_lower), 3), "\n")

if (!is.null(dgm_spline$hazard_ratios$overall)) {
  cat("Coverage probability:", round(mean(sim_study_results$covers_truth), 3), "\n")
  cat("Bias:", round(mean(sim_study_results$hr_est) -
                    dgm_spline$hazard_ratios$overall, 3), "\n")
}

cat("Power (proportion p < 0.05):", round(mean(sim_study_results$p_value < 0.05), 3), "\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n\n========================================\n")
cat("Summary\n")
cat("========================================\n")

cat("\nThe enhanced generate_aft_dgm_flex function now supports:\n")
cat("1. Standard subgroup models (alt/null) with flexible cutpoint definitions\n")
cat("2. NEW: Spline-based hazard ratio models for smooth biomarker effects\n")
cat("\nKey features of the spline model:\n")
cat("- Allows smooth variation of treatment effect with a continuous biomarker\n")
cat("- Supports specification of log hazard ratios at key points\n")
cat("- Creates realistic non-linear treatment-biomarker interactions\n")
cat("- Compatible with existing simulation and analysis framework\n")
cat("\nUse model='spline' and provide spline_params to activate this functionality.\n")


# Minimal example



# Create a simple example dataset
set.seed(42)
n <- 300

# Generate example data similar to GBSG structure
example_data <- data.frame(
  # Continuous variables
  age = rnorm(n, 50, 10),
  size = rlnorm(n, log(30), 0.5),
  nodes = rpois(n, 3),  # This will be our biomarker for spline model
  pgr = rlnorm(n, log(50), 1),
  er = rlnorm(n, log(40), 1),

  # Factor variables
  meno = factor(sample(0:1, n, replace = TRUE)),
  grade = factor(sample(1:3, n, replace = TRUE)),

  # Treatment
  hormon = rbinom(n, 1, 0.5)
)

# Generate survival times with treatment effect
lp <- -2 + 0.02 * example_data$age - 0.3 * example_data$hormon +
  0.05 * example_data$nodes - 0.1 * example_data$hormon * (example_data$nodes > 3)
example_data$rfstime <- rweibull(n, shape = 1.5, scale = exp(-lp))
example_data$status <- rbinom(n, 1, 0.7)

cat("========================================\n")
cat("Example 1: Standard Alt Model\n")
cat("========================================\n\n")

# Standard model with subgroups
dgm_standard <- generate_aft_dgm_flex(
  data = example_data,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "nodes"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.25),  # Lower quartile
    nodes = 3  # Fixed cutpoint
  ),
  model = "alt",
  n_super = 1000,
  verbose = TRUE
)

# Simulate data
sim_standard <- simulate_from_dgm(
  dgm = dgm_standard,
  n_sim = 500,
  seed = 123,
  verbose = TRUE
)

# Basic Cox analysis
cox_standard <- coxph(Surv(time, event) ~ treat, data = sim_standard)
cat("\nCox model results (standard):\n")
print(summary(cox_standard)$conf.int)

cat("\n\n========================================\n")
cat("Example 2: Spline Model\n")
cat("========================================\n\n")

# Spline model with smooth biomarker effect
dgm_spline <- generate_aft_dgm_flex(
  data = example_data,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  model = "spline",
  spline_params = list(
    knot = 3,                           # Knot at nodes = 3
    zeta = 8,                           # Second point at nodes = 8
    log_hrs = log(c(0.5, 0.7, 1.2)),   # HR changes from protective to harmful
    biomarker_var = "nodes"            # Use nodes as the biomarker
  ),
  n_super = 1000,
  verbose = TRUE
)

# Simulate data
sim_spline <- simulate_from_dgm(
  dgm = dgm_spline,
  n_sim = 500,
  seed = 456,
  verbose = TRUE
)

# Basic Cox analysis
cox_spline <- coxph(Surv(time, event) ~ treat, data = sim_spline)
cat("\nCox model results (spline):\n")
print(summary(cox_spline)$conf.int)

# Visualize the spline effect
cat("\n\n========================================\n")
cat("Visualizing Spline Effect\n")
cat("========================================\n\n")

if (dgm_spline$model_type == "spline") {
  # Extract spline parameters
  knot <- dgm_spline$spline_info$knot
  b0 <- dgm_spline$model_params$b0_weibull

  # Create biomarker range
  z_range <- seq(0, 15, by = 0.5)

  # Calculate log HR at each point
  # This shows how treatment effect varies with the biomarker
  log_hr <- b0[1] +                                          # Base treatment effect
    b0[3] * z_range +                                # Linear interaction with z
    b0[5] * (z_range - knot) * (z_range > knot)    # Change in slope after knot

  # Convert to HR
  hr <- exp(log_hr)

  # Print summary at key points
  cat("Treatment Hazard Ratio at key biomarker values:\n")
  key_points <- c(0, knot, dgm_spline$spline_info$zeta, 10)
  for (z_val in key_points) {
    log_hr_z <- b0[1] + b0[3] * z_val + b0[5] * (z_val - knot) * (z_val > knot)
    cat(sprintf("  At nodes = %2d: HR = %.3f (log HR = %.3f)\n",
                z_val, exp(log_hr_z), log_hr_z))
  }

  # Create a simple text-based visualization
  cat("\nTreatment Effect Pattern (HR vs Biomarker):\n")
  cat("HR:  0.5    1.0    1.5    2.0\n")
  cat("     |      |      |      |\n")

  for (i in seq(1, length(z_range), by = 2)) {
    z_val <- z_range[i]
    hr_val <- hr[i]

    # Create a simple bar
    n_chars <- round(hr_val * 10)
    n_chars <- pmin(40, pmax(1, n_chars))

    bar <- paste(rep("=", n_chars), collapse = "")
    marker <- ifelse(z_val == knot, " <-- KNOT", "")

    cat(sprintf("z=%2.0f %s%s\n", z_val, bar, marker))
  }

  cat("\nInterpretation:\n")
  cat("- The treatment effect varies smoothly with the biomarker (nodes)\n")
  cat("- At low biomarker values, treatment is",
      ifelse(hr[1] < 1, "protective", "harmful"), "\n")
  cat("- The effect changes slope at the knot (nodes =", knot, ")\n")
  cat("- At high biomarker values, treatment is",
      ifelse(hr[length(hr)] < 1, "protective", "harmful"), "\n")
}

cat("\n\n========================================\n")
cat("Comparison Summary\n")
cat("========================================\n\n")

cat("Model Comparison:\n")
cat("-----------------\n")
cat("Standard Alt Model:\n")
cat("  - Uses discrete subgroups based on cutpoints\n")
cat("  - Treatment effect is constant within/outside subgroup\n")
cat("  - Simpler interpretation\n\n")

cat("Spline Model:\n")
cat("  - Treatment effect varies continuously with biomarker\n")
cat("  - Allows for non-linear relationships\n")
cat("  - More flexible but complex interpretation\n\n")

cat("Both models:\n")
cat("  - Generate realistic survival data\n")
cat("  - Support the same simulation framework\n")
cat("  - Can be analyzed with standard survival methods\n")

cat("\n========================================\n")
cat("Example Complete\n")
cat("========================================\n")



## NEW Refactored




# Load example data
data(cancer)  # Loads gbsg dataset

# ==============================================================================
# EXAMPLE 1: THRESHOLD-BASED SUBGROUPS (Original Approach)
# ==============================================================================

cat("\n=== EXAMPLE 1: Threshold-based Subgroups ===\n\n")

# Method A: Using the convenience wrapper (easiest)
dgm_threshold_simple <- generate_aft_dgm_threshold(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "pgr"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.25),   # er <= 25th percentile
    pgr = list(type = "function", fun = median)   # pgr <= median
  ),
  model = "alt",
  k_treat = 0.9,
  k_inter = 1.5,
  n_super = 5000,
  verbose = TRUE
)

print(dgm_threshold_simple)

# Method B: Using explicit SubgroupDefinition object (more flexible)
subgroup_def_threshold <- ThresholdSubgroup(
  vars = c("er", "meno"),
  cuts = list(
    er = 20,  # Simple fixed cutpoint
    meno = 0  # meno == 0
  )
)

dgm_threshold_explicit <- generate_aft_dgm_refactored(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_def = subgroup_def_threshold,
  model = "alt",
  n_super = 5000,
  verbose = TRUE
)

# Simulate data from the threshold DGM
sim_threshold <- simulate_from_dgm_refactored(
  dgm = dgm_threshold_simple,
  n = 700,
  rand_ratio = 1,
  max_follow = 84,
  seed = 123
)

cat("\nSimulated data dimensions:", dim(sim_threshold), "\n")
cat("Event rate:", mean(sim_threshold$event_sim), "\n")

# ==============================================================================
# EXAMPLE 2: SPLINE-BASED SUBGROUPS (New Capability)
# ==============================================================================

cat("\n\n=== EXAMPLE 2: Spline-based Subgroups ===\n\n")

# Method A: Using the convenience wrapper with automatic knots
dgm_spline_auto <- generate_aft_dgm_spline(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  spline_var = "er",           # Use ER as spline variable
  knots = NULL,                # Automatic knot selection
  degree = 3,                  # Cubic spline
  model = "alt",
  k_treat = 0.9,
  k_inter = 1.0,
  n_super = 5000,
  verbose = TRUE
)

print(dgm_spline_auto)

# Method B: Using explicit SplineSubgroup with custom knots
subgroup_def_spline <- SplineSubgroup(
  var = "er",
  knots = c(10, 30, 100),     # Custom knot positions
  degree = 3,
  boundary_knots = c(0, 300)
)

dgm_spline_custom <- generate_aft_dgm_refactored(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_def = subgroup_def_spline,
  model = "alt",
  n_super = 5000,
  verbose = TRUE
)

# Simulate data from the spline DGM
sim_spline <- simulate_from_dgm(
  dgm = dgm_spline_auto,
  n = 700,
  rand_ratio = 1,
  max_follow = 84,
  seed = 456
)

cat("\nSimulated spline data dimensions:", dim(sim_spline), "\n")
cat("Event rate:", mean(sim_spline$event_sim), "\n")

# ==============================================================================
# EXAMPLE 3: ADVANCED THRESHOLD SPECIFICATIONS
# ==============================================================================

cat("\n\n=== EXAMPLE 3: Advanced Threshold Specifications ===\n\n")

# Demonstrate various cutpoint types
dgm_advanced <- generate_aft_dgm_threshold(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("age", "nodes", "grade"),
  subgroup_cuts = list(
    age = list(type = "range", min = 40, max = 60),      # Age between 40-60
    nodes = list(type = "greater", quantile = 0.75),     # High node count
    grade = list(type = "multiple", values = c(2, 3))    # Grade 2 or 3
  ),
  model = "alt",
  k_inter = 2.0,  # Strong interaction
  verbose = TRUE
)

sim_advanced <- simulate_from_dgm(dgm_advanced, n = 1000, seed = 789)

# ==============================================================================
# EXAMPLE 4: NULL MODEL (No Subgroup Effects)
# ==============================================================================

cat("\n\n=== EXAMPLE 4: Null Model (No Subgroup Effects) ===\n\n")

dgm_null <- generate_aft_dgm_threshold(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = NULL,  # No subgroups
  model = "null",        # Null model
  k_treat = 0.85,
  verbose = TRUE
)

sim_null <- simulate_from_dgm(dgm_null, n = 700, seed = 999)

# ==============================================================================
# EXAMPLE 5: COMPARING THRESHOLD VS SPLINE MODELS
# ==============================================================================

cat("\n\n=== EXAMPLE 5: Comparing Threshold vs Spline ===\n\n")

# Threshold model: ER <= 25th percentile
dgm_er_threshold <- generate_aft_dgm_threshold(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = "er",
  subgroup_cuts = list(er = list(type = "quantile", value = 0.25)),
  model = "alt",
  k_inter = 1.5,
  seed = 111,
  verbose = FALSE
)

# Spline model: Smooth treatment effect across ER
dgm_er_spline <- generate_aft_dgm_spline(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  spline_var = "er",
  degree = 3,
  model = "alt",
  k_inter = 1.5,
  seed = 111,
  verbose = FALSE
)

cat("Threshold model hazard ratios:\n")
print(dgm_er_threshold$hazard_ratios)

cat("\nSpline model hazard ratios:\n")
print(dgm_er_spline$hazard_ratios)

# Simulate from both
sim_threshold_er <- simulate_from_dgm(dgm_er_threshold, n = 1000, seed = 222)
sim_spline_er <- simulate_from_dgm(dgm_er_spline, n = 1000, seed = 222)

# Compare event rates
cat("\nThreshold model event rate:", mean(sim_threshold_er$event_sim), "\n")
cat("Spline model event rate:", mean(sim_spline_er$event_sim), "\n")

# ==============================================================================
# EXAMPLE 6: CUSTOM CUTPOINT FUNCTION
# ==============================================================================

cat("\n\n=== EXAMPLE 6: Custom Cutpoint Function ===\n\n")

# Define a custom function for complex subgroup definition
# Example: Patients with both low ER OR high PGR
dgm_custom <- generate_aft_dgm_threshold(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = "er",
  subgroup_cuts = list(
    er = list(
      type = "custom",
      fun = function(x) {
        # Complex condition: lowest 30% OR highest 10%
        x <= quantile(x, 0.30) | x >= quantile(x, 0.90)
      }
    )
  ),
  model = "alt",
  verbose = TRUE
)

# ==============================================================================
# EXAMPLE 7: FLEXIBLE CENSORING
# ==============================================================================

cat("\n\n=== EXAMPLE 7: Different Censoring Mechanisms ===\n\n")

# Weibull censoring (default - fitted from data)
dgm_weibull_cens <- generate_aft_dgm_threshold(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = "er",
  subgroup_cuts = list(er = 20),
  cens_type = "weibull",
  model = "alt",
  verbose = FALSE
)

# Uniform censoring
dgm_uniform_cens <- generate_aft_dgm_threshold(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = "er",
  subgroup_cuts = list(er = 20),
  cens_type = "uniform",
  cens_params = list(min = 12, max = 100),
  model = "alt",
  verbose = FALSE
)

# Simulate with different censoring adjustments
sim_light_cens <- simulate_from_dgm(dgm_weibull_cens, n = 700,
                                    cens_adjust = -0.5, seed = 333)
sim_heavy_cens <- simulate_from_dgm(dgm_weibull_cens, n = 700,
                                    cens_adjust = 0.5, seed = 333)

cat("Light censoring event rate:", mean(sim_light_cens$event_sim), "\n")
cat("Heavy censoring event rate:", mean(sim_heavy_cens$event_sim), "\n")

# ==============================================================================
# EXAMPLE 8: WORKING WITH CUSTOM DATA
# ==============================================================================

cat("\n\n=== EXAMPLE 8: Custom Dataset ===\n\n")

# Create synthetic custom dataset
set.seed(42)
custom_data <- data.frame(
  patient_id = 1:500,
  time = rexp(500, rate = 0.01),
  status = rbinom(500, 1, 0.7),
  age = rnorm(500, 50, 10),
  biomarker = rgamma(500, 2, 1),
  stage = sample(1:4, 500, replace = TRUE),
  sex = sample(c("M", "F"), 500, replace = TRUE),
  treatment = rbinom(500, 1, 0.5)
)

# Generate DGM with spline on biomarker
dgm_custom <- generate_aft_dgm_spline(
  data = custom_data,
  continuous_vars = c("age", "biomarker"),
  factor_vars = c("stage", "sex"),
  outcome_var = "time",
  event_var = "status",
  treatment_var = "treatment",
  spline_var = "biomarker",
  degree = 3,
  model = "alt",
  n_super = 2000,
  verbose = TRUE
)

sim_custom <- simulate_from_dgm(dgm_custom, n = 800, seed = 444)

# ==============================================================================
# EXAMPLE 9: EFFECT MODIFICATION
# ==============================================================================

cat("\n\n=== EXAMPLE 9: Effect Modification ===\n\n")

# Different treatment effect strengths
dgm_weak_effect <- generate_aft_dgm_threshold(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = "er",
  subgroup_cuts = list(er = 20),
  model = "alt",
  k_treat = 0.5,   # Weak main effect
  k_inter = 3.0,   # Strong interaction
  verbose = FALSE
)

dgm_strong_effect <- generate_aft_dgm_threshold(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = "er",
  subgroup_cuts = list(er = 20),
  model = "alt",
  k_treat = 2.0,   # Strong main effect
  k_inter = 0.5,   # Weak interaction
  verbose = FALSE
)

cat("Weak treatment, strong interaction HRs:\n")
print(dgm_weak_effect$hazard_ratios)

cat("\nStrong treatment, weak interaction HRs:\n")
print(dgm_strong_effect$hazard_ratios)

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n\n" ,rep("=", 70), "\n", sep = "")
cat("SUMMARY OF EXAMPLES\n")
cat(rep("=", 70), "\n", sep = "")
cat("\n")
cat("1. Threshold subgroups      - Traditional binary/categorical subgroups\n")
cat("2. Spline subgroups         - Smooth treatment effect variation\n")
cat("3. Advanced thresholds      - Multiple complex cutpoint types\n")
cat("4. Null model               - No subgroup effects\n")
cat("5. Model comparison         - Threshold vs Spline\n")
cat("6. Custom functions         - User-defined subgroup criteria\n")
cat("7. Censoring mechanisms     - Weibull vs Uniform censoring\n")
cat("8. Custom data              - Non-GBSG datasets\n")
cat("9. Effect modification      - Varying treatment and interaction effects\n")
cat("\n")
cat("Key Benefits of Refactored Code:\n")
cat("  ✓ Modular design with clear separation of concerns\n")
cat("  ✓ Support for both threshold and spline subgroups\n")
cat("  ✓ Flexible cutpoint specifications\n")
cat("  ✓ Improved code reusability and maintainability\n")
cat("  ✓ Backward compatible with original functions\n")
cat("  ✓ Extensive documentation and examples\n")
cat("\n")

# ==============================================================================
# PERFORMANCE COMPARISON (Optional)
# ==============================================================================

if (FALSE) {  # Set to TRUE to run benchmarking

  cat("\n\n=== PERFORMANCE COMPARISON ===\n\n")

  library(microbenchmark)

  mb_result <- microbenchmark(
    threshold = generate_aft_dgm_threshold(
      data = gbsg,
      continuous_vars = c("age", "size", "nodes", "pgr", "er"),
      factor_vars = c("meno", "grade"),
      outcome_var = "rfstime",
      event_var = "status",
      treatment_var = "hormon",
      subgroup_vars = "er",
      subgroup_cuts = list(er = 20),
      model = "alt",
      verbose = FALSE
    ),

    spline = generate_aft_dgm_spline(
      data = gbsg,
      continuous_vars = c("age", "size", "nodes", "pgr", "er"),
      factor_vars = c("meno", "grade"),
      outcome_var = "rfstime",
      event_var = "status",
      treatment_var = "hormon",
      spline_var = "er",
      model = "alt",
      verbose = FALSE
    ),

    times = 10
  )

  print(mb_result)
}




