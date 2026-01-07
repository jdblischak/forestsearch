# Enhanced AFT DGM with Spline Model Support

## Overview
This enhanced version of `generate_aft_dgm_flex()` incorporates the spline model functionality from `get_dgm_spline()` as an additional option, while preserving all original functionality.

## Files
- `generate_aft_dgm_flex_enhanced.R` - Main enhanced function (syntax error fixed)
- `simulate_from_dgm_flex.R` - Simulation function compatible with both models  
- `demo_enhanced_dgm.R` - Demonstration script with examples

## Key Enhancement
Added `model = "spline"` option that implements a spline-based hazard ratio model where the log hazard ratio varies smoothly with a continuous biomarker.

## Fixed Issues
- Changed reserved keyword `function` to `fun_name` in the process_cutpoint helper function (line 296)

## Usage

### Standard Model (Original)
```r
dgm <- generate_aft_dgm_flex(
  data = mydata,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  model = "alt",  # or "null"
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = 20,
    meno = 0
  )
)
```

### Spline Model (New)
```r
dgm <- generate_aft_dgm_flex(
  data = mydata,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime", 
  event_var = "status",
  treatment_var = "hormon",
  model = "spline",  # NEW option
  spline_params = list(
    knot = 5,                          # Knot location
    zeta = 10,                         # Evaluation point
    log_hrs = log(c(0.5, 0.75, 1.0)), # Log HRs at z=0, knot, zeta
    biomarker_var = "nodes"           # Which variable is the biomarker
  )
)
```

## Spline Model Parameters
- `knot`: Location of the spline change point (default: 5)
- `zeta`: Second evaluation point for log HR specification (default: 10)
- `log_hrs`: Vector of 3 log hazard ratios at z=0, z=knot, and z=zeta
- `biomarker_var`: Name of the continuous variable to use as biomarker (default: first continuous variable)

## Model Types
- `"alt"`: Alternative model with subgroup effects (original)
- `"null"`: Null model without subgroup effects (original)  
- `"spline"`: Spline-based hazard ratio model (NEW)

## Key Features of Spline Model
1. Creates smooth treatment effect modification based on continuous biomarker
2. Allows specification of hazard ratios at key biomarker values
3. Generates spline terms: z.k = (z - knot) * I(z > knot)
4. Compatible with existing simulation and analysis framework

## Integration Notes
The spline model is fully integrated while preserving the original logic:
- All original cutpoint specification methods remain unchanged
- Standard models ("alt", "null") work exactly as before
- Backward compatible with existing code
- Same output structure with additional spline_info when applicable

## Example Workflow
```r
# 1. Source the functions
source("generate_aft_dgm_flex_enhanced.R")
source("simulate_from_dgm_flex.R")

# 2. Generate DGM with spline model
dgm <- generate_aft_dgm_flex(
  data = your_data,
  model = "spline",
  spline_params = list(knot = 5, zeta = 10, log_hrs = log(c(0.7, 0.8, 1.2)))
)

# 3. Simulate data
sim_data <- simulate_from_dgm(dgm, n_sim = 1000)

# 4. Analyze
analyze_simulation(sim_data, dgm = dgm)
```

## Technical Details
The spline model creates the following terms:
- `z`: The biomarker variable
- `z.treat`: z * treatment interaction
- `z.k`: (z - knot) * I(z > knot) - spline term
- `z.k.treat`: z.k * treatment - spline-treatment interaction

These terms allow the log hazard ratio to change linearly up to the knot, then change slope after the knot, creating a piecewise linear function in log HR space.
