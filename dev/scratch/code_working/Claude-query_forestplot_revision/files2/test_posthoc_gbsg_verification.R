#' =============================================================================
#' TEST SCRIPT: Verify posthoc subgroup handling in plot_subgroup_results_forestplot
#' =============================================================================
#'
#' This script uses the GBSG dataset to test and verify:
#' 1. Correct handling of subgroup_list with type="posthoc"
#' 2. Correct labeling of identified subgroups (benefit vs harm)
#' 3. Proper alignment of bias-corrected estimates
#'
#' @author ForestSearch Development Team

# ==============================================================================
# SETUP
# ==============================================================================

library(survival)

# Load GBSG data
data(gbsg, package = "survival")

# Prepare analysis dataset
df_gbsg <- gbsg
df_gbsg$rfs_months <- df_gbsg$rfstime / 30.44

# Create binary markers for subgroup definitions
df_gbsg$er_high <- ifelse(df_gbsg$er > median(df_gbsg$er, na.rm = TRUE), 1, 0)
df_gbsg$pgr_high <- ifelse(df_gbsg$pgr > median(df_gbsg$pgr, na.rm = TRUE), 1, 0)
df_gbsg$er_low <- 1 - df_gbsg$er_high
df_gbsg$pgr_low <- 1 - df_gbsg$pgr_high

# Analysis parameters
outcome.name <- "rfs_months"
event.name <- "status"
treat.name <- "hormon"

cat("=================================================================\n")
cat("TEST: Posthoc Subgroup Handling Verification\n")
cat("=================================================================\n\n")

cat("Dataset: GBSG (German Breast Cancer Study Group)\n")
cat("N =", nrow(df_gbsg), "\n")
cat("Treatment (hormon=1):", sum(df_gbsg$hormon), "\n")
cat("Control (hormon=0):", sum(1-df_gbsg$hormon), "\n")
cat("Events:", sum(df_gbsg$status), "\n\n")

# ==============================================================================
# TEST 1: Verify subgroup_list posthoc entries are NOT processed
# ==============================================================================

cat("=================================================================\n")
cat("TEST 1: subgroup_list posthoc entries are NOT being added\n")
cat("=================================================================\n\n")

cat("CURRENT CODE (lines 267-278):\n")
cat("```r\n")
cat("posthoc_sgs <- Filter(function(x) x$type == 'posthoc', subgroup_list)\n")
cat("if (length(posthoc_sgs) > 0) {\n")
cat("  separator_row <- create_header_row('Post-hoc subgroups', E.name, C.name)\n")
cat("  dt <- rbind(dt, separator_row)\n")
cat("  row_types <- c(row_types, 'separator')\n")
cat("  \n")
cat("  # Add blank row\n")
cat("  blank_row <- create_header_row(' ', E.name, C.name)\n")
cat("  dt <- rbind(dt, blank_row)\n")
cat("  row_types <- c(row_types, 'blank')\n")
cat("}\n")
cat("```\n\n")

cat("PROBLEM: The code filters posthoc subgroups but NEVER adds them!\n")
cat("         Only separator and blank rows are added.\n")
cat("         There is no loop to process each posthoc subgroup.\n\n")

# Define test subgroups to demonstrate
test_subgroups <- list(
  ref_er = list(
    subset_expr = "er_high == 1",
    name = "ER High (Reference)",
    type = "reference"
  ),
  posthoc_er_pgr = list(
    subset_expr = "er_high == 1 & pgr_high == 1",
    name = "ER High & PgR High",
    type = "posthoc"
  )
)

posthoc_entries <- Filter(function(x) x$type == "posthoc", test_subgroups)
cat("Example subgroup_list with posthoc entries:\n")
for (nm in names(posthoc_entries)) {
  sg <- posthoc_entries[[nm]]
  df_temp <- subset(df_gbsg, eval(parse(text = sg$subset_expr)))
  cat(sprintf("  - %s: n=%d (would NOT appear in current output)\n", sg$name, nrow(df_temp)))
}

cat("\nFIX NEEDED: Add loop to process posthoc_sgs after separator:\n")
cat("```r\n")
cat("for (sg in posthoc_sgs) {\n")
cat("  df_sg <- subset(df_analysis, eval(parse(text = sg$subset_expr)))\n")
cat("  if (nrow(df_sg) > 10) {\n")
cat("    sg_row <- create_hr_row(df_sg, sg$name, ...)\n")
cat("    dt <- rbind(dt, sg_row)\n")
cat("    row_types <- c(row_types, 'posthoc_list')\n")
cat("  }\n")
cat("}\n")
cat("```\n\n")

# ==============================================================================
# TEST 2: Verify ForestSearch labeling conventions
# ==============================================================================

cat("=================================================================\n")
cat("TEST 2: ForestSearch Labeling Convention\n")
cat("=================================================================\n\n")

cat("ForestSearch naming conventions:\n")
cat("  - sg.harm: Contains the definition of the HARM/QUESTIONABLE subgroup (H)\n")
cat("  - treat.recommend == 0: Patient is IN the harm subgroup (H)\n")
cat("  - treat.recommend == 1: Patient is in the COMPLEMENT (Hc, benefit)\n\n")

cat("For est.scale = 'hr' (searching for harm subgroup):\n")
cat("  - H  (treat.recommend=0): Subgroup defined by sg.harm (elevated HR, harm)\n")
cat("  - Hc (treat.recommend=1): Complement of sg.harm (potential benefit)\n\n")

# ==============================================================================
# TEST 3: Demonstrate the labeling bug with synthetic ForestSearch output
# ==============================================================================

cat("=================================================================\n")
cat("TEST 3: Demonstrate Labeling Bug\n")
cat("=================================================================\n\n")

# Simulate ForestSearch output
# Suppose ForestSearch identified: sg.harm = c("{er<=27}", "{pgr<=41}")
# This means H = patients with er<=27 AND pgr<=41 (harm group)

cat("Simulated ForestSearch scenario:\n")
cat("  - sg.harm = c('{er<=27}', '{pgr<=41}')\n")
cat("  - Interpretation: H = patients with er<=27 AND pgr<=41\n")
cat("  - These patients have ELEVATED HR (harm/questionable benefit)\n\n")

# Create simulated treat.recommend
er_cutoff <- 27
pgr_cutoff <- 41

df_gbsg$treat.recommend <- ifelse(
  df_gbsg$er <= er_cutoff & df_gbsg$pgr <= pgr_cutoff, 
  0,  # In harm subgroup (H) - treat.recommend = 0
  1   # In complement (Hc) - treat.recommend = 1
)

cat("Subgroup distribution:\n")
n_H <- sum(df_gbsg$treat.recommend == 0)
n_Hc <- sum(df_gbsg$treat.recommend == 1)
cat(sprintf("  - H  (treat.recommend=0): n=%d - HARM group (er<=27 & pgr<=41)\n", n_H))
cat(sprintf("  - Hc (treat.recommend=1): n=%d - BENEFIT group (complement)\n", n_Hc))

# Compute HRs to verify
compute_hr <- function(df, outcome.name, event.name, treat.name) {
  formula_str <- paste0("Surv(", outcome.name, ", ", event.name, ") ~ ", treat.name)
  fit <- survival::coxph(as.formula(formula_str), data = df)
  ci <- summary(fit)$conf.int[c(1, 3, 4)]
  c(HR = ci[1], lower = ci[2], upper = ci[3])
}

hr_H <- compute_hr(subset(df_gbsg, treat.recommend == 0), outcome.name, event.name, treat.name)
hr_Hc <- compute_hr(subset(df_gbsg, treat.recommend == 1), outcome.name, event.name, treat.name)

cat(sprintf("\n  - H  HR = %.3f (%.3f, %.3f) - should be HIGHER (less benefit)\n", 
            hr_H[1], hr_H[2], hr_H[3]))
cat(sprintf("  - Hc HR = %.3f (%.3f, %.3f) - should be LOWER (more benefit)\n",
            hr_Hc[1], hr_Hc[2], hr_Hc[3]))

# ==============================================================================
# TEST 4: Show the labeling bug in current code
# ==============================================================================

cat("\n=================================================================\n")
cat("TEST 4: The Labeling Bug\n")
cat("=================================================================\n\n")

cat("CURRENT CODE (lines 293-312):\n")
cat("```r\n")
cat("sg_harm <- fs.est$sg.harm\n")
cat("if (!is.null(sg_harm)) {\n")
cat("  sg1_label <- paste(sg_harm, collapse = ' & ')  # = '{er<=27} & {pgr<=41}'\n")
cat("}\n")
cat("\n")
cat("if (est.scale == 'hr') {\n")
cat("  df_benefit <- subset(df_sg, treat.recommend == 1)  # Correct: Hc\n")
cat("  df_question <- subset(df_sg, treat.recommend == 0) # Correct: H\n")
cat("  benefit_name <- sg1_label                          # BUG!\n")
cat("  question_name <- paste0('NOT (', sg1_label, ')')   # BUG!\n")
cat("}\n")
cat("```\n\n")

cat("THE BUG:\n")
cat("  - sg1_label = sg.harm = '{er<=27} & {pgr<=41}'\n")
cat("  - sg.harm defines the HARM subgroup (H), NOT the benefit subgroup (Hc)\n")
cat("  - Current code assigns:\n")
cat("      benefit_name = '{er<=27} & {pgr<=41}'  <- WRONG! This is the harm definition\n")
cat("      question_name = 'NOT ({er<=27} & {pgr<=41})'  <- WRONG! This is actually benefit\n\n")

cat("WHAT THE FOREST PLOT WOULD SHOW (INCORRECTLY):\n")
cat("  Row: '{er<=27} & {pgr<=41}'         <- Labeled as 'benefit' but actually HARM\n")
cat(sprintf("        Full-analysis: HR = %.2f     <- This is the Hc HR, mislabeled\n", hr_Hc[1]))
cat("  Row: 'NOT ({er<=27} & {pgr<=41})'   <- Labeled as 'questionable' but actually BENEFIT\n")
cat(sprintf("        Full-analysis: HR = %.2f     <- This is the H HR, mislabeled\n\n", hr_H[1]))

cat("CORRECT LABELING SHOULD BE:\n")
cat("  For est.scale = 'hr':\n")
cat("    benefit_name = 'NOT ({er<=27} & {pgr<=41})'  <- Hc is complement of sg.harm\n")
cat("    question_name = '{er<=27} & {pgr<=41}'       <- H is defined by sg.harm\n\n")

# ==============================================================================
# TEST 5: Proposed corrected code
# ==============================================================================

cat("=================================================================\n")
cat("TEST 5: Corrected Code\n")
cat("=================================================================\n\n")

cat("CORRECTED CODE:\n")
cat("```r\n")
cat("# Get subgroup definition from sg.harm\n")
cat("# IMPORTANT: sg.harm defines the HARM subgroup (H), not benefit\n")
cat("sg_harm <- fs.est$sg.harm\n")
cat("if (!is.null(sg_harm)) {\n")
cat("  harm_label <- paste(sg_harm, collapse = ' & ')\n")
cat("} else {\n")
cat("  harm_label <- 'Identified Subgroup'\n")
cat("}\n")
cat("\n")
cat("if (est.scale == 'hr') {\n")
cat("  df_benefit <- subset(df_sg, treat.recommend == 1)   # Hc (complement)\n")
cat("  df_question <- subset(df_sg, treat.recommend == 0)  # H (harm)\n")
cat("  \n")
cat("  # FIXED: benefit is COMPLEMENT of sg.harm, question IS sg.harm\n")
cat("  benefit_name <- paste0('NOT (', harm_label, ')')    # Hc\n")
cat("  question_name <- harm_label                          # H\n")
cat("} else {\n")
cat("  # For 1/hr scale: roles are reversed\n")
cat("  df_benefit <- subset(df_sg, treat.recommend == 0)   # H is now benefit\n")
cat("  df_question <- subset(df_sg, treat.recommend == 1)  # Hc is now questionable\n")
cat("  \n")
cat("  benefit_name <- harm_label                          # H\n")
cat("  question_name <- paste0('NOT (', harm_label, ')')   # Hc\n")
cat("}\n")
cat("```\n\n")

# ==============================================================================
# TEST 6: Verify with actual HR computation
# ==============================================================================

cat("=================================================================\n")
cat("TEST 6: Verification with Actual HRs\n")
cat("=================================================================\n\n")

# Using our simulated data
cat("With simulated ForestSearch output (sg.harm = er<=27 & pgr<=41):\n\n")

# Current (BUGGY) output would show:
cat("CURRENT (BUGGY) FOREST PLOT OUTPUT:\n")
cat("------------------------------------\n")
cat(sprintf("  '{er<=27} & {pgr<=41}' (labeled 'Benefit')\n"))
cat(sprintf("    Full-analysis: HR = %.2f  <- WRONG: Hc HR labeled with H definition\n\n", hr_Hc[1]))
cat(sprintf("  'NOT ({er<=27} & {pgr<=41})' (labeled 'Questionable')\n"))
cat(sprintf("    Full-analysis: HR = %.2f  <- WRONG: H HR labeled with Hc definition\n\n", hr_H[1]))

# Correct output should show:
cat("CORRECT FOREST PLOT OUTPUT:\n")
cat("---------------------------\n")
cat(sprintf("  'NOT ({er<=27} & {pgr<=41})' (Benefit - Hc)\n"))
cat(sprintf("    Full-analysis: HR = %.2f  <- CORRECT: Hc HR with Hc definition\n\n", hr_Hc[1]))
cat(sprintf("  '{er<=27} & {pgr<=41}' (Questionable - H)\n"))
cat(sprintf("    Full-analysis: HR = %.2f  <- CORRECT: H HR with H definition\n\n", hr_H[1]))

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("=================================================================\n")
cat("SUMMARY OF ISSUES FOUND\n")
cat("=================================================================\n\n")

cat("ISSUE 1: Posthoc subgroups from subgroup_list are NOT added\n")
cat("  Location: Lines 267-278\n")
cat("  Problem:  Filter identifies posthoc entries but never adds them to plot\n")
cat("  Impact:   User-specified posthoc subgroups are silently ignored\n")
cat("  Fix:      Add loop to process each posthoc subgroup after separator\n\n")

cat("ISSUE 2: ForestSearch subgroup labels are REVERSED\n")
cat("  Location: Lines 293-312\n")
cat("  Problem:  benefit_name gets sg.harm (which defines H, not Hc)\n")
cat("            question_name gets NOT(sg.harm) (which defines Hc, not H)\n")
cat("  Impact:   Forest plot shows correct HRs but with SWAPPED labels\n")
cat("  Fix:      Swap the label assignments:\n")
cat("            - benefit_name = 'NOT (' + harm_label + ')' for est.scale='hr'\n")
cat("            - question_name = harm_label for est.scale='hr'\n\n")

cat("ISSUE 3: Variable naming is confusing\n")
cat("  Problem:  'sg1_label' suggests 'subgroup 1' but contains harm definition\n")
cat("  Fix:      Rename to 'harm_label' for clarity\n\n")

cat("=================================================================\n")
cat("END OF TEST\n")
cat("=================================================================\n")
