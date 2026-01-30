#' =============================================================================
#' GBSG Forest Plot Example
#' =============================================================================
#'
#' This script demonstrates how to create a subgroup results forest plot using
#' the German Breast Cancer Study Group (GBSG) dataset from the survival package.
#'
#' The GBSG dataset contains data from a randomized trial comparing hormonal
#' therapy vs. control in node-positive breast cancer patients.
#'
#' @author ForestSearch Development Team
#' @date 2025

# ==============================================================================
# SETUP AND DATA PREPARATION
# ==============================================================================

# Required packages
library(survival)
library(forestploter)
library(grid)
library(data.table)

# Load the GBSG dataset
data(gbsg, package = "survival")

# Examine the data structure
cat("GBSG Dataset Structure:\n")
cat("========================\n")
cat("N =", nrow(gbsg), "patients\n\n")
cat("Variables:\n")
str(gbsg)

# Key variables in gbsg:
# - pid: patient id
# - age: age in years
# - meno: menopausal status (0=premenopausal, 1=postmenopausal)
# - size: tumor size in mm
# - grade: tumor grade (1-3)
# - nodes: number of positive nodes
# - pgr: progesterone receptor (fmol/l)
# - er: estrogen receptor (fmol/l)
# - hormon: hormonal therapy (0=no, 1=yes) - THIS IS THE TREATMENT
# - rfstime: recurrence-free survival time (days)
# - status: censoring status (0=censored, 1=event)

# Prepare the analysis dataset
df_gbsg <- gbsg

# Create binary variables for subgroup definitions
df_gbsg$er_high <- ifelse(df_gbsg$er > median(df_gbsg$er, na.rm = TRUE), 1, 0)
df_gbsg$pgr_high <- ifelse(df_gbsg$pgr > median(df_gbsg$pgr, na.rm = TRUE), 1, 0)
df_gbsg$age_gt50 <- ifelse(df_gbsg$age > 50, 1, 0)
df_gbsg$size_large <- ifelse(df_gbsg$size > 25, 1, 0)
df_gbsg$nodes_high <- ifelse(df_gbsg$nodes > 3, 1, 0)
df_gbsg$grade_high <- ifelse(df_gbsg$grade == 3, 1, 0)

# Convert time to months for easier interpretation
df_gbsg$rfs_months <- df_gbsg$rfstime / 30.44

cat("\n\nPrepared dataset summary:\n")
cat("Treatment (hormon=1):", sum(df_gbsg$hormon), "patients\n
")
cat("Control (hormon=0):", sum(1 - df_gbsg$hormon), "patients\n")
cat("Events:", sum(df_gbsg$status), "\n")
cat("Median follow-up (months):", round(median(df_gbsg$rfs_months), 1), "\n")

# ==============================================================================
# HELPER FUNCTIONS (Simplified versions for this example)
# ==============================================================================

#' Compute Hazard Ratio for a Subgroup
#'
#' @param df Data frame with survival data
#' @param outcome.name Name of time variable
#' @param event.name Name of event indicator
#' @param treat.name Name of treatment variable
#' @return Named vector with HR, lower CI, upper CI, SE

compute_hr <- function(df, outcome.name, event.name, treat.name) {
  formula_str <- paste0("Surv(", outcome.name, ", ", event.name, ") ~ ", treat.name)
  fit <- survival::coxph(as.formula(formula_str), data = df)
  ci <- summary(fit)$conf.int[c(1, 3, 4)]
  se <- (ci[3] - ci[1]) / 1.96
  c(est = ci[1], low = ci[2], hi = ci[3], se = se)
}


#' Create Forest Plot Data Row
#'
#' @param df Data frame for analysis
#' @param sg_name Subgroup display name
#' @param outcome.name Name of time variable
#' @param event.name Name of event indicator
#' @param treat.name Name of treatment variable
#' @param E.name Label for treatment arm
#' @param C.name Label for control arm
#' @return Data frame row for forest plot

create_fp_row <- function(df, sg_name, outcome.name, event.name, treat.name,
                          E.name = "Hormon", C.name = "Control") {
  
  hr <- compute_hr(df, outcome.name, event.name, treat.name)
  n_treat <- sum(df[[treat.name]])
  n_control <- sum(1 - df[[treat.name]])
  
  data.frame(
    Subgroup = sg_name,
    n_E = n_treat,
    n_C = n_control,
    est = hr["est"],
    low = hr["low"],
    hi = hr["hi"],
    se = hr["se"],
    row.names = NULL
  )
}


#' Create Header/Separator Row
#'
#' @param text Display text for the row
#' @return Data frame row with NA estimates

create_header_row <- function(text) {
  data.frame(
    Subgroup = text,
    n_E = NA,
    n_C = NA,
    est = NA,
    low = NA,
    hi = NA,
    se = NA
  )
}


#' Generate Cross-Validation Sensitivity Text
#'
#' Simulates CV metrics for demonstration purposes
#'
#' @param cv_found Proportion of CV folds that found a subgroup
#' @param sens_benefit Sensitivity for benefit classification
#' @param sens_question Sensitivity for questionable classification
#' @return Formatted text string

generate_cv_text <- function(cv_found, sens_benefit, sens_question) {
  cv_text <- paste0("CV found = ", round(100 * cv_found, 0), "%")
  aa <- paste0(round(100 * sens_benefit, 0), "%,")
  bb <- paste0(round(100 * sens_question, 0), "%")
  paste(cv_text, ", Agreement(+,-) = ", aa, bb, sep = "")
}


# ==============================================================================
# ANALYSIS: Compute Hazard Ratios for Different Subgroups
# ==============================================================================

cat("\n\n")
cat("=================================================================\n")
cat("HAZARD RATIO ANALYSIS\n")
cat("=================================================================\n\n")

# Analysis parameters
outcome.name <- "rfs_months"
event.name <- "status"
treat.name <- "hormon"
E.name <- "Hormon"
C.name <- "Control"

# ITT Analysis
cat("ITT Analysis:\n")
hr_itt <- compute_hr(df_gbsg, outcome.name, event.name, treat.name)
cat(sprintf("  HR = %.2f (%.2f - %.2f)\n\n", hr_itt["est"], hr_itt["low"], hr_itt["hi"]))

# Reference Subgroups
cat("Reference Subgroups:\n")

# ER High
hr_er_high <- compute_hr(subset(df_gbsg, er_high == 1), outcome.name, event.name, treat.name)
cat(sprintf("  ER High (n=%d): HR = %.2f (%.2f - %.2f)\n", 
            sum(df_gbsg$er_high), hr_er_high["est"], hr_er_high["low"], hr_er_high["hi"]))

# ER Low
hr_er_low <- compute_hr(subset(df_gbsg, er_high == 0), outcome.name, event.name, treat.name)
cat(sprintf("  ER Low (n=%d): HR = %.2f (%.2f - %.2f)\n",
            sum(1 - df_gbsg$er_high), hr_er_low["est"], hr_er_low["low"], hr_er_low["hi"]))

# PgR High
hr_pgr_high <- compute_hr(subset(df_gbsg, pgr_high == 1), outcome.name, event.name, treat.name)
cat(sprintf("  PgR High (n=%d): HR = %.2f (%.2f - %.2f)\n",
            sum(df_gbsg$pgr_high), hr_pgr_high["est"], hr_pgr_high["low"], hr_pgr_high["hi"]))

# Age > 50
hr_age_gt50 <- compute_hr(subset(df_gbsg, age_gt50 == 1), outcome.name, event.name, treat.name)
cat(sprintf("  Age > 50 (n=%d): HR = %.2f (%.2f - %.2f)\n",
            sum(df_gbsg$age_gt50), hr_age_gt50["est"], hr_age_gt50["low"], hr_age_gt50["hi"]))

# Simulated "Identified" Subgroups (for demonstration)
# In a real analysis, these would come from ForestSearch

cat("\nPost-hoc Identified Subgroups (Simulated):\n")

# Subgroup 1: ER High & PgR High - might show enhanced benefit
df_sg1 <- subset(df_gbsg, er_high == 1 & pgr_high == 1)
if (nrow(df_sg1) > 20) {
  hr_sg1 <- compute_hr(df_sg1, outcome.name, event.name, treat.name)
  cat(sprintf("  ER High & PgR High (n=%d): HR = %.2f (%.2f - %.2f)\n",
              nrow(df_sg1), hr_sg1["est"], hr_sg1["low"], hr_sg1["hi"]))
}

# Complement of Subgroup 1
df_sg1_c <- subset(df_gbsg, !(er_high == 1 & pgr_high == 1))
if (nrow(df_sg1_c) > 20) {
  hr_sg1_c <- compute_hr(df_sg1_c, outcome.name, event.name, treat.name)
  cat(sprintf("  NOT (ER High & PgR High) (n=%d): HR = %.2f (%.2f - %.2f)\n",
              nrow(df_sg1_c), hr_sg1_c["est"], hr_sg1_c["low"], hr_sg1_c["hi"]))
}

# Subgroup 2: Postmenopausal & ER High
df_sg2 <- subset(df_gbsg, meno == 1 & er_high == 1)
if (nrow(df_sg2) > 20) {
  hr_sg2 <- compute_hr(df_sg2, outcome.name, event.name, treat.name)
  cat(sprintf("  Postmenopausal & ER High (n=%d): HR = %.2f (%.2f - %.2f)\n",
              nrow(df_sg2), hr_sg2["est"], hr_sg2["low"], hr_sg2["hi"]))
}

# Complement of Subgroup 2
df_sg2_c <- subset(df_gbsg, !(meno == 1 & er_high == 1))
if (nrow(df_sg2_c) > 20) {
  hr_sg2_c <- compute_hr(df_sg2_c, outcome.name, event.name, treat.name)
  cat(sprintf("  NOT (Postmeno & ER High) (n=%d): HR = %.2f (%.2f - %.2f)\n",
              nrow(df_sg2_c), hr_sg2_c["est"], hr_sg2_c["low"], hr_sg2_c["hi"]))
}


# ==============================================================================
# CREATE THE FOREST PLOT
# ==============================================================================

cat("\n\n")
cat("=================================================================\n")
cat("CREATING FOREST PLOT\n")
cat("=================================================================\n\n")

# Build the data frame for the forest plot
dt <- data.frame()

# Row 1: ITT
dt <- rbind(dt, create_fp_row(df_gbsg, "ITT (All Patients)", 
                               outcome.name, event.name, treat.name, E.name, C.name))

# Reference Subgroups
dt <- rbind(dt, create_fp_row(subset(df_gbsg, er_high == 1), "ER > median",
                               outcome.name, event.name, treat.name, E.name, C.name))
dt <- rbind(dt, create_fp_row(subset(df_gbsg, er_high == 0), "ER <= median",
                               outcome.name, event.name, treat.name, E.name, C.name))
dt <- rbind(dt, create_fp_row(subset(df_gbsg, pgr_high == 1), "PgR > median",
                               outcome.name, event.name, treat.name, E.name, C.name))
dt <- rbind(dt, create_fp_row(subset(df_gbsg, meno == 1), "Postmenopausal",
                               outcome.name, event.name, treat.name, E.name, C.name))

# Count rows before post-hoc section
n_ref <- nrow(dt)

# Separator for post-hoc subgroups
dt <- rbind(dt, create_header_row("Post-hoc Identified Subgroups"))
dt <- rbind(dt, create_header_row(" "))

# Post-hoc Subgroup 1: ER High & PgR High
sg1_name <- "ER High & PgR High"
dt <- rbind(dt, create_header_row(paste0("                   ", sg1_name)))
dt <- rbind(dt, create_fp_row(df_sg1, "  Full-analysis",
                               outcome.name, event.name, treat.name, E.name, C.name))
# Simulated bias-corrected (typically would shift toward null)
bc_row1 <- create_fp_row(df_sg1, "  Bias-corrected*",
                          outcome.name, event.name, treat.name, E.name, C.name)
bc_row1$est <- bc_row1$est * 1.05  # Shift slightly toward null for demo
bc_row1$low <- bc_row1$low * 1.03
bc_row1$hi <- bc_row1$hi * 1.07
dt <- rbind(dt, bc_row1)

# Complement
sg1_c_name <- paste0("NOT (", sg1_name, ")")
dt <- rbind(dt, create_header_row(sg1_c_name))
dt <- rbind(dt, create_fp_row(df_sg1_c, "  Full-analysis",
                               outcome.name, event.name, treat.name, E.name, C.name))

# Post-hoc Subgroup 2: Postmenopausal & ER High  
sg2_name <- "Postmenopausal & ER High"
dt <- rbind(dt, create_header_row(paste0("                   ", sg2_name)))
dt <- rbind(dt, create_fp_row(df_sg2, "  Full-analysis",
                               outcome.name, event.name, treat.name, E.name, C.name))
# Simulated bias-corrected
bc_row2 <- create_fp_row(df_sg2, "  Bias-corrected*",
                          outcome.name, event.name, treat.name, E.name, C.name)
bc_row2$est <- bc_row2$est * 1.05
bc_row2$low <- bc_row2$low * 1.03
bc_row2$hi <- bc_row2$hi * 1.07
dt <- rbind(dt, bc_row2)

# Complement
sg2_c_name <- paste0("NOT (", sg2_name, ")")
dt <- rbind(dt, create_header_row(sg2_c_name))
dt <- rbind(dt, create_fp_row(df_sg2_c, "  Full-analysis",
                               outcome.name, event.name, treat.name, E.name, C.name))

# Reset row names
rownames(dt) <- NULL

# Rename columns for display
names(dt)[names(dt) == "n_E"] <- E.name
names(dt)[names(dt) == "n_C"] <- C.name

# Convert n columns to character for display (handles NA in headers)
dt[[E.name]] <- ifelse(is.na(dt[[E.name]]), "", as.character(dt[[E.name]]))
dt[[C.name]] <- ifelse(is.na(dt[[C.name]]), "", as.character(dt[[C.name]]))

# Create color scheme
# ITT: yellow, Reference: powderblue, Post-hoc header: yellowgreen
# Post-hoc benefit: powderblue, Post-hoc complement: beige
n_rows <- nrow(dt)
sg_colors <- c(
  "yellow",                              # ITT
  rep("powderblue", n_ref - 1),          # Reference subgroups
  "yellowgreen",                         # Post-hoc header
  "white",                               # Blank
  "powderblue", "powderblue", "powderblue",  # SG1 benefit group
  "beige", "beige",                      # SG1 complement
  "powderblue", "powderblue", "powderblue",  # SG2 benefit group
  "beige", "beige"                       # SG2 complement
)

# Extend colors if needed
if (length(sg_colors) < n_rows) {
  sg_colors <- c(sg_colors, rep("white", n_rows - length(sg_colors)))
}
sg_colors <- sg_colors[1:n_rows]

# Create theme
tm <- forest_theme(
  core = list(
    fg_params = list(hjust = 1, x = 0.9),
    bg_params = list(fill = sg_colors)
  ),
  colhead = list(fg_params = list(hjust = 0.5, x = 0.5)),
  footnote_gp = gpar(cex = 0.65, fontface = "italic", col = "darkcyan")
)

# Add spacing column and HR text
dt$` ` <- paste(rep(" ", 20), collapse = " ")
dt$`HR (95% CI)` <- ifelse(
  is.na(dt$se), "",
  sprintf("%.2f (%.2f to %.2f)", dt$est, dt$low, dt$hi)
)

# Generate the forest plot
title_text <- "GBSG: Subgroup Treatment Effects (Hormonal Therapy vs Control)"
arrow_text <- c("Favors Hormonal Therapy", "Favors Control")
footnote_text <- "* Bias-corrected estimates adjust for post-hoc selection; Example demonstration only"

p <- forest(
  dt[, c("Subgroup", E.name, C.name, " ", "HR (95% CI)")],
  title = title_text,
  est = dt$est,
  lower = dt$low,
  upper = dt$hi,
  sizes = 0.4,
  ci_column = 4,
  ref_line = 1,
  arrow_lab = arrow_text,
  xlim = c(0.4, 2.0),
  ticks_at = c(0.5, 0.75, 1.0, 1.5, 2.0),
  footnote = footnote_text,
  theme = tm
)

# Add CV metrics text (simulated for demonstration)
cv_text_sg1 <- generate_cv_text(0.85, 0.78, 0.72)
cv_text_sg2 <- generate_cv_text(0.75, 0.70, 0.68)

# Insert CV text after each subgroup block
# Row positions need to account for the plot structure
g <- p
g <- insert_text(g, text = cv_text_sg1, row = 13, just = "left",
                 gp = gpar(cex = 0.75, col = "darkblue", fontface = "italic"))
g <- insert_text(g, text = cv_text_sg2, row = 18, just = "left",
                 gp = gpar(cex = 0.75, col = "darkblue", fontface = "italic"))

# Display the plot
cat("Generating forest plot...\n")
plot(g)

cat("\nForest plot created successfully!\n")

# ==============================================================================
# SAVE THE PLOT (Optional)
# ==============================================================================

# Save as PNG
# png("gbsg_subgroup_forestplot.png", width = 12, height = 10, units = "in", res = 300)
# plot(g)
# dev.off()

# Save as PDF
# pdf("gbsg_subgroup_forestplot.pdf", width = 12, height = 10)
# plot(g)
# dev.off()


# ==============================================================================
# SUMMARY DATA TABLE
# ==============================================================================

cat("\n\n")
cat("=================================================================\n")
cat("SUMMARY TABLE\n")
cat("=================================================================\n\n")

# Print the data used for the forest plot
summary_dt <- dt[, c("Subgroup", E.name, C.name, "est", "low", "hi")]
names(summary_dt) <- c("Subgroup", "N(Hormon)", "N(Control)", "HR", "Lower", "Upper")
summary_dt$HR <- round(summary_dt$HR, 3)
summary_dt$Lower <- round(summary_dt$Lower, 3)
summary_dt$Upper <- round(summary_dt$Upper, 3)
print(summary_dt, row.names = FALSE)


# ==============================================================================
# INTEGRATION WITH FORESTSEARCH (Template)
# ==============================================================================

cat("\n\n")
cat("=================================================================\n")
cat("TEMPLATE: Integration with ForestSearch Package\n")
cat("=================================================================\n\n")

cat("
# To use with actual ForestSearch results, replace the simulated subgroups with:

# library(forestsearch)

# Step 1: Prepare confounders
# confounders.name <- c('er', 'pgr', 'age', 'size', 'nodes', 'grade', 'meno')

# Step 2: Run ForestSearch
# fs.est <- forestsearch(
#   df.analysis = df_gbsg,
#   outcome.name = 'rfs_months',
#   event.name = 'status',
#   treat.name = 'hormon',
#   id.name = 'pid',
#   confounders.name = confounders.name,
#   est.scale = 'hr',
#   hr.threshold = 1.25,
#   hr.consistency = 1.0,
#   pconsistency.threshold = 0.80,
#   n.min = 50,
#   maxk = 2,
#   details = TRUE
# )

# Step 3: Run Bootstrap Bias Correction
# fs_bc <- forestsearch_bootstrap_dofuture(
#   fs.est = fs.est,
#   nb_boots = 500,
#   parallel_args = list(plan = 'multisession', workers = 4)
# )

# Step 4: Run K-fold Cross-Validation
# fs_kfold <- forestsearch_tenfold(
#   fs.est = fs.est,
#   sims = 5,
#   Kfolds = 10
# )

# Step 5: Create Forest Plot using wrapper
# source('plot_subgroup_results_forestplot.R')
# 
# result <- plot_subgroup_results_forestplot(
#   fs_results = list(
#     fs.est = fs.est,
#     fs_bc = fs_bc,
#     fs_kfold = fs_kfold
#   ),
#   df_analysis = df_gbsg,
#   outcome.name = 'rfs_months',
#   event.name = 'status',
#   treat.name = 'hormon',
#   E.name = 'Hormonal',
#   C.name = 'Control',
#   title_text = 'GBSG: Identified Subgroups'
# )
# 
# plot(result$plot)
")

cat("\n\nExample complete!\n")
