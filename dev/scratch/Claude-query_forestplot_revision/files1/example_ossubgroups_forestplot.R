#' Example: Creating the OSsubgroups Forest Plot
#'
#' This script demonstrates how to recreate the forest plot from the
#' "OSsubgroups" chunk in SubgroupIdentification_Summary_v0details.Rmd
#' using the ForestSearch package and the new wrapper functions.
#'
#' @author ForestSearch Development Team
#' @date 2025

# ==============================================================================
# SETUP
# ==============================================================================

# Load required packages
library(survival)
library(forestploter)
library(grid)
library(data.table)

# Source the wrapper functions
# source("plot_subgroup_results_forestplot.R")

# ==============================================================================
# LEGACY HELPER FUNCTIONS (from subgroup_results_summary_v0.R)
# These are provided for backward compatibility with existing analysis workflows
# ==============================================================================

#' Create HR Table Row for Subgroup Analysis
#'
#' Creates a data frame row with hazard ratio estimates for forest plot.
#' This is the legacy SG_HRtable2 function from the original codebase.
#'
#' @param dfa Data frame for analysis (NULL if using fs_bc estimates)
#' @param df.OOB Data frame for out-of-bag validation (optional)
#' @param fa_SG Full-analysis subgroup estimates from bootstrap (optional)
#' @param ntreat_fa Number of treated in full analysis
#' @param ncontrol_fa Number of control in full analysis
#' @param outcome.name Name of survival time variable
#' @param event.name Name of event indicator variable
#' @param treat.name Name of treatment variable
#' @param sg_name Display name for the subgroup
#' @param E.name Label for experimental arm (default: "E")
#' @param C.name Label for control arm (default: "C")
#'
#' @return Data frame row(s) with HR estimates suitable for forest plot
#' @export

SG_HRtable2 <- function(dfa = NULL, df.OOB = NULL, fa_SG = NULL, 
                        ntreat_fa = NA, ncontrol_fa = NA,
                        outcome.name, event.name, treat.name,
                        sg_name, E.name = "E", C.name = "C") {
  
  # Case 1: Direct analysis on dfa (no bootstrap)
  if (is.null(df.OOB) & is.null(fa_SG)) {
    sf <- paste0("Surv(", outcome.name, ",", event.name, ") ~ ", treat.name)
    cox.formula <- as.formula(sf)
    
    hr <- summary(survival::coxph(cox.formula, data = dfa))$conf.int[c(1, 3, 4)]
    
    ntreat <- sum(dfa[, treat.name])
    ncontrol <- sum(1 - dfa[, treat.name])
    
    est <- hr[1]
    low <- hr[2]
    hi <- hr[3]
    se <- (hi - est) / 1.96
    
    hr_vec <- c(ntreat, ncontrol, est, low, hi, se)
    names(hr_vec) <- c(E.name, C.name, "est", "low", "hi", "se")
    aa <- as.data.frame(t(hr_vec))
    aa$Subgroup <- sg_name
    aa <- aa[, c("Subgroup", E.name, C.name, "est", "low", "hi", "se")]
    return(aa)
  }
  
  # Case 2: Bias-corrected estimates only (no OOB)
  if (is.null(df.OOB) & !is.null(fa_SG)) {
    resSG <- fa_SG[, c("H2", "H2_lower", "H2_upper")]
    Hstat <- round(unlist(resSG[1, ]), 3)
    
    est <- Hstat[1]
    low <- Hstat[2]
    hi <- Hstat[3]
    se <- (hi - est) / 1.96
    
    hr_SG <- c(ntreat_fa, ncontrol_fa, est, low, hi, se)
    names(hr_SG) <- c(E.name, C.name, "est", "low", "hi", "se")
    aa <- as.data.frame(t(hr_SG))
    aa$Subgroup <- sg_name
    aa <- aa[, c("Subgroup", E.name, C.name, "est", "low", "hi", "se")]
    return(aa)
  }
  
  # Case 3: Both bootstrap and OOB estimates
  if (!is.null(df.OOB) & !is.null(fa_SG)) {
    ans <- list()
    
    # Full-analysis bias-corrected
    resSG <- fa_SG[, c("H2", "H2_lower", "H2_upper")]
    Hstat <- round(unlist(resSG[1, ]), 3)
    
    est <- Hstat[1]
    low <- Hstat[2]
    hi <- Hstat[3]
    se <- (hi - est) / 1.96
    
    hr_SG <- c(ntreat_fa, ncontrol_fa, est, low, hi, se)
    names(hr_SG) <- c(E.name, C.name, "est", "low", "hi", "se")
    aa <- as.data.frame(t(hr_SG))
    aa$Subgroup <- "full-Analysis"
    aa <- aa[, c("Subgroup", E.name, C.name, "est", "low", "hi", "se")]
    ans$fa <- aa
    
    # OOB estimates
    sf <- paste0("Surv(", outcome.name, ",", event.name, ") ~ ", treat.name)
    cox.formula <- as.formula(sf)
    
    hr <- summary(survival::coxph(cox.formula, data = df.OOB))$conf.int[c(1, 3, 4)]
    ntreat <- sum(df.OOB[, treat.name])
    ncontrol <- sum(1 - df.OOB[, treat.name])
    
    est <- hr[1]
    low <- hr[2]
    hi <- hr[3]
    se <- (hi - est) / 1.96
    
    hr_vec <- c(ntreat, ncontrol, est, low, hi, se)
    names(hr_vec) <- c(E.name, C.name, "est", "low", "hi", "se")
    bb <- as.data.frame(t(hr_vec))
    bb$Subgroup <- "N-fold"
    bb <- bb[, c("Subgroup", E.name, C.name, "est", "low", "hi", "se")]
    ans$oob <- bb
    
    return(ans)
  }
}


#' Create Data Frame for Subgroup Forest Plot
#'
#' Creates formatted data for forest plot from ForestSearch results with
#' bootstrap bias correction and out-of-bag validation.
#' This is the legacy df_sgforest function.
#'
#' @param fs_OOB Out-of-bag cross-validation results data frame
#' @param df_sg Data frame with treat.recommend assignments
#' @param fs_bc Bootstrap bias-corrected results
#' @param outcome.name Name of survival time variable
#' @param treat.name Name of treatment variable
#' @param sg1_name Name for the benefitting subgroup
#' @param sg0_name Name for the complement subgroup
#' @param E.name Label for experimental arm
#' @param C.name Label for control arm
#' @param est.scale Estimate scale: "hr" or "1/hr"
#'
#' @return List with res_sg1 and res_sg0 data frames
#' @export

df_sgforest <- function(fs_OOB, df_sg, fs_bc, outcome.name, treat.name,
                        sg1_name, sg0_name, E.name, C.name, est.scale = "hr") {
  ans <- list()
  
  # Get event.name - assume it follows naming convention
  event.name <- gsub("_time", "_event", outcome.name)
  if (!event.name %in% names(df_sg)) {
    event.name <- gsub("time", "event", outcome.name)
  }
  
  if (est.scale == "hr") {
    # Benefitting subgroup (Hc) - treat.recommend == 1
    df.OOB <- subset(fs_OOB, treat.recommend == 1)
    ntreat_fa <- nrow(subset(df_sg, treat.recommend == 1 & treat == 1))
    ncontrol_fa <- nrow(subset(df_sg, treat.recommend == 1 & treat == 0))
    
    temp <- SG_HRtable2(
      df.OOB = df.OOB, fa_SG = fs_bc$Hc_estimates,
      ntreat_fa = ntreat_fa, ncontrol_fa = ncontrol_fa,
      outcome.name = outcome.name, event.name = event.name, treat.name = treat.name,
      sg_name = sg1_name, E.name = E.name, C.name = C.name
    )
    
    bb <- temp$fa
    cc <- temp$oob
    
    # Add header row
    aa <- bb
    aa$Subgroup <- sg1_name
    aa[, c(E.name, C.name)] <- ""
    aa[, c("est", "low", "hi", "se")] <- NA
    ans$res_sg1 <- rbind(aa, bb, cc)
    
    # Complement (H) - treat.recommend == 0
    df.OOB <- subset(fs_OOB, treat.recommend == 0)
    ntreat_fa <- nrow(subset(df_sg, treat.recommend == 0 & treat == 1))
    ncontrol_fa <- nrow(subset(df_sg, treat.recommend == 0 & treat == 0))
    
    temp <- SG_HRtable2(
      df.OOB = df.OOB, fa_SG = fs_bc$H_estimates,
      ntreat_fa = ntreat_fa, ncontrol_fa = ncontrol_fa,
      outcome.name = outcome.name, event.name = event.name, treat.name = treat.name,
      sg_name = sg0_name, E.name = E.name, C.name = C.name
    )
    
    bb <- temp$fa
    cc <- temp$oob
    
    # Add header row
    aa <- bb
    aa$Subgroup <- sg0_name
    aa[, c(E.name, C.name)] <- ""
    aa[, c("est", "low", "hi", "se")] <- NA
    ans$res_sg0 <- rbind(aa, bb, cc)
    
  } else if (est.scale == "1/hr") {
    # Reverse roles: H becomes benefitting, Hc becomes complement
    
    # Benefitting subgroup (H)
    df.OOB <- subset(fs_OOB, treat.recommend == 0)
    ntreat_fa <- nrow(subset(df_sg, treat.recommend == 0 & treat == 1))
    ncontrol_fa <- nrow(subset(df_sg, treat.recommend == 0 & treat == 0))
    
    temp <- SG_HRtable2(
      df.OOB = df.OOB, fa_SG = fs_bc$H_estimates,
      ntreat_fa = ntreat_fa, ncontrol_fa = ncontrol_fa,
      outcome.name = outcome.name, event.name = event.name, treat.name = "treat2",
      sg_name = sg1_name, E.name = E.name, C.name = C.name
    )
    
    bb <- temp$fa
    cc <- temp$oob
    
    aa <- bb
    aa$Subgroup <- sg1_name
    aa[, c(E.name, C.name)] <- ""
    aa[, c("est", "low", "hi", "se")] <- NA
    ans$res_sg1 <- rbind(aa, bb, cc)
    
    # Complement (Hc)
    df.OOB <- subset(fs_OOB, treat.recommend == 1)
    ntreat_fa <- nrow(subset(df_sg, treat.recommend == 1 & treat == 1))
    ncontrol_fa <- nrow(subset(df_sg, treat.recommend == 1 & treat == 0))
    
    temp <- SG_HRtable2(
      df.OOB = df.OOB, fa_SG = fs_bc$Hc_estimates,
      ntreat_fa = ntreat_fa, ncontrol_fa = ncontrol_fa,
      outcome.name = outcome.name, event.name = event.name, treat.name = "treat2",
      sg_name = sg0_name, E.name = E.name, C.name = C.name
    )
    
    bb <- temp$fa
    cc <- temp$oob
    
    aa <- bb
    aa$Subgroup <- sg0_name
    aa[, c(E.name, C.name)] <- ""
    aa[, c("est", "low", "hi", "se")] <- NA
    ans$res_sg0 <- rbind(aa, bb, cc)
  }
  
  return(ans)
}


#' Generate CV Sensitivity Text
#'
#' Creates formatted text for cross-validation agreement metrics.
#' This is the legacy sens_text function.
#'
#' @param fs_kfold K-fold cross-validation results
#' @param est.scale Estimate scale: "hr" or "1/hr"
#'
#' @return Character string with formatted metrics
#' @export

sens_text <- function(fs_kfold, est.scale = "hr") {
  cv <- fs_kfold$find_summary["Any"]
  
  if (est.scale == "hr") {
    Q <- fs_kfold$sens_summary["sens_H"]
    B <- fs_kfold$sens_summary["sens_Hc"]
  } else {
    Q <- fs_kfold$sens_summary["sens_Hc"]
    B <- fs_kfold$sens_summary["sens_H"]
  }
  
  cv_text <- paste0("CV found = ", round(100 * cv, 0), "%")
  aa <- paste0(round(100 * B, 0), "%,")
  bb <- paste0(round(100 * Q, 0), "%")
  sense_text <- paste("Agreement(+,-) = ", aa, bb, collapse = ",")
  sg_text <- paste(cv_text, sense_text, sep = ", ")
  
  return(sg_text)
}


# ==============================================================================
# EXAMPLE: Recreating the OSsubgroups Figure
# ==============================================================================

#' Create OSsubgroups Forest Plot
#'
#' This function recreates the forest plot from the OSsubgroups chunk,
#' showing multiple subgroup analyses with cross-validation metrics.
#'
#' @param df_itt ITT data frame with all analysis variables
#' @param fs_results_list Named list of ForestSearch result sets, where each 
#'   element is a list containing:
#'   \itemize{
#'     \item \code{fs.est}: ForestSearch estimation object
#'     \item \code{fs_bc}: Bootstrap bias-corrected results
#'     \item \code{fs_OOB}: Out-of-bag results
#'     \item \code{fs_kfold}: K-fold CV results
#'     \item \code{sg1_name}: Benefitting subgroup name
#'     \item \code{sg0_name}: Complement subgroup name
#'   }
#' @param reference_subgroups Named list of reference subgroup definitions
#' @param outcome.name Name of survival time variable
#' @param event.name Name of event indicator variable
#' @param treat.name Name of treatment variable
#' @param E.name Label for experimental arm
#' @param C.name Label for control arm
#' @param title_text Plot title
#' @param arrow_text Arrow labels
#' @param footnote_text Footnote text
#' @param xlim X-axis limits
#' @param ticks_at X-axis tick positions
#'
#' @return forestploter grob object
#' @export

create_ossubgroups_forestplot <- function(
    df_itt,
    fs_results_list,
    reference_subgroups = NULL,
    outcome.name,
    event.name,
    treat.name,
    E.name = "E",
    C.name = "C",
    title_text = "OS Subgroups Identified (Post-Hoc)",
    arrow_text = c("favors Treatment", "Control"),
    footnote_text = NULL,
    xlim = c(0.25, 1.5),
    ticks_at = c(0.25, 0.70, 1.0, 1.5)
) {
  
  # ==========================================================================
  # Build ITT Row
  # ==========================================================================
  
  res_itt <- SG_HRtable2(
    dfa = df_itt, df.OOB = NULL, fa_SG = NULL,
    ntreat_fa = NA, ncontrol_fa = NA,
    outcome.name = outcome.name, event.name = event.name, treat.name = treat.name,
    sg_name = "ITT", E.name = E.name, C.name = C.name
  )
  
  dt <- res_itt
  
  # ==========================================================================
  # Add Reference Subgroups (if provided)
  # ==========================================================================
  
  if (!is.null(reference_subgroups)) {
    for (ref_name in names(reference_subgroups)) {
      ref <- reference_subgroups[[ref_name]]
      dfa_ref <- subset(df_itt, eval(parse(text = ref$subset_expr)))
      
      if (nrow(dfa_ref) > 10) {
        ref_row <- SG_HRtable2(
          dfa = dfa_ref, df.OOB = NULL, fa_SG = NULL,
          ntreat_fa = NA, ncontrol_fa = NA,
          outcome.name = outcome.name, event.name = event.name, treat.name = treat.name,
          sg_name = ref$name, E.name = E.name, C.name = C.name
        )
        dt <- rbind(dt, ref_row)
      }
    }
  }
  
  # Count reference/fixed rows for coloring
  kfix <- nrow(dt)
  
  # ==========================================================================
  # Add Post-hoc Header
  # ==========================================================================
  
  zz <- res_itt
  zz$Subgroup <- "Post-hoc subgroups"
  zz[, c(E.name, C.name)] <- ""
  zz[, c("est", "low", "hi", "se")] <- NA
  dt <- rbind(dt, zz)
  
  # Blank row
  zzz <- res_itt
  zzz$Subgroup <- " "
  zzz[, c(E.name, C.name)] <- ""
  zzz[, c("est", "low", "hi", "se")] <- NA
  
  # ==========================================================================
  # Add Post-hoc Subgroup Results
  # ==========================================================================
  
  sg_texts <- list()
  sg_rows <- list()
  kposthoc <- 0
  
  for (sg_key in names(fs_results_list)) {
    fs_res <- fs_results_list[[sg_key]]
    
    if (!is.null(fs_res$fs.est) && !is.null(fs_res$fs_bc)) {
      
      df_sg <- fs_res$fs.est$df.est
      est.scale <- fs_res$fs.est$est.scale
      
      # Handle treatment assignment based on scale
      if (est.scale == "1/hr") {
        df_sg$treat <- 1 - df_sg[, treat.name]
      } else {
        df_sg$treat <- df_sg[, treat.name]
      }
      
      # Get OOB data if available
      fs_OOB <- fs_res$fs_OOB
      if (!is.null(fs_OOB)) {
        fs_OOB <- fs_OOB$resCV
      }
      
      # Create subgroup forest data
      if (!is.null(fs_OOB)) {
        dfForest_sg <- df_sgforest(
          fs_OOB = fs_OOB, df_sg = df_sg, fs_bc = fs_res$fs_bc,
          outcome.name = outcome.name, treat.name = treat.name,
          sg1_name = fs_res$sg1_name, sg0_name = fs_res$sg0_name,
          E.name = E.name, C.name = C.name, est.scale = est.scale
        )
        
        res_sg <- rbind(dfForest_sg$res_sg1, dfForest_sg$res_sg0)
      } else {
        # Fallback without OOB
        res_sg <- create_posthoc_rows_simple(
          df_sg = df_sg, fs_bc = fs_res$fs_bc,
          sg1_name = fs_res$sg1_name, sg0_name = fs_res$sg0_name,
          outcome.name = outcome.name, event.name = event.name,
          treat.name = treat.name, E.name = E.name, C.name = C.name,
          est.scale = est.scale
        )
      }
      
      dt <- rbind(dt, res_sg)
      kposthoc <- kposthoc + 1
      
      # Get CV text if available
      if (!is.null(fs_res$fs_kfold)) {
        sg_texts[[sg_key]] <- sens_text(fs_res$fs_kfold, est.scale)
      }
    }
  }
  
  # ==========================================================================
  # Create Color Scheme
  # ==========================================================================
  
  # Color logic from original:
  # - ITT: yellow
  # - Reference subgroups: powderblue
  # - Post-hoc header: yellowgreen
  # - Post-hoc subgroups: alternating powderblue/beige groups
  
  n_rows <- nrow(dt)
  sg_colors <- rep("white", n_rows)
  
  # ITT row
  sg_colors[1] <- "yellow"
  
  # Reference rows
  if (kfix > 1) {
    sg_colors[2:kfix] <- "powderblue"
  }
  
  # Post-hoc header
  sg_colors[kfix + 1] <- "yellowgreen"
  
  # Post-hoc subgroup rows (alternating powderblue/beige by group)
  if (kposthoc > 0) {
    posthoc_start <- kfix + 2
    # Each posthoc group has 6 rows: header + FA + OOB for benefit and complement
    for (i in seq_len(kposthoc)) {
      group_start <- posthoc_start + (i - 1) * 6
      group_end <- min(group_start + 5, n_rows)
      
      if (group_start <= n_rows) {
        # First 3 rows (benefit subgroup) - powderblue
        for (j in group_start:min(group_start + 2, n_rows)) {
          sg_colors[j] <- "powderblue"
        }
        # Next 3 rows (complement) - beige
        for (j in (group_start + 3):min(group_end, n_rows)) {
          sg_colors[j] <- "beige"
        }
      }
    }
  }
  
  # ==========================================================================
  # Create Forest Plot
  # ==========================================================================
  
  tm <- forestploter::forest_theme(
    core = list(
      fg_params = list(hjust = 1, x = 0.9),
      bg_params = list(fill = sg_colors)
    ),
    colhead = list(fg_params = list(hjust = 0.5, x = 0.5)),
    footnote_gp = grid::gpar(cex = 0.65, fontface = "italic", col = "darkcyan")
  )
  
  dt$` ` <- paste(rep(" ", 20), collapse = " ")
  
  dt$`HR (95% CI)` <- ifelse(
    is.na(dt$se), "",
    sprintf("%.2f (%.2f to %.2f)", dt$est, dt$low, dt$hi)
  )
  
  p <- forestploter::forest(
    dt[, c("Subgroup", E.name, C.name, " ", "HR (95% CI)")],
    title = title_text,
    est = dt$est,
    lower = dt$low,
    upper = dt$hi,
    sizes = 0.4,
    ci_column = 4,
    ref_line = 1,
    arrow_lab = arrow_text,
    xlim = xlim,
    ticks_at = ticks_at,
    footnote = footnote_text,
    theme = tm
  )
  
  # ==========================================================================
  # Add CV Metrics Text
  # ==========================================================================
  
  g <- p
  
  if (length(sg_texts) > 0) {
    # Calculate row positions for each CV text
    # Each posthoc group occupies 6 rows, text goes after each group
    posthoc_start <- kfix + 2
    
    for (i in seq_along(sg_texts)) {
      row_loc <- posthoc_start + (i - 1) * 6 + 6 + 1
      
      if (row_loc <= n_rows + length(sg_texts)) {
        g <- forestploter::insert_text(
          g,
          text = sg_texts[[i]],
          row = row_loc,
          just = "left",
          gp = grid::gpar(cex = 0.8, col = "black", fontface = "italic")
        )
      }
    }
  }
  
  return(g)
}


#' Create Post-hoc Rows Without OOB (Simplified)
#'
#' Helper function to create post-hoc subgroup rows when OOB data is unavailable.
#'
#' @keywords internal

create_posthoc_rows_simple <- function(df_sg, fs_bc, sg1_name, sg0_name,
                                       outcome.name, event.name, treat.name,
                                       E.name, C.name, est.scale) {
  
  if (est.scale == "hr") {
    df_benefit <- subset(df_sg, treat.recommend == 1)
    df_complement <- subset(df_sg, treat.recommend == 0)
    bc_benefit <- fs_bc$Hc_estimates
    bc_complement <- fs_bc$H_estimates
  } else {
    df_benefit <- subset(df_sg, treat.recommend == 0)
    df_complement <- subset(df_sg, treat.recommend == 1)
    bc_benefit <- fs_bc$H_estimates
    bc_complement <- fs_bc$Hc_estimates
  }
  
  # Benefit header
  header_benefit <- data.frame(
    Subgroup = sg1_name,
    stringsAsFactors = FALSE
  )
  header_benefit[[E.name]] <- ""
  header_benefit[[C.name]] <- ""
  header_benefit$est <- NA
  header_benefit$low <- NA
  header_benefit$hi <- NA
  header_benefit$se <- NA
  
  # Benefit row
  row_benefit <- SG_HRtable2(
    dfa = NULL, df.OOB = NULL, fa_SG = bc_benefit,
    ntreat_fa = sum(df_benefit$treat),
    ncontrol_fa = sum(1 - df_benefit$treat),
    outcome.name = outcome.name, event.name = event.name, treat.name = treat.name,
    sg_name = "  full-Analysis", E.name = E.name, C.name = C.name
  )
  
  # Complement header
  header_complement <- data.frame(
    Subgroup = sg0_name,
    stringsAsFactors = FALSE
  )
  header_complement[[E.name]] <- ""
  header_complement[[C.name]] <- ""
  header_complement$est <- NA
  header_complement$low <- NA
  header_complement$hi <- NA
  header_complement$se <- NA
  
  # Complement row
  row_complement <- SG_HRtable2(
    dfa = NULL, df.OOB = NULL, fa_SG = bc_complement,
    ntreat_fa = sum(df_complement$treat),
    ncontrol_fa = sum(1 - df_complement$treat),
    outcome.name = outcome.name, event.name = event.name, treat.name = treat.name,
    sg_name = "  full-Analysis", E.name = E.name, C.name = C.name
  )
  
  rbind(header_benefit, row_benefit, header_complement, row_complement)
}


# ==============================================================================
# USAGE EXAMPLE (commented out - requires actual data)
# ==============================================================================

# # Load your data and ForestSearch results
# load("df_ITT.Rdata")
# 
# # Define reference subgroups for context
# reference_subgroups <- list(
#   cps1 = list(subset_expr = "cps > 1", name = "CPS > 1"),
#   cps3 = list(subset_expr = "cps > 3", name = "CPS > 3"),
#   ts19 = list(subset_expr = "tmrsize > 19", name = "Tumor Size > 19"),
#   ts62 = list(subset_expr = "tmrsize > 62", name = "Tumor Size > 62")
# )
# 
# # Load ForestSearch results for each identified subgroup
# load("sg1_v0b.Rdata")
# load("sg1_boots=1000_v0b.Rdata")
# load("sg1_OOB_v0b.Rdata")
# load("sg1_20fold_v0b.Rdata")
# 
# fs_results_list <- list(
#   sg1 = list(
#     fs.est = fs.est,
#     fs_bc = fs_bc,
#     fs_OOB = fs_OOB,
#     fs_kfold = fs_kfold,
#     sg1_name = "CPS > 3 or tumor size > 62",
#     sg0_name = "CPS <= 3 & tumor size <= 62"
#   )
# )
# 
# # Create the forest plot
# g <- create_ossubgroups_forestplot(
#   df_itt = df_itt,
#   fs_results_list = fs_results_list,
#   reference_subgroups = reference_subgroups,
#   outcome.name = "os_time",
#   event.name = "os_event",
#   treat.name = "combo",
#   E.name = "Pembro+CT",
#   C.name = "CT",
#   title_text = "OS Subgroups Identified (Post-Hoc)",
#   arrow_text = c("favors Pembro+CT", "CT"),
#   footnote_text = "CV metrics: % found, Agreement(+,-)"
# )
# 
# # Display the plot
# plot(g)
