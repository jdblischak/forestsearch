# =============================================================================
# ForestSearch S3 Summary Methods
# =============================================================================
#
# CRAN COMPLIANCE NOTE:
# All summary methods for ForestSearch classes are defined in this single file.
# This avoids duplicate method registration errors.
#
# =============================================================================

#' Summary Method for forestsearch Objects
#'
#' Provides detailed summary of ForestSearch analysis including parameters,
#' variable selection results, consistency evaluation, and selected subgroup.
#'
#' @param object A forestsearch object returned by \code{\link{forestsearch}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' \dontrun{
#' result <- forestsearch(df.analysis = trial_data)
#' summary(result)
#' }
#'
#' @export
summary.forestsearch <- function(object, ...) {
  cat("ForestSearch Summary\n")
  cat("====================\n\n")

  # ---------------------------------------------------------------------------
  # Analysis Parameters
  # ---------------------------------------------------------------------------
  cat("Analysis Parameters:\n")

  # Try different locations for parameters
  params <- object$args_call_all %||% object$parameters

  if (!is.null(params)) {
    .print_param(params, "sg_focus")
    .print_param(params, "hr.threshold")
    .print_param(params, "hr.consistency")
    .print_param(params, "pconsistency.threshold")
    .print_param(params, "n.min")
    .print_param(params, "fs.splits", alt_name = "n.splits")
    .print_param(params, "maxk")
    .print_param(params, "use_twostage")
  }
  cat("\n")

  # ---------------------------------------------------------------------------
  # Variable Selection
  # ---------------------------------------------------------------------------
  cat("Variable Selection:\n")

  if (!is.null(object$confounders.candidate)) {
    cat("  Candidate confounders:", length(object$confounders.candidate), "\n")
  }
  if (!is.null(object$confounders.evaluated)) {
    cat("  Confounders evaluated:", length(object$confounders.evaluated), "\n")
  }
  if (!is.null(object$grf_importance)) {
    cat("  GRF variables:", length(object$grf_importance), "\n")
  }
  if (!is.null(object$lasso_selected)) {
    cat("  LASSO selected:", length(object$lasso_selected), "\n")
  }
  cat("\n")

  # ---------------------------------------------------------------------------
  # Consistency Evaluation
  # ---------------------------------------------------------------------------
  if (!is.null(object$grp.consistency)) {
    cat("Consistency Evaluation:\n")
    gc <- object$grp.consistency

    # Algorithm used
    if (!is.null(gc$algorithm)) {
      cat("  Algorithm:", gc$algorithm, "\n")
    }

    # Candidates evaluated
    if (!is.null(gc$n_candidates_evaluated)) {
      cat("  Candidates evaluated:", gc$n_candidates_evaluated, "\n")
      cat("  Candidates passed:", gc$n_passed, "\n")
    } else if (!is.null(gc$all_candidates)) {
      n_cand <- nrow(gc$all_candidates)
      cat("  Candidates evaluated:", n_cand, "\n")

      # Count meeting thresholds
      if (!is.null(params)) {
        hr_thresh <- params$hr.threshold %||% 1.0
        pcons_thresh <- params$pconsistency.threshold %||% 0.9

        n_hr <- sum(gc$all_candidates$hr > hr_thresh, na.rm = TRUE)
        n_pcons <- sum(gc$all_candidates$Pcons >= pcons_thresh, na.rm = TRUE)

        cat("  Meeting hr.threshold:", n_hr, "\n")
        cat("  Meeting consistency:", n_pcons, "\n")
      }
    }
    cat("\n")
  }

  # ---------------------------------------------------------------------------
  # Results by sg_focus
  # ---------------------------------------------------------------------------
  .print_results_by_focus(object)

  # ---------------------------------------------------------------------------
  # Selected Subgroup
  # ---------------------------------------------------------------------------
  if (!is.null(object$sg.harm)) {
    cat("\nSelected Subgroup:\n")
    cat("  Definition:", paste(object$sg.harm, collapse = " & "), "\n")

    if (!is.null(object$grp.consistency$out_sg$result)) {
      top <- object$grp.consistency$out_sg$result[1L, ]
      if (!is.null(top$N)) cat("  Sample size:", top$N, "\n")
      if (!is.null(top$hr)) cat("  Hazard ratio:", round(as.numeric(top$hr), 3), "\n")
      if (!is.null(top$Pcons)) cat("  Consistency:", round(as.numeric(top$Pcons) * 100, 1), "%\n")
      if (!is.null(top$d0)) cat("  Control events:", top$d0, "\n")
      if (!is.null(top$d1)) cat("  Treatment events:", top$d1, "\n")
    }
  } else {
    cat("\nNo subgroup identified.\n")
  }

  # ---------------------------------------------------------------------------
  # Computation Time
  # ---------------------------------------------------------------------------
  if (!is.null(object$minutes_all)) {
    cat("\nComputation time:", round(object$minutes_all, 2), "minutes\n")
  }

  invisible(object)
}


#' Summary Method for aft_dgm_flex Objects
#'
#' Provides detailed summary of AFT data generating mechanism.
#'
#' @param object An aft_dgm_flex object.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
summary.aft_dgm_flex <- function(object, ...) {
  cat("AFT DGM Summary (Flexible)\n")
  cat("==========================\n\n")

  # Basic info
  cat("Configuration:\n")
  cat("  Model type:", object$model_type, "\n")
  cat("  Super population:", object$n_super, "\n")
  cat("  Seed:", object$seed, "\n\n")

  # Model parameters
  if (!is.null(object$model_params)) {
    cat("Model Parameters:\n")
    mp <- object$model_params
    if (!is.null(mp$mu)) cat("  Intercept (mu):", round(mp$mu, 4), "\n")
    if (!is.null(mp$tau)) cat("  Scale (tau):", round(mp$tau, 4), "\n")
    cat("\n")
  }

  # Subgroup information
  if (!is.null(object$subgroup_info)) {
    cat("Subgroup Information:\n")
    si <- object$subgroup_info

    # Variables used
    if (!is.null(si$vars)) {
      cat("  Variables:", paste(si$vars, collapse = ", "), "\n")
    }

    # Definitions
    if (!is.null(si$definitions)) {
      cat("  Definitions:\n")
      for (def in si$definitions) {
        cat("    -", def, "\n")
      }
    }

    cat("  Size:", si$size, sprintf("(%.1f%%)\n", si$proportion * 100))
    cat("\n")
  }

  # Hazard ratios
  if (!is.null(object$hazard_ratios)) {
    cat("True Hazard Ratios:\n")
    hr <- object$hazard_ratios

    .print_hr("Overall", hr$overall)
    .print_hr("Harm subgroup (H)", hr$harm_subgroup)
    .print_hr("No-harm subgroup (Hc)", hr$no_harm_subgroup)

    if (!is.null(hr$AHR)) {
      cat("\n  Average Hazard Ratios:\n")
      .print_hr("AHR (overall)", hr$AHR, indent = 4)
      .print_hr("AHR_H", hr$AHR_harm, indent = 4)
      .print_hr("AHR_Hc", hr$AHR_no_harm, indent = 4)
    }
  }

  # Super population summary
  if (!is.null(object$df_super)) {
    cat("\nSuper Population Summary:\n")
    df <- object$df_super

    cat("  N:", nrow(df), "\n")
    if ("treat" %in% names(df)) {
      cat("  Treated:", sum(df$treat), sprintf("(%.1f%%)\n", 100 * mean(df$treat)))
    }
    if ("flag_harm" %in% names(df)) {
      cat("  In harm subgroup:", sum(df$flag_harm),
          sprintf("(%.1f%%)\n", 100 * mean(df$flag_harm)))
    }
  }

  invisible(object)
}


#' Summary Method for gbsg_dgm Objects
#'
#' Provides detailed summary of GBSG-based data generating mechanism.
#'
#' @param object A gbsg_dgm object.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
summary.gbsg_dgm <- function(object, ...) {
  cat("GBSG DGM Summary\n")
  cat("================\n\n")

  # Basic config
  cat("Configuration:\n")
  cat("  Model type:", object$model_type, "\n")
  cat("  Super population:", object$n_super, "\n")
  cat("  Seed:", object$seed, "\n\n")

  # Effect modifiers
  if (!is.null(object$effect_modifiers)) {
    cat("Effect Modifiers:\n")
    em <- object$effect_modifiers
    cat("  k_treat:", em$k_treat, "(treatment effect multiplier)\n")
    cat("  k_inter:", em$k_inter, "(interaction multiplier)\n")
    cat("  k_z3:", em$k_z3, "(menopausal status effect)\n\n")
  }

  # Subgroup
  if (!is.null(object$subgroup_info)) {
    cat("Harm Subgroup:\n")
    si <- object$subgroup_info
    cat("  Definition:", si$definition, "\n")
    cat("  ER threshold:", si$er_threshold, "\n")
    cat("  Size:", si$size, sprintf("(%.1f%%)\n", si$proportion * 100))
    cat("\n")
  }

  # True hazard ratios
  cat("True Hazard Ratios (Cox-based):\n")
  .print_hr("Overall (ITT)", object$hr_causal)
  .print_hr("Harm subgroup (H)", object$hr_H_true)
  .print_hr("No-harm subgroup (Hc)", object$hr_Hc_true)

  # AHR metrics
  if (!is.null(object$AHR)) {
    cat("\nAverage Hazard Ratios:\n")
    .print_hr("AHR", object$AHR)
    .print_hr("AHR_H", object$AHR_H_true)
    .print_hr("AHR_Hc", object$AHR_Hc_true)
  }

  # Model parameters
  if (!is.null(object$model_params)) {
    cat("\nAFT Model Parameters:\n")
    mp <- object$model_params
    cat("  mu (intercept):", round(mp$mu, 4), "\n")
    cat("  sigma (scale):", round(mp$sigma, 4), "\n")
  }

  invisible(object)
}


#' Summary Method for fs_bootstrap Objects
#'
#' Provides detailed summary of bootstrap bias correction results.
#'
#' @param object An fs_bootstrap object.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
summary.fs_bootstrap <- function(object, ...) {
  cat("ForestSearch Bootstrap Summary\n")
  cat("==============================\n\n")

  # Bootstrap configuration
  cat("Configuration:\n")
  if (!is.null(object$n_boots)) {
    cat("  Bootstrap iterations:", object$n_boots, "\n")
  }
  if (!is.null(object$parallel_args)) {
    cat("  Workers:", object$parallel_args$workers, "\n")
    cat("  Plan:", object$parallel_args$plan, "\n")
  }
  cat("\n")

  # Success statistics
  if (!is.null(object$results)) {
    cat("Bootstrap Results:\n")
    res <- object$results
    n_total <- nrow(res)

    # Count successful iterations
    n_H <- sum(!is.na(res$H_biasadj_2))
    n_Hc <- sum(!is.na(res$Hc_biasadj_2))

    cat("  Total iterations:", n_total, "\n")
    cat("  Successful (H):", n_H, sprintf("(%.1f%%)\n", 100 * n_H / n_total))
    cat("  Successful (Hc):", n_Hc, sprintf("(%.1f%%)\n", 100 * n_Hc / n_total))
    cat("\n")
  }

  # Bias-corrected estimates
  if (!is.null(object$FSsg_tab)) {
    cat("Bias-Corrected Estimates:\n")
    print(object$FSsg_tab)
    cat("\n")
  }

  # Detailed summary table
  if (!is.null(object$summary) && !is.null(object$summary$table)) {
    cat("Summary Statistics:\n")
    print(object$summary$table)
    cat("\n")
  }

  # Timing
  if (!is.null(object$timing)) {
    cat("Timing:\n")
    cat("  Total:", round(object$timing$total_minutes, 2), "minutes\n")
    if (!is.null(object$timing$mean_per_iteration)) {
      cat("  Mean per iteration:", round(object$timing$mean_per_iteration, 2), "seconds\n")
    }
  }

  invisible(object)
}


# =============================================================================
# Internal Helper Functions
# =============================================================================

#' Null Coalescing Operator
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Print Parameter Value
#' @keywords internal
.print_param <- function(params, name, alt_name = NULL, indent = 2L) {
  val <- params[[name]]
  if (is.null(val) && !is.null(alt_name)) {
    val <- params[[alt_name]]
  }

  if (!is.null(val)) {
    spaces <- paste(rep(" ", indent), collapse = "")
    cat(spaces, name, ": ", val, "\n", sep = "")
  }
}


#' Print Hazard Ratio
#' @keywords internal
.print_hr <- function(label, value, indent = 2L) {
  if (!is.null(value) && !is.na(value)) {
    spaces <- paste(rep(" ", indent), collapse = "")
    cat(spaces, label, ": ", sprintf("%.3f", value), "\n", sep = "")
  }
}


#' Print Results by sg_focus
#' @keywords internal
.print_results_by_focus <- function(object) {
  gc <- object$grp.consistency
  if (is.null(gc)) return()

  focus_types <- c("hr", "maxSG", "minSG")
  has_results <- FALSE

  for (focus in focus_types) {
    result_name <- paste0("out_", focus)
    out <- gc[[result_name]]

    if (!is.null(out) && !is.null(out$sg.harm)) {
      if (!has_results) {
        cat("Results by sg_focus:\n")
        has_results <- TRUE
      }

      sg <- if (!is.null(out$result)) out$result[1L, ] else NULL

      if (!is.null(sg)) {
        cat(sprintf("  %-8s: N=%d, HR=%.2f, Pcons=%.1f%%\n",
                    focus,
                    as.integer(sg$N),
                    as.numeric(sg$hr),
                    as.numeric(sg$Pcons) * 100))
      }
    }
  }
}
