# =============================================================================
# ForestSearch S3 Print Methods
# =============================================================================
#
# CRAN COMPLIANCE NOTE:
# All print methods for ForestSearch classes are defined in this single file.
# This avoids duplicate method registration errors.
#
# =============================================================================

#' Print Method for forestsearch Objects
#'
#' Displays a concise summary of ForestSearch analysis results including
#' the identified subgroup, key parameters, and computation time.
#'
#' @param x A forestsearch object returned by \code{\link{forestsearch}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' \dontrun
#' result <- forestsearch(df.analysis = trial_data)
#' print(result)
#' }
#'
#' @export
print.forestsearch <- function(x, ...) {
  cat("ForestSearch Results\n")
  cat("====================\n\n")

  # Handle case where no subgroup was identified
 if (is.null(x$sg.harm)) {
    cat("No subgroup identified.\n")
    .print_computation_time(x)
    return(invisible(x))
  }

  # Print selected subgroup information
  cat("Selected Subgroup:\n")
  cat("  Definition:", paste(x$sg.harm, collapse = " & "), "\n")

  # Print consistency results if available
  if (!is.null(x$grp.consistency)) {
    gc <- x$grp.consistency

    # Print sg_focus
    if (!is.null(gc$sg_focus)) {
      cat("  sg_focus:", gc$sg_focus, "\n")
    }

    # Print detailed results from out_sg
    if (!is.null(gc$out_sg) && !is.null(gc$out_sg$result)) {
      top <- gc$out_sg$result[1L, ]
      cat("  N:", top$N, "\n")
      cat("  HR:", round(as.numeric(top$hr), 3), "\n")
      cat("  Pcons:", round(as.numeric(top$Pcons), 3), "\n")
    }

    # Print algorithm used
    if (!is.null(gc$algorithm)) {
      cat("  Algorithm:", gc$algorithm, "\n")
    }

    # Print evaluation statistics
    if (!is.null(gc$n_candidates_evaluated)) {
      cat("  Candidates evaluated:", gc$n_candidates_evaluated, "\n")
      cat("  Candidates passed:", gc$n_passed, "\n")
    }
  }

  .print_computation_time(x)
  invisible(x)
}


#' Print Method for aft_dgm_flex Objects
#'
#' Displays summary information about an AFT data generating mechanism.
#'
#' @param x An aft_dgm_flex object from \code{\link{generate_aft_dgm_flex}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.aft_dgm_flex <- function(x, ...) {
  cat("AFT Data Generating Mechanism (Flexible)\n")
  cat("=========================================\n\n")

  # Model type
  cat("Model type:", x$model_type, "\n")
  cat("Super population size:", x$n_super, "\n\n")

  # Subgroup information
  if (!is.null(x$subgroup_info)) {
    cat("Subgroup Definition:\n")

    if (!is.null(x$subgroup_info$definitions)) {
      for (def in x$subgroup_info$definitions) {
        cat("  -", def, "\n")
      }
    }

    cat("  Size:", x$subgroup_info$size,
        sprintf("(%.1f%%)", x$subgroup_info$proportion * 100), "\n\n")
  }

  # Hazard ratios
  if (!is.null(x$hazard_ratios)) {
    cat("Hazard Ratios:\n")
    hr <- x$hazard_ratios

    cat("  Overall:", sprintf("%.3f", hr$overall), "\n")

    if (!is.null(hr$harm_subgroup)) {
      cat("  Harm subgroup (H):", sprintf("%.3f", hr$harm_subgroup), "\n")
    }
    if (!is.null(hr$no_harm_subgroup)) {
      cat("  No-harm subgroup (Hc):", sprintf("%.3f", hr$no_harm_subgroup), "\n")
    }
    if (!is.null(hr$AHR)) {
      cat("  AHR (average):", sprintf("%.3f", hr$AHR), "\n")
    }
  }

  cat("\nSeed:", x$seed, "\n")

  invisible(x)
}


#' Print Method for gbsg_dgm Objects
#'
#' Displays summary information about a GBSG-based data generating mechanism.
#'
#' @param x A gbsg_dgm object from \code{\link{create_gbsg_dgm}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.gbsg_dgm <- function(x, ...) {
  cat("GBSG-Based Data Generating Mechanism\n")
  cat("=====================================\n\n")

  # Model type
  cat("Model type:", x$model_type, "\n")
  cat("Super population size:", x$n_super, "\n\n")

  # Effect modifiers
  if (!is.null(x$effect_modifiers)) {
    cat("Effect Modifiers:\n")
    em <- x$effect_modifiers
    cat("  k_treat:", em$k_treat, "\n")
    cat("  k_inter:", em$k_inter, "\n")
    cat("  k_z3:", em$k_z3, "\n\n")
  }

  # Subgroup information
  if (!is.null(x$subgroup_info)) {
    cat("Subgroup:\n")
    si <- x$subgroup_info
    cat("  Definition:", si$definition, "\n")
    cat("  Size:", si$size,
        sprintf("(%.1f%%)", si$proportion * 100), "\n\n")
  }

  # Hazard ratios
  cat("Hazard Ratios:\n")
  cat("  Overall:", sprintf("%.3f", x$hr_causal), "\n")
  cat("  Harm subgroup (H):", sprintf("%.3f", x$hr_H_true), "\n")
  cat("  No-harm subgroup (Hc):", sprintf("%.3f", x$hr_Hc_true), "\n")

  # AHR if available
  if (!is.null(x$AHR)) {
    cat("\nAverage Hazard Ratios:\n")
    cat("  AHR:", sprintf("%.3f", x$AHR), "\n")
    cat("  AHR_H:", sprintf("%.3f", x$AHR_H_true), "\n")
    cat("  AHR_Hc:", sprintf("%.3f", x$AHR_Hc_true), "\n")
  }

  cat("\nSeed:", x$seed, "\n")

  invisible(x)
}


#' Print Method for fs_bootstrap Objects
#'
#' Displays summary of bootstrap bias correction results.
#'
#' @param x An fs_bootstrap object from
#'   \code{\link{forestsearch_bootstrap_dofuture}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.fs_bootstrap <- function(x, ...) {
  cat("ForestSearch Bootstrap Results\n")
  cat("==============================\n\n")

  # Number of bootstraps
  if (!is.null(x$n_boots)) {
    cat("Bootstrap iterations:", x$n_boots, "\n")
  }

  # Success rate
  if (!is.null(x$results)) {
    n_total <- nrow(x$results)
    n_success <- sum(!is.na(x$results$H_biasadj_2))
    cat("Successful iterations:", n_success, "/", n_total,
        sprintf("(%.1f%%)", 100 * n_success / n_total), "\n\n")
  }

  # Bias-corrected estimates
  if (!is.null(x$FSsg_tab)) {
    cat("Bias-Corrected Estimates:\n")
    print(x$FSsg_tab)
  }

  # Timing
  if (!is.null(x$timing)) {
    cat("\nComputation time:", round(x$timing$total_minutes, 2), "minutes\n")
  }

  invisible(x)
}


# =============================================================================
# Internal Helper Functions
# =============================================================================

#' Print Computation Time
#' @keywords internal
.print_computation_time <- function(x) {
  if (!is.null(x$minutes_all)) {
    cat("\nComputation time:", round(x$minutes_all, 2), "minutes\n")
  } else if (!is.null(x$computation_time)) {
    cat("\nComputation time:", sprintf("%.1f", x$computation_time), "seconds\n")
  }
}
