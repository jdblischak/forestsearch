# ============================================================================
# tag_internal_exports.R
#
# Automates Tier 2 export audit: adds `@keywords internal` to exported
# functions that should remain in the NAMESPACE but be hidden from the
# user-facing pkgdown reference index.
#
# Usage:
#   source("tag_internal_exports.R")
#   # Dry run (default) — shows what would change, writes nothing:
#   tag_internal_exports(pkg_dir = ".")
#   # Apply changes:
#   tag_internal_exports(pkg_dir = ".", dry_run = FALSE)
#
# What it does:
#   For each function in the Tier 2 list, finds its roxygen block in R/*.R
#   and inserts `#' @keywords internal` immediately before the `#' @export`
#   line, if not already present.
#
# Safety:
#   - Never modifies files in dry_run = TRUE mode (the default)
#   - Skips functions that already have @keywords internal
#   - Skips functions it cannot locate (with a warning)
#   - Reports every change made
#   - Creates no backup files (use version control instead)
# ============================================================================

tag_internal_exports <- function(pkg_dir = ".", dry_run = TRUE) {


  # --------------------------------------------------------------------------
  # Tier 2 functions: keep @export, add @keywords internal
  # --------------------------------------------------------------------------
  tier2_functions <- c(
    # Parallel worker requirements
    "bootstrap_results",
    "bootstrap_ystar",
    "count_boot_id",
    "get_dfRes",
    "build_cox_formula",
    "fit_cox_models",
    "get_Cox_sg",
    "get_split_hr_fast",
    "run_single_consistency_split",
    "setup_parallel_SGcons",
    "sg_consistency_out",
    "evaluate_consistency_twostage",
    "wilson_ci",
    "early_stop_decision",
    "evaluate_comparison",
    "evaluate_cuts_once",
    "sort_subgroups",
    "select_best_subgroup",
    "analyze_subgroup",
    "assign_subgroup_membership",
    "extract_subgroup",
    "get_subgroup_membership",
    "prepare_subgroup_data",
    "evaluate_subgroup_consistency",
    "generate_bootstrap_synthetic",
    "generate_bootstrap_with_noise",
    "generate_gbsg_bootstrap_general",

    # Algorithm internals (not parallel, but used by exported functions)
    "filter_call_args",
    "filter_by_lassokeep",
    "get_combinations_info",
    "get_conf_force",
    "get_covs_in",
    "process_conf_force_expr",
    "create_grf_config",
    "validate_grf_data",
    "fit_causal_forest",
    "fit_policy_trees",
    "compute_node_metrics",
    "find_leaf_split",
    "print_grf_details"
  )

  r_dir <- file.path(pkg_dir, "R")
  if (!dir.exists(r_dir)) {
    stop("Directory not found: ", r_dir, "\n",
         "Set pkg_dir to your package root.")
  }

  r_files <- list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE)

  if (length(r_files) == 0) {
    stop("No .R files found in: ", r_dir)
  }

  # --------------------------------------------------------------------------
  # For each function, find its @export line and insert @keywords internal
  # --------------------------------------------------------------------------
  results <- data.frame(
    func        = character(),
    file        = character(),
    status      = character(),
    line_number = integer(),
    stringsAsFactors = FALSE
  )

  # Track which files need writing (file -> modified lines)
  file_contents <- list()

  for (fn in tier2_functions) {

    # Pattern to match the function definition (handles assignment styles)
    # Matches: fn <- function(  or  fn = function(
    fn_def_pattern <- paste0("^\\s*", gsub("\\.", "\\\\.", fn),
                             "\\s*(<-|=)\\s*function\\s*\\(")

    found <- FALSE

    for (rfile in r_files) {

      # Read file (cache to avoid re-reading)
      if (is.null(file_contents[[rfile]])) {
        file_contents[[rfile]] <- list(
          lines    = readLines(rfile, warn = FALSE),
          modified = FALSE
        )
      }

      lines <- file_contents[[rfile]]$lines

      # Find the function definition line
      def_lines <- grep(fn_def_pattern, lines)

      if (length(def_lines) == 0) next

      # Use the first match
      def_line <- def_lines[1]

      # Walk backward from definition to find the roxygen block
      export_line <- NULL
      already_internal <- FALSE
      roxygen_start <- NULL

      for (i in seq(def_line - 1, max(1, def_line - 80), by = -1)) {
        ln <- trimws(lines[i])

        # Stop if we hit a non-roxygen, non-blank line
        if (!grepl("^#'", ln) && ln != "") break

        # Track roxygen boundaries
        if (grepl("^#'", ln)) {
          roxygen_start <- i

          if (grepl("^#'\\s*@export\\s*$", ln)) {
            export_line <- i
          }
          if (grepl("^#'\\s*@keywords\\s+internal", ln)) {
            already_internal <- TRUE
          }
        }
      }

      if (is.null(export_line)) {
        # No @export found in roxygen block — skip
        next
      }

      found <- TRUE

      if (already_internal) {
        results <- rbind(results, data.frame(
          func = fn, file = basename(rfile),
          status = "SKIP (already has @keywords internal)",
          line_number = export_line,
          stringsAsFactors = FALSE
        ))
        break
      }

      # Insert @keywords internal on the line before @export
      new_line <- "#' @keywords internal"
      lines <- append(lines, new_line, after = export_line - 1)

      file_contents[[rfile]]$lines <- lines
      file_contents[[rfile]]$modified <- TRUE

      results <- rbind(results, data.frame(
        func = fn, file = basename(rfile),
        status = if (dry_run) "WOULD ADD" else "ADDED",
        line_number = export_line,
        stringsAsFactors = FALSE
      ))

      break
    }

    if (!found) {
      results <- rbind(results, data.frame(
        func = fn, file = "(not found)",
        status = "WARNING: function not located",
        line_number = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }
  }

  # --------------------------------------------------------------------------
  # Write modified files (unless dry run)
  # --------------------------------------------------------------------------
  files_modified <- 0
  if (!dry_run) {
    for (rfile in names(file_contents)) {
      if (file_contents[[rfile]]$modified) {
        writeLines(file_contents[[rfile]]$lines, rfile)
        files_modified <- files_modified + 1
      }
    }
  }

  # --------------------------------------------------------------------------
  # Report
  # --------------------------------------------------------------------------
  cat("\n")
  cat("========================================\n")
  cat("  Tier 2 @keywords internal Tagger\n")
  cat("========================================\n")
  cat("  Mode:", if (dry_run) "DRY RUN" else "APPLIED", "\n")
  cat("  Package:", normalizePath(pkg_dir), "\n")
  cat("  Functions checked:", length(tier2_functions), "\n")
  cat("\n")

  n_add  <- sum(grepl("ADD", results$status))
  n_skip <- sum(grepl("SKIP", results$status))
  n_warn <- sum(grepl("WARNING", results$status))

  cat("  Will add @keywords internal:", n_add, "\n")
  cat("  Already tagged (skipped):   ", n_skip, "\n")
  cat("  Not found (warnings):       ", n_warn, "\n")

  if (!dry_run) {
    cat("  Files written:              ", files_modified, "\n")
  }

  cat("\n")

  # Print detailed results
  if (nrow(results) > 0) {
    cat("Details:\n")
    cat(sprintf("  %-42s %-38s %s\n", "Function", "File", "Status"))
    cat(sprintf("  %-42s %-38s %s\n",
                strrep("-", 42), strrep("-", 38), strrep("-", 30)))
    for (i in seq_len(nrow(results))) {
      cat(sprintf("  %-42s %-38s %s\n",
                  results$func[i], results$file[i], results$status[i]))
    }
  }

  cat("\n")
  if (dry_run) {
    cat("This was a DRY RUN. To apply changes, run:\n")
    cat('  tag_internal_exports(pkg_dir = ".", dry_run = FALSE)\n')
  } else {
    cat("Done. Run devtools::document() to regenerate NAMESPACE.\n")
    cat("Then run devtools::check() to verify.\n")
  }
  cat("\n")

  invisible(results)
}
