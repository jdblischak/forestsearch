# ============================================================================
# unexport_internals.R
#
# Automates Tier 3 export audit: removes `@export` from functions that
# should become fully internal, replacing with `@keywords internal` or
# `@noRd` as specified.
#
# Usage:
#   source("unexport_internals.R")
#
#   # Step 1: Safety scan — check for parallel usage (ALWAYS run first):
#   scan_parallel_usage(pkg_dir = ".")
#
#   # Step 2: Dry run — shows what would change, writes nothing:
#   unexport_internals(pkg_dir = ".")
#
#   # Step 3: Apply changes:
#   unexport_internals(pkg_dir = ".", dry_run = FALSE)
#
#   # Step 4: Rebuild
#   devtools::document()
#   devtools::check()
#
# What it does:
#   For each Tier 3 function, finds its roxygen block and either:
#     - Replaces `#' @export` with `#' @keywords internal` (documented but
#       hidden from pkgdown index, removed from NAMESPACE)
#     - Replaces `#' @export` with `#' @noRd` (no documentation page at all,
#       removed from NAMESPACE)
#
# IMPORTANT:
#   Unlike Tier 2, this REMOVES functions from the NAMESPACE. Any function
#   called by parallel workers (foreach/doFuture/callr) will break if
#   un-exported. Always run scan_parallel_usage() first.
# ============================================================================


# ============================================================================
# Step 1: Safety scanner — detect parallel usage before un-exporting
# ============================================================================

#' Scan for Parallel Usage of Tier 3 Functions
#'
#' Searches all R source files for references to Tier 3 candidate functions
#' inside foreach/doFuture blocks, .export arguments, and other parallel
#' contexts. Run this BEFORE un-exporting anything.
#'
#' @param pkg_dir Character. Path to package root. Default: "."
#' @return Data frame of findings (invisibly).
scan_parallel_usage <- function(pkg_dir = ".") {

  tier3_names <- names(tier3_functions())

  r_dir <- file.path(pkg_dir, "R")
  r_files <- list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE)

  # Patterns that indicate parallel execution contexts
  parallel_patterns <- c(
    "foreach",
    "%dofuture%",
    "%dopar%",
    "future_lapply",
    "future_sapply",
    "future_apply",
    "\\.export\\s*=",
    "clusterExport",
    "clusterCall"
  )

  cat("\n")
  cat("==========================================\n")
  cat("  Parallel Usage Scan for Tier 3 Functions\n")
  cat("==========================================\n\n")

  findings <- data.frame(
    func    = character(),
    file    = character(),
    line    = integer(),
    context = character(),
    stringsAsFactors = FALSE
  )

  for (rfile in r_files) {
    lines <- readLines(rfile, warn = FALSE)
    fname <- basename(rfile)

    # Identify lines that are in a parallel context (rough heuristic:
    # within 50 lines after a foreach/future_lapply call)
    parallel_zones <- integer()
    for (pp in parallel_patterns) {
      hits <- grep(pp, lines)
      for (h in hits) {
        parallel_zones <- c(parallel_zones, seq(h, min(h + 50, length(lines))))
      }
    }
    parallel_zones <- unique(parallel_zones)

    # Check each Tier 3 function name
    for (fn in tier3_names) {
      # Word-boundary match to avoid partial hits
      fn_pattern <- paste0("\\b", gsub("\\.", "\\\\.", fn), "\\b")
      ref_lines <- grep(fn_pattern, lines)

      # Exclude the function's own definition or alias assignment line
      def_pattern <- paste0("^\\s*", gsub("\\.", "\\\\.", fn),
                            "\\s*(<-|=)\\s*")
      def_lines <- grep(def_pattern, lines)
      ref_lines <- setdiff(ref_lines, def_lines)

      # Also exclude roxygen comment lines
      ref_lines <- ref_lines[!grepl("^\\s*#", lines[ref_lines])]

      for (rl in ref_lines) {
        in_parallel <- rl %in% parallel_zones
        context <- if (in_parallel) {
          "** IN PARALLEL ZONE **"
        } else {
          "regular call"
        }

        findings <- rbind(findings, data.frame(
          func = fn, file = fname, line = rl,
          context = context,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Also scan for explicit .export references
  cat("Checking .export arguments for Tier 3 names...\n\n")
  for (rfile in r_files) {
    lines <- readLines(rfile, warn = FALSE)
    export_lines <- grep("\\.export\\s*=", lines)
    for (el in export_lines) {
      for (fn in tier3_names) {
        if (grepl(fn, lines[el], fixed = TRUE)) {
          findings <- rbind(findings, data.frame(
            func = fn, file = basename(rfile), line = el,
            context = "** EXPLICIT .export REFERENCE **",
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  # Report
  parallel_hits <- findings[grepl("PARALLEL|EXPLICIT", findings$context), ]

  if (nrow(parallel_hits) > 0) {
    cat("!! PARALLEL USAGE DETECTED — review before un-exporting:\n\n")
    cat(sprintf("  %-35s %-35s Line   Context\n", "Function", "File"))
    cat(sprintf("  %-35s %-35s %-6s %s\n",
                strrep("-", 35), strrep("-", 35), strrep("-", 6), strrep("-", 30)))
    for (i in seq_len(nrow(parallel_hits))) {
      cat(sprintf("  %-35s %-35s %-6d %s\n",
                  parallel_hits$func[i],
                  parallel_hits$file[i],
                  parallel_hits$line[i],
                  parallel_hits$context[i]))
    }
    cat("\n  These functions should NOT be un-exported.\n")
    cat("  Move them from Tier 3 to Tier 2 (keep @export, add @keywords internal).\n")
  } else {
    cat("  No parallel usage detected for any Tier 3 function.\n")
    cat("  Safe to proceed with un-exporting.\n")
  }

  # Show all references for context
  if (nrow(findings) > 0) {
    cat("\nAll references found (", nrow(findings), " total):\n\n")
    cat(sprintf("  %-35s %-35s Line   Context\n", "Function", "File"))
    cat(sprintf("  %-35s %-35s %-6s %s\n",
                strrep("-", 35), strrep("-", 35), strrep("-", 6), strrep("-", 30)))
    for (i in seq_len(nrow(findings))) {
      cat(sprintf("  %-35s %-35s %-6d %s\n",
                  findings$func[i],
                  findings$file[i],
                  findings$line[i],
                  findings$context[i]))
    }
  } else {
    cat("\n  No references found to any Tier 3 function outside their definitions.\n")
  }

  cat("\n")
  invisible(findings)
}


# ============================================================================
# Tier 3 function registry
# ============================================================================

#' Tier 3 Function List with Target Tag
#'
#' Returns a named list: name = function name, value = target tag.
#' "noRd"     -> replace @export with @noRd (no help page)
#' "internal" -> replace @export with @keywords internal (hidden help page)
#'
#' @return Named character vector.
#' @noRd
tier3_functions <- function() {
  c(
    # --- Bitwise / encoding helpers -> @noRd ---
    "one.zero"                          = "noRd",
    "ztrail"                            = "noRd",
    "acm.disjctif"                      = "noRd",
    "dummy"                             = "noRd",
    "dummy2"                            = "noRd",

    # --- Data preparation internals -> @keywords internal ---
    "get_cut_name"                      = "internal",
    "is_flag_continuous"                 = "internal",
    "is_flag_drop"                      = "internal",
    "is.continuous"                      = "internal",
    "cut_var"                           = "internal",
    "detect_variable_types"             = "internal",
    "add_id_column"                     = "internal",

    # --- Tree cut extraction -> @keywords internal ---
    "extract_all_tree_cuts"             = "internal",
    "extract_selected_tree_cuts"        = "internal",
    "extract_tree_cuts"                 = "internal",
    "extract_idx_flagredundancy"        = "internal",

    # --- Subgroup sorting / filtering -> @keywords internal ---
    "sort_subgroups_preview"            = "internal",
    "remove_near_duplicate_subgroups"   = "internal",
    "remove_redundant_subgroups"        = "internal",

    # --- Summary table variants -> @keywords internal ---
    "create_summary_table_compact"      = "internal",
    "create_summary_table_minimal"      = "internal",
    "create_summary_table_presentation" = "internal",
    "create_summary_table_publication"  = "internal",
    "create_sample_size_table"          = "internal",
    "create_subgroup_summary_df"        = "internal",

    # --- Bootstrap / formatting internals -> @keywords internal ---
    "summarize_bootstrap_subgroups"     = "internal",
    "summarize_bootstrap_events"        = "internal",
    "summarize_factor_presence_robust"  = "internal",
    "format_bootstrap_diagnostics_table" = "internal",
    "format_bootstrap_timing_table"     = "internal",
    "format_subgroup_summary_tables"    = "internal",
    "create_factor_summary_tables"      = "internal",
    "format_results"                    = "internal",
    "format_oc_results"                 = "internal",

    # --- Simulation internals -> @keywords internal ---
    "get_dgm_with_output"              = "internal",
    "validate_k_inter_effect"          = "internal",
    "compute_detection_probability"    = "internal",
    "create_null_result"               = "internal",
    "create_success_result"            = "internal",
    "default_fs_params"                = "internal",
    "default_grf_params"               = "internal",

    # --- Miscellaneous helpers -> @keywords internal ---
    "ci_est"                           = "internal",
    "calc_cov"                         = "internal",
    "get_targetEst"                    = "internal",
    "calculate_counts"                 = "internal",
    "calculate_potential_hr"           = "internal",
    "density_threshold_both"           = "internal",
    "find_quantile_for_proportion"     = "internal",
    "qlow"                             = "internal",
    "qhigh"                            = "internal",
    "get_best_survreg"                 = "internal",
    "cox_summary_batch"                = "internal",
    "cox_summary_legacy"               = "internal",
    "cox_summary_vectorized"           = "internal",
    "sens_text"                        = "internal",
    "figure_note"                      = "internal",
    "km_summary"                       = "internal",
    "plot_subgroup"                    = "internal"
  )
}


# ============================================================================
# Step 2 & 3: Apply the un-export changes
# ============================================================================

#' Un-export Tier 3 Internal Functions
#'
#' For each Tier 3 function, removes `@export` from the roxygen block and
#' replaces it with `@keywords internal` or `@noRd` as specified.
#'
#' @param pkg_dir Character. Path to package root. Default: "."
#' @param dry_run Logical. If TRUE (default), report but don't write files.
#' @return Data frame of results (invisibly).
unexport_internals <- function(pkg_dir = ".", dry_run = TRUE) {

  tier3 <- tier3_functions()
  fn_names <- names(tier3)

  r_dir <- file.path(pkg_dir, "R")
  if (!dir.exists(r_dir)) {
    stop("Directory not found: ", r_dir)
  }

  r_files <- list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE)

  results <- data.frame(
    func        = character(),
    file        = character(),
    action      = character(),
    tag         = character(),
    line_number = integer(),
    stringsAsFactors = FALSE
  )

  # Cache file contents
  file_contents <- list()

  for (fn in fn_names) {

    target_tag <- tier3[[fn]]

    # Build regex for function definition or alias assignment
    # Matches: fn <- function(  OR  fn <- other_name  OR  fn = function(
    fn_def_pattern <- paste0("^\\s*", gsub("\\.", "\\\\.", fn),
                             "\\s*(<-|=)\\s*")

    found <- FALSE

    for (rfile in r_files) {

      if (is.null(file_contents[[rfile]])) {
        file_contents[[rfile]] <- list(
          lines    = readLines(rfile, warn = FALSE),
          modified = FALSE
        )
      }

      lines <- file_contents[[rfile]]$lines

      def_lines <- grep(fn_def_pattern, lines)
      if (length(def_lines) == 0) next

      # Try each match — the first may be a local variable, not the
      # exported function definition. Pick the one with @export above it.
      for (def_line in def_lines) {

        # Walk backward to find roxygen block
        export_line <- NULL
        internal_line <- NULL
        nord_line <- NULL

        for (i in seq(def_line - 1, max(1, def_line - 80), by = -1)) {
          ln <- trimws(lines[i])
          if (!grepl("^#'", ln) && ln != "") break

          if (grepl("^#'\\s*@export\\s*$", ln)) export_line <- i
          if (grepl("^#'\\s*@keywords\\s+internal", ln)) internal_line <- i
          if (grepl("^#'\\s*@noRd\\s*$", ln)) nord_line <- i
        }

        if (!is.null(export_line)) break
      }

      if (is.null(export_line)) next

      found <- TRUE

      # Collect line indices to remove (applied in one pass at end)
      lines_to_remove <- integer()

      # Determine what to do
      if (target_tag == "noRd") {
        if (!is.null(nord_line)) {
          # Already has @noRd — just remove @export
          lines_to_remove <- c(lines_to_remove, export_line)
          action <- "removed @export (already has @noRd)"
        } else {
          # Replace @export with @noRd
          lines[export_line] <- "#' @noRd"
          action <- "replaced @export with @noRd"
        }
        # Also remove @keywords internal if present (redundant with @noRd)
        if (!is.null(internal_line)) {
          lines_to_remove <- c(lines_to_remove, internal_line)
          action <- paste(action, "+ removed redundant @keywords internal")
        }

      } else {
        # target_tag == "internal"
        if (!is.null(internal_line)) {
          # Already has @keywords internal — just remove @export
          lines_to_remove <- c(lines_to_remove, export_line)
          action <- "removed @export (already has @keywords internal)"
        } else {
          # Replace @export with @keywords internal
          lines[export_line] <- "#' @keywords internal"
          action <- "replaced @export with @keywords internal"
        }
      }

      # Remove collected lines in one pass (negative indexing)
      if (length(lines_to_remove) > 0) {
        lines <- lines[-lines_to_remove]
      }

      file_contents[[rfile]]$lines <- lines
      file_contents[[rfile]]$modified <- TRUE

      results <- rbind(results, data.frame(
        func = fn, file = basename(rfile),
        action = if (dry_run) paste("WOULD:", action) else action,
        tag = target_tag,
        line_number = export_line,
        stringsAsFactors = FALSE
      ))

      break
    }

    if (!found) {
      results <- rbind(results, data.frame(
        func = fn, file = "(not found)",
        action = "WARNING: function not located",
        tag = target_tag,
        line_number = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }
  }

  # --------------------------------------------------------------------------
  # Write modified files
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
  cat("==========================================\n")
  cat("  Tier 3 Un-Export Tool\n")
  cat("==========================================\n")
  cat("  Mode:", if (dry_run) "DRY RUN" else "APPLIED", "\n")
  cat("  Package:", normalizePath(pkg_dir), "\n")
  cat("  Functions checked:", length(fn_names), "\n\n")

  n_nord     <- sum(results$tag == "noRd" & !grepl("WARNING", results$action))
  n_internal <- sum(results$tag == "internal" & !grepl("WARNING", results$action))
  n_warn     <- sum(grepl("WARNING", results$action))

  cat("  -> @noRd (fully hidden):       ", n_nord, "\n")
  cat("  -> @keywords internal (hidden): ", n_internal, "\n")
  cat("  Not found (warnings):           ", n_warn, "\n")

  if (!dry_run) {
    cat("  Files written:                  ", files_modified, "\n")
  }

  cat("\nDetails:\n")
  cat(sprintf("  %-40s %-12s %-35s %s\n",
              "Function", "Target", "File", "Action"))
  cat(sprintf("  %-40s %-12s %-35s %s\n",
              strrep("-", 40), strrep("-", 12), strrep("-", 35), strrep("-", 40)))
  for (i in seq_len(nrow(results))) {
    cat(sprintf("  %-40s %-12s %-35s %s\n",
                results$func[i],
                results$tag[i],
                results$file[i],
                results$action[i]))
  }

  cat("\n")
  if (dry_run) {
    cat("This was a DRY RUN. To apply:\n")
    cat("  1. Run scan_parallel_usage() first to verify safety\n")
    cat("  2. Run unexport_internals(pkg_dir = \".\", dry_run = FALSE)\n")
    cat("  3. Run devtools::document()\n")
    cat("  4. Run devtools::check()\n")
  } else {
    cat("Done. Now run:\n")
    cat("  devtools::document()   # regenerate NAMESPACE\n")
    cat("  devtools::check()      # verify no breakage\n")
  }
  cat("\n")

  invisible(results)
}
