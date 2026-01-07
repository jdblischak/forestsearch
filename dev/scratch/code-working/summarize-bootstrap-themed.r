#' Enhanced Bootstrap Results Summary with Themed Tables
#'
#' Creates comprehensive output including themed formatted tables, diagnostic plots,
#' bootstrap quality metrics, and detailed timing analysis. This function follows
#' the same logic as summarize_bootstrap_results() but applies consistent GT themes
#' to all table outputs.
#'
#' @param sgharm Subgroup definition from original analysis (e.g., fs.est$sgharm)
#' @param boot_results List. Output from forestsearch_bootstrap_dofuture() containing:
#'   \itemize{
#'     \item FSsg_tab: Results table with HR estimates
#'     \item results: Data.table of bootstrap iterations
#'     \item H_estimates: List with H subgroup estimates
#'     \item Hc_estimates: List with Hc subgroup estimates
#'   }
#' @param create_plots Logical. Generate diagnostic plots (default: FALSE)
#' @param est.scale Character. "hr" or "1/hr" for effect scale (default: "hr")
#'
#' @return List with the following components:
#'   \itemize{
#'     \item table: Themed gt table with treatment effects
#'     \item diagnostics: List of diagnostic metrics
#'     \item diagnostics_table: Themed gt table with bootstrap quality metrics
#'     \item timing: List with timing information (if available)
#'     \item timing_table: Themed gt table with timing analysis (if available)
#'     \item plots: List of ggplot2 diagnostic plots (if create_plots=TRUE)
#'     \item subgroup_summary: Subgroup analysis summary (if available)
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run bootstrap analysis
#' boot_results <- forestsearch_bootstrap_dofuture(
#'   data = df,
#'   nb_boots = 1000,
#'   n_cores = 8
#' )
#'
#' # Create themed summary
#' summary <- summarize_bootstrap_results_themed(
#'   sgharm = fs.est$sgharm,
#'   boot_results = boot_results,
#'   create_plots = TRUE
#' )
#'
#' # View themed tables
#' summary$table              # Results with consistent styling
#' summary$diagnostics_table  # Diagnostics with theme
#' summary$timing_table       # Timing with theme
#'
#' # Save themed tables
#' gt::gtsave(summary$table, "bootstrap_results.html")
#' }
summarize_bootstrap_results_themed <- function(sgharm, 
                                               boot_results, 
                                               create_plots = FALSE,
                                               est.scale = "hr") {

  # ===========================================================================
  # SECTION 1: EXTRACT AND VALIDATE COMPONENTS
  # ===========================================================================
  
  # Extract components from boot_results
  FSsg_tab <- boot_results$FSsg_tab
  results <- boot_results$results
  H_estimates <- boot_results$H_estimates
  Hc_estimates <- boot_results$Hc_estimates
  
  # Validate required components
  if (is.null(FSsg_tab)) {
    stop("boot_results$FSsg_tab is missing", call. = FALSE)
  }
  
  if (is.null(results)) {
    stop("boot_results$results is missing", call. = FALSE)
  }

  # ===========================================================================
  # SECTION 2: CALCULATE BOOTSTRAP SUCCESS METRICS
  # ===========================================================================
  
  # Calculate bootstrap success rate
  boot_success_rate <- mean(!is.na(results$H_biasadj_2))
  nb_boots <- nrow(results)
  
  # Create diagnostics summary
  diagnostics <- list(
    n_boots = nb_boots,
    success_rate = boot_success_rate,
    n_successful = sum(!is.na(results$H_biasadj_2)),
    n_failed = sum(is.na(results$H_biasadj_2))
  )

  # ===========================================================================
  # SECTION 3: CREATE THEMED RESULTS TABLE
  # ===========================================================================
  
  formatted_table <- format_bootstrap_table_enhanced(
    FSsg_tab = FSsg_tab,
    nb_boots = nb_boots,
    est.scale = est.scale,
    boot_success_rate = boot_success_rate,
    add_color_coding = TRUE,
    add_interpretation = TRUE
  )

  # ===========================================================================
  # SECTION 4: CREATE THEMED DIAGNOSTICS TABLE
  # ===========================================================================
  
  diagnostics_table_gt <- format_diagnostics_enhanced(
    diagnostics = diagnostics,
    nb_boots = nb_boots,
    results = results,
    add_visual_indicators = TRUE
  )

  # ===========================================================================
  # SECTION 5: EXTRACT AND ANALYZE TIMING INFORMATION
  # ===========================================================================
  
  has_timing <- !is.null(attr(results, "timing"))
  timing_table_gt <- NULL
  timing_list <- NULL
  
  if (has_timing) {
    # Overall timing from attributes
    overall_timing <- attr(results, "timing")
    
    # Per-iteration statistics (if columns exist)
    has_iteration_timing <- "tmins_iteration" %in% names(results)
    has_search_timing <- "tmins_search" %in% names(results)
    
    # Calculate iteration statistics
    if (has_iteration_timing) {
      iteration_stats <- list(
        mean = mean(results$tmins_iteration, na.rm = TRUE),
        median = median(results$tmins_iteration, na.rm = TRUE),
        sd = sd(results$tmins_iteration, na.rm = TRUE),
        min = min(results$tmins_iteration, na.rm = TRUE),
        max = max(results$tmins_iteration, na.rm = TRUE),
        q25 = quantile(results$tmins_iteration, 0.25, na.rm = TRUE),
        q75 = quantile(results$tmins_iteration, 0.75, na.rm = TRUE)
      )
    } else {
      iteration_stats <- NULL
    }
    
    # ForestSearch timing (subset where forestsearch ran)
    if (has_search_timing) {
      fs_times <- results$tmins_search[!is.na(results$tmins_search)]
      if (length(fs_times) > 0) {
        fs_stats <- list(
          n_runs = length(fs_times),
          pct_runs = length(fs_times) / nb_boots * 100,
          mean = mean(fs_times),
          median = median(fs_times),
          sd = sd(fs_times),
          min = min(fs_times),
          max = max(fs_times),
          total = sum(fs_times)
        )
      } else {
        fs_stats <- NULL
      }
    } else {
      fs_stats <- NULL
    }
    
    # Overhead timing (time not in forestsearch)
    if (has_iteration_timing && has_search_timing) {
      overhead_times <- results$tmins_iteration - results$tmins_search
      overhead_times <- overhead_times[!is.na(overhead_times)]
      
      if (length(overhead_times) > 0) {
        overhead_stats <- list(
          mean = mean(overhead_times),
          median = median(overhead_times),
          total = sum(overhead_times),
          pct_of_total = sum(overhead_times) / overall_timing$total_minutes * 100
        )
      } else {
        overhead_stats <- NULL
      }
    } else {
      overhead_stats <- NULL
    }
    
    # Compile timing list
    timing_list <- list(
      overall = overall_timing,
      iteration_stats = iteration_stats,
      fs_stats = fs_stats,
      overhead_stats = overhead_stats
    )
    
    # Create themed timing table
    timing_table_gt <- format_timing_enhanced(
      timing_list = timing_list,
      nb_boots = nb_boots,
      boot_success_rate = boot_success_rate,
      add_progress_bars = TRUE
    )
  }

  # ===========================================================================
  # SECTION 6: CREATE DIAGNOSTIC PLOTS (IF REQUESTED)
  # ===========================================================================
  
  plots <- NULL
  if (create_plots && requireNamespace("ggplot2", quietly = TRUE)) {
    
    plots <- create_bootstrap_diagnostic_plots(
      results = results,
      H_estimates = H_estimates,
      Hc_estimates = Hc_estimates
    )
    
    # Add timing plots if data available
    if (has_timing && has_iteration_timing) {
      timing_plots <- create_bootstrap_timing_plots(results)
      plots$timing <- timing_plots
    }
    
    # Add combined plot if patchwork is available
    if (requireNamespace("patchwork", quietly = TRUE)) {
      plots$combined <- (plots$H_distribution | plots$Hc_distribution) +
        patchwork::plot_annotation(
          title = "Bootstrap Distributions",
          subtitle = sprintf("%d iterations, %.1f%% successful",
                             nrow(results), boot_success_rate * 100),
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(size = 16, face = "bold")
          )
        )
      
      # Add timing plots to combined if available
      if (!is.null(plots$timing) && !is.null(plots$timing$iteration_time)) {
        plots$combined_with_timing <- 
          (plots$H_distribution | plots$Hc_distribution) /
          (plots$timing$iteration_time | plots$timing$search_time) +
          patchwork::plot_annotation(
            title = "Bootstrap Analysis: Distributions and Timing",
            subtitle = sprintf(
              "%d iterations, %.1f%% successful, %.1f min total",
              nrow(results), 
              boot_success_rate * 100,
              overall_timing$total_minutes
            ),
            theme = ggplot2::theme(
              plot.title = ggplot2::element_text(size = 16, face = "bold")
            )
          )
      }
    }
  }

  # ===========================================================================
  # SECTION 7: SUMMARIZE BOOTSTRAP SUBGROUPS (IF AVAILABLE)
  # ===========================================================================
  
  subgroup_summary <- NULL
  if (!is.null(results) && "Pcons" %in% names(results)) {
    
    tryCatch({
      subgroup_summary <- summarize_bootstrap_subgroups(
        results = results,
        nb_boots = nb_boots,
        original_sg = sgharm,
        maxk = 2  # Or extract from boot_results if stored
      )
    }, error = function(e) {
      warning("Could not create subgroup summary: ", e$message)
      NULL
    })
  }

  # ===========================================================================
  # SECTION 8: PRINT COMPREHENSIVE SUMMARY
  # ===========================================================================
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("           BOOTSTRAP ANALYSIS SUMMARY (THEMED)                 \n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  # Success metrics
  cat("BOOTSTRAP SUCCESS METRICS:\n")
  cat("─────────────────────────────────────────────────────────────\n")
  cat(sprintf("  Total iterations:              %d\n", diagnostics$n_boots))
  cat(sprintf("  Successful subgroup ID:        %d (%.1f%%)\n",
              diagnostics$n_successful, diagnostics$success_rate * 100))
  cat(sprintf("  Failed to find subgroup:       %d (%.1f%%)\n",
              diagnostics$n_failed, (1 - diagnostics$success_rate) * 100))
  cat("\n")
  
  # Timing summary (if available)
  if (has_timing) {
    cat("TIMING ANALYSIS:\n")
    cat("─────────────────────────────────────────────────────────────\n")
    cat("Overall:\n")
    cat(sprintf("  Total bootstrap time:          %.2f minutes (%.2f hours)\n",
                overall_timing$total_minutes, overall_timing$total_hours))
    cat(sprintf("  Average per iteration:         %.2f min (%.1f sec)\n",
                overall_timing$avg_minutes_per_boot,
                overall_timing$avg_seconds_per_boot))
    
    # Projected times
    if (nb_boots < 1000) {
      projected_1000 <- overall_timing$avg_minutes_per_boot * 1000
      cat(sprintf("  Projected for 1000 boots:      %.2f minutes (%.2f hours)\n",
                  projected_1000, projected_1000 / 60))
    }
    cat("\n")
    
    # Per-iteration timing
    if (!is.null(iteration_stats)) {
      cat("Per-iteration timing:\n")
      cat(sprintf("  Mean:                          %.2f min (%.1f sec)\n",
                  iteration_stats$mean, iteration_stats$mean * 60))
      cat(sprintf("  Median:                        %.2f min (%.1f sec)\n",
                  iteration_stats$median, iteration_stats$median * 60))
      cat(sprintf("  Range:                         [%.2f, %.2f] minutes\n",
                  iteration_stats$min, iteration_stats$max))
      cat("\n")
    }
    
    # ForestSearch timing
    if (!is.null(fs_stats)) {
      cat("ForestSearch timing (successful iterations only):\n")
      cat(sprintf("  Iterations with FS:            %d (%.1f%%)\n",
                  fs_stats$n_runs, fs_stats$pct_runs))
      cat(sprintf("  Mean FS time:                  %.2f min (%.1f sec)\n",
                  fs_stats$mean, fs_stats$mean * 60))
      cat(sprintf("  Total FS time:                 %.2f minutes\n",
                  fs_stats$total))
      cat(sprintf("  FS as %% of total:              %.1f%%\n",
                  fs_stats$total / overall_timing$total_minutes * 100))
      cat("\n")
    }
    
    # Overhead timing
    if (!is.null(overhead_stats)) {
      cat("Overhead timing (Cox models, bias correction, etc.):\n")
      cat(sprintf("  Mean overhead:                 %.2f min (%.1f sec)\n",
                  overhead_stats$mean, overhead_stats$mean * 60))
      cat(sprintf("  Total overhead:                %.2f minutes\n",
                  overhead_stats$total))
      cat(sprintf("  Overhead as %% of total:        %.1f%%\n",
                  overhead_stats$pct_of_total))
      cat("\n")
    }
    
    # Performance assessment
    cat("PERFORMANCE ASSESSMENT:\n")
    cat("─────────────────────────────────────────────────────────────\n")
    
    if (!is.null(iteration_stats)) {
      avg_sec <- overall_timing$avg_seconds_per_boot
      
      if (avg_sec < 5) {
        performance <- "Excellent ✓✓✓"
      } else if (avg_sec < 15) {
        performance <- "Good ✓✓"
      } else if (avg_sec < 30) {
        performance <- "Acceptable ✓"
      } else if (avg_sec < 60) {
        performance <- "Slow ⚠"
      } else {
        performance <- "Very Slow ⚠⚠"
      }
      
      cat(sprintf("  Performance rating:            %s\n", performance))
      cat(sprintf("  Average iteration speed:       %.1f seconds\n", avg_sec))
      
      # Recommendations based on performance
      if (avg_sec > 30) {
        cat("\n  Recommendations for improvement:\n")
        if (!is.null(fs_stats) && fs_stats$mean > 0.5) {
          cat("    • Consider reducing max.minutes in forestsearch\n")
          cat("    • Consider reducing maxk if currently > 2\n")
        }
        if (nb_boots > 500) {
          cat("    • Consider reducing nb_boots for initial testing\n")
        }
        cat("    • Ensure sufficient parallel workers are allocated\n")
      }
    }
    cat("\n")
  }
  
  # Subgroup summary (if available)
  if (!is.null(subgroup_summary)) {
    cat("SUBGROUP ANALYSIS:\n")
    cat("─────────────────────────────────────────────────────────────\n")
    cat(sprintf("  Subgroups identified:          %d (%.1f%%)\n",
                subgroup_summary$n_found,
                subgroup_summary$pct_found))
    cat("\n")
  }
  
  # Theme information
  cat("THEMED TABLES:\n")
  cat("─────────────────────────────────────────────────────────────\n")
  cat("  ✓ Results table               (apply_results_theme)\n")
  cat("  ✓ Diagnostics table            (apply_diagnostics_theme)\n")
  if (!is.null(timing_table_gt)) {
    cat("  ✓ Timing table                 (apply_timing_theme)\n")
  }
  cat("\nAll tables feature consistent styling:\n")
  cat("  • Professional borders and spacing\n")
  cat("  • Color-coded success indicators\n")
  cat("  • Interactive hover effects\n")
  cat("  • Accessibility features (WCAG)\n")
  cat("  • Publication-ready formatting\n")
  cat("\n")
  
  cat("═══════════════════════════════════════════════════════════════\n\n")

  # ===========================================================================
  # SECTION 9: COMPILE AND RETURN OUTPUT
  # ===========================================================================
  
  output <- list(
    # Main themed tables
    table = formatted_table,
    diagnostics_table = diagnostics_table_gt,
    
    # Diagnostics data
    diagnostics = diagnostics,
    
    # Timing (if available)
    timing = if (has_timing) {
      list(
        overall = overall_timing,
        iteration_stats = iteration_stats,
        fs_stats = fs_stats,
        overhead_stats = overhead_stats
      )
    } else NULL,
    
    timing_table = timing_table_gt,
    
    # Plots (if requested)
    plots = plots,
    
    # Subgroup summary (if available)
    subgroup_summary = subgroup_summary,
    
    # Metadata
    metadata = list(
      nb_boots = nb_boots,
      boot_success_rate = boot_success_rate,
      est_scale = est.scale,
      has_timing = has_timing,
      has_plots = !is.null(plots),
      has_subgroup_summary = !is.null(subgroup_summary),
      timestamp = Sys.time()
    )
  )
  
  # Set class for print method
  class(output) <- c("bootstrap_summary_themed", "list")
  
  invisible(output)
}


# ==============================================================================
# PRINT METHOD FOR THEMED SUMMARY
# ==============================================================================

#' Print Method for Themed Bootstrap Summary
#' 
#' @param x Object of class "bootstrap_summary_themed"
#' @param ... Additional arguments (ignored)
#' @export
print.bootstrap_summary_themed <- function(x, ...) {
  
  cat("\n")
  cat("Bootstrap Analysis Summary (Themed)\n")
  cat("═══════════════════════════════════════════\n\n")
  
  # Success metrics
  cat(sprintf("Iterations:        %d (%.1f%% successful)\n",
              x$metadata$nb_boots,
              x$metadata$boot_success_rate * 100))
  
  # Timing
  if (x$metadata$has_timing) {
    cat(sprintf("Total time:        %.2f minutes\n",
                x$timing$overall$total_minutes))
    cat(sprintf("Avg per boot:      %.1f seconds\n",
                x$timing$overall$avg_seconds_per_boot))
  }
  
  # Components
  cat("\nAvailable components:\n")
  cat("  • $table                 - Themed results table\n")
  cat("  • $diagnostics_table     - Themed diagnostics\n")
  if (!is.null(x$timing_table)) {
    cat("  • $timing_table          - Themed timing analysis\n")
  }
  if (x$metadata$has_plots) {
    cat("  • $plots                 - Diagnostic plots\n")
  }
  if (x$metadata$has_subgroup_summary) {
    cat("  • $subgroup_summary      - Subgroup analysis\n")
  }
  cat("  • $diagnostics           - Diagnostic metrics\n")
  cat("  • $timing                - Timing data\n")
  cat("  • $metadata              - Analysis metadata\n")
  
  cat("\nUse gt::gtsave() to export tables\n")
  cat("═══════════════════════════════════════════\n\n")
  
  invisible(x)
}


# ==============================================================================
# CONVENIENCE FUNCTION: EXPORT ALL THEMED TABLES
# ==============================================================================

#' Export All Themed Bootstrap Tables
#' 
#' Convenience function to save all themed tables from a bootstrap summary
#' to files in a specified output directory.
#' 
#' @param summary Output from summarize_bootstrap_results_themed()
#' @param output_dir Directory path for output files (default: "bootstrap_output")
#' @param prefix Filename prefix (default: "bootstrap")
#' @param format Output format: "html", "png", or "both" (default: "html")
#' @param overwrite Logical. Overwrite existing files? (default: FALSE)
#' 
#' @return Invisible list of file paths created
#' @export
#' 
#' @examples
#' \dontrun{
#' summary <- summarize_bootstrap_results_themed(...)
#' 
#' # Export all tables as HTML
#' export_themed_tables(summary, output_dir = "results", format = "html")
#' 
#' # Export as both HTML and PNG
#' export_themed_tables(summary, format = "both")
#' }
export_themed_tables <- function(summary,
                                 output_dir = "bootstrap_output",
                                 prefix = "bootstrap",
                                 format = c("html", "png", "both"),
                                 overwrite = FALSE) {
  
  format <- match.arg(format)
  
  if (!inherits(summary, "bootstrap_summary_themed")) {
    stop("summary must be output from summarize_bootstrap_results_themed()")
  }
  
  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }
  
  files_created <- character()
  
  # Helper function to save table
  save_table <- function(table, name, format) {
    if (is.null(table)) return(NULL)
    
    base_path <- file.path(output_dir, paste0(prefix, "_", name))
    
    if (format %in% c("html", "both")) {
      html_path <- paste0(base_path, ".html")
      if (file.exists(html_path) && !overwrite) {
        warning("File exists (use overwrite=TRUE): ", html_path)
      } else {
        gt::gtsave(table, html_path)
        files_created <<- c(files_created, html_path)
        message("Created: ", html_path)
      }
    }
    
    if (format %in% c("png", "both")) {
      if (!requireNamespace("webshot2", quietly = TRUE)) {
        warning("webshot2 package required for PNG export. Skipping PNG files.")
      } else {
        png_path <- paste0(base_path, ".png")
        if (file.exists(png_path) && !overwrite) {
          warning("File exists (use overwrite=TRUE): ", png_path)
        } else {
          gt::gtsave(table, png_path)
          files_created <<- c(files_created, png_path)
          message("Created: ", png_path)
        }
      }
    }
  }
  
  # Save each table
  save_table(summary$table, "results", format)
  save_table(summary$diagnostics_table, "diagnostics", format)
  save_table(summary$timing_table, "timing", format)
  
  message("\nExport complete! ", length(files_created), " file(s) created.")
  
  invisible(files_created)
}


# ==============================================================================
# CONVENIENCE FUNCTION: COMPARE OLD VS NEW
# ==============================================================================

#' Compare Original vs Themed Bootstrap Summaries
#' 
#' Helper function to compare output from summarize_bootstrap_results()
#' and summarize_bootstrap_results_themed() side-by-side.
#' 
#' @param boot_results Output from forestsearch_bootstrap_dofuture()
#' @param sgharm Original subgroup definition
#' 
#' @return List with both summaries for comparison
#' @export
#' 
#' @examples
#' \dontrun{
#' comparison <- compare_bootstrap_summaries(boot_results, sgharm)
#' 
#' # View original
#' comparison$original$table
#' 
#' # View themed
#' comparison$themed$table
#' }
compare_bootstrap_summaries <- function(boot_results, sgharm = NULL) {
  
  message("Creating original summary...")
  original <- summarize_bootstrap_results(
    sgharm = sgharm,
    boot_results = boot_results,
    create_plots = FALSE
  )
  
  message("Creating themed summary...")
  themed <- summarize_bootstrap_results_themed(
    sgharm = sgharm,
    boot_results = boot_results,
    create_plots = FALSE
  )
  
  message("\nComparison complete!")
  message("Use $original and $themed to view each version")
  
  list(
    original = original,
    themed = themed
  )
}
