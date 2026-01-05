# Fixed Version of Bootstrap Summary Helpers
# CRAN-compliant version with proper documentation, exports, and Unicode fixes
#
# Fixes applied:
# - Unicode characters replaced with ASCII alternatives
# - Proper @importFrom declarations
# - @export tags added
# - data.table handling improved

#' Summarize Bootstrap Subgroup Results
#'
#' Analyzes bootstrap replication results to summarize the stability and
#' characteristics of identified subgroups across bootstrap iterations.
#'
#' @param results A data.table or data.frame containing bootstrap results with
#'   subgroup characteristics. Expected columns include \code{Pcons} (consistency),
#'   \code{M.1}, \code{M.2}, etc. (factor definitions), \code{hr_sg} (subgroup HR),
#'   and \code{N_sg} (subgroup size).
#' @param nb_boots Integer. Total number of bootstrap iterations attempted.
#' @param original_sg Character vector. Original subgroup definition from main
#'   analysis (e.g., \code{c("M.1", "M.2")} values from \code{fs.est}).
#' @param maxk Integer. Maximum number of factors allowed in subgroup definition.
#'   Default is 2.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{basic_stats}{data.table with summary statistics (success rate, Pcons, HR)}
#'   \item{consistency_dist}{data.table with distribution of consistency values}
#'   \item{size_dist}{data.table with distribution of subgroup sizes}
#'   \item{factor_freq}{data.table with frequency of each factor appearing}
#'   \item{agreement}{data.table with frequency of subgroup definitions}
#'   \item{factor_presence}{data.table with base factor frequencies}
#'   \item{factor_presence_specific}{data.table with specific factor definition frequencies}
#'   \item{original_agreement}{data.table comparing bootstrap results to original}
#'   \item{n_found}{Integer. Number of successful iterations}
#'   \item{pct_found}{Numeric. Percentage of successful iterations}
#' }
#'
#' @examples
#' \donttest{
#' # Create example bootstrap results
#' boot_results <- data.frame(
#'   Pcons = runif(50, 0.5, 1),
#'   M.1 = sample(c("age<=50", "biomarker<=2"), 50, replace = TRUE),
#'   M.2 = sample(c("", "stage=1"), 50, replace = TRUE),
#'   K_sg = sample(1:2, 50, replace = TRUE),
#'   hr_sg = runif(50, 0.5, 1.5),
#'   N_sg = sample(50:200, 50, replace = TRUE)
#' )
#'
#' summary <- summarize_bootstrap_subgroups(
#'   results = boot_results,
#'   nb_boots = 50,
#'   original_sg = c("age<=50", "stage=1"),
#'   maxk = 2
#' )
#'
#' print(summary$basic_stats)
#' }
#'
#' @importFrom data.table data.table as.data.table copy setnames setcolorder .N
#' @importFrom stats quantile sd
#' @export
summarize_bootstrap_subgroups <- function(results,
                                          nb_boots,
                                          original_sg = NULL,
                                          maxk = 2) {

  # ==========================================================================
  # Ensure results is a proper data.table
  # ==========================================================================

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install with: install.packages('data.table')")
  }

  # Convert to data.table if needed
 if (!inherits(results, "data.table")) {
    if (is.data.frame(results)) {
      results <- data.table::as.data.table(results)
    } else if (is.matrix(results)) {
      results <- data.table::as.data.table(as.data.frame(results))
    } else {
      stop("results must be a data.frame, data.table, or matrix")
    }
  }

  # Make a copy to avoid modifying the original
  results <- data.table::copy(results)

  # ==========================================================================
  # SECTION 1: Filter to Successful Iterations
  # ==========================================================================

  if ("Pcons" %in% names(results)) {
    sg_found <- results[!is.na(results$Pcons), ]
  } else {
    warning("Column 'Pcons' not found in results")
    sg_found <- results
  }

  n_found <- nrow(sg_found)
  pct_found <- 100 * n_found / nb_boots

  if (n_found == 0) {
    warning("No subgroups identified in any bootstrap iteration")
    return(list(
      basic_stats = NULL,
      consistency_dist = NULL,
      size_dist = NULL,
      factor_freq = NULL,
      agreement = NULL,
      factor_presence = NULL,
      factor_presence_specific = NULL,
      original_agreement = NULL,
      n_found = 0,
      pct_found = 0
    ))
  }

  # ==========================================================================
  # SECTION 2: Basic Statistics Table
  # ==========================================================================

  basic_stats_list <- list(
    Metric = character(),
    Value = character()
  )

  # Add basic counts
  basic_stats_list$Metric <- c(
    "Total bootstrap iterations",
    "Subgroups identified",
    "Success rate (%)"
  )
  basic_stats_list$Value <- c(
    as.character(nb_boots),
    as.character(n_found),
    sprintf("%.1f%%", pct_found)
  )

  # Add Pcons statistics if available
  if ("Pcons" %in% names(sg_found) && sum(!is.na(sg_found$Pcons)) > 0) {
    pcons_vals <- sg_found$Pcons[!is.na(sg_found$Pcons)]
    basic_stats_list$Metric <- c(basic_stats_list$Metric,
                                  "",
                                  "Consistency (Pcons)",
                                  "  Mean",
                                  "  Median",
                                  "  SD",
                                  "  Min",
                                  "  Max",
                                  "  Q25",
                                  "  Q75")
    basic_stats_list$Value <- c(basic_stats_list$Value,
                                 "",
                                 "",
                                 sprintf("%.3f", mean(pcons_vals)),
                                 sprintf("%.3f", median(pcons_vals)),
                                 sprintf("%.3f", sd(pcons_vals)),
                                 sprintf("%.3f", min(pcons_vals)),
                                 sprintf("%.3f", max(pcons_vals)),
                                 sprintf("%.3f", quantile(pcons_vals, 0.25)),
                                 sprintf("%.3f", quantile(pcons_vals, 0.75)))
  }

  # Add hr_sg statistics if available
  if ("hr_sg" %in% names(sg_found) && sum(!is.na(sg_found$hr_sg)) > 0) {
    hr_vals <- sg_found$hr_sg[!is.na(sg_found$hr_sg)]
    basic_stats_list$Metric <- c(basic_stats_list$Metric,
                                  "",
                                  "Hazard Ratio (hr_sg)",
                                  "  Mean",
                                  "  Median",
                                  "  SD",
                                  "  Min",
                                  "  Max",
                                  "  Q25",
                                  "  Q75")
    basic_stats_list$Value <- c(basic_stats_list$Value,
                                 "",
                                 "",
                                 sprintf("%.2f", mean(hr_vals)),
                                 sprintf("%.2f", median(hr_vals)),
                                 sprintf("%.2f", sd(hr_vals)),
                                 sprintf("%.2f", min(hr_vals)),
                                 sprintf("%.2f", max(hr_vals)),
                                 sprintf("%.2f", quantile(hr_vals, 0.25)),
                                 sprintf("%.2f", quantile(hr_vals, 0.75)))
  }

  # Convert to data.table
  basic_stats <- data.table::data.table(
    Metric = basic_stats_list$Metric,
    Value = basic_stats_list$Value
  )

  # ==========================================================================
  # SECTION 3: Factor Frequency Table
  # ==========================================================================

  factor_cols <- paste0("M.", seq_len(maxk))
  factor_freq_list <- list()

  for (i in seq_len(maxk)) {
    col <- paste0("M.", i)
    if (col %in% names(sg_found)) {
      valid_rows <- !is.na(sg_found[[col]]) & sg_found[[col]] != ""
      if (sum(valid_rows) > 0) {
        factor_vals <- sg_found[[col]][valid_rows]
        freq_table <- table(factor_vals)

        freq <- data.table::data.table(
          Factor = names(freq_table),
          N = as.integer(freq_table),
          Position = paste0("M.", i),
          Percent = 100 * as.integer(freq_table) / n_found
        )

        freq <- freq[order(-freq$N), ]
        factor_freq_list[[i]] <- freq
      }
    }
  }

  # Combine data.tables
  if (length(factor_freq_list) > 0) {
    factor_freq <- data.table::rbindlist(factor_freq_list, fill = TRUE)
    data.table::setcolorder(factor_freq, c("Position", "Factor", "N", "Percent"))
  } else {
    factor_freq <- data.table::data.table(
      Position = character(),
      Factor = character(),
      N = integer(),
      Percent = numeric()
    )
  }

  # ==========================================================================
  # SECTION 4: Subgroup Definition Agreement
  # ==========================================================================

  agreement <- NULL

  if (maxk == 1) {
    if ("M.1" %in% names(sg_found)) {
      valid_rows <- !is.na(sg_found$M.1) & sg_found$M.1 != ""
      if (sum(valid_rows) > 0) {
        freq_table <- table(sg_found$M.1[valid_rows])
        agreement <- data.table::data.table(
          Subgroup = names(freq_table),
          K_sg = 1L,
          N = as.integer(freq_table)
        )
      }
    }
  } else if (maxk == 2) {
    if (all(c("M.1", "M.2") %in% names(sg_found))) {
      valid_rows <- !is.na(sg_found$M.1) & sg_found$M.1 != ""
      if (sum(valid_rows) > 0) {
        sg_temp <- sg_found[valid_rows, ]
        sg_temp$Subgroup <- ifelse(
          is.na(sg_temp$M.2) | sg_temp$M.2 == "",
          as.character(sg_temp$M.1),
          paste(sg_temp$M.1, "&", sg_temp$M.2)
        )

        # Count combinations using data.table
        if ("K_sg" %in% names(sg_temp)) {
          agreement <- sg_temp[, list(N = .N), by = list(Subgroup, K_sg)]
        } else {
          agreement <- sg_temp[, list(N = .N), by = list(Subgroup)]
          agreement$K_sg <- NA_integer_
        }
      }
    }
  }

  # Add percentage and sort
  if (!is.null(agreement) && nrow(agreement) > 0) {
    agreement$Percent_of_successful <- 100 * agreement$N / n_found
    agreement <- agreement[order(-agreement$N), ]
    agreement$Rank <- seq_len(nrow(agreement))
  }

  # ==========================================================================
  # SECTION 5: Consistency Distribution
  # ==========================================================================

  consistency_dist <- NULL

  if ("Pcons" %in% names(sg_found) && sum(!is.na(sg_found$Pcons)) > 0) {
    # Create bins - FIXED: replaced corrupted Unicode with ASCII
    sg_found$Pcons_bin <- cut(
      sg_found$Pcons,
      breaks = c(0, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0),
      labels = c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.9", "0.9-0.95", ">=0.95"),
      include.lowest = TRUE
    )

    if (sum(!is.na(sg_found$Pcons_bin)) > 0) {
      freq_table <- table(sg_found$Pcons_bin[!is.na(sg_found$Pcons_bin)])
      consistency_dist <- data.table::data.table(
        Consistency_Range = names(freq_table),
        Count = as.integer(freq_table),
        Percent = 100 * as.integer(freq_table) / n_found
      )
    }
  }

  # ==========================================================================
  # SECTION 6: Size Distribution
  # ==========================================================================

  size_dist <- NULL

  if ("N_sg" %in% names(sg_found) && sum(!is.na(sg_found$N_sg)) > 0) {
    # Create bins - FIXED: replaced corrupted Unicode with ASCII
    size_breaks <- c(0, 50, 100, 150, 200, 300, Inf)
    size_labels <- c("<50", "50-99", "100-149", "150-199", "200-299", ">=300")

    sg_found$N_sg_bin <- cut(
      sg_found$N_sg,
      breaks = size_breaks,
      labels = size_labels,
      include.lowest = TRUE
    )

    if (sum(!is.na(sg_found$N_sg_bin)) > 0) {
      freq_table <- table(sg_found$N_sg_bin[!is.na(sg_found$N_sg_bin)])
      size_dist <- data.table::data.table(
        Size_Range = names(freq_table),
        Count = as.integer(freq_table),
        Percent = 100 * as.integer(freq_table) / n_found
      )
    }
  }

  # ==========================================================================
  # SECTION 7: Factor Presence Analysis
  # ==========================================================================

  factor_presence_results <- NULL

  if (n_found > 0) {
    tryCatch({
      factor_presence_results <- summarize_factor_presence(
        sg_found,
        maxk = maxk,
        threshold = 20
      )
    }, error = function(e) {
      warning("Could not summarize factor presence: ", e$message)
      NULL
    })
  }

  # ==========================================================================
  # SECTION 8: Original Agreement Calculation
  # ==========================================================================

  original_agreement <- NULL

  if (!is.null(original_sg) && n_found > 0) {
    # Find the subgroup column
    subgroup_col <- NULL
    for (possible_col in c("M.1", "Subgroup", "H_bootstrap", "H", "H_star")) {
      if (possible_col %in% names(sg_found)) {
        subgroup_col <- possible_col
        break
      }
    }

    if (!is.null(subgroup_col)) {
      orig_sg_char <- as.character(original_sg[1])  # Use first element
      matches <- sum(as.character(sg_found[[subgroup_col]]) == orig_sg_char, na.rm = TRUE)

      original_agreement <- data.table::data.table(
        Metric = c(
          "Total bootstrap iterations",
          "Successful iterations",
          "Failed iterations (no subgroup)",
          "Exact match with original",
          "Different from original"
        ),
        Value = c(
          as.character(nb_boots),
          as.character(n_found),
          as.character(nb_boots - n_found),
          sprintf("%d (%.1f%%)", matches, 100 * matches / n_found),
          sprintf("%d (%.1f%%)", n_found - matches, 100 * (n_found - matches) / n_found)
        )
      )
    }
  }

  # ==========================================================================
  # Return Compiled Results
  # ==========================================================================

  list(
    basic_stats = basic_stats,
    consistency_dist = consistency_dist,
    size_dist = size_dist,
    factor_freq = factor_freq,
    agreement = agreement,
    factor_presence = if (!is.null(factor_presence_results)) factor_presence_results$base_factors else NULL,
    factor_presence_specific = if (!is.null(factor_presence_results)) factor_presence_results$specific_factors else NULL,
    original_agreement = original_agreement,
    n_found = n_found,
    pct_found = pct_found
  )
}


#' Summarize Factor Presence in Bootstrap Results
#'
#' Analyzes which base factors and specific factor definitions appear
#' most frequently across bootstrap iterations.
#'
#' @param results data.table. Bootstrap results with M.1, M.2, etc. columns
#' @param maxk Integer. Maximum number of factors. Default is 2.
#' @param threshold Numeric. Minimum percentage to include in specific factors
#'   table. Default is 20.
#'
#' @return A list with:
#' \describe{
#'   \item{base_factors}{data.table of base factor frequencies}
#'   \item{specific_factors}{data.table of specific factor definition frequencies}
#'   \item{threshold}{The threshold value used}
#' }
#'
#' @importFrom data.table data.table setcolorder .N
#' @keywords internal
summarize_factor_presence <- function(results, maxk = 2, threshold = 20) {

  # Ensure data.table
  if (!inherits(results, "data.table")) {
    results <- data.table::as.data.table(results)
  }

  # Filter successful iterations
  if ("Pcons" %in% names(results)) {
    sg_found <- results[!is.na(results$Pcons), ]
  } else {
    sg_found <- results
  }

  n_found <- nrow(sg_found)
  if (n_found == 0) return(NULL)

  # Extract all factors
  all_factors <- character()
  for (i in seq_len(maxk)) {
    col <- paste0("M.", i)
    if (col %in% names(sg_found)) {
      factors <- sg_found[[col]][!is.na(sg_found[[col]]) & sg_found[[col]] != ""]
      all_factors <- c(all_factors, as.character(factors))
    }
  }

  if (length(all_factors) == 0) return(NULL)

  # Function to extract base factor name
  extract_base_factor <- function(factor_string) {
    factor_string <- trimws(factor_string)
    factor_string <- gsub("^!\\{", "{", factor_string)
    factor_string <- gsub("[{}]", "", factor_string)
    base_name <- sub("^([a-zA-Z_][a-zA-Z0-9_.]*).*", "\\1", factor_string)
    return(base_name)
  }

  # Count base factors
  base_factors <- vapply(all_factors, extract_base_factor,
                          FUN.VALUE = character(1), USE.NAMES = FALSE)
  base_factor_counts <- table(base_factors)

  base_factor_summary <- data.table::data.table(
    Factor = names(base_factor_counts),
    Count = as.integer(base_factor_counts),
    Percent = 100 * as.integer(base_factor_counts) / n_found
  )

  base_factor_summary <- base_factor_summary[order(-base_factor_summary$Count), ]
  base_factor_summary$Rank <- seq_len(nrow(base_factor_summary))
  data.table::setcolorder(base_factor_summary, c("Rank", "Factor", "Count", "Percent"))

  # Specific factors
  specific_factor_counts <- table(all_factors)

  specific_factor_summary <- data.table::data.table(
    Factor_Definition = names(specific_factor_counts),
    Count = as.integer(specific_factor_counts),
    Percent = 100 * as.integer(specific_factor_counts) / n_found
  )

  # Filter by threshold
  specific_factor_summary <- specific_factor_summary[
    specific_factor_summary$Percent >= threshold,
  ]

  if (nrow(specific_factor_summary) > 0) {
    specific_factor_summary$Base_Factor <- vapply(
      specific_factor_summary$Factor_Definition,
      extract_base_factor,
      FUN.VALUE = character(1)
    )
    specific_factor_summary <- specific_factor_summary[
      order(specific_factor_summary$Base_Factor, -specific_factor_summary$Count),
    ]
    specific_factor_summary$Rank <- seq_len(nrow(specific_factor_summary))
    data.table::setcolorder(
      specific_factor_summary,
      c("Rank", "Base_Factor", "Factor_Definition", "Count", "Percent")
    )
  }

  return(list(
    base_factors = base_factor_summary,
    specific_factors = specific_factor_summary,
    threshold = threshold
  ))
}


#' Calculate Original Subgroup Agreement
#'
#' Compares bootstrap-identified subgroups to the original main analysis
#' subgroup to assess stability.
#'
#' @param bootstrap_results data.table. Bootstrap results
#' @param original_subgroup Character. Original subgroup definition
#' @param subgroup_column Character. Name of column containing subgroup
#'   definitions. If NULL, auto-detects.
#'
#' @return data.table with agreement metrics
#'
#' @importFrom data.table data.table
#' @keywords internal
calculate_original_agreement <- function(bootstrap_results,
                                          original_subgroup,
                                          subgroup_column = NULL) {

  if (is.null(original_subgroup)) return(NULL)

  # Auto-detect subgroup column
  if (is.null(subgroup_column)) {
    possible_cols <- c("M.1", "Subgroup", "H_bootstrap", "H", "H_star")
    for (col in possible_cols) {
      if (col %in% names(bootstrap_results)) {
        subgroup_column <- col
        break
      }
    }
  }

  if (is.null(subgroup_column)) {
    warning("Could not find subgroup column in bootstrap results")
    return(NULL)
  }

  # Filter to successful iterations
  if ("Pcons" %in% names(bootstrap_results)) {
    sg_found <- bootstrap_results[!is.na(bootstrap_results$Pcons), ]
  } else {
    sg_found <- bootstrap_results
  }

  n_total <- nrow(bootstrap_results)
  n_found <- nrow(sg_found)

  if (n_found == 0) {
    return(data.table::data.table(
      Metric = "No successful iterations",
      Value = "N/A"
    ))
  }

  # Compare to original
  orig_char <- as.character(original_subgroup[1])
  bootstrap_vals <- as.character(sg_found[[subgroup_column]])

  n_exact_match <- sum(bootstrap_vals == orig_char, na.rm = TRUE)
  n_partial_match <- sum(grepl(orig_char, bootstrap_vals, fixed = TRUE), na.rm = TRUE)

  data.table::data.table(
    Metric = c(
      "Total bootstrap iterations",
      "Successful iterations",
      "Failed iterations",
      "Exact match with original",
      "Partial match with original",
      "Different from original"
    ),
    Value = c(
      as.character(n_total),
      as.character(n_found),
      as.character(n_total - n_found),
      sprintf("%d (%.1f%%)", n_exact_match, 100 * n_exact_match / n_found),
      sprintf("%d (%.1f%%)", n_partial_match, 100 * n_partial_match / n_found),
      sprintf("%d (%.1f%%)", n_found - n_exact_match, 100 * (n_found - n_exact_match) / n_found)
    )
  )
}


# Backward compatibility aliases (deprecated - use main function names)
summarize_bootstrap_subgroups_fixed <- summarize_bootstrap_subgroups
summarize_factor_presence_fixed <- summarize_factor_presence
