# Quick Fix for Missing Links in ForestSearch Documentation
# Run this from your package directory

# Function to fix all missing links at once
fix_all_missing_links <- function() {
  
  # List of files and their missing links
  fixes <- list(
    "R/find_k_inter_for_target_hr.R" = c("find_k_inter_grid_search", "find_k_inter_batch"),
    "R/forestsearch.R" = c("bootstrap_forestsearch", "summarize_forestsearch"),
    "R/generate_aft_dgm_flex.R" = c("find_quantile_for_proportion")
  )
  
  for (file in names(fixes)) {
    if (file.exists(file)) {
      content <- readLines(file)
      original <- content
      
      # Fix each missing link in this file
      for (func in fixes[[file]]) {
        # Replace \link{func} with \code{func}
        pattern1 <- paste0("\\\\link\\{", func, "\\}")
        replacement1 <- paste0("\\\\code{", func, "}")
        content <- gsub(pattern1, replacement1, content, fixed = TRUE)
        
        # Also handle cases with spaces
        pattern2 <- paste0("\\\\link\\{\\s*", func, "\\s*\\}")
        content <- gsub(pattern2, replacement1, content)
      }
      
      # Check if changes were made
      if (!identical(content, original)) {
        # Backup original
        writeLines(original, paste0(file, ".bak"))
        # Write fixed version
        writeLines(content, file)
        cat("✓ Fixed", file, "\n")
      } else {
        cat("ℹ No changes needed in", file, "\n")
      }
    } else {
      cat("✗ File not found:", file, "\n")
    }
  }
  
  cat("\n✓ Done! Now run:\n")
  cat("  devtools::document()\n")
  cat("  devtools::check()\n")
}

# Alternative: Remove the @seealso sections entirely
remove_broken_seealso <- function() {
  
  files <- c(
    "R/find_k_inter_for_target_hr.R",
    "R/forestsearch.R", 
    "R/generate_aft_dgm_flex.R"
  )
  
  for (file in files) {
    if (file.exists(file)) {
      content <- readLines(file)
      original <- content
      
      # Remove lines with @seealso containing the broken links
      problem_patterns <- c(
        "find_k_inter_grid_search", "find_k_inter_batch",
        "bootstrap_forestsearch", "summarize_forestsearch",
        "find_quantile_for_proportion"
      )
      
      # Find and remove @seealso lines with these functions
      for (i in seq_along(content)) {
        if (grepl("@seealso", content[i])) {
          if (any(sapply(problem_patterns, function(p) grepl(p, content[i])))) {
            content[i] <- ""  # Remove the line
          }
        }
      }
      
      # Remove empty lines where @seealso was removed
      content <- content[content != ""]
      
      if (!identical(content, original)) {
        writeLines(original, paste0(file, ".bak"))
        writeLines(content, file)
        cat("✓ Removed @seealso from", file, "\n")
      }
    }
  }
  
  cat("\n✓ Done! Now run:\n")
  cat("  devtools::document()\n")
  cat("  devtools::check()\n")
}

# Check which functions actually exist
check_missing_functions <- function() {
  missing <- c(
    "find_k_inter_grid_search",
    "find_k_inter_batch",
    "bootstrap_forestsearch",
    "summarize_forestsearch",
    "find_quantile_for_proportion"
  )
  
  cat("Checking if functions exist in package...\n\n")
  
  r_files <- list.files("R/", pattern = "\\.R$", full.names = TRUE)
  all_functions <- character()
  
  for (file in r_files) {
    content <- readLines(file, warn = FALSE)
    # Look for function definitions
    func_defs <- grep("^[a-zA-Z_][a-zA-Z0-9_.]*\\s*<-\\s*function", content, value = TRUE)
    if (length(func_defs) > 0) {
      func_names <- trimws(sub("\\s*<-.*", "", func_defs))
      all_functions <- c(all_functions, func_names)
    }
  }
  
  for (func in missing) {
    if (func %in% all_functions) {
      cat("✓", func, "EXISTS in package\n")
    } else {
      cat("✗", func, "NOT FOUND\n")
      
      # Look for similar names
      similar <- grep(substr(func, 1, 10), all_functions, value = TRUE)
      if (length(similar) > 0) {
        cat("  Possible matches:", paste(similar, collapse = ", "), "\n")
      }
    }
  }
}

# ================================================================================
# RUN THE FIX
# ================================================================================

cat("================================================================================\n")
cat("MISSING LINKS FIX FOR FORESTSEARCH PACKAGE\n")
cat("================================================================================\n\n")

# First, check what exists
check_missing_functions()

cat("\n--------------------------------------------------------------------------------\n")
cat("Choose your fix approach:\n")
cat("1. fix_all_missing_links()  # Changes \\link{} to \\code{} (recommended)\n")
cat("2. remove_broken_seealso()  # Removes the @seealso lines entirely\n")
cat("\nRun your choice, then devtools::document()\n")
cat("--------------------------------------------------------------------------------\n")
