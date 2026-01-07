# Complete Diagnostic and Fix for Documentation Issues
# This will find the actual files regardless of their names

fix_documentation_complete <- function() {
  cat("================================================================================\n")
  cat("COMPLETE DOCUMENTATION FIX FOR FORESTSEARCH\n") 
  cat("================================================================================\n\n")
  
  # Check we're in the right place
  if (!file.exists("DESCRIPTION")) {
    stop("ERROR: DESCRIPTION file not found. Please run from package root directory.")
  }
  
  # Step 1: Delete old .Rd files that might have cached problems
  cat("Step 1: Cleaning up old documentation files...\n")
  
  problem_rd_files <- c(
    "man/find_k_inter_for_target_hr.Rd",
    "man/forestsearch.Rd",
    "man/generate_aft_dgm_flex.Rd"
  )
  
  for (rd_file in problem_rd_files) {
    if (file.exists(rd_file)) {
      file.remove(rd_file)
      cat("  ✓ Deleted:", rd_file, "\n")
    }
  }
  
  # Step 2: Find ALL R files with these function references
  cat("\nStep 2: Searching all R files for documentation issues...\n")
  
  r_files <- list.files("R/", pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
  
  if (length(r_files) == 0) {
    stop("No R files found in R/ directory")
  }
  
  files_with_issues <- character()
  all_issues <- list()
  
  for (file in r_files) {
    content <- readLines(file, warn = FALSE)
    file_has_issue <- FALSE
    
    # Check for each problematic link
    problem_links <- c(
      "find_k_inter_grid_search",
      "find_k_inter_batch", 
      "bootstrap_forestsearch",
      "summarize_forestsearch",
      "find_quantile_for_proportion"
    )
    
    for (link in problem_links) {
      # Check various patterns
      patterns <- c(
        paste0("\\\\link\\{", link, "\\}"),
        paste0("\\\\link\\[.*\\]\\{", link, "\\}"),
        paste0("@seealso.*\\\\link\\{", link, "\\}")
      )
      
      for (pattern in patterns) {
        if (any(grepl(pattern, content))) {
          cat("  ✓ Found '", link, "' in: ", basename(file), "\n", sep = "")
          file_has_issue <- TRUE
          
          if (is.null(all_issues[[file]])) {
            all_issues[[file]] <- character()
          }
          all_issues[[file]] <- c(all_issues[[file]], link)
        }
      }
    }
    
    if (file_has_issue) {
      files_with_issues <- c(files_with_issues, file)
    }
  }
  
  # Step 3: Fix the issues
  if (length(files_with_issues) > 0) {
    cat("\nStep 3: Fixing", length(files_with_issues), "file(s)...\n\n")
    
    for (file in unique(files_with_issues)) {
      cat("Fixing:", basename(file), "\n")
      
      content <- readLines(file, warn = FALSE)
      original <- content
      
      # Fix each problematic link
      for (link in all_issues[[file]]) {
        # Replace \link{} with \code{}
        content <- gsub(
          paste0("\\\\link\\{", link, "\\}"),
          paste0("\\\\code{", link, "}"),
          content,
          fixed = TRUE
        )
        
        # Replace \link[anything]{} with \code{}
        content <- gsub(
          paste0("\\\\link\\[[^]]*\\]\\{", link, "\\}"),
          paste0("\\\\code{", link, "}"),
          content
        )
      }
      
      # Save if changed
      if (!identical(content, original)) {
        writeLines(original, paste0(file, ".bak"))
        writeLines(content, file)
        cat("  ✅ Fixed and saved\n")
      }
    }
    
    cat("\n✅ All files fixed!\n\n")
    cat("Next steps:\n")
    cat("1. devtools::document()  # Regenerate documentation\n")
    cat("2. devtools::check()     # Verify all issues resolved\n")
    
  } else {
    cat("\n✓ No documentation issues found in R files!\n\n")
    cat("The warnings might be from cached .Rd files.\n")
    cat("Run: devtools::document() to regenerate them.\n")
  }
}

# Alternative: Just regenerate everything
clean_and_rebuild <- function() {
  cat("Clean rebuild of all documentation...\n\n")
  
  # Remove all .Rd files
  rd_files <- list.files("man/", pattern = "\\.Rd$", full.names = TRUE)
  if (length(rd_files) > 0) {
    cat("Removing", length(rd_files), ".Rd files...\n")
    file.remove(rd_files)
  }
  
  # Regenerate
  cat("\nRegenerating documentation...\n")
  devtools::document()
  
  cat("\n✅ Documentation regenerated!\n")
  cat("Run devtools::check() to verify issues are resolved.\n")
}

# Quick diagnostic
diagnose_documentation <- function() {
  cat("DOCUMENTATION DIAGNOSTIC\n")
  cat("========================\n\n")
  
  # Check for .Rd files
  cat("Checking man/ directory:\n")
  rd_files <- list.files("man/", pattern = "\\.Rd$")
  cat("  Found", length(rd_files), ".Rd files\n")
  
  # Check for problematic .Rd files
  problem_rd <- c("find_k_inter_for_target_hr.Rd", "forestsearch.Rd", "generate_aft_dgm_flex.Rd")
  for (rd in problem_rd) {
    if (rd %in% rd_files) {
      cat("  ✓", rd, "exists\n")
      
      # Check content
      content <- readLines(file.path("man", rd), warn = FALSE)
      bad_links <- grep("\\\\link\\{(find_k_inter_grid_search|find_k_inter_batch|bootstrap_forestsearch|summarize_forestsearch|find_quantile_for_proportion)\\}", 
                       content, value = TRUE)
      if (length(bad_links) > 0) {
        cat("    ⚠️ Contains problematic links\n")
      }
    }
  }
  
  # Check R files
  cat("\nChecking R/ directory:\n")
  r_files <- list.files("R/", pattern = "\\.R$")
  cat("  Found", length(r_files), ".R files\n")
  
  # List them
  cat("\nR files present:\n")
  for (f in r_files) {
    cat("  -", f, "\n")
  }
}

# ================================================================================
# RUN DIAGNOSTIC AND FIX
# ================================================================================

cat("\n🔧 ForestSearch Documentation Complete Fix\n")
cat("==========================================\n\n")

# First run diagnostic
diagnose_documentation()

cat("\n--------------------------------------------------------------------------------\n")
cat("Choose your approach:\n")
cat("1. fix_documentation_complete()  # Smart fix (recommended)\n")
cat("2. clean_and_rebuild()           # Delete all .Rd and rebuild\n")
cat("3. diagnose_documentation()      # Just check status\n\n")

# Auto-run the fix
cat("Auto-running smart fix...\n\n")
fix_documentation_complete()
