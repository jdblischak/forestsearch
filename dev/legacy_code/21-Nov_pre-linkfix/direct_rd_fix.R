#!/usr/bin/env Rscript

# Direct fix for .Rd files with missing links
# This edits the .Rd files directly to fix the warnings

fix_rd_files_directly <- function() {
  cat("================================================================================\n")
  cat("DIRECTLY FIXING .Rd FILES WITH MISSING LINKS\n")
  cat("================================================================================\n\n")
  
  # The problematic .Rd files and their missing links
  files_to_fix <- list(
    "man/find_k_inter_for_target_hr.Rd" = c("find_k_inter_grid_search", "find_k_inter_batch"),
    "man/forestsearch.Rd" = c("bootstrap_forestsearch", "summarize_forestsearch"),
    "man/generate_aft_dgm_flex.Rd" = c("find_quantile_for_proportion")
  )
  
  files_fixed <- 0
  
  for (rd_file in names(files_to_fix)) {
    if (!file.exists(rd_file)) {
      cat("✗ File not found:", rd_file, "\n")
      next
    }
    
    cat("Processing:", rd_file, "\n")
    
    # Read the .Rd file
    content <- readLines(rd_file, warn = FALSE)
    original_content <- content
    
    # Fix each problematic link
    for (func in files_to_fix[[rd_file]]) {
      # Pattern 1: \link{function}
      pattern1 <- paste0("\\\\link\\{", func, "\\}")
      replacement <- paste0("\\\\code{", func, "}")
      
      # Count occurrences before replacement
      occurrences <- sum(grepl(pattern1, content, fixed = TRUE))
      
      if (occurrences > 0) {
        content <- gsub(pattern1, replacement, content, fixed = TRUE)
        cat("  ✓ Replaced", occurrences, "occurrence(s) of \\link{", func, "}\n", sep = " ")
      }
      
      # Pattern 2: \link[package]{function}
      pattern2 <- paste0("\\\\link\\[[^]]*\\]\\{", func, "\\}")
      if (any(grepl(pattern2, content))) {
        content <- gsub(pattern2, replacement, content)
        cat("  ✓ Replaced \\link[...]{", func, "}\n", sep = "")
      }
    }
    
    # Write the file if changed
    if (!identical(content, original_content)) {
      # Create backup
      backup_file <- paste0(rd_file, ".bak")
      writeLines(original_content, backup_file)
      
      # Write fixed content
      writeLines(content, rd_file)
      cat("  ✅ File updated successfully\n")
      cat("  📁 Backup saved to:", backup_file, "\n")
      files_fixed <- files_fixed + 1
    } else {
      cat("  ℹ️ No changes needed\n")
    }
    cat("\n")
  }
  
  cat("================================================================================\n")
  cat("SUMMARY\n")
  cat("================================================================================\n\n")
  
  if (files_fixed > 0) {
    cat("✅ Fixed", files_fixed, ".Rd file(s)\n\n")
    cat("⚠️ IMPORTANT: These changes are to generated files!\n")
    cat("   The next time you run devtools::document(), these changes may be lost.\n")
    cat("   To make permanent fixes, you need to update the source .R files.\n\n")
    cat("Next steps:\n")
    cat("1. Run: devtools::check()  # Verify warnings are gone\n")
    cat("2. Find and fix the source .R files that generate these .Rd files\n")
  } else {
    cat("ℹ️ No .Rd files were modified\n")
  }
  
  return(files_fixed)
}

# Find the R source files that generate these .Rd files
find_source_files <- function() {
  cat("\n================================================================================\n")
  cat("FINDING SOURCE R FILES\n")
  cat("================================================================================\n\n")
  
  # Function names from .Rd files (without .Rd extension)
  target_functions <- c(
    "find_k_inter_for_target_hr",
    "forestsearch",
    "generate_aft_dgm_flex"
  )
  
  r_files <- list.files("R/", pattern = "\\.R$", full.names = TRUE)
  
  cat("Looking for source files that define these functions...\n\n")
  
  for (func_name in target_functions) {
    found <- FALSE
    
    for (r_file in r_files) {
      content <- readLines(r_file, warn = FALSE)
      
      # Look for function definition
      pattern1 <- paste0("^", func_name, "\\s*<-\\s*function")
      pattern2 <- paste0("^", func_name, "\\s*=\\s*function")
      pattern3 <- paste0("^#'\\s*@name\\s+", func_name)
      pattern4 <- paste0("^#'\\s*@rdname\\s+", func_name)
      
      if (any(grepl(pattern1, content)) || 
          any(grepl(pattern2, content)) ||
          any(grepl(pattern3, content)) ||
          any(grepl(pattern4, content))) {
        cat("✓ Found '", func_name, "' in: ", basename(r_file), "\n", sep = "")
        found <- TRUE
        
        # Check for @seealso with problematic links
        seealso_lines <- grep("^#'\\s*@seealso", content)
        if (length(seealso_lines) > 0) {
          cat("    Contains @seealso at line(s):", paste(seealso_lines, collapse = ", "), "\n")
          
          # Show the actual @seealso content
          for (line_num in seealso_lines) {
            cat("    Line", line_num, ":", trimws(content[line_num]), "\n")
          }
        }
        break
      }
    }
    
    if (!found) {
      cat("✗ Could not find source for '", func_name, "'\n", sep = "")
      
      # Look for similar names
      for (r_file in r_files) {
        if (grepl(func_name, basename(r_file), ignore.case = TRUE)) {
          cat("    Possible file:", basename(r_file), "\n")
        }
      }
    }
  }
}

# Fix source R files permanently
fix_source_r_files <- function() {
  cat("\n================================================================================\n")
  cat("PERMANENTLY FIXING SOURCE R FILES\n")
  cat("================================================================================\n\n")
  
  r_files <- list.files("R/", pattern = "\\.R$", full.names = TRUE)
  
  # Links to fix
  problem_links <- c(
    "find_k_inter_grid_search",
    "find_k_inter_batch",
    "bootstrap_forestsearch",
    "summarize_forestsearch",
    "find_quantile_for_proportion"
  )
  
  files_modified <- 0
  
  for (r_file in r_files) {
    content <- readLines(r_file, warn = FALSE)
    original_content <- content
    file_modified <- FALSE
    
    # Check each line for problematic links in comments
    for (i in seq_along(content)) {
      line <- content[i]
      
      # Only process roxygen2 comment lines
      if (grepl("^#'", line)) {
        for (link in problem_links) {
          # Replace \link{} with \code{}
          if (grepl(paste0("\\\\link\\{", link, "\\}"), line)) {
            content[i] <- gsub(
              paste0("\\\\link\\{", link, "\\}"),
              paste0("\\\\code{", link, "}"),
              content[i],
              fixed = TRUE
            )
            cat("✓ Fixed \\link{", link, "} in ", basename(r_file), " at line ", i, "\n", sep = "")
            file_modified <- TRUE
          }
        }
      }
    }
    
    # Save if modified
    if (file_modified) {
      writeLines(original_content, paste0(r_file, ".bak"))
      writeLines(content, r_file)
      files_modified <- files_modified + 1
      cat("  ✅ Saved:", basename(r_file), "\n\n")
    }
  }
  
  if (files_modified > 0) {
    cat("✅ Modified", files_modified, "source file(s)\n\n")
    cat("Next steps:\n")
    cat("1. Run: devtools::document()  # Regenerate .Rd files\n")
    cat("2. Run: devtools::check()     # Verify fixes\n")
  } else {
    cat("ℹ️ No source files needed modification\n")
    cat("The links might be generated programmatically.\n")
  }
  
  return(files_modified)
}

# ================================================================================
# MAIN EXECUTION
# ================================================================================

cat("\n🔧 ForestSearch Documentation Direct Fix\n")
cat("=========================================\n\n")

# Step 1: Fix .Rd files directly (immediate relief)
cat("STEP 1: Direct .Rd file fix (immediate relief from warnings)\n")
cat("--------------------------------------------------------------\n")
rd_fixed <- fix_rd_files_directly()

# Step 2: Find source files
cat("\n\nSTEP 2: Finding source R files\n")
cat("-------------------------------\n")
find_source_files()

# Step 3: Fix source files
cat("\n\nSTEP 3: Fixing source R files (permanent fix)\n")
cat("----------------------------------------------\n")
r_fixed <- fix_source_r_files()

# Summary
cat("\n\n================================================================================\n")
cat("FINAL SUMMARY\n")
cat("================================================================================\n\n")

if (rd_fixed > 0) {
  cat("✅ Immediate fix applied to", rd_fixed, ".Rd file(s)\n")
  cat("   R CMD check warnings should be gone now!\n\n")
}

if (r_fixed > 0) {
  cat("✅ Permanent fix applied to", r_fixed, "R source file(s)\n")
  cat("   Run devtools::document() to regenerate documentation\n\n")
} else {
  cat("⚠️ No R source files were modified\n")
  cat("   The .Rd fixes are temporary and may be overwritten by devtools::document()\n")
  cat("   You may need to manually find and fix the source files\n\n")
}

cat("Recommended actions:\n")
cat("1. Run: devtools::check()     # Verify warnings are gone\n")
if (r_fixed > 0) {
  cat("2. Run: devtools::document()  # Regenerate from fixed sources\n")
  cat("3. Run: devtools::check()     # Verify again\n")
}
