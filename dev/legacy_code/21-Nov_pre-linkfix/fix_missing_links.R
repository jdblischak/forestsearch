# Fix Missing Documentation Links in ForestSearch Package
# Run this from your package directory

# ================================================================================
# STEP 1: Identify where the links are
# ================================================================================

cat("Searching for documentation with broken links...\n\n")

# Function to search for link references in R files
search_for_links <- function() {
  r_files <- list.files("R/", pattern = "\\.R$", full.names = TRUE)
  
  links_to_find <- c(
    "find_k_inter_grid_search",
    "find_k_inter_batch", 
    "bootstrap_forestsearch",
    "summarize_forestsearch",
    "find_quantile_for_proportion"
  )
  
  for (file in r_files) {
    content <- readLines(file)
    for (link in links_to_find) {
      if (any(grepl(paste0("\\\\link\\{", link, "\\}"), content))) {
        cat("Found \\link{", link, "} in: ", file, "\n")
        line_nums <- which(grepl(paste0("\\\\link\\{", link, "\\}"), content))
        cat("  Line(s):", paste(line_nums, collapse = ", "), "\n")
      }
    }
  }
}

# Run the search
search_for_links()

# ================================================================================
# STEP 2: Check if these functions actually exist
# ================================================================================

cat("\n\nChecking if functions exist in package...\n")

check_function_exists <- function(func_name) {
  # Check in R/ directory
  r_files <- list.files("R/", pattern = "\\.R$", full.names = TRUE)
  found <- FALSE
  
  for (file in r_files) {
    content <- readLines(file)
    # Look for function definition
    if (any(grepl(paste0("^", func_name, "\\s*<-\\s*function"), content))) {
      cat("✓ Found", func_name, "in", file, "\n")
      found <- TRUE
      
      # Check if it's exported
      if (any(grepl(paste0("@export.*", func_name), content)) || 
          any(grepl("@export", content))) {
        cat("  Status: Exported\n")
      } else {
        cat("  Status: Not exported (internal)\n")
      }
      break
    }
  }
  
  if (!found) {
    cat("✗", func_name, "NOT FOUND\n")
  }
  
  return(found)
}

# Check each function
missing_functions <- c(
  "find_k_inter_grid_search",
  "find_k_inter_batch",
  "bootstrap_forestsearch",
  "summarize_forestsearch",
  "find_quantile_for_proportion"
)

for (func in missing_functions) {
  check_function_exists(func)
}

# ================================================================================
# STEP 3: Automated fix suggestions
# ================================================================================

cat("\n\n================================================================================\n")
cat("RECOMMENDED FIXES:\n")
cat("================================================================================\n\n")

cat("Option 1: Quick fix - Replace \\link{} with \\code{}\n")
cat("------------------------------------------------\n")
cat("This prevents R CMD check warnings while still showing the function names.\n\n")

fix_links_to_code <- function() {
  r_files <- list.files("R/", pattern = "\\.R$", full.names = TRUE)
  
  links_to_fix <- c(
    "find_k_inter_grid_search",
    "find_k_inter_batch",
    "bootstrap_forestsearch",
    "summarize_forestsearch", 
    "find_quantile_for_proportion"
  )
  
  for (file in r_files) {
    content <- readLines(file)
    original_content <- content
    modified <- FALSE
    
    for (link in links_to_fix) {
      pattern <- paste0("\\\\link\\{", link, "\\}")
      replacement <- paste0("\\\\code{", link, "}")
      
      if (any(grepl(pattern, content))) {
        content <- gsub(pattern, replacement, content, fixed = TRUE)
        modified <- TRUE
        cat("Fixed \\link{", link, "} -> \\code{", link, "} in ", basename(file), "\n")
      }
    }
    
    if (modified) {
      # Create backup
      writeLines(original_content, paste0(file, ".bak"))
      # Write fixed version
      writeLines(content, file)
    }
  }
}

cat("\nTo apply this fix, run:\n")
cat("  fix_links_to_code()\n\n")

# ================================================================================
# STEP 4: Check for possible typos
# ================================================================================

cat("\nOption 2: Check for possible typos\n")
cat("----------------------------------\n")

# Look for similar function names
find_similar_functions <- function(target) {
  r_files <- list.files("R/", pattern = "\\.R$", full.names = TRUE)
  all_functions <- character()
  
  for (file in r_files) {
    content <- readLines(file)
    # Extract function names
    func_lines <- grep("^[a-zA-Z_][a-zA-Z0-9_.]*\\s*<-\\s*function", content, value = TRUE)
    if (length(func_lines) > 0) {
      func_names <- sub("\\s*<-.*", "", func_lines)
      func_names <- trimws(func_names)
      all_functions <- c(all_functions, func_names)
    }
  }
  
  # Find similar names
  all_functions <- unique(all_functions)
  
  # Calculate similarity
  similarities <- sapply(all_functions, function(f) {
    # Simple similarity based on common substrings
    common_chars <- sum(strsplit(target, "")[[1]] %in% strsplit(f, "")[[1]])
    return(common_chars / max(nchar(target), nchar(f)))
  })
  
  # Get top matches
  top_matches <- head(all_functions[order(similarities, decreasing = TRUE)], 3)
  
  if (length(top_matches) > 0) {
    cat("\nPossible matches for '", target, "':\n")
    for (match in top_matches) {
      cat("  -", match, "\n")
    }
  }
}

cat("\nSearching for similar function names...\n")

# These might be typos:
possibly_wrong <- list(
  "bootstrap_forestsearch" = "forestsearch_bootstrap",
  "summarize_forestsearch" = "summarize_bootstrap_subgroups"
)

for (wrong in names(possibly_wrong)) {
  cat("\n'", wrong, "' might be '", possibly_wrong[[wrong]], "'\n")
  find_similar_functions(wrong)
}

# ================================================================================
# STEP 5: Manual fix instructions
# ================================================================================

cat("\n\n================================================================================\n")
cat("MANUAL FIX INSTRUCTIONS:\n")
cat("================================================================================\n\n")

cat("1. For each missing link, decide:\n")
cat("   a) If the function doesn't exist -> Change \\link{} to \\code{}\n")
cat("   b) If the function exists but isn't exported -> Add #' @export\n")
cat("   c) If it's a typo -> Fix the function name\n\n")

cat("2. Files to check:\n")
cat("   - R/find_k_inter_for_target_hr.R\n")
cat("   - R/forestsearch.R\n") 
cat("   - R/generate_aft_dgm_flex.R\n\n")

cat("3. After fixing, run:\n")
cat("   devtools::document()\n")
cat("   devtools::check()\n\n")

# ================================================================================
# FINAL RECOMMENDATION
# ================================================================================

cat("================================================================================\n")
cat("RECOMMENDED ACTION:\n")
cat("================================================================================\n\n")

cat("Run this to apply the quick fix (changes \\link to \\code):\n\n")

cat("# Define the function\n")
cat(deparse(fix_links_to_code), sep = "\n")
cat("\n\n# Run it\n")
cat("fix_links_to_code()\n")
cat("devtools::document()\n")
cat("devtools::check()\n")
