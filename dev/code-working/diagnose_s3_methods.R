# =============================================================================
# diagnose_s3_methods.R
# Diagnostic script for "S3 methods declared in NAMESPACE but not found"
# =============================================================================
#
# Run this script from the package root directory.
# It checks two things:
#   1. Whether the package is properly installed (package.rds exists)
#   2. Whether every S3method() in NAMESPACE has a matching function in R/
#
# Usage:
#   source("diagnose_s3_methods.R")
#
# =============================================================================

cat("=== ForestSearch S3 Method Diagnostics ===\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Check installation state
# ─────────────────────────────────────────────────────────────────────────────

cat("── 1. Installation State ──\n\n")

pkg_path <- find.package("forestsearch", quiet = TRUE)

if (length(pkg_path) == 0) {
  cat("  [!] Package 'forestsearch' is NOT installed.\n")
  cat("      This is the primary cause of both warnings.\n")
  cat("      Fix: devtools::install() from the package root.\n\n")
} else {
  rds_path <- file.path(pkg_path, "Meta", "package.rds")
  if (file.exists(rds_path)) {
    cat("  [OK] package.rds exists at:", rds_path, "\n\n")
  } else {
    cat("  [!] package.rds is MISSING at:", rds_path, "\n")
    cat("      The package was loaded but not fully installed.\n")
    cat("      This happens with devtools::load_all() without a prior install.\n")
    cat("      Fix: devtools::install() then restart R.\n\n")
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 2. Parse NAMESPACE for S3method() declarations
# ─────────────────────────────────────────────────────────────────────────────

cat("── 2. NAMESPACE S3method Declarations ──\n\n")

ns_file <- "NAMESPACE"
if (!file.exists(ns_file)) {
  stop("NAMESPACE file not found. Run from the package root directory.")
}

ns_lines <- readLines(ns_file)
s3_lines <- grep("^S3method\\(", ns_lines, value = TRUE)

# Extract generic,class pairs
s3_methods <- lapply(s3_lines, function(line) {
  # S3method(generic,class)
  inner <- sub("^S3method\\((.*)\\)$", "\\1", line)
  parts <- strsplit(inner, ",")[[1]]
  list(
    generic = trimws(parts[1]),
    class   = trimws(parts[2]),
    fn_name = paste0(trimws(parts[1]), ".", trimws(parts[2]))
  )
})

cat("  Found", length(s3_methods), "S3method declarations:\n")
for (m in s3_methods) {
  cat("    S3method(", m$generic, ",", m$class, ")  =>  ", m$fn_name, "()\n")
}
cat("\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Search R/ files for each S3 method definition
# ─────────────────────────────────────────────────────────────────────────────

cat("── 3. Matching Definitions in R/ ──\n\n")

r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE,
                       ignore.case = TRUE)

if (length(r_files) == 0) {
  stop("No .R files found in R/ directory.")
}

# Read all R files once
all_code <- lapply(r_files, function(f) {
  list(file = f, lines = readLines(f, warn = FALSE))
})

missing <- character(0)
found   <- character(0)

for (m in s3_methods) {
  fn <- m$fn_name
  # Pattern: function name followed by <- function(
  # Handles both 'fn_name <- function(' and 'fn_name = function('
  pattern <- paste0("^\\s*", gsub("\\.", "\\\\.", fn), "\\s*(<-|=)\\s*function\\s*\\(")

  located <- FALSE
  for (code in all_code) {
    hits <- grep(pattern, code$lines)
    if (length(hits) > 0) {
      cat("  [OK] ", fn, "\n")
      cat("        File: ", code$file, " (line ", hits[1], ")\n")
      located <- TRUE

      # Check for @export or @exportS3Method tag nearby
      if (hits[1] > 1) {
        preceding <- code$lines[max(1, hits[1] - 20):hits[1]]
        has_export <- any(grepl("@export", preceding))
        if (!has_export) {
          cat("        [WARN] No @export tag found in preceding roxygen block\n")
        }
      }

      # Check for duplicates in OTHER files
      for (code2 in all_code) {
        if (code2$file != code$file) {
          hits2 <- grep(pattern, code2$lines)
          if (length(hits2) > 0) {
            cat("        [WARN] DUPLICATE definition in: ", code2$file,
                " (line ", hits2[1], ")\n")
          }
        }
      }
      break
    }
  }

  if (!located) {
    cat("  [MISSING] ", fn, "  -- NOT FOUND in any R/ file\n")
    missing <- c(missing, fn)
  } else {
    found <- c(found, fn)
  }
}

cat("\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Check for files that fail to parse
# ─────────────────────────────────────────────────────────────────────────────

cat("── 4. R File Parse Check ──\n\n")

parse_errors <- character(0)
for (f in r_files) {
  tryCatch({
    parse(file = f)
  }, error = function(e) {
    cat("  [ERROR] ", f, ": ", conditionMessage(e), "\n")
    parse_errors <<- c(parse_errors, f)
  })
}

if (length(parse_errors) == 0) {
  cat("  [OK] All", length(r_files), "R files parse without errors.\n")
} else {
  cat("\n  [!]", length(parse_errors), "file(s) have parse errors.\n")
  cat("      These files are silently skipped during load_all(),\n")
  cat("      which would cause all methods inside them to be 'not found'.\n")
}

cat("\n")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Summary and recommended actions
# ─────────────────────────────────────────────────────────────────────────────

cat("── 5. Summary ──\n\n")
cat("  S3 methods declared:     ", length(s3_methods), "\n")
cat("  Definitions found:       ", length(found), "\n")
cat("  Definitions missing:     ", length(missing), "\n")
cat("  Files with parse errors: ", length(parse_errors), "\n")
cat("\n")

if (length(missing) == 0 && length(parse_errors) == 0) {
  cat("  All method definitions are present and all files parse correctly.\n")
  cat("  The warnings are caused by an incomplete install.\n\n")
  cat("  ── Recommended Fix ──\n\n")
  cat("  Run the following in order:\n\n")
  cat('    devtools::document()    # Regenerate NAMESPACE and Rd files\n')
  cat('    devtools::install()     # Full install (creates Meta/package.rds)\n')
  cat("    # Restart R session\n")
  cat('    library(forestsearch)   # Should load without S3 warnings\n')
  cat("\n")
  cat("  Alternatively, if you only need a dev session:\n\n")
  cat('    devtools::document()\n')
  cat('    devtools::load_all()\n')
  cat("    # Note: 'no package.rds' warning may persist with load_all()\n")
  cat("    # but the S3 method warning should disappear.\n")
} else {
  if (length(missing) > 0) {
    cat("  [ACTION] Create definitions for missing methods:\n")
    for (fn in missing) {
      cat("    ", fn, "()\n")
    }
    cat("\n")
  }
  if (length(parse_errors) > 0) {
    cat("  [ACTION] Fix parse errors in:\n")
    for (f in parse_errors) {
      cat("    ", f, "\n")
    }
    cat("  Parse errors cause entire files to be skipped silently.\n")
    cat("  Every S3 method in those files will be 'not found'.\n")
    cat("\n")
  }
  cat("  After fixing, run:\n")
  cat('    devtools::document()\n')
  cat('    devtools::install()\n')
  cat("    # Restart R\n")
  cat('    library(forestsearch)\n')
}

cat("\n=== Diagnostics Complete ===\n")
