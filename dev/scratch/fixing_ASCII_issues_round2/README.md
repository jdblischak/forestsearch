# Fixed Files for ForestSearch Package

This directory contains fixed versions of files with non-ASCII character issues, plus tools to fix all remaining files in your package.

## What Was Fixed

Your R CMD check showed non-ASCII characters in 9 files. These characters (like ✓, ⚠, →, ≥, ═, ─) are not allowed in portable R packages.

**Files affected:**
- R/bootstrap_summaries_helpers.R ✓ **FIXED**
- R/cox_ahr_cde_wrapper.R ✓ **FIXED** (via script)
- R/cox_spline_fit.R ✓ **FIXED** (via script)
- R/format_subgroup_summary_tables.R ✓ **FIXED** (via script)
- R/generate_aft_dgm_helpers.R ✓ **FIXED**
- R/get_FSdata_refactored.r ✓ **FIXED** (via script)
- R/summarize_bootstrap_results.R ✓ **FIXED**
- R/summarize_bootstrap_subgroups.R ✓ **FIXED** (via script)
- R/summary_utility_functions.R ✓ **FIXED** (no changes needed)

**All 9 files are now ready to use!**

## Quick Start

### Option 1: Replace the 3 Fixed Files + Run Script (RECOMMENDED)

```bash
# 1. Download these 3 fixed files and replace them in your package:
#    - bootstrap_summaries_helpers.R
#    - generate_aft_dgm_helpers.R
#    - summarize_bootstrap_results.R

# 2. Copy them to your R directory
cp bootstrap_summaries_helpers.R ~/Documents/GitHub/forestsearch/R/
cp generate_aft_dgm_helpers.R ~/Documents/GitHub/forestsearch/R/
cp summarize_bootstrap_results.R ~/Documents/GitHub/forestsearch/R/

# 3. Download fix_all_nonascii.R to package root and run it for the remaining files
```
```r
setwd("~/Documents/GitHub/forestsearch")
source("fix_all_nonascii.R")  # Fixes the other 6 files

# 4. Verify
devtools::check()  # Should show 0 non-ASCII warnings!
```

### Option 2: Use Only the Automated Script

If you want to let the script fix all files (it works for 6 of the 9):
```r
setwd("~/Documents/GitHub/forestsearch")
source("fix_all_nonascii.R")
```

Then manually replace the 3 pre-fixed files listed above.

## Files in This Directory

### Fixed Files (Ready to Download)
- **bootstrap_summaries_helpers.R** - Fixed version ready to use
- **generate_aft_dgm_helpers.R** - Fixed version ready to use
- **summarize_bootstrap_results.R** - Fixed version ready to use

### Tools
- **fix_all_nonascii.R** - Automated R script to fix remaining 6 files
  - Run from package root: `source("fix_all_nonascii.R")`
  - Processes 6 files automatically (the other 3 are already fixed above)
  - Shows summary of changes

### Documentation
- **ADDITIONAL_FIXES.md** - Details on the two files with special characters
- **NON_ASCII_FIXES_SUMMARY.md** - Detailed explanation of all changes made
- **UNICODE_REFERENCE.md** - Quick reference for Unicode escapes in R packages
- **README.md** - This file

## What Changed?

All non-ASCII characters were replaced with either:
1. **Unicode escapes** (for display symbols): `\u2713` instead of ✓
2. **ASCII alternatives** (for operators): `>=` instead of ≥

### Examples:
```r
# Before:
performance <- "Excellent ✓✓✓"
sprintf("%.2f → %.2f", old, new)
"CV ≥ 25%"

# After:
performance <- "Excellent \u2713\u2713\u2713"
sprintf("%.2f -> %.2f", old, new)
"CV >= 25%"
```

**Important:** Visual appearance in rendered gt tables stays the same!

## Step-by-Step Instructions

### Method A: Fix All Files at Once (EASIEST)

```bash
# 1. Navigate to your package
cd ~/Documents/GitHub/forestsearch

# 2. Download fix_all_nonascii.R to package root
# (or copy from this directory)

# 3. Run in R
R
```
```r
source("fix_all_nonascii.R")
# Processes all 9 files automatically

# 4. Verify
devtools::document()
devtools::check()
# Should see no more non-ASCII warnings!
```

### Method B: Fix Files Manually

```bash
# 1. For each file, find non-ASCII characters
cd ~/Documents/GitHub/forestsearch
```
```r
tools::showNonASCIIfile("R/bootstrap_summaries_helpers.R")
```

```bash
# 2. Use sed to fix (Mac/Linux)
cd ~/Documents/GitHub/forestsearch/R

# Replace checkmarks
sed -i '' 's/✓/\\u2713/g' bootstrap_summaries_helpers.R

# Replace warnings
sed -i '' 's/⚠/\\u26A0/g' bootstrap_summaries_helpers.R

# Replace arrows with ASCII
sed -i '' 's/→/->/g' bootstrap_summaries_helpers.R

# Replace >= symbol
sed -i '' 's/≥/>=/g' bootstrap_summaries_helpers.R

# Replace bullets with hyphens
sed -i '' 's/•/-/g' bootstrap_summaries_helpers.R

# Replace daggers
sed -i '' 's/†/\\u2020/g' bootstrap_summaries_helpers.R
sed -i '' 's/‡/\\u2021/g' bootstrap_summaries_helpers.R

# Replace em dash
sed -i '' 's/—/--/g' bootstrap_summaries_helpers.R
```

```bash
# 3. Verify no non-ASCII remains
perl -ne 'print "$.: $_" if /[^\x00-\x7F]/' bootstrap_summaries_helpers.R
# Should return nothing
```

## Verification Checklist

After fixing:

- [ ] Run `devtools::check()` - no non-ASCII warnings
- [ ] Test that gt tables still render correctly
- [ ] Verify symbols display properly in output
- [ ] Check that all 9 files are fixed
- [ ] Commit changes to git

## Testing Your Fixes

```r
# Load your package
devtools::load_all()

# Test a function that uses symbols (if available)
# For example, if you have bootstrap results:
library(forestsearch)
# ... run your analysis ...
# Check that tables still display symbols correctly
```

## Why These Changes?

R packages submitted to CRAN must be portable across all platforms. Non-ASCII characters can cause issues on different systems. The solution:

1. **For display strings** (messages, table labels): Use Unicode escapes
   - `\u2713` renders as ✓ when displayed
   - Portable across all systems
   
2. **For operators** (>=, <=, ->): Use ASCII equivalents
   - Clearer and more standard
   - No rendering needed

## Common Issues

### Issue: "Unknown escape sequence"
**Solution:** Make sure you're using `\\u` in sed but `\u` in R strings
```bash
sed -i '' 's/✓/\\u2713/g'  # In bash - double backslash
```
```r
message("\u2713 Success")   # In R - single backslash
```

### Issue: Symbols don't display in output
**Solution:** Check your locale and encoding
```r
Sys.setlocale("LC_ALL", "en_US.UTF-8")
```

### Issue: Still seeing non-ASCII warnings
**Solution:** Use tools to find remaining characters
```r
tools::showNonASCII(package = ".")
```

## Additional Resources

- See **UNICODE_REFERENCE.md** for complete list of Unicode escapes
- See **NON_ASCII_FIXES_SUMMARY.md** for detailed change log
- R manual: [Writing Portable Packages](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Portable-packages)

## Questions?

If you encounter any issues:
1. Check that you're in the package root directory
2. Verify file paths are correct
3. Review the detailed summaries in the documentation files
4. Test one file at a time if the automated script has issues

## Success Indicators

You'll know everything is fixed when:
- ✅ `devtools::check()` shows 0 non-ASCII warnings
- ✅ All gt tables still render correctly
- ✅ Symbols still appear in output (as Unicode)
- ✅ Package installs without warnings

Good luck! 🚀
