# Unicode Escape Reference Guide

When you need to use special symbols in your R package strings (for display purposes), use these Unicode escapes instead of the actual characters.

## Common Symbols

| Symbol | Unicode Escape | Name | Usage |
|--------|---------------|------|-------|
| † | `\u2020` | Dagger | Footnote marker (1st) |
| ‡ | `\u2021` | Double dagger | Footnote marker (2nd) |
| ✓ | `\u2713` | Checkmark | Success indicator |
| ✗ | `\u2717` | X mark | Failure indicator |
| ⚠ | `\u26A0` | Warning sign | Warning/caution |
| ℹ | `\u2139` | Information | Info symbol |
| ★ | `\u2605` | Star (filled) | Rating/importance |
| ☆ | `\u2606` | Star (empty) | Rating/importance |
| → | `->` (ASCII) | Right arrow | Use ASCII instead |
| ← | `<-` (ASCII) | Left arrow | Use ASCII instead |
| ≥ | `>=` (ASCII) | Greater than or equal | Use ASCII instead |
| ≤ | `<=` (ASCII) | Less than or equal | Use ASCII instead |

## Subscripts and Superscripts

| Symbol | Unicode Escape | Example |
|--------|---------------|---------|
| ₀ | `\u2080` | X₀ = `X\u2080` |
| ₁ | `\u2081` | X₁ = `X\u2081` |
| ₂ | `\u2082` | X₂ = `X\u2082` |
| ⁰ | `\u2070` | X⁰ = `X\u2070` |
| ¹ | `\u00B9` | X¹ = `X\u00B9` |
| ² | `\u00B2` | X² = `X\u00B2` |
| ³ | `\u00B3` | X³ = `X\u00B3` |

**Note:** For gt tables, use HTML instead: `<sub>1</sub>` or `<sup>2</sup>`

## Greek Letters

| Symbol | Unicode Escape | Name |
|--------|---------------|------|
| α | `\u03B1` | alpha |
| β | `\u03B2` | beta |
| γ | `\u03B3` | gamma |
| δ | `\u03B4` | delta |
| θ | `\u03B8` | theta |
| λ | `\u03BB` | lambda |
| μ | `\u03BC` | mu |
| π | `\u03C0` | pi |
| σ | `\u03C3` | sigma |
| τ | `\u03C4` | tau |
| φ | `\u03C6` | phi |
| χ | `\u03C7` | chi |

## Usage Examples

### In Regular Strings
```r
# Performance ratings
performance <- "Excellent \u2713\u2713\u2713"
warning_msg <- "\u26A0 High variability detected"

# Messages with symbols
message("Success \u2713")
cat("Processing", "\u2192", "Complete\n")  # Better: use "->" instead
```

### In gt Tables
```r
# Use Unicode escapes in gt::md()
gt::cols_label(
  hr = "HR",
  pvalue = gt::md("P<sup>\u2020</sup>")  # Footnote marker
)

# Or use HTML (preferred for subscripts/superscripts)
gt::cols_label(
  med_t = gt::md("Med<sub>T</sub>"),
  hr_squared = gt::md("HR<sup>2</sup>")
)
```

### In Plot Labels (ggplot2)
```r
# ggplot2 supports Unicode escapes
ggplot(data, aes(x, y)) +
  labs(
    title = "Performance Rating: \u2713\u2713\u2713",
    subtitle = "\u26A0 Use with caution"
  )

# Better: use expression() for math symbols
labs(
  x = expression(mu ~ "±" ~ sigma),
  y = expression(alpha[1])
)
```

## When to Use ASCII Instead

Some symbols have perfectly good ASCII alternatives. Use these instead of Unicode:

| Instead of | Use ASCII | Reason |
|-----------|-----------|--------|
| → | `->` | Clearer, more standard |
| ← | `<-` | R assignment operator |
| × | `*` | Standard multiplication |
| ÷ | `/` | Standard division |
| − | `-` | Standard minus/hyphen |
| ≥ | `>=` | Standard comparison |
| ≤ | `<=` | Standard comparison |
| • | `-` or `*` | Markdown bullets |
| — | `--` or `-` | Em/en dash |

## Tips for Package Development

### 1. Always Use Unicode Escapes in Strings
```r
# BAD - Direct Unicode character
message("✓ Success")

# GOOD - Unicode escape
message("\u2713 Success")
```

### 2. Comments Can Use Direct Unicode
```r
# This is OK in comments: ✓ ⚠ →
# But code must use: \u2713 \u26A0 ->
```

### 3. Check for Non-ASCII
```r
# From your package root
tools::showNonASCII(package = ".")

# Or check specific file
tools::showNonASCIIfile("R/your_file.R")
```

### 4. Prefer HTML in gt Tables
```r
# Unicode escapes work
gt::md("HR<sup>\u2020</sup>")

# But HTML is often clearer
gt::md("HR<sup>&dagger;</sup>")
gt::md("Med<sub>T</sub>")
```

## Complete List of Replacements Made

From your package:
- `†` → `\u2020`
- `‡` → `\u2021`
- `✓` → `\u2713`
- `⚠` → `\u26A0`
- `•` → `-`
- `→` → `->`
- `≥` → `>=`
- `—` → `--`

These maintain visual appearance while being CRAN-compliant!
