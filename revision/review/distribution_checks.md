# Repository package checks — September 10, 2026

Checks run from the new repository branch:

- All 53 R source files parsed successfully.
- `make_tables.R` regenerated the tables from the included summary CSVs.
- `R CMD INSTALL` compiled and installed `lrcbart` into a temporary library.
- The single-leaf and two-leaf tests passed all 102 assertions, with no test failures, skips, or test warnings. R reported that the installed testthat package was built under R 4.4.3.
- A clean `latexmk` build produced an 82-page manuscript. The final LaTeX log contained no warnings, undefined references, or overfull boxes.

The complete simulations were not rerun for packaging. The committed PDF is the previously visually inspected manuscript; the repository rebuild checks the distributed source and bibliography.
