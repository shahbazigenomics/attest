# Runs the whole attest test suite without installing the package, and prints
# one table: a row per test file and a single total at the bottom.
# Usage, from the attest project:  source("run_tests.R")
if (!requireNamespace("testthat", quietly = TRUE)) stop('install.packages("testthat") first')
for (f in list.files("R", full.names = TRUE)) source(f)
testthat::test_dir("tests/testthat", load_package = "none", stop_on_failure = FALSE)
