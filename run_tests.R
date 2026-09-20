# Runs the attest test suite without installing the package.
# Usage, from the attest project:  source("run_tests.R")
if (!requireNamespace("testthat", quietly = TRUE)) stop('install.packages("testthat") first')
source("R/counts.R")
testthat::test_file("tests/testthat/test-counts.R")
