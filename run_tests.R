# Runs the attest test suite without installing the package.
# Usage, from the attest project:  source("run_tests.R")
if (!requireNamespace("testthat", quietly = TRUE)) stop('install.packages("testthat") first')
for (f in list.files("R", full.names = TRUE)) source(f)
for (f in list.files("tests/testthat", pattern = "^test-.*\\.R$", full.names = TRUE))
  testthat::test_file(f)
