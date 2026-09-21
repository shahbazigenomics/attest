## Submission

First submission of attest 0.1.0.

## Test environments

* local: Ubuntu 24.04, R version 4.3.3 (2024-02-29) (R CMD check --as-cran: Status OK)
* local: macOS (arm64), R 4.x
* win-builder: R-devel  (to run before submission)

## R CMD check results

0 errors | 0 warnings | 1 note (expected on CRAN only)

* This is a new submission.

## Suggested packages

SummarizedExperiment, DESeq2 and edgeR are Bioconductor packages in Suggests.
They are used only to accept their objects as input; every use is guarded by
`requireNamespace()` or `methods::is()`, and the tests that construct those
objects are skipped when the packages are not installed.
