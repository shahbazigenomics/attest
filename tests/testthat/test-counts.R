test_that("verdicts are correct on real count matrices and their transforms", {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"            # package not installed yet
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds" # run from tests/testthat/
  fx <- readRDS(f)
  set.seed(1)
  for (nm in names(fx)) {
    m <- fx[[nm]]
    len <- round(runif(nrow(m), 500, 8000))
    expect_equal(attest_counts(m)$verdict, "PERMITTED", info = nm)
    expect_equal(attest_counts(round(m / 40))$verdict, "PERMITTED", info = nm)
    expect_equal(attest_counts(m[rowSums(m) >= 10, ])$verdict, "PERMITTED", info = nm)

    cpm <- t(t(m) / colSums(m)) * 1e6
    expect_equal(attest_counts(cpm)$verdict, "NOT PERMITTED", info = nm)
    expect_equal(attest_counts(round(cpm))$verdict, "NOT PERMITTED", info = nm)
    expect_equal(attest_counts(log2(cpm + 1))$verdict, "NOT PERMITTED", info = nm)

    fpkm <- (m / len) / rep(colSums(m) / 1e9, each = nrow(m))
    expect_equal(attest_counts(fpkm)$verdict, "NOT PERMITTED", info = nm)

    r <- m / len
    expect_equal(attest_counts(round(t(t(r) / colSums(r)) * 1e6))$verdict, "NOT PERMITTED", info = nm)

    expect_equal(attest_counts(round(t(t(m) / at_size_factors(m))))$verdict, "CAUTION", info = nm)

    est <- m; est[m > 0] <- m[m > 0] * runif(sum(m > 0), 0.9, 1.1)
    expect_equal(attest_counts(est)$verdict, "CAUTION", info = nm)
  }
})

test_that("bad input returns UNKNOWN instead of stopping", {
  expect_equal(attest_counts("not a matrix")$verdict, "UNKNOWN")
  expect_equal(attest_counts(matrix(1:4, 2))$verdict, "UNKNOWN")          # too small
  expect_equal(attest_counts(matrix(c(NA, 1:99), 10))$verdict, "UNKNOWN") # NA present
  expect_equal(attest_counts(NULL)$verdict, "UNKNOWN")
  expect_equal(attest_counts(data.frame(a = letters[1:5]))$verdict, "UNKNOWN")
})

test_that("small panels are judged, not refused, and unrun checks are reported", {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  m <- readRDS(f)$airway

  panel <- m[1:30, ]                       # 30 genes: below the old floor of 50
  res <- attest_counts(panel)
  expect_equal(res$verdict, "PERMITTED")
  expect_true(length(res$not_assessed) >= 1)          # says which checks could not run
  expect_true(any(grepl("size factors", res$not_assessed)))

  expect_equal(attest_counts(m[1:5, ])$verdict, "UNKNOWN")   # 5 genes: still too small

  full <- attest_counts(m)
  expect_equal(length(full$not_assessed), 0)          # nothing unrun on a normal matrix
  expect_true(nzchar(full$headline))
})

test_that("a scale factor is claimed only when samples agree", {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  m <- readRDS(f)$airway
  small <- round(m / 60)
  cpm <- round(t(t(small) / colSums(small)) * 1e6)
  lat <- at_lattice(cpm)
  expect_true(!is.null(lat$factor_estimate))
  expect_true(lat$n_agree >= 2)
  expect_null(at_lattice(m)$factor_estimate)          # raw counts: no lattice
})
