# "Already normalised" needs positive evidence that each sample was divided by
# a factor and rounded. Flat size factors alone are not that evidence: raw
# libraries of near-equal depth whose totals differ through a few genes look
# the same (whole blood without globin depletion is the common case).

at_fixture_airway_r <- function() {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  readRDS(f)$airway
}
at_thin <- function(m, prob) {
  p <- rep(prob, each = nrow(m))
  matrix(stats::rbinom(length(m), as.vector(m), p), nrow(m), dimnames = dimnames(m))
}

test_that("raw libraries of equal depth with variable globin stay raw", {
  set.seed(11)
  m  <- at_fixture_airway_r()
  mb <- at_thin(m, 0.9 * min(colSums(m)) / colSums(m))          # depth-balanced
  # alternating 20% / 80% globin: totals differ ~4x at equal depth, which is
  # guaranteed to trip the size-factor rule this test exists to outlive
  share <- rep(c(0.2, 0.8), length.out = ncol(mb))
  hb <- round(colSums(mb) * share / (1 - share))
  x  <- rbind(mb, HBB = round(hb * 0.6), HBA1 = round(hb * 0.4))

  r <- attest_counts(x)
  expect_equal(r$verdict, "PERMITTED")
  expect_lt(r$measurements$sf_vs_colsum, 0.3)          # the old rule would have fired
  expect_lt(r$measurements$aliasing_max, 1.5)
  expect_true(any(grepl("globin", r$evidence)))
})

test_that("normalised and rounded counts are caught at depths the size factors miss", {
  set.seed(12)
  m <- at_fixture_airway_r()
  for (keep in c(1, 0.3, 0.1, 0.03)) {                  # down to ~30k reads per sample
    x  <- at_thin(m, keep)
    nr <- round(t(t(x) / at_size_factors(x)))
    r  <- attest_counts(nr)
    expect_equal(r$verdict, "CAUTION", info = paste("keep", keep))
    expect_gt(r$measurements$aliasing_max, 1.5)
    expect_equal(attest_counts(x)$verdict, "PERMITTED", info = paste("raw, keep", keep))
  }
})

test_that("the rescaling measure only runs where it means something", {
  m <- at_fixture_airway_r()
  expect_true(all(is.na(at_aliasing(m + 0.25))))       # non-integer: not applicable
  a <- at_aliasing(m)
  expect_equal(length(a), ncol(m))
  expect_true(all(a < 1.5, na.rm = TRUE))
})
