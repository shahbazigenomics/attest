# Detectability answers "what could DESeq2, run on this design, have detected?"
# Calibrated against DESeq2 itself in validation/detectability_calibration_deseq2.R.

at_fx_det <- function() {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  m <- readRDS(f)$airway
  sh <- data.frame(cell = factor(rep(c("N61311", "N052611", "N080611", "N061011"), each = 2)),
                   dex  = factor(rep(c("untrt", "trt"), times = 4), levels = c("untrt", "trt")))
  list(m = m, sh = sh)
}

test_that("with DESeq2 installed, its own noise estimate is used and said so", {
  skip_if_not_installed("DESeq2")
  o <- at_fx_det()
  d <- attest_detectability(o$m, group = o$sh$dex, metadata = o$sh, design = ~ cell + dex)
  expect_true(d$measurements$deseq2_dispersion)
  expect_match(d$measurements$dispersion_source, "^DESeq2's own dispersion estimates for ~cell \\+ dex")
  expect_true(any(grepl("calibrated against DESeq2 itself", d$not_assessed)))
})

test_that("the moment estimate says what it costs", {
  o <- at_fx_det()
  d <- attest_detectability(o$m, group = o$sh$dex, metadata = o$sh, design = ~ cell + dex,
                            dispersion = "moments")
  expect_false(d$measurements$deseq2_dispersion)
  expect_match(d$measurements$dispersion_source, "^attest's moment estimate after fitting")
  expect_true(any(grepl("65-78%", d$not_assessed)))
})

test_that("a paired design's blocks are not counted as noise", {
  o <- at_fx_det()
  for (disp in c("moments", if (requireNamespace("DESeq2", quietly = TRUE)) "deseq2")) {
    unpaired <- attest_detectability(o$m, group = o$sh$dex, dispersion = disp)
    paired   <- attest_detectability(o$m, group = o$sh$dex, metadata = o$sh,
                                     design = ~ cell + dex, dispersion = disp)
    expect_lt(paired$measurements$median_min_detectable_fc,
              unpaired$measurements$median_min_detectable_fc, label = disp)
  }
})

test_that("the detectable change is the one that holds in both directions", {
  z <- stats::qnorm(1 - 0.005 / 2) + stats::qnorm(0.8)
  mu <- c(10, 100, 1000, 10000); phi <- rep(0.05, 4)
  both <- at_mdfc(mu, phi, 4, z)
  symmetric <- exp(z * sqrt((2 / 4) * (1 / mu + phi)))
  expect_true(all(both >= symmetric))
  # the gap closes as counts grow: the decrease only matters where counting noise does
  expect_lt((both / symmetric)[4], (both / symmetric)[1])
  # a gene averaging a fraction of a read cannot show a decrease of any size
  expect_true(is.infinite(at_mdfc(0.3, 0.05, 4, z)))
})

test_that("genes that cannot show a change are described in words, not as Inf", {
  o <- at_fx_det()
  d <- attest_detectability(o$m, group = o$sh$dex, dispersion = "moments")
  expect_false(any(grepl("Inf", d$evidence)))
  expect_true(is.finite(d$measurements$median_min_detectable_fc))
})
