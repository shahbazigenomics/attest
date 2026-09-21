# Real Bioconductor objects, not look-alikes. A hand-built stand-in class is
# what let attest(dds) ship unable to read a single real DESeqDataSet: the
# stand-in exercised the sample-sheet extraction and never the counts.

at_fixture_objects <- function() {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  m  <- readRDS(f)$airway
  md <- data.frame(cell = factor(rep(c("N61311", "N052611", "N080611", "N061011"), each = 2)),
                   dex  = factor(rep(c("untrt", "trt"), times = 4), levels = c("untrt", "trt")),
                   row.names = colnames(m))
  list(m = m, md = md)
}

test_that("a real SummarizedExperiment is read", {
  skip_if_not_installed("SummarizedExperiment")
  o <- at_fixture_objects()
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = o$m), colData = o$md)

  expect_equal(attest_counts(se)$verdict, "PERMITTED")
  expect_equal(dim(at_as_matrix(se)), dim(o$m))
  r <- attest(se)
  expect_equal(r$checks[["value scale"]]$verdict, "PERMITTED")
  expect_true("identifiers" %in% names(r$checks))
  expect_true(grepl("sample sheet", r$source_note))

  # the counts assay is found by name, not by position
  se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(cpm = t(t(o$m) / colSums(o$m)) * 1e6, counts = o$m), colData = o$md)
  expect_equal(attest_counts(se2)$verdict, "PERMITTED")
})

test_that("a real DESeqDataSet runs every check with nothing else supplied", {
  skip_if_not_installed("DESeq2")
  o <- at_fixture_objects()
  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(o$m, o$md, ~ cell + dex))

  r <- attest(dds)
  expect_equal(r$verdict, "PERMITTED")
  expect_true(all(c("value scale", "identifiers", "completeness", "identity", "design",
                    "detectability") %in% names(r$checks)))
  expect_equal(length(r$not_run), 0L)
  expect_true(grepl("DESeqDataSet", r$source_note))
  expect_equal(r$checks$design$measurements$of_interest, "dex")

  # and a fault inside it is still found
  cpm <- round(t(t(o$m) / colSums(o$m)) * 1e6)
  storage.mode(cpm) <- "integer"
  bad <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(cpm, o$md, ~ cell + dex))
  expect_equal(attest(bad)$verdict, "NOT PERMITTED")
})

test_that("a real DGEList gives up its counts and its grouping", {
  skip_if_not_installed("edgeR")
  o <- at_fixture_objects()
  y <- edgeR::DGEList(counts = o$m, group = o$md$dex)

  expect_equal(attest_counts(y)$verdict, "PERMITTED")
  got <- at_from_object(y)
  expect_equal(got$source, "DGEList")
  expect_equal(nlevels(droplevels(as.factor(got$group))), 2L)

  r <- attest(y)
  expect_true("detectability" %in% names(r$checks))
  expect_true(any(grepl("^design adequacy", r$not_run)))   # a DGEList has no formula
})

test_that("an object that cannot be read says so, and says why", {
  r <- attest(list(a = 1))
  expect_equal(r$checks[["value scale"]]$verdict, "UNKNOWN")
  expect_true(any(grepl("could not be read", r$not_run)))
  expect_false(any(grepl("no row names", r$not_run)))
})
