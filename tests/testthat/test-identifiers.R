at_fixture <- function() {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  readRDS(f)
}

test_that("real identifiers from three species pass, including the ones that look wrong", {
  fx <- at_fixture()
  for (nm in names(fx))
    expect_equal(attest_identifiers(fx[[nm]])$verdict, "PERMITTED", info = nm)

  # fission yeast systematic names end in .NN (SPAC212.09c, SPAC977.03). A rule
  # that read a trailing .N as an Ensembl version would fire on every one.
  expect_false(any(grepl("version", attest_identifiers(fx$fission)$evidence)))
  expect_true(grepl("symbol or systematic", attest_identifiers(fx$fission)$evidence[1]))
  expect_true(grepl("Ensembl",  attest_identifiers(fx$airway)$evidence[1]))
  expect_true(grepl("FlyBase",  attest_identifiers(fx$pasilla)$evidence[1]))
})

test_that("counting summary rows are found and priced", {
  m <- at_fixture()$airway
  h <- rbind(m,
             "__no_feature"           = round(colSums(m) * 0.12),
             "__ambiguous"            = round(colSums(m) * 0.03),
             "__too_low_aQual"        = 0,
             "__not_aligned"          = 0,
             "__alignment_not_unique" = round(colSums(m) * 0.05))
  r <- attest_identifiers(h)
  expect_equal(r$verdict, "NOT PERMITTED")
  expect_equal(length(r$measurements$summary_rows), 5L)
  expect_true(all(r$measurements$summary_row_share > 0.15))
  expect_true(grepl("library", r$evidence[1]))

  s <- rbind(m, "N_unmapped" = round(colSums(m) * 0.08),
             "N_multimapping" = round(colSums(m) * 0.10),
             "N_noFeature"    = round(colSums(m) * 0.05),
             "N_ambiguous"    = round(colSums(m) * 0.02))
  expect_equal(attest_identifiers(s)$verdict, "NOT PERMITTED")

  # and it sets the whole report, since those rows are in every total
  expect_equal(attest(h)$verdict, "NOT PERMITTED")
})

test_that("Excel damage is reported, and Entrez IDs are not mistaken for it", {
  m <- at_fixture()$airway
  e <- m; rownames(e)[c(3, 10, 77)] <- c("1-Sep", "1-Mar", "Sep-2")
  r <- attest_identifiers(e)
  expect_equal(r$verdict, "CAUTION")
  expect_equal(length(r$measurements$excel_mangled), 3L)
  expect_true(any(grepl("spreadsheet", r$evidence)))

  # Entrez gene IDs are all digits; an Excel-serial rule would flag every row
  ent <- m; rownames(ent) <- as.character(seq_len(nrow(m)) + 40000)
  r2 <- attest_identifiers(ent)
  expect_equal(r2$verdict, "PERMITTED")
  expect_true(grepl("Entrez", r2$evidence[1]))
})

test_that("two sources of identifiers are separated from one", {
  m <- at_fixture()$airway

  v <- m; rownames(v)[1:1500] <- paste0(rownames(v)[1:1500], ".3")
  expect_equal(attest_identifiers(v)$verdict, "CAUTION")

  v2 <- m; rownames(v2) <- paste0(rownames(v2), ".3")   # all versioned: consistent
  expect_equal(attest_identifiers(v2)$verdict, "PERMITTED")

  x <- m; rownames(x)[1:300] <- paste0("GENE", 1:300)
  expect_equal(attest_identifiers(x)$verdict, "CAUTION")

  # spike-ins are a legitimate second naming system and must not fire
  sp <- rbind(m, matrix(stats::rpois(8 * 20, 50), 20, 8,
                        dimnames = list(sprintf("ERCC-%05d", 1:20), colnames(m))))
  r <- attest_identifiers(sp)
  expect_equal(r$verdict, "PERMITTED")
  expect_equal(r$measurements$n_spike_in, 20L)
})

test_that("duplicated and empty names are caught", {
  m <- at_fixture()$airway
  d <- m; rownames(d)[5] <- rownames(d)[1]
  expect_equal(attest_identifiers(d)$verdict, "CAUTION")
  expect_equal(attest_identifiers(d)$measurements$n_duplicate, 1L)

  b <- m; rownames(b)[7] <- ""
  expect_equal(attest_identifiers(b)$verdict, "CAUTION")
  expect_equal(attest_identifiers(b)$measurements$n_blank, 1L)
})

test_that("featureCounts' own Unassigned_* summary rows are recognised", {
  # previously absent from at_summary_rows() entirely - none of these start
  # with "__" or "N_", so they read as ordinary genes and were never priced
  m <- at_fixture()$airway
  h <- rbind(m,
             "Unassigned_Ambiguity"    = round(colSums(m) * 0.06),
             "Unassigned_MultiMapping" = round(colSums(m) * 0.04))
  r <- attest_identifiers(h)
  expect_equal(r$verdict, "NOT PERMITTED")
  expect_equal(length(r$measurements$summary_rows), 2L)
})

test_that("a totals row with no recognised name is flagged as CAUTION, not silently kept as a gene", {
  m <- at_fixture()$airway
  # a name no pipeline in at_summary_rows() knows, at the bottom of the table,
  # holding far more than any real gene ever does
  h <- rbind(m, "grand_total_reads" = round(colSums(m) * 1.3))
  r <- attest_identifiers(h)
  expect_equal(r$verdict, "CAUTION")
  expect_equal(r$measurements$suspect_summary_rows, "grand_total_reads")
  expect_true(any(grepl("bottom of the table", r$evidence)))

  # a real gene near the edge that just happens to be highly expressed must
  # not be flagged - only rows far beyond any real gene's share trigger this
  ordinary <- attest_identifiers(m)
  expect_equal(length(ordinary$measurements$suspect_summary_rows), 0L)
  expect_equal(ordinary$verdict, "PERMITTED")
})

test_that("numbered rows are a missing input, not a fault", {
  m <- at_fixture()$airway
  nr <- m; rownames(nr) <- NULL
  expect_equal(attest_identifiers(nr)$verdict, "UNKNOWN")

  r <- attest(nr)
  expect_equal(r$verdict, "PERMITTED")
  expect_false("identifiers" %in% names(r$checks))
  expect_true(any(grepl("^identifiers", r$not_run)))

  expect_true("identifiers" %in% names(attest(m)$checks))
})
