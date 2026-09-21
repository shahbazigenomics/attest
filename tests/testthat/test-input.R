at_fixture <- function() {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  readRDS(f)
}

test_that("a plain count file is read and audited", {
  m <- at_fixture()$airway
  p <- tempfile(fileext = ".tsv")
  out <- data.frame(Geneid = rownames(m), m, check.names = FALSE)
  utils::write.table(out, p, sep = "\t", quote = FALSE, row.names = FALSE)

  r <- attest_file(p)
  expect_s3_class(r, "attest_report")
  expect_equal(r$checks$input$verdict, "PERMITTED")
  expect_equal(r$checks$input$measurements$n_samples, ncol(m))
  expect_equal(r$checks$input$measurements$n_genes, nrow(m))
  expect_equal(r$checks$input$measurements$id_column, "Geneid")
  expect_equal(r$checks[["value scale"]]$verdict, "PERMITTED")
  # the identifiers survived, so the human sex markers are still findable
  expect_false(any(grepl("^sex check", r$checks$identity$not_assessed)))
})

test_that("featureCounts annotation columns are removed, not counted as samples", {
  m <- at_fixture()$airway[1:500, ]
  p <- tempfile(fileext = ".txt")
  out <- data.frame(Geneid = rownames(m),
                    Chr = "1", Start = 1:500, End = 501:1000, Strand = "+",
                    Length = round(runif(500, 500, 8000)),
                    m, check.names = FALSE)
  writeLines(c("# Program:featureCounts v2.0.6; Command:\"featureCounts\"",
               paste(names(out), collapse = "\t")), p)
  utils::write.table(out, p, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = FALSE, append = TRUE)

  r <- attest_file(p)
  expect_equal(r$checks$input$measurements$n_comment_lines, 1L)
  expect_equal(r$checks$input$measurements$n_samples, ncol(m))
  expect_true(all(c("Chr", "Start", "End", "Strand", "Length") %in%
                    r$checks$input$measurements$annotation_columns))
  expect_true(any(grepl("featureCounts output", r$checks$input$evidence)))
  # and the matrix itself is still recognised as raw counts
  expect_equal(r$checks[["value scale"]]$verdict, "PERMITTED")

  # why this has to happen at the file: with Start, End and Length left in as
  # three extra libraries, every check on the matrix still says PERMITTED
  naive <- as.matrix(out[, vapply(out, is.numeric, logical(1))])
  expect_equal(attest_counts(naive)$verdict, "PERMITTED")
})

test_that("comma separation, gzip and R's short-header row names all work", {
  m <- at_fixture()$pasilla[1:800, ]

  p1 <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(gene_id = rownames(m), m, check.names = FALSE),
                   p1, quote = FALSE, row.names = FALSE)
  expect_equal(attest_file(p1)$checks$input$measurements$sep, ",")
  expect_equal(attest_file(p1)$checks$input$measurements$n_samples, ncol(m))

  p2 <- tempfile(fileext = ".csv.gz")
  con <- gzfile(p2, "wt"); utils::write.csv(m, con, quote = FALSE); close(con)
  r2 <- attest_file(p2)                     # the unnamed first column write.csv() leaves
  expect_equal(r2$checks$input$measurements$n_samples, ncol(m))
  expect_equal(r2$checks$input$measurements$id_column, "(first column)")
  expect_equal(r2$checks[["value scale"]]$verdict, "PERMITTED")
})

test_that("duplicated gene identifiers are reported, not silently merged", {
  m <- at_fixture()$fission[1:400, ]
  ids <- rownames(m); ids[c(5, 9)] <- ids[1]
  p <- tempfile(fileext = ".tsv")
  utils::write.table(data.frame(gene = ids, m, check.names = FALSE), p,
                     sep = "\t", quote = FALSE, row.names = FALSE)

  r <- attest_file(p)
  expect_equal(r$checks$input$verdict, "CAUTION")
  expect_equal(r$checks$input$measurements$n_duplicate_ids, 2L)
  expect_equal(r$verdict, "CAUTION")        # an input fault does set the report
})

test_that("unreadable input fails typed, never raising", {
  expect_equal(attest_file(tempfile())$verdict, "UNKNOWN")

  p <- tempfile(); writeLines(c("just", "some", "prose"), p)
  expect_equal(attest_file(p)$verdict, "UNKNOWN")

  p2 <- tempfile(); writeLines(c("gene\tlen", "a\t100", "b\t200"), p2)
  r <- attest_file(p2)                      # one numeric column after annotation
  expect_equal(r$verdict, "UNKNOWN")
  expect_true(length(r$not_run) > 0)
})

test_that("the design and the sample sheet are taken from the object", {
  m  <- at_fixture()$airway
  md <- data.frame(cell = factor(rep(c("A","B","C","D"), each = 2)),
                   dex  = factor(rep(c("untrt","trt"), times = 4)),
                   row.names = colnames(m))

  # a DESeqDataSet look-alike: colData plus a design slot
  setClass("atFakeSE", representation(assays = "matrix", cd = "data.frame",
                                      design = "formula"))
  dds <- new("atFakeSE", assays = m, cd = md, design = ~ cell + dex)
  got <- at_from_object(dds)
  expect_true(inherits(got$design, "formula"))
  expect_equal(got$source, "DESeqDataSet")

  # a DGEList look-alike carries its grouping in $samples
  y <- list(counts = m, samples = data.frame(group = md$dex, lib.size = colSums(m)))
  class(y) <- "DGEList"
  got2 <- at_from_object(y)
  expect_equal(got2$source, "DGEList")
  expect_equal(nlevels(droplevels(as.factor(got2$group))), 2L)

  r <- attest(y)
  expect_true("detectability" %in% names(r$checks))
  expect_equal(length(r$not_run), 1L)       # design still not run: no formula
  expect_true(grepl("^design adequacy", r$not_run))
  expect_true(grepl("grouping taken from the DGEList", r$source_note))
})

test_that("a check with no inputs is named, not dropped from the report", {
  m <- at_fixture()$airway
  r <- attest(m)
  expect_equal(r$verdict, "PERMITTED")
  expect_equal(length(r$not_run), 2L)
  expect_true(any(grepl("^design adequacy", r$not_run)))
  expect_true(any(grepl("^detectability", r$not_run)))
  expect_equal(length(attest_as_list(r)$not_run), 2L)

  txt <- paste(utils::capture.output(print(r)), collapse = "\n")
  expect_true(grepl("not run", txt))

  md <- data.frame(cell = factor(rep(c("A","B","C","D"), each = 2)),
                   dex  = factor(rep(c("untrt","trt"), times = 4)))
  expect_equal(length(attest(m, metadata = md, design = ~ cell + dex)$not_run), 0L)
})

test_that("a spreadsheet's totals row and trailing blank rows are found (GEO GSE185245 shape)", {
  m <- at_fixture()$airway[1:500, ]
  f <- tempfile(fileext = ".tsv")
  body <- cbind(gene = rownames(m), as.data.frame(m))
  blank <- as.data.frame(matrix(NA, 2, ncol(body), dimnames = list(NULL, names(body))))
  blank$gene <- ""
  total <- body[1, ]; total$gene <- ""; total[-1] <- colSums(m)
  utils::write.table(rbind(body, blank, total), f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  r <- attest_file(f)
  expect_equal(r$checks$input$verdict, "NOT PERMITTED")
  expect_equal(r$checks$input$measurements$n_empty_rows, 2L)
  expect_equal(length(r$checks$input$measurements$no_id_rows), 1L)
  expect_true(grepl("no gene identifier", r$checks$input$headline))
  expect_equal(r$verdict, "NOT PERMITTED")

  # a row with a gap in it: the value checks decline instead of failing
  body2 <- body; body2[3, 2] <- NA
  utils::write.table(body2, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  r2 <- attest_file(f)
  expect_equal(r2$checks$input$verdict, "CAUTION")
  expect_equal(r2$checks$identity$verdict, "UNKNOWN")
  expect_equal(r2$checks[["value scale"]]$verdict, "UNKNOWN")

  # a genomic Position column is annotation, not a library
  body3 <- cbind(body[, 1, drop = FALSE], Position = seq_len(nrow(body)), body[, -1])
  utils::write.table(body3, f, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_equal(ncol(at_read_counts(f)$counts), ncol(m))
})

test_that("tab wins over a delimiter-like character inside an annotation column (GEO GSE115255 shape)", {
  m <- at_fixture()$airway[1:300, ]
  # featureCounts' Chr/Start/End join multiple exons with ';' - a DIFFERENT
  # number per gene (real genes have different exon counts), so semicolons
  # are inconsistent line to line even though tabs stay fixed at 6 + ncol(m)
  set.seed(1)
  nexon <- sample(1:40, nrow(m), replace = TRUE)
  joined <- function(x) vapply(nexon, function(k) paste(rep(x, k), collapse = ";"), "")
  body <- data.frame(Geneid = rownames(m), Chr = joined("chrX"), Start = joined("100"),
                     End = joined("200"), Strand = joined("+"), Length = 1000, m, check.names = FALSE)
  f <- tempfile(fileext = ".tsv")
  utils::write.table(body, f, sep = "\t", quote = FALSE, row.names = FALSE)
  rd <- at_read_counts(f)
  expect_equal(rd$check$measurements$sep, "\t")
  expect_equal(ncol(rd$counts), ncol(m))
  expect_equal(nrow(rd$counts), nrow(m))
})

test_that("a title and merged-header rows above the real header are skipped (GEO GSE162669 shape)", {
  m <- at_fixture()$airway[1:400, ]
  f <- tempfile(fileext = ".tsv")
  ncol_out <- ncol(m) + 2
  blank_row <- function(first) c(first, rep("", ncol_out - 1))
  lines <- c(
    paste(blank_row("Raw counts of sequencing reads (TPM)"), collapse = "\t"),
    paste(c("", "", "Sample Name", rep("", ncol_out - 3)), collapse = "\t"),
    paste(c("", "", colnames(m)), collapse = "\t"),
    paste(blank_row("Sample No."), collapse = "\t"),
    paste(c("ID", "symbol", colnames(m)), collapse = "\t"))
  body <- apply(cbind(rownames(m), "SYM", m), 1, paste, collapse = "\t")
  writeLines(c(lines, body), f)
  rd <- at_read_counts(f)
  expect_equal(rd$check$measurements$n_banner_rows, 4L)
  expect_equal(dim(rd$counts), dim(m))
  expect_equal(colnames(rd$counts), colnames(m))
  expect_true(any(grepl("rows above the header", rd$check$evidence)))

  # an ordinary file with no banner rows is read exactly as before (n_banner = 0)
  f2 <- tempfile(fileext = ".tsv")
  utils::write.table(cbind(gene = rownames(m), as.data.frame(m)), f2, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_equal(at_read_counts(f2)$check$measurements$n_banner_rows, 0L)
})
