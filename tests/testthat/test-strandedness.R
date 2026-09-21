# Files written by real STAR and real featureCounts on simulated libraries of
# known protocol (validation/strand_sim/): reverse-stranded (dUTP), forward,
# unstranded, and reverse with 60% of reads outside every gene.

strand_file <- function(...) {
  d <- system.file("extdata", "strand", package = "attest")
  if (!nzchar(d)) d <- "inst/extdata/strand"
  if (!dir.exists(d)) d <- "../../inst/extdata/strand"
  file.path(d, c(...))
}
rpg <- function(...) strand_file(paste0(c(...), ".ReadsPerGene.out.tab"))

test_that("STAR's three columns give the protocol", {
  r <- attest_strandedness(rpg("reverse_s1", "reverse_s2", "reverse_s3"))
  expect_equal(r$verdict, "PERMITTED")
  expect_equal(r$measurements$recommended_column, "reverse")
  expect_true(all(r$measurements$forward_share < 0.2))

  expect_equal(attest_strandedness(rpg("forward_s1"))$measurements$recommended_column, "forward")
  expect_equal(attest_strandedness(rpg("unstranded_s1", "unstranded_s2"))$measurements$recommended_column, "unstranded")
  # background reads do not move the call
  expect_equal(attest_strandedness(rpg("reverse_bg_s1"))$measurements$protocol[[1]], "reverse")
  # libraries made two ways
  expect_equal(attest_strandedness(rpg("reverse_s1", "forward_s1"))$verdict, "CAUTION")
})

test_that("the matrix is matched to the STAR column it came from", {
  files <- rpg("reverse_s1", "reverse_s2", "reverse_s3")
  cols  <- lapply(files, function(p) at_read_star(p)$counts)
  take  <- function(k) sapply(cols, function(m) m[, k])

  expect_equal(attest_strandedness(files, counts = take("reverse"))$verdict, "PERMITTED")
  wrong <- attest_strandedness(files, counts = take("forward"))
  expect_equal(wrong$verdict, "NOT PERMITTED")
  expect_true(all(wrong$measurements$matrix_column == "forward"))
  expect_equal(attest_strandedness(files, counts = take("unstranded"))$verdict, "CAUTION")

  un <- rpg("unstranded_s1", "unstranded_s2")
  half <- sapply(lapply(un, function(p) at_read_star(p)$counts), function(m) m[, "reverse"])
  expect_equal(attest_strandedness(un, counts = half)$verdict, "CAUTION")

  # a matrix that matches no column is reported as unassessed, not guessed
  r <- attest_strandedness(files, counts = take("reverse") + 1L)
  expect_true(any(grepl("no matrix column is identical", r$not_assessed)))
})

test_that("a featureCounts summary raises a collapse and confirms nothing else", {
  s1 <- attest_strandedness(strand_file("fc_s1.txt", "fc_s1.txt.summary"))
  expect_equal(s1$verdict, "CAUTION")
  expect_equal(s1$measurements$strand_setting, 1L)
  expect_true(all(s1$measurements$assigned_share < 0.25))

  # the right setting cannot be confirmed from one run...
  expect_equal(attest_strandedness(strand_file("fc_s2.txt", "fc_s2.txt.summary"))$verdict, "UNKNOWN")
  # ...and neither can the trap: an unstranded library counted as stranded
  trap <- attest_strandedness(strand_file("fc_unstranded_s2.txt", "fc_unstranded_s2.txt.summary"))
  expect_equal(trap$verdict, "UNKNOWN")
  expect_true(all(trap$measurements$assigned_share > 0.5))
  # unstranded counting is valid for any library
  expect_equal(attest_strandedness(strand_file("fc_s0.txt", "fc_s0.txt.summary"))$verdict, "PERMITTED")
})

test_that("attest_file() finds the summary beside featureCounts output", {
  r1 <- attest_file(strand_file("fc_s1.txt"))
  expect_equal(r1$checks$strandedness$verdict, "CAUTION")
  expect_equal(r1$verdict, "CAUTION")

  # a summary that cannot answer is "not run", and leaves the verdict alone
  r2 <- attest_file(strand_file("fc_s2.txt"))
  expect_false("strandedness" %in% names(r2$checks))
  expect_true(any(grepl("^strandedness", r2$not_run)))
  expect_equal(r2$verdict, "PERMITTED")
})

test_that("bad input fails typed", {
  expect_equal(attest_strandedness(tempfile())$verdict, "UNKNOWN")
  p <- tempfile(); writeLines(c("a\tb", "1\t2"), p)
  expect_equal(attest_strandedness(p)$verdict, "UNKNOWN")
})
