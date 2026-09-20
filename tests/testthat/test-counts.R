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

test_that("attest() merges checks into one report", {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  m <- readRDS(f)$airway

  rep_ok <- attest(m)
  expect_s3_class(rep_ok, "attest_report")
  expect_equal(rep_ok$verdict, "PERMITTED")
  expect_equal(rep_ok$checks[["value scale"]]$verdict, "PERMITTED")

  rep_bad <- attest(round(t(t(m) / colSums(m)) * 1e6))
  expect_equal(rep_bad$verdict, "NOT PERMITTED")

  # the report is as strong as its weakest check
  expect_equal(at_worst(c("PERMITTED", "CAUTION")), "CAUTION")
  expect_equal(at_worst(c("PERMITTED", "UNKNOWN", "NOT PERMITTED")), "NOT PERMITTED")
  expect_equal(at_worst(c("PERMITTED", "PERMITTED")), "PERMITTED")

  l <- attest_as_list(rep_bad)
  expect_equal(l$schema, "attest/report/1")
  expect_equal(l$verdict, "NOT PERMITTED")
  expect_true(nzchar(l$checks[["value scale"]]$headline))
})

test_that("completeness distinguishes complete matrices from filtered ones", {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  fx <- readRDS(f)

  # every complete matrix keeps genes that are zero everywhere, whatever the species
  for (nm in names(fx)) {
    res <- attest_completeness(fx[[nm]])
    expect_equal(res$verdict, "PERMITTED", info = nm)
    expect_true(res$measurements$n_all_zero > 0, info = nm)
  }

  m <- fx$airway
  drop_zero <- attest_completeness(m[rowSums(m) > 0, ])
  expect_equal(drop_zero$verdict, "CAUTION")
  expect_equal(drop_zero$measurements$min_row_total, 1)

  thresholded <- attest_completeness(m[rowSums(m) >= 10, ])
  expect_equal(thresholded$verdict, "CAUTION")
  expect_equal(thresholded$measurements$min_row_total, 10)   # the floor recovers the filter

  # gene selection is only assessable when the annotation size is supplied
  expect_true(any(grepl("n_expected", attest_completeness(m)$not_assessed)))
  sel <- attest_completeness(m[1:1000, ], n_expected = nrow(m))   # a tenth of the fixture
  expect_equal(sel$verdict, "CAUTION")

  # and the merged report now carries both checks
  rep <- attest(m)
  expect_equal(length(rep$checks), 2)
  expect_true("completeness" %in% names(rep$checks))
  expect_equal(attest(m[rowSums(m) >= 10, ])$verdict, "CAUTION")  # weakest check wins
})

test_that("identity infers sex, catches mislabels and duplicate libraries", {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  fx <- readRDS(f)
  m <- fx$airway

  res <- attest_identity(m)
  expect_equal(res$verdict, "PERMITTED")
  expect_equal(unname(table(res$measurements$sex_called)[c("female", "male")]), c(2L, 6L))

  md <- data.frame(sex = c(rep("male", 6), "female", "female"))
  expect_equal(attest_identity(m, metadata = md)$verdict, "PERMITTED")

  md_wrong <- md; md_wrong$sex[7] <- "male"
  bad <- attest_identity(m, metadata = md_wrong)
  expect_equal(bad$verdict, "NOT PERMITTED")
  expect_equal(bad$measurements$sex_mismatches, 1)

  # duplicates: none in any clean cohort, caught when planted, in both species
  for (nm in names(fx)) expect_equal(nrow(at_duplicate_pairs(fx[[nm]])$pairs), 0, info = nm)
  expect_equal(nrow(at_duplicate_pairs(cbind(m, copy = m[, 1]))$pairs), 1)
  set.seed(1)
  reseq <- rbinom(nrow(fx$fission), size = as.integer(fx$fission[, 1]), prob = 0.9)
  expect_equal(nrow(at_duplicate_pairs(cbind(fx$fission, reseq = reseq))$pairs), 1)

  # non-human data says so instead of guessing
  expect_true(any(grepl("sex check", attest_identity(fx$fission)$not_assessed)))

  expect_equal(length(attest(m)$checks), 3)
  expect_equal(attest(m, metadata = md_wrong)$verdict, "NOT PERMITTED")
})
