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
  # the Poisson-floor check needs 200 expressed genes, so it can never run on 30
  expect_true(any(grepl("variance/mean", res$not_assessed)))

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
    expect_true(as.integer(res$measurements$n_all_zero) > 0L, info = nm)
  }

  m <- fx$airway
  drop_zero <- attest_completeness(m[rowSums(m) > 0, ])
  expect_equal(drop_zero$verdict, "CAUTION")
  expect_equal(as.numeric(drop_zero$measurements$min_row_total), 1)

  thresholded <- attest_completeness(m[rowSums(m) >= 10, ])
  expect_equal(thresholded$verdict, "CAUTION")
  expect_equal(as.numeric(thresholded$measurements$min_row_total), 10)   # the floor recovers the filter

  # gene selection is only assessable when the annotation size is supplied
  expect_true(any(grepl("n_expected", attest_completeness(m)$not_assessed)))
  sel <- attest_completeness(m[1:1000, ], n_expected = nrow(m))   # a tenth of the fixture
  expect_equal(sel$verdict, "CAUTION")

  # and the merged report now carries both checks
  rep <- attest(m)
  expect_true(all(c("value scale", "completeness", "identity") %in% names(rep$checks)))
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
  called <- res$measurements$sex_called
  expect_equal(sum(called == "female"), 2L)
  expect_equal(sum(called == "male"), 6L)

  md <- data.frame(sex = c(rep("male", 6), "female", "female"))
  expect_equal(attest_identity(m, metadata = md)$verdict, "PERMITTED")

  md_wrong <- md; md_wrong$sex[7] <- "male"
  bad <- attest_identity(m, metadata = md_wrong)
  expect_equal(bad$verdict, "NOT PERMITTED")
  expect_equal(as.integer(bad$measurements$sex_mismatches), 1L)

  # duplicates: none in any clean cohort, caught when planted, in both species
  for (nm in names(fx)) expect_equal(nrow(at_duplicate_pairs(fx[[nm]])$pairs), 0L, info = nm)
  expect_equal(nrow(at_duplicate_pairs(cbind(m, copy = m[, 1]))$pairs), 1L)
  set.seed(1)
  reseq <- rbinom(nrow(fx$fission), size = as.integer(fx$fission[, 1]), prob = 0.9)
  expect_equal(nrow(at_duplicate_pairs(cbind(fx$fission, reseq = reseq))$pairs), 1L)

  # non-human data says so instead of guessing
  expect_true(any(grepl("sex check", attest_identity(fx$fission)$not_assessed)))

  expect_true(all(c("value scale", "completeness", "identity") %in% names(attest(m)$checks)))
  expect_equal(attest(m, metadata = md_wrong)$verdict, "NOT PERMITTED")
})

test_that("normalised matrices are caught at any depth, and flat totals are flagged", {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  fx <- readRDS(f)

  for (nm in names(fx)) {
    m <- fx[[nm]]
    # a fixed size-factor threshold is depth-dependent (1.0057 at 22M reads,
    # 1.0154 at 1M); the size-factor-to-total ratio is not
    expect_equal(attest_counts(round(t(t(m) / at_size_factors(m))))$verdict, "CAUTION", info = nm)
    # scaled to a common total that is not 1e6: not counts, or rarefied
    expect_equal(attest_counts(round(t(t(m) / colSums(m)) * 5e5))$verdict, "CAUTION", info = nm)
    # and the same matrix shallower is still raw counts
    expect_equal(attest_counts(round(m / 40))$verdict, "PERMITTED", info = nm)
  }
})

test_that("design adequacy: estimability, partial confounding, replication", {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  m <- readRDS(f)$airway

  balanced  <- data.frame(cond = factor(rep(c("ctrl", "trt"), each = 4)),
                          batch = factor(rep(c("A", "B"), times = 4)))
  confounded <- data.frame(cond = factor(rep(c("ctrl", "trt"), each = 4)),
                           batch = factor(rep(c("A", "B"), each = 4)))
  partial    <- data.frame(cond = factor(rep(c("ctrl", "trt"), each = 4)),
                           batch = factor(c("A", "A", "A", "A", "B", "B", "B", "A")))

  ok <- attest_design(m, balanced, ~ batch + cond)
  expect_equal(ok$verdict, "PERMITTED")
  expect_equal(unname(ok$measurements$max_vif), 1)
  expect_equal(unname(ok$measurements$effective_n), 8)

  # complete confounding: caught before DESeq2 refuses and before limma quietly
  # returns NA coefficients
  conf <- attest_design(m, confounded, ~ batch + cond)
  expect_equal(conf$verdict, "NOT PERMITTED")
  expect_true(length(conf$measurements$not_estimable) >= 1)

  # partial confounding: nothing else reports this at all
  part <- attest_design(m, partial, ~ batch + cond)
  expect_equal(unname(round(part$measurements$max_vif, 2)), 2.5)
  expect_equal(unname(round(part$measurements$effective_n, 1)), 3.2)
  expect_equal(part$verdict, "CAUTION")

  # replication
  k <- c(1, 2, 5, 6)
  expect_equal(attest_design(m[, k], balanced[k, ], ~ batch + cond)$verdict, "CAUTION")
  k1 <- c(1, 2, 3, 5)
  expect_equal(attest_design(m[, k1], balanced[k1, ], ~ cond)$verdict, "NOT PERMITTED")

  # malformed input never stops
  expect_equal(attest_design(m, balanced, ~ batch + missing_col)$verdict, "UNKNOWN")
  expect_equal(attest_design(m, balanced[1:6, ], ~ batch + cond)$verdict, "UNKNOWN")

  # runs inside the report only when a design is supplied
  expect_false("design" %in% names(attest(m)$checks))
  expect_true("design" %in% names(attest(m, metadata = balanced, design = ~ batch + cond)$checks))
  expect_equal(attest(m, metadata = confounded, design = ~ batch + cond)$verdict, "NOT PERMITTED")
})

test_that("detectability reports per-gene thresholds and needs the groups", {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  m <- readRDS(f)$airway
  dex <- factor(rep(c("untrt", "trt"), times = 4))

  d <- attest_detectability(m, group = dex)
  mdfc <- d$measurements$min_detectable_fc
  expect_equal(length(mdfc), nrow(m))
  expect_true(all(mdfc >= 1))                       # a threshold below 1 is nonsense
  expect_true(all(is.infinite(mdfc[rowSums(m) == 0])))  # unexpressed genes: never

  # more samples must never need a larger change
  d4 <- attest_detectability(m, group = dex)$measurements$median_min_detectable_fc
  k  <- c(1, 2, 5, 6)
  d2 <- attest_detectability(m[, k], group = dex[k])$measurements$median_min_detectable_fc
  expect_true(d4 < d2)

  # highly expressed genes must be easier to detect than quiet ones
  mu <- d$measurements$mean_count
  hi <- mu > stats::quantile(mu[mu > 0], 0.9); lo <- mu > 0 & mu < stats::quantile(mu[mu > 0], 0.25)
  expect_true(stats::median(mdfc[hi]) < stats::median(mdfc[lo]))

  # inputs it cannot answer
  expect_equal(attest_detectability(m)$verdict, "UNKNOWN")                       # no groups
  expect_equal(attest_detectability(m, group = rep("a", 8))$verdict, "UNKNOWN")  # one level
  expect_equal(attest_detectability(m[, 1:3], group = factor(c("a","a","b")))$verdict, "UNKNOWN")

  # runs in the report when a two-level factor is available
  md <- data.frame(cell = factor(rep(c("A","B","C","D"), each = 2)), dex = dex)
  expect_true("detectability" %in% names(attest(m, metadata = md, design = ~ cell + dex)$checks))
  expect_false("detectability" %in% names(attest(m)$checks))
})

test_that("scope checks report but do not set the overall verdict", {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  m <- readRDS(f)$airway
  md <- data.frame(cell = factor(rep(c("A","B","C","D"), each = 2)),
                   dex  = factor(rep(c("untrt","trt"), times = 4)))

  r <- attest(m, metadata = md, design = ~ cell + dex)
  expect_equal(r$checks[["detectability"]]$kind, "scope")
  expect_equal(r$checks[["value scale"]]$kind, "fault")

  # detectability is a limitation, not a defect: a good dataset stays PERMITTED
  # even when most genes cannot show a 2-fold change
  expect_equal(r$checks[["detectability"]]$verdict, "CAUTION")
  expect_equal(r$verdict, "PERMITTED")

  # a real fault still sets the verdict
  cpm <- round(t(t(m) / colSums(m)) * 1e6)
  expect_equal(attest(cpm, metadata = md, design = ~ cell + dex)$verdict, "NOT PERMITTED")

  # and the distinction survives into the machine-readable form
  kinds <- vapply(attest_as_list(r)$checks, function(c) c$kind, character(1))
  expect_equal(sum(kinds == "scope"), 1L)
})
