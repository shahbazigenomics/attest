at_fixture_airway <- function() {
  f <- system.file("extdata", "fixtures.rds", package = "attest")
  if (!nzchar(f)) f <- "inst/extdata/fixtures.rds"
  if (!file.exists(f)) f <- "../../inst/extdata/fixtures.rds"
  readRDS(f)$airway
}
cond8 <- factor(rep(c("ctrl", "trt"), each = 4))

test_that("confounding costs precision, and detectability is where that is charged", {
  m <- at_fixture_airway()
  partial  <- data.frame(cond = cond8, batch = factor(c("A","A","A","A","B","B","B","A")))
  balanced <- data.frame(cond = cond8, batch = factor(rep(c("A","B"), times = 4)))
  confounded <- data.frame(cond = cond8, batch = factor(rep(c("A","B"), each = 4)))

  d0 <- attest_detectability(m, group = cond8)
  d1 <- attest_detectability(m, group = cond8, metadata = partial, design = ~ batch + cond)
  expect_gt(d1$measurements$median_min_detectable_fc, d0$measurements$median_min_detectable_fc)
  expect_true(any(grepl("treated as", d1$evidence)))
  expect_equal(round(d1$measurements$design_vif, 2), 2.5)

  # a balanced design costs no precision to confounding: no inflation charged
  # (its noise can still differ - fitting the batch term is the point)
  db <- attest_detectability(m, group = cond8, metadata = balanced, design = ~ batch + cond)
  expect_equal(db$measurements$design_vif, 1)
  expect_false(any(grepl("treated as", db$evidence)))

  # an effect the design cannot estimate has no precision to report
  expect_equal(attest_detectability(m, group = cond8, metadata = confounded,
                                    design = ~ batch + cond)$verdict, "UNKNOWN")

  # in the report: the verdict stands, the cost shows up under scope
  r <- attest(m, metadata = partial, design = ~ batch + cond)
  expect_equal(r$verdict, "PERMITTED")
  expect_equal(round(r$checks$detectability$measurements$design_vif, 2), 2.5)
})

test_that("a column the design leaves out is caught when it lines up with the effect", {
  m <- at_fixture_airway()

  # batch nested in condition, not in the formula
  nest <- data.frame(cond = cond8, batch = factor(c("A","A","B","B","C","C","D","D")))
  r <- attest_design(m, nest, ~ cond)
  expect_equal(r$verdict, "CAUTION")
  expect_equal(names(r$measurements$omitted_flags), "batch")

  # two libraries per donor, analysed as four independent replicates
  donor <- data.frame(cond = cond8, donor = factor(rep(paste0("d", 1:4), each = 2)))
  r <- attest_design(m, donor, ~ cond)
  expect_equal(r$verdict, "CAUTION")
  expect_true(any(grepl("from 2 donor values", r$evidence)))
  # and putting it in the formula is what DESeq2 would refuse
  expect_equal(attest_design(m, donor, ~ donor + cond)$verdict, "NOT PERMITTED")

  # a covariate that separates the groups, with the chance of that stated
  rin <- data.frame(cond = cond8, RIN = c(6.1, 6.5, 6.8, 7.0, 8.2, 8.5, 9.0, 9.3))
  r <- attest_design(m, rin, ~ cond)
  expect_equal(r$verdict, "CAUTION")
  expect_equal(r$measurements$omitted_flags$RIN$chance, 2 / choose(8, 4))
  expect_true(any(grepl("1 time in 35", r$evidence)))
})

test_that("columns that line up for harmless reasons are left alone", {
  m <- at_fixture_airway()

  # crossed with the condition: no problem, whether in the formula or not
  crossed <- data.frame(cond = cond8, batch = factor(rep(c("A","B"), times = 4)))
  expect_equal(attest_design(m, crossed, ~ cond)$verdict, "PERMITTED")

  # identifiers, constants, a relabelling of the condition, a count-derived
  # column, and a row index on a sheet sorted by condition
  sheet <- data.frame(cond = cond8, sample = paste0("S", 1:8), constant = "x",
                      treatment = ifelse(cond8 == "ctrl", "none", "drug"),
                      sizeFactor = c(.5, .6, .7, .8, 1.2, 1.3, 1.4, 1.5), row = 1:8)
  r <- attest_design(m, sheet, ~ cond)
  expect_equal(r$verdict, "PERMITTED")
  expect_true(all(c("sample", "constant", "treatment", "sizeFactor", "row") %in%
                    r$measurements$omitted_skipped))

  # overlapping numeric covariate
  rin <- data.frame(cond = cond8, RIN = c(6.1, 8.5, 6.8, 9.0, 8.2, 6.5, 9.3, 7.0))
  expect_equal(attest_design(m, rin, ~ cond)$verdict, "PERMITTED")

  # an airway-shaped sample sheet
  aw <- data.frame(SampleName = paste0("GSM", 1:8),
                   cell = factor(rep(c("N61311", "N052611", "N080611", "N061011"), each = 2)),
                   dex = factor(rep(c("untrt", "trt"), times = 4)), albut = "untrt",
                   Run = paste0("SRR", 1:8), avgLength = c(126, 126, 126, 87, 120, 126, 101, 98))
  expect_equal(attest_design(m, aw, ~ cell + dex)$verdict, "PERMITTED")
})

test_that("a variable's coefficient columns are not confused with another variable's", {
  # grep("^var", colnames(mm)) matched any OTHER term whose name happened to
  # start with var's name too (~ cell + celltype attributed celltype's
  # coefficients, and its VIF, to cell as well)
  md <- data.frame(cell = factor(rep(c("A", "B", "C", "D"), each = 2)),
                   celltype = factor(rep(c("x", "y"), 4)))
  mm <- stats::model.matrix(~ cell + celltype, md)

  cell_cols <- at_coefficient_columns(mm, "cell", md)
  expect_equal(sort(colnames(mm)[cell_cols]), c("cellB", "cellC", "cellD"))
  expect_false(any(grepl("^celltype", colnames(mm)[cell_cols])))

  celltype_cols <- at_coefficient_columns(mm, "celltype", md)
  expect_equal(colnames(mm)[celltype_cols], "celltypey")

  # a continuous variable's own column, not confused with a longer name either
  md2 <- data.frame(rin = c(6.1, 6.5, 6.8, 7.0, 8.2, 8.5, 9.0, 9.3), rinbatch = factor(rep(1:2, 4)))
  mm2 <- stats::model.matrix(~ rin + rinbatch, md2)
  rin_cols <- at_coefficient_columns(mm2, "rin", md2)
  expect_equal(colnames(mm2)[rin_cols], "rin")
})

test_that("a design term with an NA silently dropping rows is caught, not miscounted", {
  # model.matrix()'s default na.action is na.omit: an NA in a design term
  # (batch here) drops that sample's row without erroring, so rank, group
  # sizes and VIF were previously computed on fewer rows than n_samples claimed
  m <- at_fixture_airway()
  md <- data.frame(cond = cond8, batch = factor(c("A", "A", "B", "B", "A", "A", NA, "B")))
  r <- attest_design(m, md, ~ batch + cond)
  expect_equal(r$verdict, "UNKNOWN")
  expect_match(r$headline, "dropped")
})

test_that("the chance the message quotes is the chance that happens", {
  # an unrelated covariate separates 4 + 4 samples 2 / choose(8, 4) of the time;
  # 2,000 draws puts the realised rate within ~0.008 of that with high probability
  set.seed(42)
  hits <- vapply(seq_len(2000), function(i)
    length(at_omitted_variables(data.frame(cond = cond8, x = stats::rnorm(8)),
                                "cond", "cond")$flags) > 0, logical(1))
  expect_lt(abs(mean(hits) - 2 / 70), 0.012)
})
