# Does "detectable at 80% power" hold when the test is DESeq2 itself?
#
# Templates come from DESeq2's own fit to real data - per-gene baseline, block
# effects (cell line / donor) and dispersions, and the real size factors - so
# the simulated experiments have real structure. For each template:
#   1. a pilot experiment: 10% of genes change 2-fold (random direction), so
#      the number of truly changing genes, which sets the BH level, is realistic;
#   2. attest_detectability() on the pilot states each gene's minimum
#      detectable fold change and the per-gene level it will be judged at;
#   3. a test experiment: the same 10% of genes now change by exactly their
#      claimed amount; DESeq2 with the design tests every gene;
#   4. realised power = the share of those genes DESeq2 detects, at the
#      claimed per-gene level (checks the noise model) and at padj < 0.05
#      (end to end, with DESeq2's own BH).
# Both of attest's noise estimates are scored on the same simulated data:
# DESeq2's own dispersions (the default when DESeq2 is installed) and attest's
# moment estimate (the fallback). Both are given the design.
#
# History: the first version of detectability pooled noise within groups and
# plugged in each gene's average count. Against DESeq2 it gave 0.66 on an
# unpaired design, and 0.63-0.88 on paired ones only because an overestimate
# of noise happened to offset an optimistic formula.
#
# Needs DESeq2, validation/raw_matrices.rds and validation/umi_pseudobulk.rds.
# Writes validation/detectability_calibration_deseq2.csv.

suppressPackageStartupMessages(library(DESeq2))
for (f in list.files("R", full.names = TRUE)) source(f)
set.seed(20260921)

template <- function(m, sh, design, n_genes = 4000) {
  dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(m, sh, design), quiet = TRUE))
  keep <- which(mcols(dds)$baseMean >= 10 & is.finite(dispersions(dds)))
  keep <- sort(sample(keep, min(n_genes, length(keep))))
  # DESeq2's coefficients come in model-matrix order, under different names
  # ("Intercept", "cell_B_vs_A"); the tested term is the last column of both
  cf <- coef(dds)[keep, , drop = FALSE]
  mm <- stats::model.matrix(design, sh)
  stopifnot(ncol(cf) == ncol(mm))
  base <- cf[, -ncol(cf), drop = FALSE]                           # intercept + blocks, log2
  mm_base <- mm[, -ncol(mm), drop = FALSE]
  list(log2mu = base %*% t(mm_base), disp = dispersions(dds)[keep],
       sf = sizeFactors(dds), sheet = sh, design = design,
       cond = droplevels(factor(sh[[utils::tail(all.vars(design), 1)]])))
}

simulate <- function(tp, lfc) {                                   # lfc: log2, per gene
  trt <- as.integer(tp$cond == levels(tp$cond)[2])
  mu <- 2^(tp$log2mu + outer(lfc, trt)) * rep(tp$sf, each = nrow(tp$log2mu))
  y <- matrix(stats::rnbinom(length(mu), mu = mu, size = 1 / tp$disp), nrow(mu))
  dimnames(y) <- list(paste0("g", seq_len(nrow(y))), rownames(tp$sheet))
  y
}

score <- function(tp, label, reps = 3) {
  out <- list()
  G <- nrow(tp$log2mu)
  for (r in seq_len(reps)) {
    de <- sample(G, round(0.1 * G))
    dir <- sample(c(-1, 1), length(de), replace = TRUE)
    lfc0 <- numeric(G); lfc0[de] <- dir
    pilot <- simulate(tp, lfc0)
    for (method in c("deseq2", "moments")) {
      d <- attest_detectability(pilot, group = tp$cond, metadata = tp$sheet, design = tp$design,
                                dispersion = method)
      fc <- d$measurements$min_detectable_fc
      use <- de[is.finite(fc[de]) & fc[de] < 64]
      lfc <- numeric(G); lfc[use] <- dir[match(use, de)] * log2(fc[use])
      y <- simulate(tp, lfc)
      dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(y, tp$sheet, tp$design), quiet = TRUE))
      res <- results(dds, alpha = 0.05)
      out[[length(out) + 1]] <- data.frame(
        template = label, n_per_group = min(table(tp$cond)), method = method, rep = r,
        median_claimed_fc = stats::median(fc[use]),
        power_at_claimed_level = mean(res$pvalue[use] < d$measurements$per_gene_alpha, na.rm = TRUE),
        power_padj = mean(res$padj[use] < 0.05, na.rm = TRUE))
    }
  }
  do.call(rbind, out)
}

raw <- readRDS("validation/raw_matrices.rds")
m <- raw$airway[rowSums(raw$airway) > 0, ]
sa <- data.frame(cell = factor(rep(c("N61311", "N052611", "N080611", "N061011"), each = 2)),
                 dex  = factor(rep(c("untrt", "trt"), times = 4), levels = c("untrt", "trt")),
                 row.names = colnames(m))
k <- readRDS("validation/umi_pseudobulk.rds")
kc <- k$counts[rowSums(k$counts) > 0, ]; ks <- k$sheet

res <- list()
res[[1]] <- score(template(m, sa, ~ cell + dex), "airway, paired (4 cell lines)")
res[[2]] <- score(template(m, sa, ~ dex), "airway, unpaired")
for (nd in c(3, 4, 6, 8)) {
  donors <- levels(ks$donor)[seq_len(nd)]
  keep <- ks$donor %in% donors
  sh <- droplevels(ks[keep, ]); y <- kc[, keep]
  res[[length(res) + 1]] <- score(template(y, sh, ~ donor + cond), sprintf("Kang, paired (%d donors)", nd))
}
res <- do.call(rbind, res)
utils::write.csv(res, "validation/detectability_calibration_deseq2.csv", row.names = FALSE)

agg <- stats::aggregate(cbind(median_claimed_fc, power_at_claimed_level, power_padj) ~ template + method,
                        res, mean)
agg <- agg[order(agg$template, agg$method), ]
cat("\nclaimed power 0.80; realised = share DESeq2 detects (mean of 3 simulated experiments)\n\n")
print(format(agg, digits = 3), row.names = FALSE)
