# How far is attest_detectability() from DESeq2 on paired designs, and why?
#
# Two routes by which a paired design (~ block + condition) makes the current
# detectability conservative:
#   1. the number of genes that really differ (R), which sets the BH level each
#      gene is judged at - estimated by a quick two-group Wald test;
#   2. the per-gene dispersion, pooled within condition groups, which leaves the
#      block (donor, cell line) variance in.
# This measures both against DESeq2 on airway (cell lines) and Kang 2018
# pseudobulk (donors). Needs DESeq2 and validation/raw_matrices.rds,
# validation/umi_pseudobulk.rds.

suppressPackageStartupMessages(library(DESeq2))
for (f in list.files("R", full.names = TRUE)) source(f)

phi_design <- function(norm, mm) {
  r <- stats::lm.fit(mm, t(norm))$residuals
  v <- colSums(r^2) / (nrow(mm) - ncol(mm)); mu <- rowMeans(norm)
  phi <- (v - mu) / pmax(mu^2, 1e-12); phi[!is.finite(phi)] <- NA; phi <- pmax(phi, 0)
  fill <- stats::median(phi[phi > 0], na.rm = TRUE); phi[is.na(phi)] <- fill
  0.5 * phi + 0.5 * fill
}
mdfc <- function(mu, phi, n_eff, a) exp((qnorm(1 - a/2) + qnorm(0.8)) * sqrt((2/n_eff) * (1/mu + phi)))

compare <- function(label, m, sh, block, cond) {
  g <- droplevels(factor(sh[[cond]]))
  R <- vapply(list(stats::as.formula(paste("~", cond)),
                   stats::as.formula(paste("~", block, "+", cond))), function(des) {
    dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(m, sh, des), quiet = TRUE))
    sum(results(dds, alpha = 0.05)$padj < 0.05, na.rm = TRUE) }, numeric(1))
  quick <- attest_detectability(m, group = g)$measurements$n_rejected_estimate
  sf <- at_size_factors(m); norm <- t(t(m) / sf); mu <- rowMeans(norm); ok <- mu >= 10
  n <- table(g); n_eff <- 2 / (1/n[1] + 1/n[2]); n_test <- sum(mu > 0)
  phi_g <- at_dispersion_moments(norm, g)
  phi_d <- phi_design(norm, stats::model.matrix(stats::as.formula(paste("~", block, "+", cond)), sh))
  a_q <- 0.05 * quick / n_test; a_d <- 0.05 * R[2] / n_test
  med <- function(phi, a) stats::median(mdfc(mu[ok], phi[ok], n_eff, a))
  data.frame(dataset = label, R_quick = quick, R_deseq2_unpaired = R[1], R_deseq2_paired = R[2],
             mdfc_now = med(phi_g, a_q), mdfc_R_fixed = med(phi_g, a_d),
             mdfc_noise_fixed = med(phi_d, a_q), mdfc_both = med(phi_d, a_d),
             disp_within_group = stats::median(phi_g[ok]), disp_design = stats::median(phi_d[ok]))
}

m  <- readRDS("validation/raw_matrices.rds")$airway
sa <- data.frame(cell = factor(rep(c("N61311", "N052611", "N080611", "N061011"), each = 2)),
                 dex  = factor(rep(c("untrt", "trt"), times = 4), levels = c("untrt", "trt")),
                 row.names = colnames(m))
k  <- readRDS("validation/umi_pseudobulk.rds")
res <- rbind(compare("airway", m, sa, "cell", "dex"),
             compare("Kang 2018 pseudobulk", k$counts, k$sheet, "donor", "cond"))
print(res, digits = 3, row.names = FALSE)
utils::write.csv(res, "validation/detectability_vs_deseq2.csv", row.names = FALSE)
