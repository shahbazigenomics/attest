# How often does attest_counts mis-flag genuine raw counts?
# Simulated regimes only: real-data validation is harness.R (needs Bioconductor).
for (f in list.files("R", full.names = TRUE)) source(f)
set.seed(42)

sim_counts <- function(n_genes, n_samples, mean_lib, lib_cv, size = 4, umi = FALSE) {
  mu <- if (umi) rgamma(n_genes, 0.15, 0.6) else rgamma(n_genes, 0.6, 0.02)
  mu <- mu * (mean_lib / sum(mu))
  depth <- exp(rnorm(n_samples, 0, lib_cv))
  m <- vapply(seq_len(n_samples),
              function(j) rnbinom(n_genes, mu = mu * depth[j], size = size),
              numeric(n_genes))
  dimnames(m) <- list(paste0("g", seq_len(n_genes)), paste0("s", seq_len(n_samples)))
  m
}

regimes <- list(
  list(id = "bulk, typical (20M, CV 0.35)",      n = 20000, s = 8,  lib = 20e6, cv = 0.35),
  list(id = "bulk, balanced (20M, CV 0.05)",     n = 20000, s = 8,  lib = 20e6, cv = 0.05),
  list(id = "bulk, very balanced (20M, CV 0.02)",n = 20000, s = 8,  lib = 20e6, cv = 0.02),
  list(id = "bulk, n=4 (20M, CV 0.35)",          n = 20000, s = 4,  lib = 20e6, cv = 0.35),
  list(id = "bulk, large cohort (n=50)",         n = 20000, s = 50, lib = 20e6, cv = 0.35),
  list(id = "shallow 3' (1M, CV 0.35)",          n = 20000, s = 8,  lib = 1e6,  cv = 0.35),
  list(id = "very shallow (0.3M, CV 0.35)",      n = 20000, s = 8,  lib = 3e5,  cv = 0.35),
  list(id = "pseudobulk UMI (0.5M, sparse)",     n = 20000, s = 8,  lib = 5e5,  cv = 0.35, umi = TRUE),
  list(id = "full annotation (60k genes)",       n = 60000, s = 8,  lib = 20e6, cv = 0.35)
)

reps <- 20
rows <- list()
for (r in regimes) {
  v <- replicate(reps, {
    m <- sim_counts(r$n, r$s, r$lib, r$cv, umi = isTRUE(r$umi))
    attest_counts(m)$verdict
  })
  rows[[length(rows) + 1]] <- data.frame(
    regime = r$id,
    PERMITTED = sum(v == "PERMITTED"),
    CAUTION   = sum(v == "CAUTION"),
    NOT       = sum(v == "NOT PERMITTED"),
    UNKNOWN   = sum(v == "UNKNOWN"),
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, rows)
res$false_alarm_pct <- round(100 * (res$CAUTION + res$NOT + res$UNKNOWN) / reps, 1)
print(res, row.names = FALSE)
write.csv(res, "validation/simulation_specificity.csv", row.names = FALSE)
