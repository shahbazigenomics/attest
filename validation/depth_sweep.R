# Does attest_counts() still get it right at shallow depth?
#
# Binomial thinning of real raw matrices is the same thing, statistically, as
# sequencing the same libraries less deeply: every read is kept independently
# with probability f. One common f per depth keeps each dataset's real
# library-size spread and its real biological variance. Each thinned matrix is
# then transformed the way the harness transforms full-depth data, and every
# verdict is compared with the one it should get.
#
# Usage, from the attest project:
#   Rscript validation/depth_sweep.R validation/raw_matrices.rds
# Writes validation/depth_sweep.csv.

args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1] else "validation/raw_matrices.rds"
for (f in list.files("R", full.names = TRUE)) source(f)
raw <- readRDS(path)

depths <- c(10e6, 3e6, 1e6, 3e5, 1e5, 3e4)       # mean reads per sample
reps   <- 5
expect <- c(raw = "PERMITTED", raw_filtered = "PERMITTED",
            norm_rounded = "CAUTION", scaled_5e5 = "CAUTION",
            cpm = "NOT PERMITTED", cpm_rounded = "NOT PERMITTED",
            tpm_rounded = "NOT PERMITTED", fpkm = "NOT PERMITTED",
            log2cpm = "NOT PERMITTED")

transforms <- function(m, len) {
  cs  <- colSums(m)
  cpm <- t(t(m) / cs) * 1e6
  rpk <- m / len
  tpm <- t(t(rpk) / colSums(rpk)) * 1e6
  sf  <- at_size_factors(m)
  list(raw          = m,
       raw_filtered = m[rowSums(m) >= 10, , drop = FALSE],
       norm_rounded = if (all(is.finite(sf))) round(t(t(m) / sf)) else NULL,
       scaled_5e5   = round(t(t(m) / cs) * 5e5),
       cpm          = cpm,
       cpm_rounded  = round(cpm),
       tpm_rounded  = round(tpm),
       fpkm         = t(t(m / len * 1e3) / cs) * 1e6,
       log2cpm      = log2(cpm + 1))
}

set.seed(20260921)
rows <- list()
for (ds in names(raw)) {
  m0  <- raw[[ds]]
  len <- round(stats::runif(nrow(m0), 500, 8000))          # gene lengths for TPM/FPKM
  mean_depth <- mean(colSums(m0))
  for (d in depths) {
    if (d >= mean_depth) next
    f <- d / mean_depth
    for (r in seq_len(reps)) {
      m <- matrix(stats::rbinom(length(m0), as.vector(m0), f), nrow(m0),
                  dimnames = dimnames(m0))
      tr <- transforms(m, len)
      for (t in names(expect)) {
        x <- tr[[t]]
        got <- if (is.null(x)) "NOT BUILT" else attest_counts(x)$verdict
        rows[[length(rows) + 1]] <- data.frame(
          dataset = ds, depth = d, rep = r, transform = t,
          expected = expect[[t]], got = got,
          min_reads = min(colSums(m)), stringsAsFactors = FALSE)
      }
    }
  }
}
res <- do.call(rbind, rows)
res$correct <- res$got == res$expected
utils::write.csv(res, "validation/depth_sweep.csv", row.names = FALSE)

cat("\nwrong verdicts, by dataset x depth x transform (of 5 replicates):\n")
wrong <- res[!res$correct, ]
if (!nrow(wrong)) cat("  none\n") else
  print(stats::aggregate(rep ~ dataset + depth + transform + expected + got, wrong, length))

cat("\naccuracy by depth (all datasets and transforms):\n")
print(stats::aggregate(correct ~ depth, res, function(x) sprintf("%d/%d", sum(x), length(x))))
