# Does the claimed minimum detectable fold change actually hold?
# For each gene we take attest's claimed MDFC, simulate data with exactly that
# effect, test it, and count how often it is called. Claimed power is 80%.
for (f in list.files("R", full.names = TRUE)) source(f)
suppressPackageStartupMessages(library(MASS))
set.seed(11)

test_gene <- function(y, g) {                       # NB GLM Wald test, one gene
  fit <- try(suppressWarnings(MASS::glm.nb(y ~ g)), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  cf <- summary(fit)$coefficients
  if (nrow(cf) < 2) return(NA_real_)
  cf[2, 4]
}

calibrate <- function(n_per_group = 4, n_genes = 300, reps = 1, mu_range = c(5, 5000)) {
  mu  <- exp(runif(n_genes, log(mu_range[1]), log(mu_range[2])))
  phi <- pmax(0.02, 0.4 / sqrt(mu) + rnorm(n_genes, 0.05, 0.02))
  g   <- factor(rep(c("a", "b"), each = n_per_group))

  # a pilot dataset under the null, from which attest estimates detectability
  pilot <- t(vapply(seq_len(n_genes), function(i)
    rnbinom(2 * n_per_group, mu = mu[i], size = 1/phi[i]), numeric(2 * n_per_group)))
  rownames(pilot) <- paste0("g", seq_len(n_genes))
  claimed <- attest_detectability(pilot, group = g)$measurements$min_detectable_fc

  # now simulate each gene with exactly its claimed effect and test it
  hits <- numeric(0)
  for (r in seq_len(reps)) {
    y <- t(vapply(seq_len(n_genes), function(i) {
      fc <- claimed[i]
      if (!is.finite(fc)) return(rep(NA_real_, 2 * n_per_group))
      c(rnbinom(n_per_group, mu = mu[i],      size = 1/phi[i]),
        rnbinom(n_per_group, mu = mu[i] * fc, size = 1/phi[i]))
    }, numeric(2 * n_per_group)))
    a <- attest_detectability(pilot, group = g)$measurements$per_gene_alpha
    p <- vapply(seq_len(n_genes), function(i)
      if (anyNA(y[i, ])) NA_real_ else test_gene(y[i, ], g), numeric(1))
    hits <- c(hits, mean(p < a, na.rm = TRUE))
  }
  data.frame(n_per_group = n_per_group,
             median_claimed_fc = round(median(claimed[is.finite(claimed)]), 2),
             realised_power = round(mean(hits), 3))
}

res <- do.call(rbind, lapply(c(3, 4, 6, 10), calibrate))
print(res, row.names = FALSE)
cat("\nclaimed power was 0.80; realised power is the fraction actually detected\n")
write.csv(res, "validation/detectability_calibration.csv", row.names = FALSE)
