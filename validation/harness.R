# Real-data validation harness for attest_counts().
# Run in the attest project. Needs the data packages below:
#   BiocManager::install(c("airway","pasilla","fission","parathyroidSE"))
# Writes validation/real_data_results.csv

source("R/check_scale.R")

# ---- collect raw count matrices -------------------------------------------
raw_sets <- list()

add <- function(name, expr) {
  m <- tryCatch(suppressMessages(expr), error = function(e) NULL)
  if (!is.null(m) && is.matrix(m) && nrow(m) > 100) raw_sets[[name]] <<- m
  else message("skipped: ", name)
}

add("airway (human, polyA)", {
  library(airway); data(airway)
  as.matrix(SummarizedExperiment::assay(airway))
})
add("pasilla (fly)", {
  fn <- system.file("extdata", "pasilla_gene_counts.tsv", package = "pasilla")
  as.matrix(read.table(fn, header = TRUE, row.names = 1))
})
add("fission (yeast)", {
  library(fission); data(fission)
  as.matrix(SummarizedExperiment::assay(fission))
})
add("parathyroid (human)", {
  library(parathyroidSE); data(parathyroidGenesSE)
  as.matrix(SummarizedExperiment::assay(parathyroidGenesSE))
})

message("datasets loaded: ", paste(names(raw_sets), collapse = ", "))

# ---- transformations, with the verdict each one should get -----------------
gene_len <- function(m) round(runif(nrow(m), 500, 8000))   # stand-in lengths

variants <- list(
  "raw counts"            = list(f = function(m) m,                                          want = "PERMITTED"),
  "CPM"                   = list(f = function(m) t(t(m) / colSums(m)) * 1e6,                  want = "NOT PERMITTED"),
  "CPM rounded"           = list(f = function(m) round(t(t(m) / colSums(m)) * 1e6),           want = "NOT PERMITTED"),
  "TPM rounded"           = list(f = function(m) { r <- m / gene_len(m); round(t(t(r) / colSums(r)) * 1e6) },
                                 want = "NOT PERMITTED"),
  "FPKM"                  = list(f = function(m) (m / gene_len(m)) / rep(colSums(m) / 1e9, each = nrow(m)),
                                 want = "NOT PERMITTED"),
  "size-factor normalised"= list(f = function(m) round(t(t(m) / at_size_factors(m))),         want = "CAUTION"),
  "log2(CPM+1)"           = list(f = function(m) log2(t(t(m) / colSums(m)) * 1e6 + 1),        want = "NOT PERMITTED"),
  # salmon/kallisto estimates: fractional where a gene has reads, exactly 0 where it has none.
  # (An earlier fixture added noise to the zeros too, which made it look like FPKM.)
  "estimated counts"      = list(f = function(m) { e <- m; e[m > 0] <- m[m > 0] * runif(sum(m > 0), 0.9, 1.1); e },
                                 want = "CAUTION"),
  "downsampled 1/40"      = list(f = function(m) round(m / 40),                               want = "PERMITTED"),
  "expressed genes only"  = list(f = function(m) m[rowSums(m) >= 10, , drop = FALSE],         want = "PERMITTED")
)

rows <- list()
for (ds in names(raw_sets)) {
  m0 <- raw_sets[[ds]]
  for (v in names(variants)) {
    got <- tryCatch(attest_counts(variants[[v]]$f(m0))$verdict,
                    error = function(e) paste("ERROR:", conditionMessage(e)))
    rows[[length(rows) + 1]] <- data.frame(
      dataset = ds, variant = v,
      expected = variants[[v]]$want, got = got,
      agree = identical(got, variants[[v]]$want),
      stringsAsFactors = FALSE)
  }
}

res <- do.call(rbind, rows)
print(res, row.names = FALSE)
cat("\nagreement:", sum(res$agree), "/", nrow(res), "\n")
cat("\nmismatches:\n"); print(res[!res$agree, ], row.names = FALSE)

dir.create("validation", showWarnings = FALSE)
write.csv(res, "validation/real_data_results.csv", row.names = FALSE)
cat("\nwritten: validation/real_data_results.csv\n")
