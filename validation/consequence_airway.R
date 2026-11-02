# What happens to a DESeq2 analysis when the airway matrix is not raw counts?
# Design ~ cell + dex, adjusted p-value < 0.05 (results(..., alpha = 0.05)).
# Run from the attest project root:  source("validation/consequence_airway.R")
# Needs: DESeq2, airway  (BiocManager::install(c("DESeq2", "airway")))
# Writes validation/consequence_airway.csv and the session info beside it.
suppressPackageStartupMessages({ library(DESeq2); library(airway) })
data("airway", package = "airway")
m  <- assay(airway, "counts")
cd <- as.data.frame(colData(airway))[, c("cell", "dex")]
cd$dex <- relevel(factor(cd$dex), ref = "untrt")

# gene length = total length of the union of the gene's exons
len <- sum(width(GenomicRanges::reduce(rowRanges(airway))))
stopifnot(all(len > 0), length(len) == nrow(m))

dds0 <- estimateSizeFactors(DESeqDataSetFromMatrix(m, cd, ~ cell + dex))
forms <- list(
  "raw counts"                 = m,
  "DESeq2-normalised, rounded" = round(counts(dds0, normalized = TRUE)),
  "rounded CPM"                = round(t(t(m) / colSums(m)) * 1e6),
  "rounded TPM"                = { r <- m / len * 1e3; round(t(t(r) / colSums(r)) * 1e6) }
)
de <- function(x) {
  storage.mode(x) <- "integer"
  dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(x, cd, ~ cell + dex), quiet = TRUE))
  res <- results(dds, alpha = 0.05)
  rownames(res)[which(res$padj < 0.05)]
}
sets <- lapply(forms, de)
out <- data.frame(input = names(sets),
                  de_genes = vapply(sets, length, integer(1)),
                  shared_with_raw = vapply(sets, function(s) length(intersect(s, sets[["raw counts"]])), integer(1)),
                  row.names = NULL)
out$lost_pct <- round(100 * (1 - out$shared_with_raw / out$de_genes[1]), 1)
print(out)
write.csv(out, "validation/consequence_airway.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "validation/consequence_airway_sessionInfo.txt")
cat("written: validation/consequence_airway.csv and _sessionInfo.txt\n")
