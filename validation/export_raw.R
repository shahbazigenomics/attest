# Exports the raw count matrices of the example datasets so they can be
# validated outside this machine. Writes to ~/Documents/attest/validation/.
out <- path.expand("~/Documents/attest/validation")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

sets <- list()
try({ library(airway); data(airway)
      sets[["airway"]] <- as.matrix(SummarizedExperiment::assay(airway)) }, silent = TRUE)
try({ fn <- system.file("extdata", "pasilla_gene_counts.tsv", package = "pasilla")
      sets[["pasilla"]] <- as.matrix(read.table(fn, header = TRUE, row.names = 1)) }, silent = TRUE)
try({ library(fission); data(fission)
      sets[["fission"]] <- as.matrix(SummarizedExperiment::assay(fission)) }, silent = TRUE)
try({ library(parathyroidSE); data(parathyroidGenesSE)
      sets[["parathyroid"]] <- as.matrix(SummarizedExperiment::assay(parathyroidGenesSE)) }, silent = TRUE)

saveRDS(sets, file.path(out, "raw_matrices.rds"))
cat("exported:", paste(names(sets), collapse = ", "), "\n")
cat("file:", file.path(out, "raw_matrices.rds"), "\n")
