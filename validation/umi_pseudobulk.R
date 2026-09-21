# Real UMI data: Kang et al. 2018, 10x Chromium PBMCs from 8 lupus patients,
# each sampled with and without interferon-beta (muscData::Kang18_8vs8).
# Summed per donor x condition it gives 16 pseudobulk samples - the way
# single-cell data reaches DESeq2/edgeR, and so the way it reaches attest.
#
# Run on a machine with Bioconductor, from the attest project:
#   source("validation/umi_pseudobulk.R")
# First run downloads the dataset through ExperimentHub (about 100 MB, cached).
# Writes validation/umi_pseudobulk.rds (the 16-sample matrix and its sample
# sheet, a few MB) and prints every verdict.

for (p in c("muscData", "SingleCellExperiment", "ExperimentHub"))
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
for (f in list.files("R", full.names = TRUE)) source(f)

sce <- muscData::Kang18_8vs8()
cd  <- as.data.frame(SummarizedExperiment::colData(sce))
cat("cells:", ncol(sce), " genes:", nrow(sce), "\n")
cat("sample-sheet columns:", paste(names(cd), collapse = ", "), "\n")

pick <- function(cands) { hit <- intersect(cands, names(cd)); if (!length(hit)) stop("none of ", paste(cands, collapse = "/"), " in colData"); hit[1] }
donor_col <- pick(c("ind", "sample_id", "donor", "patient"))
cond_col  <- pick(c("stim", "group_id", "condition", "treatment"))
keep <- rep(TRUE, ncol(sce))
if ("multiplets" %in% names(cd)) keep <- cd$multiplets %in% c("singlet", "Singlet")
cat("cells kept (singlets):", sum(keep), "of", length(keep), "\n")

donor <- factor(cd[[donor_col]][keep])
cond  <- factor(cd[[cond_col]][keep])
grp   <- interaction(donor, cond, sep = "_", drop = TRUE)

counts <- SummarizedExperiment::assay(sce, "counts")[, keep]
mm <- Matrix::sparse.model.matrix(~ 0 + grp)
colnames(mm) <- levels(grp)
pb <- as.matrix(counts %*% mm)                       # genes x (donor x condition)
storage.mode(pb) <- "integer"

sheet <- data.frame(
  donor = factor(sub("_.*$", "", colnames(pb))),
  cond  = factor(sub("^[^_]*_", "", colnames(pb))),
  cells = as.integer(table(grp)[colnames(pb)]),
  row.names = colnames(pb))
cat("\npseudobulk:", nrow(pb), "genes x", ncol(pb), "samples;",
    "UMIs per sample", format(min(colSums(pb)), big.mark = ","), "-",
    format(max(colSums(pb)), big.mark = ","), "\n")
print(table(sheet$donor, sheet$cond))

saveRDS(list(counts = pb, sheet = sheet, source = "muscData::Kang18_8vs8, singlets, summed per donor x condition"),
        "validation/umi_pseudobulk.rds", compress = "xz")
cat("\nsaved validation/umi_pseudobulk.rds (",
    round(file.size("validation/umi_pseudobulk.rds") / 1e6, 1), "MB )\n", sep = "")

cat("\n==================== full report, raw pseudobulk ====================\n")
print(attest(pb, metadata = sheet, design = ~ donor + cond))

cat("\n==================== transforms: what attest_counts() says ====================\n")
sf  <- at_size_factors(pb)
cpm <- t(t(pb) / colSums(pb)) * 1e6
forms <- list(
  "raw pseudobulk (expect PERMITTED)"       = pb,
  "rounded normalised (expect CAUTION)"     = round(t(t(pb) / sf)),
  "CPM (expect NOT PERMITTED)"              = cpm,
  "rounded CPM (expect NOT PERMITTED)"      = round(cpm),
  "log2(CPM + 1) (expect NOT PERMITTED)"    = log2(cpm + 1))
for (nm in names(forms))
  cat(sprintf("%-40s -> %s\n", nm, attest_counts(forms[[nm]])$verdict))
