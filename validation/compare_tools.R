# What do the standard tools say when handed an invalid matrix?
# Run in the countsworthy project.  Needs: DESeq2, edgeR, limma (RNAseqQC optional).
#   BiocManager::install(c("edgeR","limma","RNAseqQC"))
source("R/check_scale.R")

fx  <- readRDS("inst/extdata/fixtures.rds")
m   <- fx$airway
cd  <- data.frame(grp = factor(rep(c("a","b"), length.out = ncol(m))))
set.seed(1); len <- round(runif(nrow(m), 500, 8000))

mats <- list(
  "raw counts"  = m,
  "CPM rounded" = round(t(t(m) / colSums(m)) * 1e6),
  "TPM rounded" = { r <- m / len; round(t(t(r) / colSums(r)) * 1e6) },
  "FPKM"        = (m / len) / rep(colSums(m) / 1e9, each = nrow(m))
)

# capture everything a call emits, without letting it stop the loop
say <- function(expr) {
  msgs <- character(0)
  res <- withCallingHandlers(
    tryCatch(force(expr), error = function(e) { msgs <<- c(msgs, paste("ERROR:", conditionMessage(e))); NULL }),
    warning = function(w) { msgs <<- c(msgs, paste("WARNING:", conditionMessage(w))); invokeRestart("muffleWarning") },
    message = function(m2) { msgs <<- c(msgs, paste("MESSAGE:", trimws(conditionMessage(m2)))); invokeRestart("muffleMessage") })
  if (!length(msgs)) "(silent - accepted the input)" else paste(unique(msgs), collapse = " | ")
}

tools <- list(
  "DESeq2 (DESeqDataSetFromMatrix)" = function(x)
    DESeq2::DESeqDataSetFromMatrix(x, cd, ~ grp),
  "edgeR (DGEList + calcNormFactors)" = function(x)
    edgeR::calcNormFactors(edgeR::DGEList(counts = x, group = cd$grp)),
  "limma-voom" = function(x)
    limma::voom(x, model.matrix(~ grp, cd)),
  # RNAseqQC::make_dds omitted: it prompts to create an AnnotationHub cache,
  # which blocks a scripted run and is unrelated to input validation.
  "countsworthy" = function(x) {
    v <- cw_check_scale(x); stop(paste0(v$verdict, ": ", v$reasons[1]))
  }
)

rows <- list()
for (tn in names(tools)) for (mn in names(mats)) {
  out <- say(tools[[tn]](mats[[mn]]))
  rows[[length(rows) + 1]] <- data.frame(tool = tn, input = mn,
                                         response = substr(out, 1, 160),
                                         stringsAsFactors = FALSE)
}
res <- do.call(rbind, rows)
for (tn in unique(res$tool)) {
  cat("\n=====", tn, "=====\n")
  sub <- res[res$tool == tn, ]
  for (i in seq_len(nrow(sub))) cat(sprintf("  %-12s %s\n", sub$input[i], sub$response[i]))
}
write.csv(res, "validation/tool_comparison.csv", row.names = FALSE)
cat("\nwritten: validation/tool_comparison.csv\n")
