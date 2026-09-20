# What do the standard tools say about a bad DESIGN?
# Run in the attest project. Needs DESeq2, edgeR, limma; sva optional.
for (f in list.files("R", full.names = TRUE)) source(f)
suppressPackageStartupMessages({library(DESeq2); library(edgeR); library(limma)})

fx <- readRDS("inst/extdata/fixtures.rds")
m  <- fx$airway                      # 8 samples

designs <- list(
  "clean (batch balanced)" = data.frame(
    cond  = factor(rep(c("ctrl","trt"), each = 4)),
    batch = factor(rep(c("A","B"), times = 4))),
  "partly confounded (3 of 4)" = data.frame(
    cond  = factor(rep(c("ctrl","trt"), each = 4)),
    batch = factor(c("A","A","A","B", "B","B","B","A"))),
  "completely confounded" = data.frame(
    cond  = factor(rep(c("ctrl","trt"), each = 4)),
    batch = factor(rep(c("A","B"), each = 4))),
  "n = 2 per group" = data.frame(
    cond  = factor(rep(c("ctrl","trt"), each = 4)),
    batch = factor(rep(c("A","B"), times = 4)))     # only first 4 samples used below
)

say <- function(expr) {
  msgs <- character(0)
  withCallingHandlers(
    tryCatch(force(expr), error = function(e) { msgs <<- c(msgs, paste("ERROR:", conditionMessage(e))); NULL }),
    warning = function(w) { msgs <<- c(msgs, paste("WARNING:", conditionMessage(w))); invokeRestart("muffleWarning") },
    message = function(x) { msgs <<- c(msgs, paste("MESSAGE:", trimws(conditionMessage(x)))); invokeRestart("muffleMessage") })
  if (!length(msgs)) "(silent - accepted the design)" else paste(unique(msgs), collapse = " | ")
}

for (nm in names(designs)) {
  cd <- designs[[nm]]
  mm <- m
  if (nm == "n = 2 per group") { keep <- c(1, 2, 5, 6); mm <- m[, keep]; cd <- cd[keep, ] }
  cat("\n=====", nm, "=====\n")
  cat(sprintf("  %-22s %s\n", "DESeq2 (~batch+cond)",
      say({ d <- DESeqDataSetFromMatrix(mm, cd, ~ batch + cond); DESeq(d, quiet = TRUE) })))
  cat(sprintf("  %-22s %s\n", "edgeR (glmQLFit)",
      say({ y <- calcNormFactors(DGEList(mm, group = cd$cond))
            des <- model.matrix(~ batch + cond, cd)
            glmQLFit(estimateDisp(y, des), des) })))
  cat(sprintf("  %-22s %s\n", "limma-voom",
      say({ des <- model.matrix(~ batch + cond, cd); eBayes(lmFit(voom(mm, des), des)) })))
  if (requireNamespace("sva", quietly = TRUE)) {
    cat(sprintf("  %-22s %s\n", "ComBat-seq",
        say(sva::ComBat_seq(as.matrix(mm), batch = cd$batch, group = cd$cond))))
  }
}
cat("\nNote: 'n = 2 per group' uses 4 samples, so it also asks what each tool says about replication.\n")
