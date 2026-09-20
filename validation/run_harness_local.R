source("R/check_scale.R")
raw_sets <- readRDS("/mnt/user-data/uploads/Documents/attest/validation/raw_matrices.rds")
names(raw_sets) <- paste0(names(raw_sets), " (", vapply(raw_sets, function(m) paste0(nrow(m), "g x ", ncol(m), "s"), ""), ")")
set.seed(11)
gene_len <- function(m) round(runif(nrow(m), 500, 8000))
variants <- list(
  "raw counts"             = list(f = function(m) m, want="PERMITTED"),
  "CPM"                    = list(f = function(m) t(t(m)/colSums(m))*1e6, want="NOT PERMITTED"),
  "CPM rounded"            = list(f = function(m) round(t(t(m)/colSums(m))*1e6), want="NOT PERMITTED"),
  "TPM rounded"            = list(f = function(m){r<-m/gene_len(m); round(t(t(r)/colSums(r))*1e6)}, want="NOT PERMITTED"),
  "FPKM"                   = list(f = function(m) (m/gene_len(m))/rep(colSums(m)/1e9, each=nrow(m)), want="NOT PERMITTED"),
  "FPKM rounded"           = list(f = function(m) round((m/gene_len(m))/rep(colSums(m)/1e9, each=nrow(m)),2), want="NOT PERMITTED"),
  "size-factor normalised" = list(f = function(m) round(t(t(m)/at_size_factors(m))), want="CAUTION"),
  "log2(CPM+1)"            = list(f = function(m) log2(t(t(m)/colSums(m))*1e6+1), want="NOT PERMITTED"),
  "estimated counts"       = list(f = function(m){e<-m; e[m>0]<-m[m>0]*runif(sum(m>0),.9,1.1); e}, want="CAUTION"),
  "downsampled 1/40"       = list(f = function(m) round(m/40), want="PERMITTED"),
  "downsampled 1/200"      = list(f = function(m) round(m/200), want="PERMITTED"),
  "expressed genes only"   = list(f = function(m) m[rowSums(m)>=10,,drop=FALSE], want="PERMITTED"),
  "protein-coding subset"  = list(f = function(m) m[sample(nrow(m), min(20000, round(nrow(m)*0.4))),,drop=FALSE], want="PERMITTED")
)
rows <- list()
for (ds in names(raw_sets)) for (v in names(variants)) {
  got <- tryCatch(attest_counts(variants[[v]]$f(raw_sets[[ds]]))$verdict,
                  error=function(e) paste("ERROR:", conditionMessage(e)))
  rows[[length(rows)+1]] <- data.frame(dataset=ds, variant=v, expected=variants[[v]]$want, got=got,
                                       agree=identical(got, variants[[v]]$want), stringsAsFactors=FALSE)
}
res <- do.call(rbind, rows)
cat("agreement:", sum(res$agree), "/", nrow(res), "\n\n")
print(res[!res$agree, ], row.names=FALSE)
write.csv(res, "validation/real_data_results.csv", row.names=FALSE)
