# Runs every stage of the study against the mock and checks each planted case.
# From the attest project root:  Rscript validation/geo_study/mock/run_mock.R

for (f in list.files("R", full.names = TRUE)) source(f)
for (f in c("config.R", "lib_net.R", "lib_parse.R", "lib_truth.R", "stages.R", "summary.R"))
  source(file.path("validation/geo_study", f))
source("validation/geo_study/mock/build_mock.R")

tmp <- file.path(Sys.getenv("GEO_MOCK_DIR", tempdir()), "geo_mock"); unlink(tmp, recursive = TRUE)
study$cache <- file.path(tmp, "cache"); study$results <- file.path(tmp, "results")
study$sleep <- 0; study$n_candidates <- 50; study$target_usable <- 50
dir.create(study$results, recursive = TRUE); dir.create(study$cache, recursive = TRUE)
net$mode <- "mock"; net$map <- build_mock(file.path(tmp, "web")); net$annot <- NULL

stopifnot(isTRUE(stage_check()))
stage_run(); stage_consequence(); stage_summary()

res <- setNames(lapply(list.files(file.path(study$results, "series"), full.names = TRUE), readRDS),
                sub("\\.rds$", "", list.files(file.path(study$results, "series"))))
f <- rbind_fill(lapply(res, function(r) r$files))
row <- function(gse, file = NULL) { x <- f[f$gse == gse, ]; if (!is.null(file)) x <- x[x$file == file, ]; x[1, ] }
fails <- 0
ok <- function(cond, what) { cat(if (isTRUE(cond)) "ok  " else { fails <<- fails + 1; "FAIL" }, what, "\n") }

ok(res$GSE900001$status == "ok" && row("GSE900001")$v_value == "PERMITTED" && row("GSE900001")$truth_value == "raw counts" &&
   row("GSE900001")$match_how == "accession", "1 raw counts: attest PERMITTED, truth raw, matched by accession")
ok(row("GSE900002")$v_value == "NOT PERMITTED" && row("GSE900002")$truth_value == "depth removed" &&
   row("GSE900002")$labelled_as == "counts" && row("GSE900002")$match_how == "title", "2 CPM named raw_counts: caught, truth depth removed, mislabel, matched by title")
ok(row("GSE900003", "GSE900003_TPM.txt.gz")$v_value == "NOT PERMITTED" && row("GSE900003", "GSE900003_TPM.txt.gz")$truth_value == "depth removed" &&
   row("GSE900003", "GSE900003_counts.txt.gz")$truth_value == "raw counts", "3 TPM and counts in one series, titles with other punctuation")
ok(row("GSE900004")$v_completeness == "CAUTION" && row("GSE900004")$truth_completeness == "filtered", "4 filtered: attest CAUTION, truth filtered")
ok(sum(res$GSE900004$sex$mismatch) == 1, "4 the flipped sex label is found")
ok(row("GSE900001")$truth_completeness == "complete", "1 complete matrix: truth complete")
ok(row("GSE900005")$summary_rows == 5 && row("GSE900005")$v_identifiers == "NOT PERMITTED" && row("GSE900005")$match_how == "correlation" &&
   row("GSE900005")$n_matched == 8, "5 htseq summary rows found; samples matched by correlation, all 8")
ok(isTRUE(row("GSE900006")$from_excel) && row("GSE900006")$excel_dates == 3 && row("GSE900006")$truth_value == "raw counts" &&
   row("GSE900006")$id_system == "symbol or systematic", "6 Excel file: converted, 3 dates found, symbols mapped, truth raw")
ok(sum(res$GSE900006$sex$mismatch) == 0 && nrow(res$GSE900006$sex) == 16, "6 all 16 female, no mismatch")
ok(res$GSE900007$status == "per-sample files only (RAW.tar)", "7 raw tar, filelist.txt and bigwig only: excluded as per-sample files only")
ok(res$GSE900008$status == "single-cell only", "8 single-cell only: excluded")
ok(res$GSE900009$status == "no NCBI counts", "9 no NCBI counts: excluded")
ok(row("GSE900010")$truth_value == "depth removed" && row("GSE900010")$v_value == "CAUTION", "10 DESeq2-normalised and rounded: truth depth removed (slope beats whole numbers), attest CAUTION")
ok(row("GSE900011")$truth_value == "log-transformed" && row("GSE900011")$v_value == "NOT PERMITTED", "11 log2 CPM: truth log, attest NOT PERMITTED")
ok(row("GSE900012")$truth_value == "count scale, not whole numbers", "12 normalised, equal depths: count scale, which kind not guessed")
ok(row("GSE900013")$match_how == "title" && row("GSE900013")$n_matched == 4 && row("GSE900013")$truth_value == "raw counts",
   "13 NCBI has 4 of 8 samples: those 4 compared, truth raw counts")
ok(row("GSE900014")$match_how == "correlation" && row("GSE900014")$n_matched == 8 && row("GSE900014")$truth_value == "depth removed" &&
   row("GSE900014")$v_value == "NOT PERMITTED", "14 FPKM, unknown column names: matched by gene-centred correlation, truth depth removed, caught")
ok(row("GSE900015")$v_input == "NOT PERMITTED" && row("GSE900015")$truth_value == "raw counts", "15 spreadsheet totals row: attest NOT PERMITTED (input)")

# NCBI accepts the counts filter but ignores it (as on the first real run):
# the frame falls back to all human expression-by-sequencing series
fr <- resolve_frame(); ok(fr$term == study$term, "counts filter applied by NCBI: frame is the filtered search")
k <- "re:esearch\\.fcgi\\?db=gds&term=%22Homo"; keep <- net$map[[k]]
net$map[[k]] <- sub("es_all", "es_ids", keep)
fr <- resolve_frame(); ok(fr$term == study$term_frame && grepl("per series", fr$basis), "counts filter ignored by NCBI: frame falls back, counts checked per series")
net$map[[k]] <- keep

# NCBI answers a series' download page with a bot check: not "no NCBI counts",
# and not saved, so the next run tries it again
pk <- sprintf("https://www.ncbi.nlm.nih.gov/geo/download/?acc=%s", "GSE900001"); keep <- net$map[[pk]]
net$map[[pk]] <- file.path(dirname(net$map[[grep("annot\\.tsv\\.gz$", names(net$map), value = TRUE)[1]]]), "captcha.html")
ok(startsWith(process_series("GSE900001")$status, "unavailable"), "bot check on the download page: 'unavailable' (retried), not 'no NCBI counts'")
net$map[[pk]] <- keep

# NCBI answers the annotation request with a bot-check page, already cached
# from an earlier run: the page is not used, and NCBI Gene's table stands in
good <- net$annot; net$annot <- NULL
real <- net$map[[grep("annot\\.tsv\\.gz$", names(net$map), value = TRUE)[1]]]
for (k in grep("annot\\.tsv\\.gz$", names(net$map), value = TRUE)) net$map[[k]] <- file.path(dirname(real), "captcha.html")
invisible(file.copy(file.path(dirname(real), "captcha.html"), file.path(study$cache, "ncbi_annotation.tsv.gz"), overwrite = TRUE))
fb <- get_annot(ncbi_counts_urls("GSE900001")$annot)
ok(!file.exists(file.path(study$cache, "ncbi_annotation.tsv.gz")) && identical(fb$GeneID, good$GeneID) &&
   identical(fb$Symbol, good$Symbol) && identical(ifelse(is.na(fb$Ensembl), NA, fb$Ensembl), ifelse(is.na(good$Ensembl) | good$Ensembl == "-", NA, good$Ensembl)),
   "bot-check page instead of the annotation: not cached, NCBI Gene table used, same genes")
ok(looks_like_page(file.path(dirname(real), "captcha.html")) && !looks_like_page(real), "a web page is told apart from a data file")
net$annot <- good

cons <- utils::read.csv(file.path(study$results, "consequence.csv"))
c2 <- cons[cons$gse == "GSE900002", ]
ok(nrow(c2) == 1 && c2$de_ncbi_raw > c2$de_author_file, "consequence: CPM file finds fewer DE genes than NCBI raw counts")
cat(sprintf("     GSE900002: NCBI raw %d DE genes, author CPM file %d, both %d\n", c2$de_ncbi_raw, c2$de_author_file, c2$de_both))
sm <- readLines(file.path(study$results, "summary.md"))
ok(any(grepl("Sensitivity", sm)) && any(grepl("Specificity", sm)) && any(grepl("Named as counts but not counts", sm)), "summary.md written with the key tables")
cat("\n", if (fails) paste(fails, "FAILED") else "ALL PLANTED CASES CORRECT", "\n")
cat("\n---- summary.md ----\n"); cat(sm, sep = "\n")
