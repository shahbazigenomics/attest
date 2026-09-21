# A mock of GEO and NCBI built from real counts, for testing the pipeline
# offline. Returns the URL -> file map for net$map. Every series carries a
# planted, known situation; run_mock.R checks each comes out as planted.
#
# NCBI's counts are simulated as a *different pipeline* over the same reads:
# each gene gets its own fixed efficiency (log-normal, sd 0.25) and the counts
# are redrawn - so NCBI and the author never agree exactly, as in real data.

build_mock <- function(dir, raw_rds = "validation/raw_matrices.rds",
                       kang_rds = "validation/umi_pseudobulk.rds") {
  set.seed(1)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  map <- list(); put <- function(url, file) map[[url]] <<- file
  P <- function(...) file.path(dir, ...)
  gzw <- function(lines, path) { con <- gzfile(path, "w"); writeLines(lines, con); close(con); path }
  tsvgz <- function(df, path) { con <- gzfile(path, "w"); utils::write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con); path }
  ncbi <- "https://www.ncbi.nlm.nih.gov"

  airway <- readRDS(raw_rds)$airway
  kang <- readRDS(kang_rds)
  markers <- c(ENSG00000229807 = "XIST", ENSG00000129824 = "RPS4Y1", ENSG00000012817 = "KDM5D",
               ENSG00000067048 = "DDX3Y", ENSG00000183878 = "UTY", ENSG00000114374 = "USP9Y",
               ENSG00000198692 = "EIF1AY", ENSG00000165246 = "NLGN4Y")

  # --- one annotation for both gene universes -----------------------------------
  ens <- rownames(airway)
  sym_a <- ifelse(ens %in% names(markers), markers[ens], paste0("AW", seq_along(ens)))
  ksym <- setdiff(rownames(kang$counts), markers)
  annot <- data.frame(GeneID = c(100000 + seq_along(ens), 500000 + seq_along(ksym)),
                      Symbol = c(sym_a, ksym),
                      Description = "mock", EnsemblGeneID = c(ens, rep("", length(ksym))),
                      stringsAsFactors = FALSE)
  tsvgz(annot, P("annot.tsv.gz"))
  put(paste0(ncbi, "/geo/download/?type=rnaseq_counts&format=file&file=Human.GRCh38.p13.annot.tsv.gz"), P("annot.tsv.gz"))
  gid_of_sym <- setNames(as.character(annot$GeneID), annot$Symbol)

  # NCBI's raw counts: same reads, another pipeline
  eff_a <- exp(stats::rnorm(nrow(airway), 0, 0.25))
  ncbi_from <- function(m, gids, eff) {
    x <- matrix(stats::rpois(length(m), m * eff * 0.9), nrow(m), dimnames = list(gids, colnames(m)))
    x
  }

  sex_airway <- ifelse(at_sex_from_expression(airway)$call == "female", "female", "male")
  esearch_ids <- character(0)

  series <- function(gse, gsm, titles, author, ncbi_counts = TRUE, R = NULL, chars = list(),
                     extra_files = character(0), matrix_files = TRUE) {
    num <- as.numeric(sub("GSE", "", gse)); esearch_ids <<- c(esearch_ids, sprintf("2%08d", num))
    d <- P(gse); dir.create(d, showWarnings = FALSE)
    stub <- paste0("GSE", substr(num, 1, nchar(num) - 3), "nnn")
    ftp <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/", stub, gse)
    # download page
    page <- if (ncbi_counts) c("<html><body>",
      sprintf('<a href="/geo/download/?type=rnaseq_counts&amp;acc=%s&amp;format=file&amp;file=%s_raw_counts_GRCh38.p13_NCBI.tsv.gz">raw counts</a>', gse, gse),
      '<a href="/geo/download/?type=rnaseq_counts&amp;format=file&amp;file=Human.GRCh38.p13.annot.tsv.gz">annotation</a>',
      "</body></html>") else sprintf("<html><body>Download %s: %s_RAW.tar (http)</body></html>", gse, gse)
    writeLines(page, file.path(d, "page.html")); put(sprintf("%s/geo/download/?acc=%s", ncbi, gse), file.path(d, "page.html"))
    if (ncbi_counts) {
      colnames(R) <- gsm
      tsvgz(data.frame(GeneID = rownames(R), R, check.names = FALSE), file.path(d, "ncbi.tsv.gz"))
      put(sprintf("%s/geo/download/?type=rnaseq_counts&acc=%s&format=file&file=%s_raw_counts_GRCh38.p13_NCBI.tsv.gz", ncbi, gse, gse),
          file.path(d, "ncbi.tsv.gz"))
    }
    # series matrix
    q <- function(x) paste0('"', x, '"', collapse = "\t")
    sm <- c(sprintf('!Series_title\t"%s"', gse), paste0("!Sample_title\t", q(titles)),
            paste0("!Sample_geo_accession\t", q(gsm)), paste0("!Sample_source_name_ch1\t", q(rep("cells", length(gsm)))))
    for (k in names(chars)) sm <- c(sm, paste0("!Sample_characteristics_ch1\t", q(paste0(k, ": ", chars[[k]]))))
    sm <- c(sm, "!series_matrix_table_begin", "!series_matrix_table_end")
    mf <- paste0(gse, "_series_matrix.txt.gz")
    gzw(sm, file.path(d, mf))
    writeLines(sprintf('<a href="/geo/series/%s/">Parent Directory</a>\n<a href="%s">%s</a>', stub, mf, mf), file.path(d, "matrix_index.html"))
    put(paste0(ftp, "matrix/"), file.path(d, "matrix_index.html")); put(paste0(ftp, "matrix/", mf), file.path(d, mf))
    # supplementary files
    files <- c(names(author), extra_files)
    for (f in names(author)) put(paste0(ftp, "suppl/", f), author[[f]])
    for (f in extra_files) { writeLines("x", file.path(d, f)); put(paste0(ftp, "suppl/", f), file.path(d, f)) }
    writeLines(c(sprintf('<a href="/geo/series/%s/">Parent Directory</a>', stub),
                 sprintf('<a href="%s">%s</a>', files, files)), file.path(d, "suppl_index.html"))
    put(paste0(ftp, "suppl/"), file.path(d, "suppl_index.html"))
  }

  gsm_a <- function(k) sprintf("GSM%d%02d", k, 1:8)
  titles_a <- c("N61311_untreated", "N61311_dex", "N052611_untreated", "N052611_dex",
                "N080611_untreated", "N080611_dex", "N061011_untreated", "N061011_dex")
  cond_a <- rep(c("untreated", "dexamethasone"), 4)
  cell_a <- rep(c("N61311", "N052611", "N080611", "N061011"), each = 2)
  gid_a <- as.character(100000 + seq_along(ens))
  R_a <- ncbi_from(airway, gid_a, eff_a)
  author_tsv <- function(m, path, idname = "gene_id") tsvgz(data.frame(setNames(list(rownames(m)), idname), m, check.names = FALSE), path)
  cpm <- function(m) t(t(m) / colSums(m)) * 1e6

  # 1. raw counts, columns are GSM accessions -> raw; PERMITTED
  m <- airway; colnames(m) <- gsm_a(1)
  series("GSE900001", gsm_a(1), titles_a, list(GSE900001_raw_counts.txt.gz = author_tsv(m, P("a1.txt.gz"))),
         R = R_a, chars = list(`cell line` = cell_a, treatment = cond_a, Sex = sex_airway))
  # 2. rounded CPM *named* raw counts, columns are titles -> mislabel; NOT PERMITTED; consequence
  m <- round(cpm(airway)); colnames(m) <- titles_a
  series("GSE900002", gsm_a(2), titles_a, list(GSE900002_raw_counts.txt.gz = author_tsv(m, P("a2.txt.gz"))),
         R = R_a, chars = list(`cell line` = cell_a, treatment = cond_a))
  # 3. TPM (labelled) and raw counts; titles with other punctuation
  len <- round(stats::runif(nrow(airway), 500, 8000)); rpk <- airway / len
  tpm <- t(t(rpk) / colSums(rpk)) * 1e6
  m1 <- airway; colnames(m1) <- gsub("_", "-", titles_a); m2 <- tpm; colnames(m2) <- gsub("_", " ", titles_a)
  series("GSE900003", gsm_a(3), titles_a,
         list(GSE900003_counts.txt.gz = author_tsv(m1, P("a3a.txt.gz")), GSE900003_TPM.txt.gz = author_tsv(m2, P("a3b.txt.gz"))),
         R = R_a, chars = list(treatment = cond_a))
  # 4. filtered raw counts; one sample's sex label flipped
  m <- airway[rowSums(airway) >= 10, ]; colnames(m) <- gsm_a(4)
  sx <- sex_airway; sx[1] <- if (sx[1] == "male") "female" else "male"
  series("GSE900004", gsm_a(4), titles_a, list(GSE900004_featurecounts.tsv.gz = author_tsv(m, P("a4.txt.gz"))),
         R = R_a, chars = list(sex = sx, treatment = cond_a))
  # 5. htseq output with summary rows; column names match nothing -> correlation matching
  m <- rbind(airway, `__no_feature` = round(colSums(airway) * 0.12), `__ambiguous` = round(colSums(airway) * 0.03),
             `__too_low_aQual` = 0, `__not_aligned` = 0, `__alignment_not_unique` = round(colSums(airway) * 0.05))
  colnames(m) <- paste0("lib", sample(8))
  series("GSE900005", gsm_a(5), titles_a, list(GSE900005_htseq_counts.txt.gz = author_tsv(m, P("a5.txt.gz"))), R = R_a)
  # 6. Kang pseudobulk, symbols, as an Excel file with three names turned into dates
  kc <- kang$counts; kg <- rownames(kc) %in% names(gid_of_sym)
  kc <- kc[kg, ]; gsm_k <- sprintf("GSM6%03d", 1:16)
  R_k <- ncbi_from(kc, gid_of_sym[rownames(kc)], exp(stats::rnorm(nrow(kc), 0, 0.25)))
  ex <- data.frame(Gene = rownames(kc), kc, check.names = FALSE); colnames(ex)[-1] <- gsm_k
  hit <- which(rowSums(kc) > 0)[1:3]; ex$Gene[hit] <- c("1-Sep", "2-Sep", "1-Mar")
  writexl::write_xlsx(ex, P("a6.xlsx"))
  series("GSE900006", gsm_k, colnames(kang$counts), list(GSE900006_gene_counts.xlsx = P("a6.xlsx")),
         R = R_k, chars = list(stimulation = as.character(kang$sheet$cond), sex = rep("female", 16)))
  # 7. only raw tar and bigwigs
  series("GSE900007", gsm_a(7), titles_a, list(), R = R_a, extra_files = c("GSE900007_RAW.tar", "GSE900007_coverage.bw"))
  # 8. single-cell only
  series("GSE900008", gsm_a(8), titles_a, list(), R = R_a,
         extra_files = c("GSE900008_matrix.mtx.gz", "GSE900008_barcodes.tsv.gz", "GSE900008_features.tsv.gz"))
  # 9. no NCBI counts
  m <- airway; colnames(m) <- gsm_a(9)
  series("GSE900009", gsm_a(9), titles_a, list(GSE900009_counts.txt.gz = author_tsv(m, P("a9.txt.gz"))), ncbi_counts = FALSE)
  # 10. DESeq2-normalised and rounded, unhelpfully named -> depth removed; attest CAUTION
  m <- round(t(t(airway) / at_size_factors(airway))); colnames(m) <- gsm_a(10)
  series("GSE900010", gsm_a(10), titles_a, list(GSE900010_data.txt.gz = author_tsv(m, P("a10.txt.gz"))),
         R = R_a, chars = list(treatment = cond_a))
  # 11. log2(CPM + 1) named expression
  m <- log2(cpm(airway) + 1); colnames(m) <- gsm_a(11)
  series("GSE900011", gsm_a(11), titles_a, list(GSE900011_expression.txt.gz = author_tsv(m, P("a11.txt.gz"))), R = R_a)
  # 12. equal depths: truth cannot be decided
  eq <- matrix(stats::rbinom(length(airway), airway, rep(0.95 * min(colSums(airway)) / colSums(airway), each = nrow(airway))),
               nrow(airway), dimnames = dimnames(airway))
  m <- eq; colnames(m) <- gsm_a(12)
  series("GSE900012", gsm_a(12), titles_a, list(GSE900012_counts.txt.gz = author_tsv(m, P("a12.txt.gz"))),
         R = ncbi_from(eq, gid_a, eff_a))

  # esearch: count query and paged query
  writeLines(sprintf("<eSearchResult><Count>%d</Count><RetMax>0</RetMax><IdList></IdList></eSearchResult>", length(esearch_ids)), P("es_count.xml"))
  writeLines(paste0("<eSearchResult><Count>", length(esearch_ids), "</Count><IdList>",
                    paste0("<Id>", esearch_ids, "</Id>", collapse = ""), "</IdList>",
                    "<QueryTranslation>\"rnaseq counts\"[Filter] AND \"Homo sapiens\"[Organism]</QueryTranslation></eSearchResult>"), P("es_ids.xml"))
  # the same search without the counts filter: more series
  writeLines("<eSearchResult><Count>5000</Count><IdList></IdList></eSearchResult>", P("es_all.xml"))
  map <- c(list("re:esearch\\.fcgi\\?db=gds&term=%22Homo" = P("es_all.xml")), map)
  # NCBI Gene's table, for the annotation fallback
  gi <- data.frame("#tax_id" = 9606, GeneID = annot$GeneID, Symbol = annot$Symbol, LocusTag = "-",
                   Synonyms = "-", dbXrefs = ifelse(is.na(annot$Ensembl), "-", paste0("MIM:1|HGNC:HGNC:1|Ensembl:", annot$Ensembl)),
                   check.names = FALSE)
  tsvgz(gi, P("gene_info.gz"))
  map[["https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"]] <- P("gene_info.gz")
  # what NCBI sometimes sends instead of a file: a bot-check page
  writeLines("<!doctype html><html><head><base href=\"https://www.google.com/recaptcha/challengepage/\"></head></html>", P("captcha.html"))
  map[["re:esearch\\.fcgi.*retmax=0&"]] <- P("es_count.xml")
  map[["re:esearch\\.fcgi.*retmax=(5|200|10000)"]] <- P("es_ids.xml")
  map
}
