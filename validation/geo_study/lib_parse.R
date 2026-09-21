`%||%` <- function(a, b) if (is.null(a)) b else a

# Sample annotations, author file selection, reading, and mapping author
# genes and samples onto NCBI's.

# --- GEO series matrix: one row per GSM ------------------------------------
parse_series_matrix <- function(paths) {
  one <- function(p) {
    x <- readLines(gzfile(p), warn = FALSE)
    get <- function(tag) {
      l <- x[startsWith(x, tag)]
      lapply(l, function(s) gsub('^"|"$', "", strsplit(s, "\t")[[1]][-1]))
    }
    first <- function(l) if (length(l)) l[[1]] else NULL
    gsm <- first(get("!Sample_geo_accession"))
    if (is.null(gsm)) return(NULL)
    df <- data.frame(gsm = gsm, stringsAsFactors = FALSE)
    df$title  <- first(get("!Sample_title")) %||% rep(NA_character_, length(gsm))
    df$source <- first(get("!Sample_source_name_ch1")) %||% rep(NA_character_, length(gsm))
    for (ch in get("!Sample_characteristics_ch1")) {      # "key: value", one line per key
      lab <- ch[nzchar(ch) & grepl(":", ch)]
      if (!length(lab)) next
      key <- make.names(paste0("ch_", tolower(trimws(sub(":.*$", "", lab[1])))))
      val <- trimws(sub("^[^:]*:\\s*", "", ch))
      val[!grepl(":", ch)] <- NA
      k <- key; i <- 1
      while (k %in% names(df)) { i <- i + 1; k <- paste0(key, "_", i) }
      df[[k]] <- val
    }
    df
  }
  parts <- Filter(Negate(is.null), lapply(paths, one))
  if (!length(parts)) return(NULL)
  cols <- unique(unlist(lapply(parts, names)))
  do.call(rbind, lapply(parts, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
}

# --- which supplementary files are bulk count-like matrices ----------------
not_a_matrix <- "(_RAW\\.tar$|\\.tar(\\.gz)?$|\\.zip$|\\.bam$|\\.bai$|\\.bw$|\\.bigwig$|\\.bed(\\.gz)?$|\\.h5$|\\.h5ad$|\\.loom$|\\.mtx(\\.gz)?$|barcodes|features\\.tsv|genes\\.tsv|\\.rds$|\\.rda$|\\.rdata$|\\.pdf$|\\.fa(sta)?(\\.gz)?$|\\.gtf(\\.gz)?$|\\.gff|\\.vcf|\\.narrowpeak|\\.broadpeak|\\.wig)"
label_counts <- "count|raw|reads|featurecount|htseq|readspergene|umi"
label_norm   <- "tpm|fpkm|rpkm|cpm|norm|log|vst|rlog|expression_values|abundance"

classify_suppl <- function(files, sizes_mb) {
  f <- tolower(files)
  keep <- grepl("\\.(txt|tsv|csv|tab)(\\.gz)?$|\\.xlsx?$", f) & !grepl(not_a_matrix, f)
  single_cell <- grepl("single.?cell|scrna|10x|barcode|cellranger|_sc_", f)
  label <- ifelse(grepl(label_norm, f), "normalised",
           ifelse(grepl(label_counts, f), "counts", "unlabelled"))
  data.frame(file = files, size_mb = sizes_mb, candidate = keep & !single_cell,
             too_large = keep & !single_cell & !is.na(sizes_mb) & sizes_mb > study$max_file_mb,
             single_cell = single_cell, labelled_as = label, stringsAsFactors = FALSE)
}

# --- read an author file the way a user would -------------------------------
# text files go to attest_file() unchanged; spreadsheets are converted first,
# as a user would, and that is recorded
to_text <- function(path) {
  if (!grepl("\\.xlsx?$", tolower(path))) return(list(path = path, converted = FALSE))
  x <- readxl::read_excel(path, sheet = 1, .name_repair = "minimal")
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  out <- sub("\\.xlsx?$", ".from_excel.tsv", path, ignore.case = TRUE)
  utils::write.table(x, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  list(path = out, converted = TRUE)
}

# --- NCBI annotation: GeneID -> Symbol, Ensembl --------------------------------
# base R reads .gz directly; data.table::fread would need R.utils for that
read_tsv_gz <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  utils::read.delim(con, quote = "", comment.char = "", check.names = FALSE,
                    stringsAsFactors = FALSE, na.strings = c("", "NA"))
}

read_annot <- function(path) {
  a <- read_tsv_gz(path)
  pick <- function(pat) { h <- grep(pat, names(a), ignore.case = TRUE, value = TRUE); if (length(h)) h[1] else NA }
  gid <- pick("^GeneID$"); sym <- pick("^Symbol$"); ens <- pick("Ensembl")
  if (is.na(gid)) stop("annotation has no GeneID column")
  data.frame(GeneID = as.character(a[[gid]]),
             Symbol = if (!is.na(sym)) as.character(a[[sym]]) else NA_character_,
             Ensembl = if (!is.na(ens)) sub("\\.[0-9]+$", "", as.character(a[[ens]])) else NA_character_,
             stringsAsFactors = FALSE)
}

read_ncbi_counts <- function(path) {
  x <- read_tsv_gz(path)
  m <- as.matrix(x[, -1, drop = FALSE]); rownames(m) <- as.character(x[[1]])
  storage.mode(m) <- "double"
  m
}

# author row IDs -> NCBI GeneIDs, by whichever system the author used
map_ids <- function(ids, annot) {
  sys <- at_id_system(ids)
  main <- names(sort(table(sys), decreasing = TRUE))[1]
  key <- toupper(sub("\\.[0-9]+$", "", ids))
  gid <- rep(NA_character_, length(ids))
  if (main == "Ensembl")       gid <- annot$GeneID[match(key, toupper(annot$Ensembl))]
  else if (main == "Entrez")   gid <- ifelse(ids %in% annot$GeneID, ids, NA)
  else                         gid <- annot$GeneID[match(toupper(ids), toupper(annot$Symbol))]
  list(gene_id = gid, system = main, mapped_frac = mean(!is.na(gid)))
}

# which NCBI genes the author's ID system can represent at all
representable <- function(annot, system) {
  switch(system,
         Ensembl = annot$GeneID[!is.na(annot$Ensembl) & nzchar(annot$Ensembl) & annot$Ensembl != "-"],
         Entrez  = annot$GeneID,
         annot$GeneID[!is.na(annot$Symbol) & nzchar(annot$Symbol)])
}

# --- author sample columns -> GSMs ---------------------------------------------
norm_name <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

match_columns <- function(A_gene, R, meta) {
  # A_gene: author matrix with rows already NCBI GeneIDs; R: NCBI raw counts
  cn <- colnames(A_gene); gsm <- colnames(R)
  out <- setNames(rep(NA_character_, length(cn)), cn)
  how <- "none"
  # 1. columns are GSM accessions
  hit <- toupper(sub(".*(GSM[0-9]+).*", "\\1", cn, ignore.case = TRUE))
  if (mean(hit %in% gsm) >= 0.8) { out[hit %in% gsm] <- hit[hit %in% gsm]; how <- "accession" }
  # 2. columns are the sample titles
  if (how == "none" && !is.null(meta)) {
    t <- norm_name(meta$title); names(t) <- meta$gsm
    m <- match(norm_name(cn), t)
    if (mean(!is.na(m)) >= 0.8 && !anyDuplicated(m[!is.na(m)])) {
      out[!is.na(m)] <- names(t)[m[!is.na(m)]]; how <- "title"
    }
  }
  # 3. by correlation across shared genes: the same library counted by two
  #    pipelines agrees better than two libraries do - accepted only with a margin
  if (how == "none") {
    g <- intersect(rownames(A_gene), rownames(R)[rowSums(R) > 0])
    if (length(g) >= 500 && ncol(A_gene) >= 2) {
      la <- log1p(pmax(A_gene[g, , drop = FALSE], 0)); lr <- log1p(R[g, , drop = FALSE])
      cc <- suppressWarnings(stats::cor(la, lr))
      # each author sample against every NCBI sample: its best must clear 0.9,
      # beat its runner-up by 0.005, and not be claimed by another sample
      best <- apply(cc, 1, function(v) names(v)[which.max(v)])
      top  <- apply(cc, 1, max)
      gap  <- apply(cc, 1, function(v) { s <- sort(v, decreasing = TRUE); if (length(s) > 1) s[1] - s[2] else Inf })
      good <- top >= 0.9 & gap >= 0.005
      good <- good & !(best %in% best[good][duplicated(best[good])])
      out[good] <- best[good]
      if (sum(good) >= 2) how <- "correlation"
    }
  }
  list(gsm = out, how = how)
}

# collapse author rows that map to the same GeneID (sum), drop unmapped
by_gene <- function(A, gene_id) {
  keep <- !is.na(gene_id)
  if (!any(keep)) return(NULL)
  rowsum(A[keep, , drop = FALSE], gene_id[keep], reorder = FALSE)
}
