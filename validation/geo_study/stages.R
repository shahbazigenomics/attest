# The study's stages. run_study.R calls these; each is safe to re-run.

logmsg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = file.path(study$results, "log.txt"), append = TRUE)
}
series_file <- function(gse) file.path(study$results, "series", paste0(gse, ".rds"))
cache_dir   <- function(gse) file.path(study$cache, gse)

gene_info_url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"

# GEO's annotation for its counts; NCBI Gene's table if that cannot be fetched
get_annot <- function(url) {
  if (!is.null(net$annot)) return(net$annot)
  a <- tryCatch(read_annot(geo_download(url, file.path(study$cache, "ncbi_annotation.tsv.gz"))),
                error = function(e) { logmsg("GEO annotation unavailable (", conditionMessage(e),
                                             ") - using NCBI Gene's gene_info instead"); NULL })
  if (is.null(a)) a <- read_gene_info(geo_download(gene_info_url, file.path(study$cache, "Homo_sapiens.gene_info.gz")))
  net$annot <- a
  a
}

# --- check: every endpoint proves itself on one real series before the run --
stage_check <- function() {
  ok <- TRUE
  for (p in c("readxl", "DESeq2")) {
    has <- requireNamespace(p, quietly = TRUE)
    logmsg(sprintf("package %-10s %s", p, if (has) "ok" else "MISSING"))
    ok <- ok && has
  }
  fr <- resolve_frame()
  logmsg("counts filter: ", fr$filtered$count, " series; without it: ", fr$unfiltered$count,
         "  (NCBI read the search as: ", fr$filtered$translation, ")")
  if (nzchar(fr$filtered$warnings)) logmsg("  NCBI warnings: ", fr$filtered$warnings)
  for (v in c('"rnaseq counts"[Filter]', 'rnaseq_counts[Filter]', '"rnaseq counts"[All Fields]', 'gse[ETYP]', '"Homo sapiens"[Organism]'))
    logmsg(sprintf("  probe %-30s %s", v, esearch_info(v)$count))
  logmsg("sampling frame: ", fr$basis)
  logmsg("  search: ", fr$term)
  n <- fr$n
  logmsg("  entries in the frame: ", n)
  if (is.na(n) || n == 0) {
    logmsg("STOP: the frame search found nothing - try it at https://www.ncbi.nlm.nih.gov/gds and edit config.R")
    return(invisible(FALSE))
  }
  sh <- series_share(fr$term, n)
  logmsg(sprintf("  share of entries that are series (600 UIDs from start, middle, end): %.1f%%", 100 * sh))
  if (is.na(sh) || sh < 0.95) { logmsg("  WARNING: the frame is not series only - the entry-type restriction is not applied"); ok <- FALSE }
  ids <- esearch_info(fr$term, retmax = 50)$ids
  cand <- paste0("GSE", as.numeric(ids[startsWith(ids, "2")]) - 2e8)
  gse <- NA; tried <- 0
  for (g in utils::head(cand, 15)) {                  # first series in the frame that NCBI has counts for
    tried <- tried + 1; u <- ncbi_counts_urls(g)
    if (u$from_page) { gse <- g; break }
  }
  logmsg(sprintf("  series with NCBI counts among the first %d tried: %s", tried, if (is.na(gse)) "none" else "yes"))
  if (is.na(gse)) { logmsg("CHECK FAILED - no NCBI counts found for any tried series; paste this log to Claude"); return(invisible(FALSE)) }
  logmsg("test series: ", gse)
  logmsg("NCBI raw counts: ", u$raw, "  (found on the download page)")
  logmsg("  download page recognised as served: ", if (u$page_ok) "yes" else "NO - series would all be retried, never classified")
  ok <- ok && u$page_ok
  sz <- geo_size(u$raw)
  logmsg("  size: ", if (is.na(sz)) "unknown" else sprintf("%.1f MB", sz / 1e6))
  logmsg("annotation URL: ", u$annot)
  a <- tryCatch(get_annot(u$annot), error = function(e) { logmsg("annotation FAILED: ", conditionMessage(e)); NULL })
  if (!is.null(a)) logmsg(sprintf("annotation: %d genes; %d with Ensembl IDs; %d with symbols", nrow(a),
                                  sum(!is.na(a$Ensembl) & nzchar(a$Ensembl)), sum(!is.na(a$Symbol) & nzchar(a$Symbol))))
  sm <- list_dir(ftp_series(gse, "matrix")); logmsg("series matrix files: ", paste(sm, collapse = ", "))
  sp <- list_dir(ftp_series(gse, "suppl"));  logmsg("supplementary files: ", paste(utils::head(sp, 10), collapse = ", "))
  raw_ok <- tryCatch({ p <- geo_download(u$raw, file.path(cache_dir(gse), "ncbi_raw.tsv.gz"))
                       R <- read_ncbi_counts(p); logmsg(sprintf("NCBI raw counts read: %d genes x %d samples", nrow(R), ncol(R))); TRUE },
                     error = function(e) { logmsg("NCBI raw counts FAILED: ", conditionMessage(e)); FALSE })
  ok <- ok && !is.null(a) && raw_ok && length(sm) > 0
  logmsg(if (ok) "CHECK PASSED - run the pilot next" else "CHECK FAILED - paste this log to Claude")
  invisible(ok)
}

# --- sample: a fixed random draw of series -----------------------------------
stage_sample <- function() {
  f <- file.path(study$results, "sample.csv")
  if (file.exists(f)) return(utils::read.csv(f, stringsAsFactors = FALSE))
  fr <- resolve_frame()
  logmsg("sampling frame: ", fr$basis, " - ", fr$term)
  all <- esearch_gse(fr$term)
  logmsg("series available: ", length(all))
  set.seed(study$seed)
  s <- data.frame(order = seq_len(min(study$n_candidates, length(all))),
                  gse = sample(all, min(study$n_candidates, length(all))), stringsAsFactors = FALSE)
  utils::write.csv(s, f, row.names = FALSE)
  utils::write.csv(data.frame(n_available = length(all), term = fr$term, basis = fr$basis, seed = study$seed,
                              drawn = nrow(s), date = as.character(Sys.Date())),
                   file.path(study$results, "sample_frame.csv"), row.names = FALSE)
  s
}

# --- one series: fetch, attest, truth --------------------------------------------
process_series <- function(gse) {
  out <- list(gse = gse, status = "ok", files = NULL, sex = NULL, note = character(0))
  u <- ncbi_counts_urls(gse)
  # "no NCBI counts" only when NCBI's own download page says so; a page or file
  # that could not be fetched is a transient failure, not saved, retried next run
  if (!u$page_ok)   { out$status <- "unavailable: download page not served"; return(out) }
  if (!u$from_page) { out$status <- "no NCBI counts"; return(out) }
  R <- tryCatch(read_ncbi_counts(geo_download(u$raw, file.path(cache_dir(gse), "ncbi_raw.tsv.gz"))),
                error = function(e) conditionMessage(e))
  if (is.character(R)) { out$status <- paste("unavailable:", R); return(out) }
  annot <- get_annot(u$annot)
  out$n_gsm_ncbi <- ncol(R)

  mfiles <- grep("series_matrix\\.txt\\.gz$", list_dir(ftp_series(gse, "matrix")), value = TRUE)
  meta <- tryCatch(parse_series_matrix(vapply(mfiles, function(f)
    geo_download(paste0(ftp_series(gse, "matrix"), f), file.path(cache_dir(gse), f)), "")),
    error = function(e) NULL)

  sp <- list_dir(ftp_series(gse, "suppl"))
  if (!length(sp)) { out$status <- "no supplementary files"; return(out) }
  pre <- classify_suppl(sp, rep(NA_real_, length(sp)))
  sizes <- rep(NA_real_, length(sp))
  sizes[pre$candidate] <- vapply(sp[pre$candidate], function(f) geo_size(paste0(ftp_series(gse, "suppl"), f)), 0) / 1e6
  cl <- classify_suppl(sp, sizes)
  out$suppl <- cl
  use <- cl[cl$candidate & !cl$too_large, ]
  if (!nrow(use)) {
    out$status <- if (any(cl$too_large)) "matrix too large" else if (any(cl$single_cell)) "single-cell only"
                  else if (any(grepl("_RAW\\.tar$", cl$file))) "per-sample files only (RAW.tar)" else "no matrix file"
    return(out)
  }
  use <- use[order(match(use$labelled_as, c("counts", "unlabelled", "normalised"))), ]
  use <- utils::head(use, study$max_files_per_gse)

  rows <- list()
  for (i in seq_len(nrow(use))) {
    f <- use$file[i]
    row <- list(gse = gse, file = f, size_mb = use$size_mb[i], labelled_as = use$labelled_as[i])
    path <- tryCatch(geo_download(paste0(ftp_series(gse, "suppl"), f), file.path(cache_dir(gse), f)),
                     error = function(e) NULL)
    if (is.null(path)) { row$status <- "download failed"; rows[[i]] <- row; next }
    tx <- tryCatch(to_text(path), error = function(e) NULL)
    if (is.null(tx)) { row$status <- "unreadable spreadsheet"; rows[[i]] <- row; next }
    row$from_excel <- tx$converted

    rep <- tryCatch(attest_file(tx$path), error = function(e) conditionMessage(e))
    if (is.character(rep)) {                 # a crash in attest: a bug to fix, kept for inspection
      row$status <- "attest error"; row$attest_message <- rep; row$kept <- TRUE; rows[[i]] <- row; next
    }
    v <- function(k) if (is.null(rep$checks[[k]])) NA_character_ else rep$checks[[k]]$verdict
    row$input_headline <- rep$checks$input$headline
    if (is.null(rep$checks[["value scale"]])) {  # attest could not read it as a count table: kept to look at
      row$status <- "not read by attest"; row$kept <- TRUE; rows[[i]] <- row; next
    }
    row$status <- "audited"
    row$attest_overall <- rep$verdict
    row$v_input <- v("input"); row$v_value <- v("value scale"); row$v_identifiers <- v("identifiers")
    row$v_completeness <- v("completeness"); row$v_identity <- v("identity")
    row$value_headline <- if (is.null(rep$checks[["value scale"]])) NA else rep$checks[["value scale"]]$headline
    idm <- rep$checks$identifiers$measurements
    row$summary_rows <- length(idm$summary_rows); row$excel_dates <- length(idm$excel_mangled)
    row$duplicate_ids <- if (is.null(idm$n_duplicate)) NA else idm$n_duplicate
    row$annotation_cols <- length(rep$checks$input$measurements$annotation_columns)

    rd <- tryCatch(at_read_counts(tx$path), error = function(e) NULL)
    A <- if (is.null(rd)) NULL else rd$counts
    if (!is.null(A) && !is.null(rownames(A))) {
      row$n_rows <- nrow(A); row$n_cols <- ncol(A)
      mp <- map_ids(rownames(A), annot)
      row$id_system <- mp$system; row$mapped_frac <- mp$mapped_frac
      Ag <- by_gene(A, mp$gene_id)
      if (!is.null(Ag)) {
        mc <- match_columns(Ag, R, meta)
        ok <- !is.na(mc$gsm)
        row$match_how <- mc$how; row$n_matched <- sum(ok)
        if (sum(ok) >= 2) {
          tv <- truth_value_scale(A[, ok, drop = FALSE], Ag[, ok, drop = FALSE], R[, mc$gsm[ok], drop = FALSE])
          row$truth_value <- tv$label; row$truth_level <- tv$level
          row$truth_level_range <- sprintf("%.3g-%.3g", tv$level_min, tv$level_max)
          row$depth_slope <- tv$depth_slope; row$depth_slope_se <- tv$depth_slope_se
          row$depth_spread <- tv$depth_spread; row$truth_why <- tv$why
          pc <- if (!is.null(annot$Type) && any(!is.na(annot$Type))) annot$GeneID[annot$Type %in% "protein-coding"] else NULL
          tc <- truth_completeness(unique(mp$gene_id[!is.na(mp$gene_id)]), R[, mc$gsm[ok], drop = FALSE],
                                   representable(annot, mp$system), protein_coding = pc)
          row$truth_completeness <- tc$label; row$zero_present_frac <- tc$present_frac
          row$completeness_restricted <- tc$restricted
        }
      }
    }
    flagged <- isTRUE(row$v_value %in% c("NOT PERMITTED", "CAUTION")) ||
               isTRUE(row$truth_value == "depth removed") || isTRUE(row$summary_rows > 0)
    row$kept <- flagged || study$keep_author_files
    if (!row$kept) { unlink(path); if (tx$converted) unlink(tx$path) }
    rows[[i]] <- row
  }
  out$files <- rbind_fill(lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  out$sex <- tryCatch(if (!is.null(meta)) truth_sex(meta, R, annot) else NULL, error = function(e) NULL)
  out$meta <- meta
  out
}

stage_run <- function(limit = Inf) {
  dir.create(file.path(study$results, "series"), recursive = TRUE, showWarnings = FALSE)
  s <- stage_sample()
  usable <- 0; done <- 0
  for (gse in s$gse) {
    f <- series_file(gse)
    res <- if (file.exists(f)) readRDS(f) else {
      logmsg("processing ", gse)
      r <- tryCatch(process_series(gse), error = function(e)
        list(gse = gse, status = paste("error:", conditionMessage(e))))
      if (!startsWith(r$status, "unavailable")) saveRDS(r, f)          # transient: retried on the next run
      logmsg(sprintf("  %s: %s%s", gse, r$status,
                     if (!is.null(r$files)) sprintf(", %d file(s) audited", sum(r$files$status == "audited", na.rm = TRUE)) else ""))
      r
    }
    done <- done + 1
    if (!is.null(res$files) && any(res$files$status == "audited", na.rm = TRUE)) usable <- usable + 1
    if (usable >= study$target_usable || done >= limit) break
  }
  logmsg(sprintf("series processed: %d, with an auditable author matrix: %d", done, usable))
  invisible(list(done = done, usable = usable))
}

# --- consequence: DESeq2 on the author's file vs NCBI's raw counts -------------
pick_condition <- function(meta, gsms) {
  m <- meta[match(gsms, meta$gsm), , drop = FALSE]
  cand <- grep("^ch_", names(m), value = TRUE)
  cand <- cand[!grepl("sex|gender|age|individual|donor|patient|subject|replicate|batch|id$|time|lane|run", cand)]
  pref <- grepl("treat|condition|group|disease|status|genotype|agent|infect|stimul|knock|drug|diagnos", cand)
  for (k in c(cand[pref], cand[!pref])) {
    x <- m[[k]]; tab <- table(x)
    if (length(tab) == 2 && all(tab >= 3) && !anyNA(x)) return(list(key = k, cond = factor(x)))
  }
  NULL
}

stage_consequence <- function() {
  if (!requireNamespace("DESeq2", quietly = TRUE)) { logmsg("DESeq2 not installed; skipping"); return(invisible(NULL)) }
  out <- list()
  for (f in list.files(file.path(study$results, "series"), full.names = TRUE)) {
    r <- readRDS(f)
    if (is.null(r$files) || is.null(r$meta)) next
    for (i in seq_len(nrow(r$files))) {
      x <- r$files[i, ]
      if (!isTRUE(x$kept) || !isTRUE(x$truth_value == "depth removed")) next
      path <- file.path(cache_dir(r$gse), x$file)
      if (isTRUE(x$from_excel)) path <- sub("\\.xlsx?$", ".from_excel.tsv", path, ignore.case = TRUE)
      if (!file.exists(path)) next
      R <- read_ncbi_counts(file.path(cache_dir(r$gse), "ncbi_raw.tsv.gz"))
      A <- at_read_counts(path)$counts
      annot <- get_annot(ncbi_counts_urls(r$gse)$annot)
      mp <- map_ids(rownames(A), annot); Ag <- by_gene(A, mp$gene_id)
      mc <- match_columns(Ag, R, r$meta); ok <- !is.na(mc$gsm)
      pc <- pick_condition(r$meta, mc$gsm[ok])
      if (is.null(pc) || any(Ag < 0)) next
      sh <- data.frame(cond = pc$cond, row.names = mc$gsm[ok])
      de <- function(counts) {
        counts <- round(counts); storage.mode(counts) <- "integer"
        counts <- counts[rowSums(counts) > 0, , drop = FALSE]
        dds <- suppressMessages(DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(counts, sh, ~ cond), quiet = TRUE))
        res <- DESeq2::results(dds, alpha = 0.05)
        rownames(res)[!is.na(res$padj) & res$padj < 0.05]
      }
      d_ncbi <- tryCatch(de(R[, mc$gsm[ok], drop = FALSE]), error = function(e) NULL)
      Aa <- Ag[, ok, drop = FALSE]; colnames(Aa) <- mc$gsm[ok]
      d_auth <- tryCatch(de(Aa), error = function(e) NULL)
      if (is.null(d_ncbi) || is.null(d_auth)) next
      shared <- intersect(rownames(R), rownames(Ag))
      out[[length(out) + 1]] <- data.frame(
        gse = r$gse, file = x$file, labelled_as = x$labelled_as, attest_value = x$v_value,
        condition = pc$key, n_per_group = paste(table(pc$cond), collapse = "+"),
        de_ncbi_raw = length(d_ncbi), de_author_file = length(d_auth),
        de_both = length(intersect(d_ncbi, d_auth)),
        de_ncbi_on_shared_genes = length(intersect(d_ncbi, shared)), stringsAsFactors = FALSE)
      logmsg(sprintf("consequence %s %s: NCBI raw %d DE genes, author file %d, both %d",
                     r$gse, x$file, length(d_ncbi), length(d_auth), length(intersect(d_ncbi, d_auth))))
    }
  }
  res <- if (length(out)) do.call(rbind, out) else data.frame()
  utils::write.csv(res, file.path(study$results, "consequence.csv"), row.names = FALSE)
  invisible(res)
}
