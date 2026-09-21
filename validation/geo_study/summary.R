# Turns results/series/*.rds into the study's tables: results/summary.md and CSVs.

rbind_fill <- function(l) {
  l <- Filter(function(d) !is.null(d) && nrow(d) > 0, l)
  if (!length(l)) return(data.frame())
  cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
}

# Wilson 95% interval for k of n
wilson <- function(k, n) {
  if (n == 0) return(c(NA, NA))
  p <- k / n; z <- 1.96
  pmin(1, pmax(0, c((p + z^2/(2*n) - z * sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n),
                    (p + z^2/(2*n) + z * sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n))))
}
pct <- function(k, n) {
  if (n == 0) return("-")
  ci <- wilson(k, n); sprintf("%d / %d (%.1f%%, 95%% CI %.1f-%.1f)", k, n, 100*k/n, 100*ci[1], 100*ci[2])
}
md_table <- function(d) {
  if (!nrow(d)) return("(none)\n")
  paste0("| ", paste(names(d), collapse = " | "), " |\n|", paste(rep("---", ncol(d)), collapse = "|"), "|\n",
         paste(apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), collapse = "\n"), "\n")
}

stage_summary <- function() {
  res <- lapply(list.files(file.path(study$results, "series"), full.names = TRUE), readRDS)
  if (!length(res)) { logmsg("no results yet"); return(invisible(NULL)) }
  frame <- tryCatch(utils::read.csv(file.path(study$results, "sample_frame.csv")), error = function(e) NULL)

  status <- vapply(res, function(r) r$status, "")
  files <- rbind_fill(lapply(res, function(r) r$files))
  aud <- if (nrow(files)) files[files$status %in% "audited", , drop = FALSE] else files
  series_usable <- unique(aud$gse)
  sex <- rbind_fill(lapply(res, function(r) if (!is.null(r$sex)) cbind(gse = r$gse, r$sex) else NULL))
  cons <- tryCatch(utils::read.csv(file.path(study$results, "consequence.csv"), stringsAsFactors = FALSE),
                   error = function(e) data.frame())

  L <- character(0); add <- function(...) L <<- c(L, paste0(...))
  add("# attest on published GEO data\n")
  add(sprintf("Generated %s from %d processed series.\n", Sys.Date(), length(res)))
  if (!is.null(frame)) add(sprintf("Sampling frame: %s series matching `%s`; %d drawn at random (seed %s).\n",
                                   format(frame$n_available, big.mark = ","), frame$term, frame$drawn, frame$seed))

  add("## Flow\n")
  st <- as.data.frame(table(status = ifelse(status == "ok" & !vapply(res, function(r) r$gse %in% series_usable, TRUE),
                                            "ok, but no file could be audited", status)), stringsAsFactors = FALSE)
  add(md_table(st[order(-st$Freq), ]))
  add(sprintf("\nSeries with at least one audited author matrix: **%d**. Files audited: **%d**.\n",
              length(series_usable), nrow(aud)))

  if (nrow(aud)) {
    n <- nrow(aud); ns <- length(series_usable)
    by_series <- function(flag) length(unique(aud$gse[flag]))
    pos <- function(x) !is.na(x) & x > 0
    add("\n## What attest found (as a user would run it: `attest_file()` on the downloaded file)\n")
    tab <- data.frame(
      finding = c("value scale NOT PERMITTED (not raw counts)", "value scale CAUTION",
                  "counting summary rows left in (htseq/STAR)", "gene names turned into dates (Excel)",
                  "duplicated gene identifiers", "annotation columns in the table (featureCounts)",
                  "matrix filtered upstream (completeness CAUTION)"),
      files = c(pct(sum(aud$v_value %in% "NOT PERMITTED"), n), pct(sum(aud$v_value %in% "CAUTION"), n),
                pct(sum(aud$summary_rows > 0, na.rm = TRUE), n), pct(sum(aud$excel_dates > 0, na.rm = TRUE), n),
                pct(sum(aud$duplicate_ids > 0, na.rm = TRUE), n), pct(sum(aud$annotation_cols > 0, na.rm = TRUE), n),
                pct(sum(aud$v_completeness %in% "CAUTION"), n)),
      series = c(pct(by_series(aud$v_value %in% "NOT PERMITTED"), ns), pct(by_series(aud$v_value %in% "CAUTION"), ns),
                 pct(by_series(pos(aud$summary_rows)), ns), pct(by_series(pos(aud$excel_dates)), ns),
                 pct(by_series(pos(aud$duplicate_ids)), ns), pct(by_series(pos(aud$annotation_cols)), ns),
                 pct(by_series(aud$v_completeness %in% "CAUTION"), ns)))
    add(md_table(tab))

    add("\n## attest against the NCBI-based truth (value scale)\n")
    add("Truth comes from NCBI's own raw counts for the same samples: the slope of log(author total) on log(NCBI total) across samples - near 1 when sequencing depth is still in the values, near 0 when it has been divided out. It uses none of attest's rules.\n")
    det <- aud[!is.na(aud$truth_value) & !aud$truth_value %in% c("undetermined", "ambiguous"), , drop = FALSE]
    add(sprintf("Files with a determined truth: %d of %d audited (the rest: too few matched samples, too-even depths, or samples not matched).\n",
                nrow(det), n))
    if (nrow(det)) {
      ct <- as.data.frame.matrix(table(truth = det$truth_value, attest = det$v_value))
      add(md_table(cbind(truth = rownames(ct), ct)))
      notraw <- det$truth_value %in% c("depth removed", "log-transformed")
      add(sprintf("\n- **Sensitivity** - non-raw files attest did not call raw: %s\n",
                  pct(sum(notraw & det$v_value != "PERMITTED"), sum(notraw))))
      add(sprintf("- **Specificity** - raw-count files attest called raw: %s\n",
                  pct(sum(det$truth_value == "raw counts" & det$v_value == "PERMITTED"), sum(det$truth_value == "raw counts"))))
      mis <- det[det$labelled_as == "counts" & notraw, , drop = FALSE]
      add(sprintf("- **Named as counts but not counts** (file name says count/raw/reads; NCBI comparison says depth removed or log): %s of files named as counts; attest flagged %d of them\n",
                  pct(nrow(mis), sum(det$labelled_as == "counts")), sum(mis$v_value != "PERMITTED")))
      utils::write.csv(det, file.path(study$results, "truth_value_scale.csv"), row.names = FALSE)
    }

    add("\n## Completeness against the NCBI-based truth\n")
    dc <- aud[!is.na(aud$truth_completeness) & aud$truth_completeness %in% c("filtered", "complete"), , drop = FALSE]
    if (nrow(dc)) {
      ct <- as.data.frame.matrix(table(truth = dc$truth_completeness, attest = dc$v_completeness))
      add(md_table(cbind(truth = rownames(ct), ct)))
    } else add("(no file with a determined truth)\n")
  }

  add("\n## Sex: GEO annotation against XIST / Y-gene expression in NCBI's counts\n")
  if (nrow(sex)) {
    sser <- unique(sex$gse)
    add(sprintf("- series with sex annotated for >= 2 samples: %d\n", length(sser)))
    add(sprintf("- samples whose expression contradicts their label: %s\n", pct(sum(sex$mismatch), nrow(sex))))
    add(sprintf("- series with at least one such sample: %s\n", pct(length(unique(sex$gse[sex$mismatch])), length(sser))))
    utils::write.csv(sex[sex$mismatch, ], file.path(study$results, "sex_mismatches.csv"), row.names = FALSE)
  } else add("(no series with usable sex annotation)\n")

  add("\n## Consequence: DESeq2 on the author's file vs on NCBI's raw counts\n")
  add("For files NCBI shows to be non-raw, with a two-level condition (>= 3 samples each) in the GEO annotation. Same samples, same design (`~ cond`), padj < 0.05.\n")
  add(md_table(cons))

  add("\n## Exclusions\n")
  if (nrow(files)) add(md_table(as.data.frame(table(file_status = files$status), stringsAsFactors = FALSE)))
  writeLines(L, file.path(study$results, "summary.md"))
  if (nrow(files)) utils::write.csv(files, file.path(study$results, "files.csv"), row.names = FALSE)
  logmsg("wrote ", file.path(study$results, "summary.md"))
  invisible(L)
}
