#' Were the reads counted on the right strand?
#'
#' Adjudicates the claim "these counts were made with the strand setting that
#' matches the library". A stranded library counted on the wrong strand keeps
#' only the reads that happen to fall on antisense-overlapping genes - a few
#' per cent of the data, silently. Nothing in the count matrix shows it; the
#' evidence is in the files the counting step writes beside it.
#'
#' **STAR `ReadsPerGene.out.tab`** is conclusive. It holds three count columns -
#' unstranded, forward and reverse - so the library's protocol follows from the
#' forward share, forward / (forward + reverse): near 1 forward-stranded, near 0
#' reverse-stranded (dUTP, TruSeq Stranded), near 0.5 unstranded. The cut-offs
#' (0.8, 0.2, and 0.4-0.6) are the ones conventionally applied to RSeQC's
#' `infer_experiment.py`; libraries simulated and counted with real STAR gave
#' 0.91, 0.07-0.10 and 0.50, and 60% background reads left the reverse call intact
#' (`validation/strand_sim/`). Given `counts`, each of its columns is matched
#' exactly against the three STAR columns, which says which one was used.
#'
#' **A featureCounts `.summary`** is not. A stranded library counted on the
#' wrong strand collapses the assigned share (10-11% in the simulation), but so
#' can a library whose reads mostly fall outside the annotation; and an
#' unstranded library counted as stranded keeps about half its reads assigned,
#' which looks like an ordinary run while discarding half the data. So a
#' summary can raise a collapse, naming every explanation, and never confirm
#' the setting. featureCounts records its command line in the first line of
#' the counts file; given that file, the `-s` setting is read from it.
#'
#' @param files paths to STAR `ReadsPerGene.out.tab` files (one per sample), or
#'   to featureCounts output: the `.summary` file and, optionally, the counts
#'   file whose first line holds the command.
#' @param counts optional count matrix (or any object [attest()] accepts) to
#'   match against the STAR columns.
#' @return object of class "attest_check".
#' @examples
#' d <- system.file("extdata", "strand", package = "attest")
#' # reverse-stranded (dUTP) libraries, from STAR
#' attest_strandedness(list.files(d, "^reverse_s.*ReadsPerGene", full.names = TRUE))
#' # the same libraries counted by featureCounts with -s 1: the wrong strand
#' attest_strandedness(file.path(d, c("fc_s1.txt", "fc_s1.txt.summary")))
#' @export
attest_strandedness <- function(files, counts = NULL) {
  files <- as.character(files)
  if (!length(files) || !all(file.exists(files))) {
    miss <- files[!file.exists(files)]
    return(at_result("UNKNOWN",
                     sprintf("No file at %s.", paste(utils::head(miss, 3), collapse = ", ")),
                     character(0), NULL, character(0), list(files = files)))
  }
  kind <- vapply(files, at_strand_file_kind, character(1))
  if (any(kind == "star")) return(at_strand_star(files[kind == "star"], counts))
  if (any(kind == "fc_summary")) return(at_strand_fc(files[kind == "fc_summary"][1],
                                                     files[kind == "fc_counts"]))
  at_result("UNKNOWN",
            "None of these files is a STAR ReadsPerGene.out.tab or a featureCounts summary.",
            character(0), NULL, character(0), list(files = files, kind = kind))
}

at_strand_file_kind <- function(path) {
  first <- tryCatch(readLines(path, n = 5, warn = FALSE), error = function(e) character(0))
  if (!length(first)) return("unknown")
  if (grepl("^N_unmapped\t", first[1])) return("star")
  if (grepl("^Status\t", first[1])) return("fc_summary")
  if (grepl("^# Program:featureCounts", first[1])) return("fc_counts")
  "unknown"
}

at_read_star <- function(path) {
  x <- utils::read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                         col.names = c("id", "unstranded", "forward", "reverse"))
  genes <- x[!grepl("^N_", x$id), ]
  m <- as.matrix(genes[, 2:4]); rownames(m) <- genes$id
  list(sample = sub("\\.?ReadsPerGene\\.out\\.tab$", "", basename(path)), counts = m)
}

at_strand_call <- function(share)
  ifelse(is.na(share), "unclear",
  ifelse(share >= 0.8, "forward",
  ifelse(share <= 0.2, "reverse",
  ifelse(share >= 0.4 & share <= 0.6, "unstranded", "unclear"))))

at_strand_column <- c(unstranded = "unstranded", forward = "forward", reverse = "reverse")

at_strand_star <- function(files, counts) {
  st <- lapply(files, at_read_star)
  names(st) <- make.unique(vapply(st, `[[`, character(1), "sample"))
  share <- vapply(st, function(s) {
    fw <- sum(s$counts[, "forward"]); rv <- sum(s$counts[, "reverse"])
    if (fw + rv == 0) NA_real_ else fw / (fw + rv) }, numeric(1))
  call <- at_strand_call(share)
  ev <- list(files = files, forward_share = share, protocol = call)
  lines <- character(0)
  na <- character(0)

  tab <- table(call)
  lines <- c(lines, sprintf("forward share per sample %.2f-%.2f: %s",
                            min(share, na.rm = TRUE), max(share, na.rm = TRUE),
                            paste(sprintf("%d %s", tab, names(tab)), collapse = ", ")))
  protocols <- setdiff(unique(call), "unclear")
  verdict <- "PERMITTED"
  consequence <- NULL
  if (any(call == "unclear")) {
    verdict  <- "CAUTION"
    headline <- "The strandedness of some libraries is unclear."
    lines <- c(lines, sprintf("%s: forward share between the stranded and unstranded ranges",
                              paste(names(call)[call == "unclear"], collapse = ", ")))
    consequence <- "A share between 0.2 and 0.4 (or 0.6 and 0.8) usually means a stranded kit that worked poorly, or a mixture. Check the kit; if in doubt the unstranded column is the safe choice."
  } else if (length(protocols) > 1) {
    verdict  <- "CAUTION"
    headline <- "The libraries were not all made the same way."
    consequence <- "Samples prepared with different strand protocols need different STAR columns, and the protocol itself becomes a batch variable. Take each sample's own column, and check it is not confounded with the comparison."
  } else {
    headline <- sprintf("The libraries are %s; STAR's '%s' column is the one to use.",
                        switch(protocols, forward = "forward-stranded",
                               reverse = "reverse-stranded (dUTP / TruSeq Stranded)",
                               unstranded = "unstranded"),
                        at_strand_column[[protocols]])
  }
  ev$recommended_column <- if (length(protocols) == 1) at_strand_column[[protocols]] else NA_character_

  # --- which column did the matrix come from? --------------------------------
  if (!is.null(counts)) {
    m <- at_as_matrix(counts)
    if (is.null(m) || is.null(rownames(m))) {
      na <- c(na, "which STAR column the count matrix came from: the matrix could not be read, or has no gene names to align on")
    } else {
      used <- at_match_star_columns(m, st)
      ev$matrix_column <- used
      if (all(is.na(used))) {
        na <- c(na, "which STAR column the count matrix came from: no matrix column is identical to any STAR column (filtered, summed or from another run)")
      } else {
        u <- used[!is.na(used)]
        lines <- c(lines, sprintf("%d of %d matrix columns are identical to a STAR column: %s",
                                  length(u), ncol(m),
                                  paste(sprintf("%d from '%s'", as.integer(table(u)), names(table(u))), collapse = ", ")))
        wrong <- if (length(protocols) == 1) u[u != at_strand_column[[protocols]]] else character(0)
        if (length(wrong)) {
          p <- protocols
          opposite <- (p == "forward" && any(wrong == "reverse")) || (p == "reverse" && any(wrong == "forward"))
          if (opposite) {
            verdict  <- "NOT PERMITTED"
            headline <- "The matrix was counted on the wrong strand."
            sm <- stats::median(share, na.rm = TRUE)
            kept <- if (p == "reverse") sm / (1 - sm) else (1 - sm) / sm   # wrong column / right column
            consequence <- sprintf("For a %s library the '%s' column holds only reads that fall on the opposite strand of a gene - %.0f%% as many as the right column here, and those are antisense reads, not the gene's own. Rebuild the matrix from the '%s' column.",
                                   if (p == "reverse") "reverse-stranded" else "forward-stranded",
                                   wrong[1], 100 * kept, at_strand_column[[p]])
          } else if (p == "unstranded") {
            if (verdict == "PERMITTED") verdict <- "CAUTION"
            headline <- "An unstranded library was counted as if it were stranded."
            consequence <- "Each gene keeps only the reads that happened to come from one strand - about half - so every comparison has half the data it should. Rebuild the matrix from the 'unstranded' column."
          } else {
            if (verdict == "PERMITTED") verdict <- "CAUTION"
            headline <- "A stranded library was counted without strand."
            consequence <- "Reads from a gene's antisense partner are counted as its own where genes overlap on opposite strands. Most genes are unaffected; overlapping ones are not. The stranded column is the better choice."
          }
        }
      }
    }
  } else {
    na <- c(na, "which column the count matrix was built from: pass counts = <matrix> to check it against the three STAR columns")
  }

  at_result(verdict, headline, lines, consequence, na, ev)
}

# For each matrix column, the STAR column it is identical to, over shared genes
at_match_star_columns <- function(m, st) {
  out <- rep(NA_character_, ncol(m)); names(out) <- colnames(m)
  for (j in seq_len(ncol(m))) {
    for (s in st) {
      g <- intersect(rownames(m), rownames(s$counts))
      if (length(g) < 10) next
      for (k in colnames(s$counts)) {
        if (isTRUE(all(m[g, j] == s$counts[g, k]))) { out[j] <- k; break }
      }
      if (!is.na(out[j])) break
    }
  }
  out
}

at_read_fc_summary <- function(path) {
  x <- utils::read.table(path, sep = "\t", header = TRUE, check.names = FALSE,
                         stringsAsFactors = FALSE, comment.char = "")
  m <- as.matrix(x[, -1, drop = FALSE]); rownames(m) <- x[[1]]
  colnames(m) <- sub("\\.(Aligned\\.)?(sortedByCoord\\.)?(out\\.)?bam$", "", basename(colnames(m)))
  m
}

# featureCounts writes its command into the counts file's first line
at_fc_strand_setting <- function(path) {
  if (!length(path)) return(NA_integer_)
  l <- readLines(path[1], n = 1, warn = FALSE)
  s <- regmatches(l, regexec('"-s"\\s+"([0-9,]+)"', l))[[1]]
  if (length(s) < 2) s <- regmatches(l, regexec("-s\\s+([0-9,]+)", l))[[1]]
  if (length(s) < 2) return(if (grepl("featureCounts", l)) 0L else NA_integer_)   # default is -s 0
  v <- as.integer(strsplit(s[2], ",")[[1]])
  if (length(unique(v)) == 1) v[1] else NA_integer_
}

at_strand_fc <- function(summary_path, counts_path) {
  s <- tryCatch(at_read_fc_summary(summary_path), error = function(e) NULL)
  if (is.null(s) || !"Assigned" %in% rownames(s))
    return(at_result("UNKNOWN", "The featureCounts summary could not be read.",
                     character(0), NULL, character(0), list(file = summary_path)))
  get <- function(r) if (r %in% rownames(s)) s[r, ] else rep(0, ncol(s))
  decided <- get("Assigned") + get("Unassigned_NoFeatures") + get("Unassigned_Ambiguity")
  share <- ifelse(decided > 0, get("Assigned") / decided, NA_real_)
  setting <- at_fc_strand_setting(counts_path)
  ev <- list(file = summary_path, assigned_share = share, strand_setting = setting)
  med <- stats::median(share, na.rm = TRUE)
  lines <- sprintf("%.0f%% of reads that reached a feature decision were assigned (per sample %.0f%%-%.0f%%); the rest found no feature%s",
                   100 * med, 100 * min(share, na.rm = TRUE), 100 * max(share, na.rm = TRUE),
                   if (any(get("Unassigned_Ambiguity") > 0)) " or overlapped more than one" else "")
  lines <- c(lines, if (is.na(setting)) "the strand setting used is not known: pass the counts file too, whose first line records the featureCounts command"
                    else sprintf("counted with -s %d (%s)", setting,
                                 c("0" = "unstranded", "1" = "forward-stranded", "2" = "reverse-stranded")[[as.character(setting)]]))
  alt <- "or most reads fall outside the annotation - intronic or intergenic signal, rRNA, or an annotation from another genome build"

  if (isTRUE(med < 0.25)) {
    if (isTRUE(setting == 0)) {
      return(at_result("CAUTION", "Most reads were not assigned to any gene.", lines,
        paste("The counting was unstranded, so the strand setting is not the cause; the reads fall outside the annotation - intronic or intergenic signal, rRNA, or an annotation from another genome build. Worth finding out before trusting the counts."),
        character(0), ev))
    }
    return(at_result("CAUTION", "Few reads were assigned - compatible with counting on the wrong strand.", lines,
      paste("A stranded library counted on the wrong strand keeps only reads from antisense-overlapping genes (10-11% assigned in simulation),", alt,
            "- a single summary cannot tell these apart. STAR's ReadsPerGene.out.tab, a second featureCounts run with the other -s, or RSeQC infer_experiment.py on one BAM can."),
      character(0), ev))
  }
  if (isTRUE(setting == 0)) {
    return(at_result("PERMITTED", "Counted without strand, which is valid for any library.", lines,
      NULL, "whether the library was stranded: an unstranded count is valid either way, but a stranded library would have kept reads from antisense-overlapping genes apart",
      ev))
  }
  at_result("UNKNOWN", "A single featureCounts run cannot confirm the strand setting.", lines,
    NULL, paste("whether -s matches the library: an unstranded library counted as stranded still assigns about half its reads (55-57% in simulation) - an ordinary-looking run that has discarded half the data. STAR's ReadsPerGene.out.tab, a featureCounts run with another -s, or RSeQC infer_experiment.py can settle it"),
    ev)
}
