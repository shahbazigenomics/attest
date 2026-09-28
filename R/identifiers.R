#' Do the row names identify genes, and each gene once?
#'
#' Adjudicates the claim "every row of this matrix is a gene, named in a way
#' that can be joined to an annotation". Three things go wrong, all of them
#' invisible in the numbers:
#'
#' * **Counting summary rows.** htseq-count appends `__no_feature`,
#'   `__ambiguous`, `__too_low_aQual`, `__not_aligned` and
#'   `__alignment_not_unique`; STAR's `ReadsPerGene.out.tab` starts with four
#'   `N_*` rows. They are not genes, they hold a large share of every library,
#'   and left in they inflate library sizes, the multiple-testing denominator
#'   and the results table itself.
#' * **Excel-mangled symbols.** `SEPT1` becomes `1-Sep`, `MARCH1` becomes
#'   `1-Mar`. HGNC renamed those families in 2020 (`SEPT` to `SEPTIN`, `MARCH`
#'   to `MARCHF`) for this reason. A mangled symbol joins to nothing.
#' * **Identifiers from two sources.** Ensembl IDs with version suffixes mixed
#'   with bare ones, or Ensembl IDs mixed with symbols: half the rows drop out
#'   of any join, silently.
#'
#' Two rules were narrowed by real data. A trailing `.N` is *not* evidence of an
#' Ensembl version: every fission yeast systematic identifier ends that way
#' (`SPAC212.09c`, `SPAC977.03`), so the rule requires an Ensembl stem. And
#' all-digit identifiers are *not* Excel serial numbers, because Entrez gene IDs
#' are all digits; only the textual date forms are reported.
#'
#' @param x counts: matrix, data.frame, SummarizedExperiment, DESeqDataSet or
#'   DGEList.
#' @return object of class "attest_check".
#' @examples
#' fx <- readRDS(system.file("extdata", "fixtures.rds", package = "attest"))
#' m <- fx$airway
#' attest_identifiers(m)
#' # an htseq-count summary row left in the table
#' attest_identifiers(rbind(m, "__no_feature" = round(colSums(m) * 0.1)))
#' @export
attest_identifiers <- function(x) {

  m <- at_as_matrix(x)
  if (is.null(m)) {
    return(at_result("UNKNOWN", "The input could not be read as a numeric matrix.",
                     character(0), NULL, character(0), list(input_class = class(x)[1])))
  }
  id <- rownames(m)
  if (is.null(id)) {
    return(at_result("UNKNOWN",
                     "The matrix has no row names, so the identifiers cannot be checked.",
                     character(0), NULL,
                     "identifier checks: the rows are numbered, not named - read the counts with attest_file(), or set rownames(x)",
                     list(dim = dim(m))))
  }

  ev <- list(dim = dim(m))
  lines <- character(0)
  na <- character(0)
  verdict <- "PERMITTED"
  headline <- "The row names identify genes, one row each."

  # --- rows that are not genes ---------------------------------------------
  summ <- at_summary_rows(id)
  ev$summary_rows <- id[summ]
  if (any(summ)) {
    share <- colSums(m[summ, , drop = FALSE]) / pmax(colSums(m), 1)
    ev$summary_row_share <- share
    verdict  <- "NOT PERMITTED"
    headline <- "Some rows are not genes."
    lines <- c(lines, sprintf(
      "%s %s counting summary rows, not genes, and they hold %s of each library",
      paste(sprintf("'%s'", utils::head(id[summ], 6)), collapse = ", "),
      if (sum(summ) == 1) "is a" else "are",
      if (diff(range(share)) < 0.005) sprintf("about %.1f%%", 100 * mean(share))
      else sprintf("%.1f%% to %.1f%%", 100 * min(share), 100 * max(share))))
    lines <- c(lines, "they were written by htseq-count or STAR at the end (or start) of the table and never removed")
  }

  # --- rows that LOOK like an unrecognised summary row ----------------------
  # at_summary_rows() only knows htseq-count and STAR's own names; a totals
  # row from any other pipeline, or a summary row from a tool that names it
  # differently, has no name in that list and would otherwise be scored as an
  # ordinary gene with no signal at all. A totals/unassigned row is not a
  # naming pattern so much as a *count* pattern: it sits at the very top or
  # bottom of the table (where a pipeline appends or prepends it) and holds a
  # share of the library many times larger than any real gene ever does.
  # Neither signal alone is safe (the first or last gene in an arbitrarily
  # sorted table is unremarkable; one gene legitimately dominates a library in
  # some tissues), so both are required, and the result is a CAUTION, not an
  # automatic removal - this package reports what it cannot rule out rather
  # than silently deciding for the person running it.
  susp <- at_suspect_summary_rows(id, m, summ)
  ev$suspect_summary_rows <- id[susp]
  if (any(susp)) {
    share <- colSums(m[susp, , drop = FALSE]) / pmax(colSums(m), 1)
    if (verdict == "PERMITTED") {
      verdict  <- "CAUTION"
      headline <- "A row at the edge of the table may be an unrecognised summary row, not a gene."
    }
    lines <- c(lines, sprintf(
      "%s %s at the %s of the table and %s of each library on its own - many times any other single row; not a name this package recognises as a known pipeline's summary row, but worth checking it is really a gene",
      paste(sprintf("'%s'", utils::head(id[susp], 3)), collapse = ", "),
      if (sum(susp) == 1) "sits" else "sit",
      paste(unique(ifelse(which(susp) <= nrow(m) / 2, "top", "bottom")), collapse = "/"),
      if (diff(range(share)) < 0.005) sprintf("holds about %.1f%%", 100 * mean(share))
      else sprintf("holds %.1f%% to %.1f%%", 100 * min(share), 100 * max(share))))
  }

  # --- Excel ---------------------------------------------------------------
  mangled <- at_date_like(id)
  ev$excel_mangled <- id[mangled]
  if (any(mangled)) {
    if (verdict == "PERMITTED") {
      verdict  <- "CAUTION"
      headline <- "Some gene symbols have been converted to dates."
    }
    lines <- c(lines, sprintf(
      "%s %s %s: a spreadsheet turned the symbol into a date (1-Sep was SEPT1, now SEPTIN1; 1-Mar was MARCH1, now MARCHF1)",
      format(sum(mangled), big.mark = ","),
      if (sum(mangled) == 1) "row is named" else "rows are named",
      paste(sprintf("'%s'", utils::head(id[mangled], 6)), collapse = ", ")))
  }

  # --- empty, blank, duplicated -------------------------------------------
  blank <- is.na(id) | !nzchar(trimws(id))
  ev$n_blank <- sum(blank)
  if (any(blank)) {
    if (verdict == "PERMITTED") { verdict <- "CAUTION"; headline <- "Some rows have no identifier." }
    lines <- c(lines, sprintf("%s row(s) have an empty or missing name", format(sum(blank), big.mark = ",")))
  }
  dup <- duplicated(id)
  ev$n_duplicate <- sum(dup)
  if (any(dup)) {
    if (verdict == "PERMITTED") { verdict <- "CAUTION"; headline <- "Some genes appear more than once." }
    lines <- c(lines, sprintf(
      "%s identifier(s) appear on more than one row (first: %s), so that gene's reads are split between them",
      format(sum(dup), big.mark = ","),
      paste(utils::head(unique(id[dup]), 3), collapse = ", ")))
  }

  # --- which naming system, and is it one system? --------------------------
  sys <- at_id_system(id)
  real <- !summ & !blank & sys != "spike-in"
  tab <- sort(table(sys[real]), decreasing = TRUE)
  ev$systems <- tab
  ev$n_spike_in <- sum(sys == "spike-in")
  main <- names(tab)[1]
  minor <- tab[-1][tab[-1] >= 10 & tab[-1] / sum(tab) >= 0.01]
  if (length(minor)) {
    if (verdict == "PERMITTED") {
      verdict  <- "CAUTION"
      headline <- "The identifiers come from more than one source."
    }
    lines <- c(lines, sprintf(
      "%s of %s rows are %s identifiers, but %s: a join to any annotation will match one group and drop the other",
      format(as.integer(tab[1]), big.mark = ","), format(sum(tab), big.mark = ","), main,
      paste(sprintf("%s are %s", format(as.integer(minor), big.mark = ","), names(minor)), collapse = ", ")))
  } else {
    lines <- c(lines, sprintf("all %s gene rows carry %s identifiers",
                              format(sum(tab), big.mark = ","), main))
  }

  # --- Ensembl version suffixes -------------------------------------------
  ens <- grepl("^ENS[A-Z]{0,4}[GTPR][0-9]{6,}", id)
  if (any(ens)) {
    v <- grepl("^ENS[A-Z]{0,4}[GTPR][0-9]{6,}\\.[0-9]+$", id)
    ev$n_versioned <- sum(v)
    ev$n_unversioned <- sum(ens & !v)
    if (any(v) && any(ens & !v)) {
      if (verdict == "PERMITTED") {
        verdict  <- "CAUTION"
        headline <- "Some Ensembl identifiers carry a version suffix and some do not."
      }
      lines <- c(lines, sprintf(
        "%s Ensembl IDs end in a version (e.g. %s) and %s do not: the two halves join to different things",
        format(sum(v), big.mark = ","), utils::head(id[v], 1),
        format(sum(ens & !v), big.mark = ",")))
    } else if (any(v)) {
      lines <- c(lines, "every Ensembl ID carries its version suffix; strip it before joining to an annotation built without versions")
    }
  }

  if (ev$n_spike_in > 0)
    lines <- c(lines, sprintf("%d spike-in row(s) present and left out of the counts above", ev$n_spike_in))

  na <- c(na, "whether the identifiers match the annotation you will use: that needs the annotation, which this package does not bundle")

  at_result(verdict, headline, lines,
            if (any(summ))
              paste("These rows are in every total the analysis computes. Library sizes and size factors",
                    "are taken from them, the multiple-testing denominator counts them, and they will appear",
                    "in the results table, usually near the top. Drop them:",
                    "counts <- counts[!grepl('^(__|N_)', rownames(counts)), ].")
            else if (verdict == "CAUTION")
              paste("Nothing here changes the counts. It changes what they can be joined to: rows whose name",
                    "does not match the annotation are dropped by the join, silently, and a gene split across",
                    "two rows is below every filter that either row would have passed alone.")
            else NULL,
            na, ev)
}

# htseq-count, STAR and featureCounts write their totals into the table
# itself. featureCounts' own "Unassigned_*" categories (from its .summary
# file, sometimes pasted into the counts table by hand) were missing here
# entirely - none of them start with "__" or "N_", so they read as ordinary
# genes.
at_summary_rows <- function(id) {
  id <- trimws(as.character(id))
  grepl("^__", id) |
    grepl("^N_(unmapped|multimapping|noFeature|ambiguous)$", id, ignore.case = TRUE) |
    grepl("^Unassigned_(Ambiguity|Ambiguous|MultiMapping|NoFeatures|Unmapped|Secondary|Duplicate|MappingQuality|FragmentLength|Chimera|ReadType|Overlapping_Length)$",
          id, ignore.case = TRUE) |
    tolower(id) %in% c("total", "totals", "total_reads", "sum") |
    id %in% c("no_feature", "ambiguous", "too_low_aQual", "not_aligned",
              "alignment_not_unique", "not_aligned_total", "__no_feature")
}

# A row that is not on the known-name list above but still walks and talks
# like one: at the very top or bottom of the table (where every pipeline that
# appends totals puts them), and its share of the library is far outside what
# any real gene reaches. The threshold (10x the 99th percentile of every
# OTHER row's share, and at least 15% of the library on its own) is
# deliberately conservative - a hit here is reported as a CAUTION to look at,
# not treated as proven.
at_suspect_summary_rows <- function(id, m, known) {
  n <- nrow(m)
  out <- rep(FALSE, n)
  if (n < 20 || ncol(m) < 1) return(out)
  edge <- rep(FALSE, n)
  edge[seq_len(min(3L, n))] <- TRUE
  edge[seq.int(max(1L, n - 2L), n)] <- TRUE
  edge <- edge & !known
  if (!any(edge)) return(out)
  tot <- pmax(colSums(m), 1)
  row_share <- sweep(m, 2, tot, "/")
  own <- apply(row_share, 1, max)                     # each row's largest single-sample share
  rest <- own[!known]
  if (sum(!known) < 10) return(out)
  ref <- stats::quantile(rest[!edge[!known]], 0.99, na.rm = TRUE, names = FALSE)
  if (!is.finite(ref) || ref <= 0) return(out)
  out[edge] <- own[edge] >= max(10 * ref, 0.15)
  out
}

# Excel's damage. Deliberately NOT all-digit strings: Entrez gene IDs are digits.
at_date_like <- function(id) {
  id <- trimws(as.character(id))
  mon <- "(jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)"
  grepl(sprintf("^[0-9]{1,2}[-/ ]%s([-/ ][0-9]{2,4})?$", mon), id, ignore.case = TRUE) |
    grepl(sprintf("^%s[-/ ][0-9]{1,2}([-/ ][0-9]{2,4})?$", mon), id, ignore.case = TRUE) |
    grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", id)
}

# One label per identifier. The catch-all "symbol or systematic" bucket is what
# keeps fission yeast (SPAC212.09c) and other species' systematic names from
# looking like a second source.
at_id_system <- function(id) {
  id <- trimws(as.character(id))
  out <- rep("symbol or systematic", length(id))
  out[grepl("^ENS[A-Z]{0,4}[GTPR][0-9]{6,}(\\.[0-9]+)?$", id)]          <- "Ensembl"
  out[grepl("^FBgn[0-9]+$", id)]                                        <- "FlyBase"
  out[grepl("^WBGene[0-9]+$", id)]                                      <- "WormBase"
  out[grepl("^(NM|NR|XM|XR|NP|XP)_[0-9]+(\\.[0-9]+)?$", id)]            <- "RefSeq"
  out[grepl("^[0-9]+$", id)]                                            <- "Entrez"
  out[grepl("^(ERCC-[0-9]+|SIRV[0-9A-Za-z]*|gSpikein.*)$", id, ignore.case = TRUE)] <- "spike-in"
  out
}
