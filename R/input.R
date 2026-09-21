#' Audit a count file
#'
#' The same report as [attest()], for counts that are still in a file: a
#' featureCounts table, a GEO supplementary file, a CSV from a collaborator.
#' Most transformed matrices reach the person who did not produce them, and they
#' arrive as files, so reading the file is part of the audit rather than
#' something to get right before it starts.
#'
#' What it handles, each one reported rather than done silently:
#' * `#`-prefixed header lines, as featureCounts writes
#' * tab, comma or semicolon separation, detected from the file
#' * gzip (`.gz`)
#' * the gene identifier column, whether it is named (`Geneid`, `gene_id`, ...),
#'   unnamed, or already the row names
#' * featureCounts' annotation columns - `Chr`, `Start`, `End`, `Strand`,
#'   `Length`. Three of them are numeric and would otherwise be read as three
#'   extra samples, a mistake no check on the matrix can see afterwards
#' * duplicated gene identifiers
#' * for featureCounts output, the `.summary` file beside it, which goes to
#'   [attest_strandedness()] together with the `-s` setting recorded in the
#'   file's first line
#'
#' @param path path to the file. `.gz` is read directly.
#' @param sep field separator. `NULL` (default) detects it from the file.
#' @param ... passed to [attest()]: `metadata`, `design`, `group`, `n_expected`
#'   and the rest.
#' @return object of class "attest_report", whose first check, `input`, says how
#'   the file was read.
#' @examples
#' fx <- readRDS(system.file("extdata", "fixtures.rds", package = "attest"))
#' p <- tempfile(fileext = ".tsv")
#' utils::write.table(data.frame(Geneid = rownames(fx$pasilla), fx$pasilla),
#'                    p, sep = "\t", quote = FALSE, row.names = FALSE)
#' attest_file(p)
#' @export
attest_file <- function(path, sep = NULL, ...) {
  read <- at_read_counts(path, sep)
  if (is.null(read$counts)) {
    return(structure(list(checks = list(input = read$check),
                          verdict = read$check$verdict,
                          not_run = "every check on the matrix: the file could not be read as counts",
                          source_note = NULL),
                     class = "attest_report"))
  }
  args <- list(...)
  # featureCounts writes <file>.summary beside its counts; use it unless told otherwise
  if (is.null(args$strandedness) && isTRUE(read$check$measurements$featurecounts)) {
    summ <- paste0(path, ".summary")
    if (file.exists(summ)) args$strandedness <- c(path, summ)
  }
  rep <- do.call(attest, c(list(read$counts), args))
  rep$checks <- c(list(input = read$check), rep$checks)
  kinds <- vapply(rep$checks, at_kind, character(1))
  rep$verdict <- at_worst(vapply(rep$checks[kinds == "fault"], function(c) c$verdict, character(1)))
  rep
}

# featureCounts' annotation block, plus the names other pipelines use for the
# same thing. Start/End/Length are numeric, so nothing downstream can tell them
# from samples once they are in the matrix.
at_annotation_names <- function()
  c("chr", "chrs", "start", "end", "strand", "length", "genelength",
    "gene_name", "genename", "symbol", "gene_symbol", "biotype",
    "gene_biotype", "description", "seqnames", "width", "position", "pos")

# the names R invents when a column has no header of its own
at_unnamed <- function(x)
  is.na(x) || !nzchar(trimws(x)) || grepl("^(V1|X|X1|\\.\\.\\.1|NA)$", trimws(x))

at_id_names <- function()
  c("geneid", "gene_id", "gene", "genes", "id", "ids", "name", "target_id",
    "transcript_id", "ensembl_gene_id", "feature", "featureid", "x", "")

at_read_counts <- function(path, sep = NULL) {
  fail <- function(msg, ev = list())
    list(counts = NULL,
         check  = at_result("UNKNOWN", msg, character(0), NULL, character(0), ev))

  if (!length(path) || !is.character(path) || !file.exists(path))
    return(fail(sprintf("No file at %s.", paste(path, collapse = ", ")),
                list(path = path)))

  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  head_lines <- tryCatch(readLines(con, n = 30, warn = FALSE), error = function(e) NULL)
  close(con)
  if (is.null(head_lines) || !length(head_lines))
    return(fail("The file is empty or could not be read as text.", list(path = path)))

  n_comment <- sum(cumsum(!grepl("^#", head_lines)) == 0)
  body <- head_lines[seq.int(n_comment + 1, length(head_lines))]
  if (!length(body)) return(fail("The file holds only comment lines.", list(path = path)))

  if (is.null(sep)) {
    # the right separator gives the SAME field count line after line; a count
    # like ";" can win on raw occurrences while being useless as a delimiter,
    # e.g. featureCounts' Chr/Start/End columns join multiple exons with ";"
    # (hundreds per line) while tab stays at a constant field count
    consistency <- function(s) {
      w <- vapply(utils::head(body, 8), function(l) length(strsplit(l, s, fixed = TRUE)[[1]]), integer(1))
      m <- as.integer(names(sort(table(w), decreasing = TRUE))[1])
      c(fields = m, consistent = mean(w == m))
    }
    cand <- c("\t", ",", ";")
    cs <- sapply(cand, consistency)
    good <- cs["fields", ] >= 2 & cs["consistent", ] >= 0.75
    if (any(good)) {
      sep <- cand[good][which.max(cs["fields", good])]
    } else {
      per_line <- function(s) stats::median(vapply(utils::head(body, 5), function(l)
        as.numeric(length(gregexpr(s, l, fixed = TRUE)[[1]][gregexpr(s, l, fixed = TRUE)[[1]] > 0])),
        numeric(1)))
      n <- vapply(cand, per_line, numeric(1))
      if (max(n) >= 1) {
        sep <- cand[which.max(n)]
      } else if (per_line(" ") >= 1) {
        sep <- " "
      } else {
        return(fail("No tab, comma or semicolon in the first data line, so this does not look like a count table.",
                    list(first_line = substr(body[1], 1, 120))))
      }
    }
  }

  # rows above the real header: a title or a merged-cell super-header,
  # exported as mostly-blank lines padded to the file's column count. The
  # real header is the first line that is (almost) entirely filled and is
  # itself followed by a line that is also filled and mostly numeric.
  filled <- function(l) { f <- strsplit(l, sep, fixed = TRUE)[[1]]; if (!length(f)) return(0)
                          mean(nzchar(trimws(f))) }
  numeric_frac <- function(l) { f <- strsplit(l, sep, fixed = TRUE)[[1]]; f <- trimws(f)
                                f <- f[nzchar(f)]; if (!length(f)) return(0)
                                mean(!is.na(suppressWarnings(as.numeric(f)))) }
  n_banner <- 0L
  window <- seq_len(min(15L, length(body) - 1L))
  for (i in window) {
    if (filled(body[i]) >= 0.8 && filled(body[i + 1]) >= 0.8 && numeric_frac(body[i + 1]) >= 0.5) {
      n_banner <- i - 1L; break
    }
  }
  skip <- n_comment + n_banner

  df <- tryCatch(
    utils::read.table(path, sep = sep, header = TRUE, skip = skip,
                      check.names = FALSE, quote = "\"", comment.char = "",
                      stringsAsFactors = FALSE),
    error = function(e) NULL)
  if (is.null(df) || !ncol(df))
    return(fail(sprintf("The file could not be parsed as a '%s'-separated table.",
                        if (sep == "\t") "tab" else sep), list(sep = sep)))

  ev <- list(path = path, sep = sep, n_comment_lines = n_comment, n_banner_rows = n_banner,
             columns_read = ncol(df), rows_read = nrow(df),
             featurecounts = grepl("^# Program:featureCounts", head_lines[1]))
  lines <- character(0)
  na <- character(0)
  verdict <- "PERMITTED"
  if (n_banner > 0)
    lines <- c(lines, sprintf("%d row%s above the header %s skipped (blank or a title, not part of the table)",
                              n_banner, if (n_banner > 1) "s" else "", if (n_banner > 1) "were" else "was"))

  # --- the gene identifier ------------------------------------------------
  ids <- NULL
  id_index <- integer(0)
  nm  <- tolower(trimws(names(df)))
  non_numeric <- which(!vapply(df, is.numeric, logical(1)))
  if (length(non_numeric)) {
    j <- non_numeric[1]
    id_index <- j
    ids <- as.character(df[[j]])
    ev$id_column <- if (at_unnamed(names(df)[j])) "(first column)" else names(df)[j]
    lines <- c(lines, if (at_unnamed(names(df)[j]))
      "gene identifiers read from the first column, which has no header (the shape write.csv() leaves)"
      else sprintf("gene identifiers read from column '%s'", names(df)[j]))
  } else if (!is.null(rownames(df)) && !identical(rownames(df), as.character(seq_len(nrow(df))))) {
    ids <- rownames(df)                       # read.table's one-short-header case
    ev$id_column <- "(row names)"
    lines <- c(lines, "gene identifiers read from the first column (the header row is one field short, which is how R marks row names)")
  } else if (nm[1] %in% at_id_names()) {
    ids <- as.character(df[[1]])
    id_index <- 1L
    ev$id_column <- names(df)[1]
    lines <- c(lines, sprintf("gene identifiers read from column '%s'", names(df)[1]))
  } else {
    na <- c(na, "gene identifiers: every column is numeric and none is named like an identifier, so rows are numbered; the sex and marker checks need real gene names")
  }

  # --- annotation columns masquerading as samples --------------------------
  drop <- id_index
  ann <- which(nm %in% at_annotation_names())
  ann <- setdiff(ann, drop)
  if (length(ann)) {
    ev$annotation_columns <- names(df)[ann]
    fc <- all(c("chr", "start", "end", "strand", "length") %in% nm)
    lines <- c(lines, sprintf(
      "%s dropped as annotation, not samples%s",
      paste(sprintf("'%s'", names(df)[ann]), collapse = ", "),
      if (fc) " - the file is featureCounts output, whose Start, End and Length columns are numeric and are otherwise counted as three extra libraries" else ""))
    drop <- c(drop, ann)
  }

  keep <- setdiff(seq_len(ncol(df)), drop)
  keep <- keep[vapply(df[keep], is.numeric, logical(1))]
  dropped_text <- setdiff(setdiff(seq_len(ncol(df)), drop), keep)
  if (length(dropped_text)) {
    ev$non_numeric_columns <- names(df)[dropped_text]
    lines <- c(lines, sprintf("%s dropped: not numeric",
                              paste(sprintf("'%s'", names(df)[dropped_text]), collapse = ", ")))
  }
  if (length(keep) < 2)
    return(fail(sprintf("Only %d numeric column%s left after removing identifiers and annotation, so there is nothing to compare.",
                        length(keep), if (length(keep) == 1) "" else "s"),
                c(ev, list(kept = names(df)[keep]))))

  m <- as.matrix(df[, keep, drop = FALSE])
  storage.mode(m) <- "double"
  ev$n_samples <- ncol(m)
  ev$n_genes   <- nrow(m)
  lines <- c(lines, sprintf("%s genes x %d samples read from a %s-separated file%s",
                            format(nrow(m), big.mark = ","), ncol(m),
                            if (sep == "\t") "tab" else sprintf("'%s'", sep),
                            if (n_comment) sprintf(", after %d comment line%s", n_comment,
                                                   if (n_comment > 1) "s" else "") else ""))

  # --- rows that are not genes: blank rows, and values with no identifier ----
  # Spreadsheets often end in empty rows and a totals row whose label sits in
  # another column. A totals row has numbers, so it would be analysed as a gene.
  empty <- rowSums(!is.na(m)) == 0
  if (any(empty)) {
    ev$n_empty_rows <- sum(empty)
    lines <- c(lines, sprintf("%d empty row%s dropped", sum(empty), if (sum(empty) > 1) "s" else ""))
    m <- m[!empty, , drop = FALSE]
    if (!is.null(ids)) ids <- ids[!empty]
    ev$n_genes <- nrow(m)
  }
  if (!is.null(ids)) {
    no_id <- is.na(ids) | !nzchar(trimws(ids))
    if (any(no_id)) {
      ev$no_id_rows <- which(no_id)
      verdict <- "NOT PERMITTED"
      lines <- c(lines, sprintf("%d row%s %s numbers but no gene identifier (rows %s) - typically a totals row at the foot of a spreadsheet; it would be analysed as a gene. Remove %s.",
                                sum(no_id), if (sum(no_id) > 1) "s" else "", if (sum(no_id) > 1) "have" else "has",
                                paste(utils::head(which(no_id), 5), collapse = ", "),
                                if (sum(no_id) > 1) "them" else "it"))
    }
  }
  if (anyNA(m)) {
    ev$n_na_rows <- sum(rowSums(is.na(m)) > 0)
    verdict <- at_worst(c(verdict, "CAUTION"))
    lines <- c(lines, sprintf("%d row%s missing values; checks on the values cannot run until they are removed or filled",
                              ev$n_na_rows, if (ev$n_na_rows > 1) "s contain" else " contains"))
  }

  # --- duplicated identifiers ---------------------------------------------
  if (!is.null(ids)) {
    dup <- duplicated(ids)
    ev$n_duplicate_ids <- sum(dup)
    if (any(dup)) {
      verdict <- at_worst(c(verdict, "CAUTION"))
      lines <- c(lines, sprintf("%s identifier%s appear more than once (first: %s); they are kept as separate rows, made unique",
                                format(sum(dup), big.mark = ","),
                                if (sum(dup) > 1) "s" else "",
                                paste(utils::head(unique(ids[dup]), 3), collapse = ", ")))
      ids <- make.unique(ids)
    }
    rownames(m) <- ids
  }

  n_dup <- if (is.null(ev$n_duplicate_ids)) 0 else ev$n_duplicate_ids
  n_noid <- length(ev$no_id_rows)
  headline <- if (n_noid)
    "The file was read, but it has rows of numbers with no gene identifier."
  else if (n_dup)
    "The file was read, but some gene identifiers are duplicated."
  else if (!is.null(ev$n_na_rows))
    "The file was read, but some values are missing."
  else if (length(ann))
    "The file was read as counts, with its annotation columns removed."
  else "The file was read as counts."
  cost <- c(
    if (n_noid) paste("A row without an identifier is still a row of numbers: every tool downstream treats it as",
                      "a gene. A totals row is thousands of times larger than any gene, so it dominates",
                      "normalisation and the variance of every sample."),
    if (n_dup) paste("Duplicated identifiers usually mean the file was written from a transcript- or exon-level",
                     "table without aggregating, or that symbols were used where two Ensembl genes share one.",
                     "Counts for the same gene end up split across rows, which lowers every one of them below",
                     "the filters and the tests."))

  list(counts = m,
       check  = at_result(
         verdict, headline, lines,
         if (length(cost)) paste(cost, collapse = " ") else NULL,
         na, ev))
}
