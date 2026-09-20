#' Has this matrix been filtered?
#'
#' Adjudicates the claim "this matrix contains every gene the pipeline
#' quantified". A complete matrix carries genes that are zero in every sample;
#' filtering removes them and leaves a floor on the row totals, which also
#' reveals the threshold that was used.
#'
#' Fractions of all-zero genes are strongly species-dependent (airway 47.4%,
#' pasilla 15.3%, fission 4.0%), so the test is their presence, never their
#' proportion.
#'
#' @param x counts: matrix, data.frame, SummarizedExperiment, DESeqDataSet or
#'   DGEList.
#' @param n_expected optional: how many genes the annotation contains. When
#'   given, a much smaller matrix is reported as evidence of gene selection
#'   (e.g. protein-coding only).
#' @return object of class "attest_check".
#' @export
attest_completeness <- function(x, n_expected = NULL) {

  m <- at_as_matrix(x)
  if (is.null(m)) {
    return(at_result("UNKNOWN", "The input could not be read as a numeric matrix.",
                     character(0), NULL, character(0), list(input_class = class(x)[1])))
  }
  if (nrow(m) < 10 || ncol(m) < 2) {
    return(at_result("UNKNOWN",
                     sprintf("Too small to judge: %d genes x %d samples.", nrow(m), ncol(m)),
                     character(0), NULL, character(0), list(dim = dim(m))))
  }
  if (anyNA(m)) {
    return(at_result("UNKNOWN", "The matrix contains missing values (NA).",
                     character(0), NULL, character(0), list(n_na = sum(is.na(m)))))
  }

  rs <- rowSums(m)
  ev <- list(
    dim                = dim(m),
    n_all_zero         = sum(rs == 0),
    frac_all_zero_rows = mean(rs == 0),
    min_row_total      = min(rs),
    n_expected         = if (is.null(n_expected)) NA_integer_ else as.integer(n_expected)
  )

  na <- character(0)
  if (is.null(n_expected)) {
    na <- c(na, "gene-selection check (protein-coding only, a panel, or another subset): needs the annotation size, so pass n_expected = <genes in your GTF>")
  }

  # --- a complete matrix keeps its dead genes ------------------------------
  if (ev$n_all_zero > 0) {
    lines <- sprintf("%s of %s genes (%.1f%%) are zero in every sample; filtering removes exactly those, so they are the sign that nothing was removed",
                     format(ev$n_all_zero, big.mark = ","), format(nrow(m), big.mark = ","),
                     100 * ev$frac_all_zero_rows)
    verdict <- "PERMITTED"
    headline <- "Complete: the matrix still contains genes with no reads at all."
    if (!is.null(n_expected) && nrow(m) < 0.9 * n_expected) {
      lines <- c(lines, sprintf("but it holds %s genes against the %s in the annotation you gave, so a subset of genes was selected",
                                format(nrow(m), big.mark = ","), format(as.integer(n_expected), big.mark = ",")))
      verdict <- "CAUTION"
      headline <- "Some genes were selected, though no expression filter was applied."
    }
    return(at_result(verdict, headline, lines,
                     if (verdict == "CAUTION") at_filter_consequence() else NULL, na, ev))
  }

  # --- no all-zero genes: filtered ----------------------------------------
  lines <- "no gene is zero in every sample, and a complete annotation nearly always contains some"
  if (ev$min_row_total <= 1) {
    lines <- c(lines, "the lowest row total is 1, so the filter was most likely 'drop genes with no reads'")
  } else {
    lines <- c(lines, sprintf("the lowest row total is %s, so a threshold of about that size was applied (a per-sample rule, such as 'at least k reads in at least n samples', leaves a higher and less exact floor)",
                              format(ev$min_row_total, big.mark = ",")))
  }
  if (!is.null(n_expected)) {
    lines <- c(lines, sprintf("%s genes remain of the %s in the annotation you gave (%.0f%%)",
                              format(nrow(m), big.mark = ","), format(as.integer(n_expected), big.mark = ","),
                              100 * nrow(m) / n_expected))
  }
  at_result("CAUTION", "This matrix has been filtered before you received it.",
            lines, at_filter_consequence(), na, ev)
}

at_filter_consequence <- function() {
  paste("Filtering itself is normal and DESeq2 recommends a light version of it.",
        "What matters is that it happened upstream of you: the multiple-testing denominator is",
        "already smaller than the experiment produced, DESeq2's own independent filtering has less",
        "to work with, and a gene absent from the matrix is indistinguishable from a gene that was",
        "never expressed. Ask for the unfiltered matrix if you have to report how many genes were tested.")
}
