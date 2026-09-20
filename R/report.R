#' Audit expression data before analysing it
#'
#' Runs every available check on a count matrix and returns one report: an
#' overall verdict, each check's verdict with its evidence, and everything the
#' inputs could not answer. Individual checks stay callable on their own
#' (`attest_counts()`); this is the one call that runs them all.
#'
#' Checks currently implemented:
#' * value scale - is this a raw count matrix? (`attest_counts()`)
#' * completeness - has the matrix been filtered before you got it?
#'   (`attest_completeness()`)
#' * identity - do the samples match their labels, and is any library present
#'   twice? (`attest_identity()`)
#' * design - can the design estimate the effect, and at what effective sample
#'   size? (`attest_design()`; runs only when `metadata` and `design` are given)
#'
#' Planned, not yet implemented: detectability of negative claims (for which genes
#' could a change have been seen at all).
#'
#' @param x counts: matrix, data.frame, SummarizedExperiment, DESeqDataSet or
#'   DGEList.
#' @param ... passed to individual checks: `n_expected` (completeness),
#'   `metadata`, `sex_col`, `dup_floor` (identity), `design`, `of_interest`
#'   (design).
#' @return object of class "attest_report".
#' @export
attest <- function(x, ...) {
  args   <- list(...)
  checks <- list(
    "value scale"  = attest_counts(x),
    "completeness" = do.call(attest_completeness, c(list(x), args[names(args) %in% "n_expected"])),
    "identity"     = do.call(attest_identity,
                             c(list(x), args[names(args) %in% c("metadata", "sex_col", "dup_floor")])))
  if (!is.null(args$design) && !is.null(args$metadata)) {
    checks[["design"]] <- do.call(
      attest_design,
      c(list(x), args[names(args) %in% c("metadata", "design", "of_interest", "min_effective_frac")]))
  }
  structure(list(checks = checks,
                 verdict = at_worst(vapply(checks, function(c) c$verdict, character(1)))),
            class = "attest_report")
}

# severity order: the report is as strong as its weakest check
at_worst <- function(v) {
  order_ <- c("NOT PERMITTED", "CAUTION", "UNKNOWN", "PERMITTED")
  for (o in order_) if (any(v == o)) return(o)
  "UNKNOWN"
}

#' @export
print.attest_report <- function(x, ...) {
  cat("attest - can this data support the analysis?\n")
  cat("overall:", x$verdict, "\n")
  for (nm in names(x$checks)) {
    ch <- x$checks[[nm]]
    cat("\n[", nm, "] ", ch$verdict, " - ", ch$headline, "\n", sep = "")
    for (e in ch$evidence)     cat("  -", e, "\n")
    if (!is.null(ch$consequence)) cat("  cost:", ch$consequence, "\n")
    for (u in ch$not_assessed) cat("  not assessed:", u, "\n")
  }
  invisible(x)
}

#' Report as a machine-readable list
#'
#' The same content as the printed report, shaped for `jsonlite::toJSON()` or
#' for a pipeline step that acts on the verdict.
#'
#' @param x an "attest_report" or a single check.
#' @export
attest_as_list <- function(x) {
  if (inherits(x, "attest_check")) x <- structure(
    list(checks = list("value scale" = x), verdict = x$verdict), class = "attest_report")
  list(
    schema  = "attest/report/1",
    verdict = x$verdict,
    checks  = lapply(x$checks, function(ch) list(
      verdict      = ch$verdict,
      headline     = ch$headline,
      evidence     = as.character(ch$evidence),
      consequence  = if (is.null(ch$consequence)) NA_character_ else ch$consequence,
      not_assessed = as.character(ch$not_assessed),
      measurements = ch$measurements[c("dim", "col_sum_spread", "rel_dev_1e6",
                                       "frac_noninteger", "frac_all_zero_rows",
                                       "sf_spread", "dispersion_index")]
    ))
  )
}
