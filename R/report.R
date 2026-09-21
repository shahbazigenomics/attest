#' Audit expression data before analysing it
#'
#' Runs every available check on a count matrix and returns one report: an
#' overall verdict, each check's verdict with its evidence, and everything the
#' inputs could not answer. The overall verdict comes from the checks that look
#' for faults; `detectability` reports scope (what the data can support) and is
#' shown separately without changing it. Individual checks stay callable on their own
#' (`attest_counts()`); this is the one call that runs them all.
#'
#' If `x` is a `DESeqDataSet`, a `SummarizedExperiment` or a `DGEList`, the
#' sample sheet, the design formula and the grouping are taken from the object,
#' so `attest(dds)` runs all five checks with nothing else supplied. Anything
#' passed explicitly wins over what the object carries. A check that still has
#' no inputs is named under "not run" rather than quietly left out of the report.
#'
#' Checks currently implemented:
#' * value scale - is this a raw count matrix? (`attest_counts()`)
#' * identifiers - is every row a gene, named once, from one annotation?
#'   (`attest_identifiers()`; needs row names)
#' * completeness - has the matrix been filtered before you got it?
#'   (`attest_completeness()`)
#' * identity - do the samples match their labels, and is any library present
#'   twice? (`attest_identity()`)
#' * design - can the design estimate the effect, and at what effective sample
#'   size? (`attest_design()`; needs a sample sheet and a design formula)
#' * detectability - for which genes could a change have been seen at all?
#'   (`attest_detectability()`; needs the two groups being compared)
#'
#' @param x counts: matrix, data.frame, SummarizedExperiment, DESeqDataSet or
#'   DGEList. For counts in a file, see [attest_file()].
#' @param ... passed to individual checks: `n_expected` (completeness),
#'   `metadata`, `sex_col`, `dup_floor` (identity), `design`, `of_interest`
#'   (design), `group`, `target_fc`, `power` (detectability).
#' @return object of class "attest_report".
#' @export
attest <- function(x, ...) {
  args <- list(...)
  from <- at_from_object(x)

  metadata <- if (!is.null(args$metadata)) as.data.frame(args$metadata) else from$metadata
  design   <- if (!is.null(args$design))   args$design                   else from$design
  of       <- at_of_interest(design, metadata, args$of_interest)
  group    <- args$group
  if (is.null(group)) group <- at_two_level_group(of, metadata)
  if (is.null(group)) group <- from$group

  taken <- character(0)
  if (is.null(args$metadata) && !is.null(from$metadata))
    taken <- c(taken, "sample sheet")
  if (is.null(args$design) && !is.null(from$design))
    taken <- c(taken, "design formula")
  if (is.null(args$group) && is.null(args$design) && !is.null(from$group))
    taken <- c(taken, "grouping")
  source_note <- if (length(taken) && !is.null(from$source))
    sprintf("%s taken from the %s", paste(taken, collapse = " and "), from$source) else NULL

  id_args <- args[names(args) %in% c("sex_col", "dup_floor")]
  checks <- list("value scale" = attest_counts(x))
  not_run <- character(0)

  # named rows are what the identifier check works on; numbered rows are a
  # missing input, not a fault, so they belong in "not run"
  mat <- at_as_matrix(x)
  if (!is.null(rownames(mat))) {
    checks[["identifiers"]] <- attest_identifiers(x)
  } else {
    not_run <- c(not_run, paste(
      "identifiers (is every row a gene, named once, from one annotation?):",
      if (is.null(mat)) sprintf("the counts could not be read from this %s", class(x)[1])
      else "the matrix has no row names - read the counts with attest_file(), or set rownames(x)"))
  }

  checks[["completeness"]] <- do.call(attest_completeness,
                                      c(list(x), args[names(args) %in% "n_expected"]))
  checks[["identity"]] <- do.call(attest_identity,
                                  c(list(x), list(metadata = metadata), id_args))

  if (!is.null(metadata) && !is.null(design)) {
    checks[["design"]] <- do.call(
      attest_design,
      c(list(x), list(metadata = metadata, design = design),
        args[names(args) %in% "of_interest"]))
  } else {
    not_run <- c(not_run, sprintf(
      "design adequacy (is the effect estimable, and what is it worth after confounding?): %s",
      if (is.null(metadata))
        "no sample sheet - pass metadata = <data.frame> and design = ~ <terms>, or hand attest() a DESeqDataSet"
      else
        "no design formula - pass design = ~ <terms>"))
  }

  if (!is.null(group)) {
    checks[["detectability"]] <- do.call(
      attest_detectability,
      c(list(x), list(group = group, metadata = metadata, design = design,
                      of_interest = args$of_interest),
        args[names(args) %in% c("power", "alpha", "target_fc")]))
  } else {
    not_run <- c(not_run, paste(
      "detectability (for which genes could a change have been seen at all?):",
      "the two groups being compared are not known - pass group = <factor>,",
      "or a design whose variable of interest has two levels"))
  }

  kinds  <- vapply(checks, at_kind, character(1))
  faults <- vapply(checks[kinds == "fault"], function(c) c$verdict, character(1))
  structure(list(checks = checks,
                 verdict = at_worst(faults),
                 not_run = not_run,
                 source_note = source_note),
            class = "attest_report")
}

at_kind <- function(ch) if (is.null(ch$kind)) "fault" else ch$kind

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
  if (!is.null(x$source_note)) cat("inputs:", x$source_note, "\n")
  kinds <- vapply(x$checks, at_kind, character(1))
  show <- function(nm) {
    ch <- x$checks[[nm]]
    cat("\n[", nm, "] ", ch$verdict, " - ", ch$headline, "\n", sep = "")
    for (e in ch$evidence)     cat("  -", e, "\n")
    if (!is.null(ch$consequence)) cat("  cost:", ch$consequence, "\n")
    for (u in ch$not_assessed) cat("  not assessed:", u, "\n")
  }
  for (nm in names(x$checks)[kinds == "fault"]) show(nm)
  if (any(kinds == "scope")) {
    cat("\n-- what this data can support (does not change the verdict above) --\n")
    for (nm in names(x$checks)[kinds == "scope"]) show(nm)
  }
  if (length(x$not_run)) {
    cat("\n-- not run, so the verdict above does not cover it --\n")
    for (u in x$not_run) cat("  -", u, "\n")
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
    not_run = as.character(x$not_run),
    checks  = lapply(x$checks, function(ch) list(
      verdict      = ch$verdict,
      kind         = at_kind(ch),
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
