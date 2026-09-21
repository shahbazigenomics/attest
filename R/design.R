#' Can this design support the comparison you want?
#'
#' Adjudicates the claim "this design can estimate the effect I am asking for,
#' and nothing in the sample sheet stands in its way". Four questions:
#'
#' 1. Is every coefficient estimable? (fault)
#' 2. Is every group replicated? (fault)
#' 3. Does a column the formula leaves out line up with the effect? (fault) -
#'    every treated sample processed in one batch, or several libraries from
#'    each donor counted as independent replicates. No differential-expression
#'    tool sees this, because each looks only at the terms it was given.
#' 4. How much of the effect is shared with the other terms? Reported as an
#'    effective sample size. This is **scope**, not a fault: the effect is still
#'    estimated without bias, only less precisely, so it does not lower the
#'    verdict - it is carried into [attest_detectability()], which is where
#'    precision is judged.
#'
#' DESeq2 and edgeR already refuse a design whose own terms are completely
#' confounded; limma-voom only warns. None of them looks at the sample sheet
#' beyond the formula. BatchQC measures confounding between a batch and the
#' condition, but only for a batch column the user has already named - which is
#' the step that fails when nobody noticed there was one. This scans every column.
#'
#' @param x counts: matrix, data.frame, SummarizedExperiment, DESeqDataSet or
#'   DGEList. Used for the sample count and names.
#' @param metadata data.frame of sample annotation, one row per column of `x`.
#' @param design a formula, e.g. `~ batch + condition`.
#' @param of_interest name of the variable whose effect you care about. Defaults
#'   to the last term in `design`.
#' @return object of class "attest_check".
#' @examples
#' fx <- readRDS(system.file("extdata", "fixtures.rds", package = "attest"))
#' cond <- factor(rep(c("ctrl", "trt"), each = 4))
#' attest_design(fx$airway, data.frame(cond = cond), ~ cond)
#' # two libraries per donor, analysed as independent replicates
#' donor <- data.frame(cond = cond, donor = rep(paste0("d", 1:4), each = 2))
#' attest_design(fx$airway, donor, ~ cond)
#' @export
attest_design <- function(x, metadata, design, of_interest = NULL) {

  m <- at_as_matrix(x)
  if (is.null(m)) {
    return(at_result("UNKNOWN", "The input could not be read as a numeric matrix.",
                     character(0), NULL, character(0), list(input_class = class(x)[1])))
  }
  if (is.null(metadata) || missing(design) || !inherits(design, "formula")) {
    return(at_result("UNKNOWN",
                     "Design check needs a sample sheet and a design formula.",
                     character(0), NULL, character(0), list()))
  }
  metadata <- as.data.frame(metadata)
  if (nrow(metadata) != ncol(m)) {
    return(at_result("UNKNOWN",
                     sprintf("The sample sheet has %d rows but the matrix has %d samples.",
                             nrow(metadata), ncol(m)),
                     character(0), NULL, character(0),
                     list(n_meta = nrow(metadata), n_samples = ncol(m))))
  }

  terms_ <- all.vars(design)
  missing_vars <- setdiff(terms_, names(metadata))
  if (length(missing_vars)) {
    return(at_result("UNKNOWN",
                     sprintf("The design mentions %s, which %s not in the sample sheet.",
                             paste(missing_vars, collapse = ", "),
                             if (length(missing_vars) > 1) "are" else "is"),
                     character(0), NULL, character(0), list(missing = missing_vars)))
  }
  if (is.null(of_interest)) of_interest <- terms_[length(terms_)]

  mm <- tryCatch(stats::model.matrix(design, metadata), error = function(e) NULL)
  if (is.null(mm)) {
    return(at_result("UNKNOWN", "The design formula could not be turned into a model matrix.",
                     character(0), NULL, character(0), list()))
  }

  ev <- list(n_samples = ncol(m), design = deparse(design), of_interest = of_interest,
             rank = qr(mm)$rank, n_coefficients = ncol(mm))
  na <- character(0)
  lines <- character(0)

  # --- 1. is every coefficient estimable? ---------------------------------
  if (ev$rank < ev$n_coefficients) {
    keep <- qr(mm)$pivot[seq_len(ev$rank)]
    lost <- colnames(mm)[-keep]
    ev$not_estimable <- lost
    return(at_result(
      "NOT PERMITTED",
      "This design cannot estimate the effect you asked for.",
      c(sprintf("the model matrix has rank %d for %d coefficients, so %s cannot be estimated",
                ev$rank, ev$n_coefficients, paste(lost, collapse = ", ")),
        sprintf("that happens when two terms carry the same information - typically every sample of one condition sits in one batch")),
      "DESeq2 and edgeR will refuse this design outright; limma-voom fits it anyway and returns NA coefficients with only a warning. No analysis can separate the two terms: the experiment has to be redesigned or a term dropped.",
      na, ev))
  }

  # --- 4 (reported first). what the confounding costs in precision -------
  cols <- at_coefficient_columns(mm, of_interest, metadata)
  if (!length(cols)) {
    na <- c(na, sprintf("confounding among the design's own terms: no coefficient in the model matrix comes from '%s'", of_interest))
  } else {
    vifs <- vapply(cols, function(j) at_vif(mm, j), numeric(1))
    ev$vif <- vifs
    ev$max_vif <- max(vifs)
    ev$effective_n <- ev$n_samples / ev$max_vif
    if (ev$max_vif > 1.05) {
      lines <- c(lines, sprintf(
        "the effect of '%s' is partly shared with the other terms: standard errors are inflated %.2fx, so %d samples buy the precision of %.1f - still estimated without bias, only less precisely, and detectability is computed at that size",
        of_interest, sqrt(ev$max_vif), ev$n_samples, ev$effective_n))
    } else {
      lines <- c(lines, sprintf("'%s' is balanced against the other terms (variance inflation %.2fx, no effective loss of samples)",
                                of_interest, ev$max_vif))
    }
  }

  verdict  <- "PERMITTED"
  headline <- sprintf("The design can estimate the effect of '%s'.", of_interest)
  consequence <- character(0)

  # --- 2. replication ------------------------------------------------------
  grp <- metadata[[of_interest]]
  categorical <- is.factor(grp) || is.character(grp) || is.logical(grp)
  if (categorical) {
    tab <- table(grp)
    ev$group_sizes <- tab
    ev$min_group <- min(tab)
    if (ev$min_group < 2) {
      return(at_result(
        "NOT PERMITTED",
        sprintf("One level of '%s' has no replication.", of_interest),
        c(lines, sprintf("group sizes: %s", paste(sprintf("%s = %d", names(tab), tab), collapse = ", "))),
        "With a single sample in a group, within-group variability cannot be estimated for it, so a difference cannot be told from noise. The tools will still return a gene list.",
        na, ev))
    }
    lines <- c(lines, sprintf("group sizes: %s", paste(sprintf("%s = %d", names(tab), tab), collapse = ", ")))
    if (ev$min_group < 3) {
      verdict  <- "CAUTION"
      headline <- sprintf("The design estimates '%s', but on thin replication.", of_interest)
      lines    <- c(lines, "fewer than 3 samples in a group: dispersion is estimated from very little")
      consequence <- c(consequence,
        "With 2 samples per group the tools run silently and report a gene list; its dispersion estimates rest almost entirely on the other genes.")
    }
  }

  # --- 3. what the formula leaves out --------------------------------------
  om <- at_omitted_variables(metadata, of_interest, terms_)
  ev$omitted_checked <- om$checked
  ev$omitted_skipped <- om$skipped
  ev$omitted_flags   <- om$flags
  na <- c(na, om$not_assessed)
  if (length(om$flags)) {
    if (verdict == "PERMITTED") {
      verdict  <- "CAUTION"
      headline <- sprintf("A column the design leaves out lines up with '%s'.", of_interest)
    }
    lines <- c(lines, vapply(om$flags, `[[`, character(1), "line"))
    consequence <- c(consequence, paste(
      "The design does not include these columns, so the analysis cannot separate their effect from the",
      "one you are testing, and nothing in the counts says whether they have one. If a column is a",
      "processing variable (batch, run, date, operator), the comparison is confounded with it. If it is a",
      "subject (donor, patient, animal) with several samples each, those samples are not independent",
      "replicates, and the tests will be too confident. Neither reading can be ruled out from here."))
  } else if (length(om$checked)) {
    lines <- c(lines, sprintf("sample-sheet columns outside the design (%s) do not line up with '%s'",
                              paste(om$checked, collapse = ", "), of_interest))
  }

  na <- c(na, "replicability of the resulting gene list: not estimated here; it needs the differential-expression run itself (bootstrap resampling of the samples)")

  at_result(verdict, headline, lines,
            if (length(consequence)) paste(consequence, collapse = " ") else NULL,
            na, ev)
}

# Columns of the sample sheet the design leaves out, tested against the variable
# of interest. Two threshold-free tests, because at n = 8 any association
# measure is mostly noise and a cut-off on it would fire on half of all studies:
#  - categorical: every level of the column falls inside one level of the
#    variable of interest (nested), with some level repeated;
#  - numeric (two-group comparisons only): its ranges in the two groups do not
#    overlap at all.
# Identifiers (all values distinct), constants, exact relabelings of the
# variable of interest, and columns computed from the counts are skipped.
at_omitted_variables <- function(metadata, of_interest, design_vars) {
  out <- list(checked = character(0), skipped = character(0), flags = list(),
              not_assessed = character(0))
  f <- metadata[[of_interest]]
  if (!(is.factor(f) || is.character(f) || is.logical(f))) {
    out$not_assessed <- sprintf("columns outside the design: '%s' is continuous, so there are no groups for them to line up with", of_interest)
    return(out)
  }
  f <- droplevels(as.factor(f))
  n <- length(f)
  derived <- c("sizefactor", "lib.size", "norm.factors", "replaceable", "libsize")

  for (col in setdiff(names(metadata), c(design_vars, of_interest))) {
    v <- metadata[[col]]
    ok <- !is.na(v)
    if (tolower(col) %in% derived)           { out$skipped <- c(out$skipped, col); next }
    if (length(unique(v[ok])) <= 1)           { out$skipped <- c(out$skipped, col); next }

    if (is.numeric(v)) {
      # a row index (1..n in any order) only says the sheet is sorted by
      # condition, which it usually is. Any other numeric column is expected to
      # have all-distinct values, so distinctness is not a reason to skip it.
      u <- sort(unique(v[ok]))
      if (length(u) == sum(ok) && all(u == round(u)) && max(u) - min(u) == length(u) - 1) {
        out$skipped <- c(out$skipped, col); next
      }
      if (nlevels(f) != 2) {
        out$not_assessed <- c(out$not_assessed, sprintf(
          "'%s' against '%s': a numeric column is only tested against a two-group comparison", col, of_interest))
        next
      }
      out$checked <- c(out$checked, col)
      a <- v[ok & f == levels(f)[1]]; b <- v[ok & f == levels(f)[2]]
      if (length(a) && length(b) && (max(a) < min(b) || max(b) < min(a))) {
        chance <- 2 / choose(length(a) + length(b), length(a))
        lo <- if (max(a) < min(b)) levels(f)[1] else levels(f)[2]
        out$flags[[col]] <- list(column = col, kind = "separates", chance = chance, line = sprintf(
          "'%s' separates the groups completely: every %s sample is lower than every %s sample (%s); with these group sizes an unrelated column does that about 1 time in %.0f",
          col, lo, setdiff(levels(f), lo)[1],
          paste(sprintf("%s %s-%s", levels(f),
                        signif(c(min(v[ok & f == levels(f)[1]]), min(v[ok & f == levels(f)[2]])), 3),
                        signif(c(max(v[ok & f == levels(f)[1]]), max(v[ok & f == levels(f)[2]])), 3)),
                collapse = ", "),
          1 / chance))
      }
      next
    }

    # a non-numeric column with a different value on every row is an identifier
    if (length(unique(v[ok])) == sum(ok)) { out$skipped <- c(out$skipped, col); next }
    vv <- droplevels(as.factor(v))
    tab <- table(vv[ok], f[ok])
    nested   <- all(rowSums(tab > 0) == 1)
    relabel  <- nested && all(colSums(tab > 0) == 1)
    if (relabel) { out$skipped <- c(out$skipped, col); next }
    out$checked <- c(out$checked, col)
    if (nested && any(rowSums(tab) >= 2)) {
      per <- vapply(levels(f), function(l) sum(tab[, l] > 0), numeric(1))
      out$flags[[col]] <- list(column = col, kind = "nested", line = sprintf(
        "every level of '%s' falls inside a single level of '%s' (%s)",
        col, of_interest,
        paste(sprintf("%s: %d sample%s from %d %s value%s", levels(f), as.integer(table(f[ok])),
                      ifelse(as.integer(table(f[ok])) == 1, "", "s"), per, col,
                      ifelse(per == 1, "", "s")), collapse = "; ")))
    }
  }
  out
}

# which model-matrix columns came from this variable
at_coefficient_columns <- function(mm, var, metadata) {
  setdiff(grep(paste0("^", var), colnames(mm)), 1L)
}

# variance inflation for one coefficient: how much of it the other terms explain
at_vif <- function(mm, j) {
  others <- mm[, setdiff(seq_len(ncol(mm)), c(1L, j)), drop = FALSE]
  if (!ncol(others)) return(1)
  fit <- stats::lm.fit(cbind(1, others), mm[, j])
  r2  <- 1 - sum(fit$residuals^2) / sum((mm[, j] - mean(mm[, j]))^2)
  if (!is.finite(r2) || r2 >= 1) return(Inf)
  1 / (1 - r2)
}

# the design's variance inflation on the effect of interest, for detectability;
# NULL when it cannot be computed, Inf when the effect is not estimable
at_design_vif <- function(design, metadata, of_interest = NULL) {
  if (is.null(design) || is.null(metadata) || !inherits(design, "formula")) return(NULL)
  md <- as.data.frame(metadata)
  if (!all(all.vars(design) %in% names(md))) return(NULL)
  mm <- tryCatch(stats::model.matrix(design, md), error = function(e) NULL)
  if (is.null(mm) || nrow(mm) != nrow(md)) return(NULL)
  if (qr(mm)$rank < ncol(mm)) return(Inf)
  of <- if (is.null(of_interest)) utils::tail(all.vars(design), 1) else of_interest
  cols <- at_coefficient_columns(mm, of, md)
  if (!length(cols)) return(NULL)
  max(vapply(cols, function(j) at_vif(mm, j), numeric(1)))
}
