#' Can this design support the comparison you want?
#'
#' Adjudicates the claim "this design can estimate the effect I am asking for".
#' Three questions: is every coefficient estimable, how much of the effect of
#' interest is shared with the other terms, and are there enough replicates.
#'
#' DESeq2 and edgeR already refuse a design that is *completely* confounded (the
#' model matrix is not full rank) and limma-voom only warns. None of them says
#' anything about partial confounding or thin replication, which is what this
#' check is for.
#'
#' @param x counts: matrix, data.frame, SummarizedExperiment, DESeqDataSet or
#'   DGEList. Used for the sample count and names.
#' @param metadata data.frame of sample annotation, one row per column of `x`.
#' @param design a formula, e.g. `~ batch + condition`.
#' @param of_interest name of the variable whose effect you care about. Defaults
#'   to the last term in `design`.
#' @param min_effective_frac fraction of the samples the effect of interest must
#'   still be worth after confounding. Default 0.5: below that, half the
#'   experiment has been spent on separating the effect from the other terms.
#' @return object of class "attest_check".
#' @export
attest_design <- function(x, metadata, design, of_interest = NULL,
                          min_effective_frac = 0.5) {

  m <- at_as_matrix(x)
  if (is.null(m) || is.null(metadata) || missing(design)) {
    return(at_result("UNKNOWN",
                     "Design check needs counts, a sample sheet and a design formula.",
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

  # --- 2. how much of the effect is shared with the other terms? ----------
  cols <- at_coefficient_columns(mm, of_interest, metadata)
  if (!length(cols)) {
    na <- c(na, sprintf("confounding check: no coefficient in the model matrix comes from '%s'", of_interest))
  } else {
    vifs <- vapply(cols, function(j) at_vif(mm, j), numeric(1))
    ev$vif <- vifs
    ev$max_vif <- max(vifs)
    ev$effective_n <- ev$n_samples / ev$max_vif
    if (ev$max_vif > 1.05) {
      lines <- c(lines, sprintf("the effect of '%s' is partly shared with the other terms: standard errors are inflated %.2fx, so %d samples buy the precision of %.1f",
                                of_interest, sqrt(ev$max_vif), ev$n_samples, ev$effective_n))
    } else {
      lines <- c(lines, sprintf("'%s' is balanced against the other terms (variance inflation %.2fx, no effective loss of samples)",
                                of_interest, ev$max_vif))
    }
  }

  # --- 3. replication ------------------------------------------------------
  grp <- metadata[[of_interest]]
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

  verdict  <- "PERMITTED"
  headline <- sprintf("The design can estimate the effect of '%s'.", of_interest)
  consequence <- NULL

  if (isTRUE(ev$effective_n < min_effective_frac * ev$n_samples)) {
    verdict  <- "CAUTION"
    headline <- sprintf("The design estimates '%s', but most of the experiment is spent separating it from the other terms.",
                        of_interest)
    consequence <- "Nothing here is wrong, but the comparison is weaker than the sample count suggests, and the reported effect carries some of the other term with it. Worth stating in the methods."
  }
  if (ev$min_group < 3) {
    verdict  <- if (verdict == "PERMITTED") "CAUTION" else verdict
    lines    <- c(lines, "fewer than 3 samples in a group: dispersion is estimated from very little, and results from designs this small often do not replicate")
    consequence <- paste(c(consequence,
      "With 2 samples per group the tools run silently and report a gene list; published work on replicability finds such lists rarely hold up."), collapse = " ")
    if (headline == sprintf("The design can estimate the effect of '%s'.", of_interest))
      headline <- sprintf("The design estimates '%s', but on thin replication.", of_interest)
  }

  na <- c(na, "replicability of the resulting gene list: not estimated here; it needs the differential-expression run itself (bootstrap resampling of the samples)")

  at_result(verdict, headline, lines, consequence, na, ev)
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
