#' Is this a raw count matrix?
#'
#' Adjudicates one claim: "this matrix is valid input for count-based
#' differential expression (DESeq2 / edgeR / limma-voom)". It reports a verdict,
#' the evidence behind it, and any check that could not be run on these inputs.
#'
#' The function never stops on bad input; it returns a typed failure instead.
#'
#' @param x numeric matrix (genes x samples), data.frame, or an object carrying
#'   an assay: SummarizedExperiment / DESeqDataSet / DGEList.
#' @param tol_1e6 relative tolerance for "column sums pinned at 1e6".
#' @param sf_spread_max size-factor spread (max/min) below which a matrix looks
#'   already normalised. Default 1.01, set from simulation: normalised matrices
#'   (>=1M reads) land at 1.002-1.016, raw counts at >=1.012 even when library
#'   sizes differ by only 0.5%. Shallow (<1M) normalised matrices can exceed it
#'   and are missed.
#' @return object of class "attest_check": verdict, headline, evidence (text),
#'   consequence, not_assessed, measurements.
#' @export
attest_counts <- function(x, tol_1e6 = 0.01, sf_spread_max = 1.01) {

  m <- at_as_matrix(x)
  if (is.null(m)) {
    return(at_result("UNKNOWN", "The input could not be read as a numeric matrix.",
                     character(0), NULL, character(0), list(input_class = class(x)[1])))
  }
  if (nrow(m) < 10 || ncol(m) < 2) {
    return(at_result("UNKNOWN",
                     sprintf("Too small to judge: %d genes x %d samples (need at least 10 genes and 2 samples).",
                             nrow(m), ncol(m)),
                     character(0), NULL, character(0), list(dim = dim(m))))
  }
  if (anyNA(m)) {
    return(at_result("UNKNOWN", "The matrix contains missing values (NA).",
                     character(0), NULL, character(0), list(n_na = sum(is.na(m)))))
  }

  ev <- list()
  ev$dim                <- dim(m)
  ev$frac_negative      <- mean(m < 0)
  ev$max_value          <- max(m)
  ev$frac_noninteger    <- mean(abs(m - round(m)) > 1e-8)
  ev$col_sums           <- colSums(m)
  ev$col_sum_spread     <- max(ev$col_sums) / max(min(ev$col_sums), 1)
  ev$rel_dev_1e6        <- max(abs(ev$col_sums / 1e6 - 1))
  ev$frac_all_zero_rows <- mean(rowSums(m) == 0)
  ev$size_factors       <- at_size_factors(m)
  ev$sf_spread          <- if (all(is.finite(ev$size_factors)))
    max(ev$size_factors) / max(min(ev$size_factors), .Machine$double.eps) else NA_real_
  ev$lattice            <- at_lattice(m)
  ev$sf_vs_colsum       <- if (is.finite(ev$sf_spread) && ev$col_sum_spread > 1.001)
    log(ev$sf_spread) / log(ev$col_sum_spread) else NA_real_
  ev$dispersion_index   <- at_dispersion_index(m)

  # checks that these inputs cannot answer - reported, never silently skipped
  na <- character(0)
  if (!is.finite(ev$sf_spread)) {
    na <- c(na, sprintf(
      "prior-normalisation check: only %d genes are non-zero in every sample (need 10), so size factors could not be estimated",
      sum(rowSums(m > 0) == ncol(m))))
  }
  if (!is.finite(ev$dispersion_index)) {
    na <- c(na, sprintf(
      "variance/mean (Poisson floor) check: needs at least 3 samples and 200 expressed genes; this matrix has %d samples and %d expressed genes",
      ncol(m), sum(rowMeans(m) > 0)))
  }

  # --- log scale -----------------------------------------------------------
  if (ev$frac_negative > 0 || (ev$frac_noninteger > 0 && ev$max_value < 30)) {
    return(at_result(
      "NOT PERMITTED",
      "These are log-transformed values (vst, rlog or log-CPM), not counts.",
      sprintf("the largest value is %.1f and %.1f%% of values are negative; counts are whole numbers and never negative",
              ev$max_value, 100 * ev$frac_negative),
      "Count models assume counts. Fed log values, they are meaningless.",
      na, ev))
  }

  # --- scaled to a fixed total (CPM / TPM) ---------------------------------
  if (ev$rel_dev_1e6 <= tol_1e6) {
    lines <- sprintf("every sample totals about 1,000,000 (%s), which happens when each column is divided by its own library size",
                     if (all(ev$col_sums <= 1e6)) "all slightly below, so the values were scaled down then rounded"
                     else if (all(ev$col_sums >= 1e6)) "all slightly above, so the values were scaled up then rounded"
                     else "all at 1e6")
    if (ev$frac_all_zero_rows > 0.6) {
      lines <- c(lines, sprintf("%.0f%% of genes are zero in every sample, so rounding has already destroyed the low-expression genes",
                                100 * ev$frac_all_zero_rows))
    }
    if (!is.null(ev$lattice$factor_estimate)) {
      lines <- c(lines, sprintf("small values are spaced about %.2f apart in %d of %d samples, implying original libraries of roughly %s reads",
                                ev$lattice$factor_estimate, ev$lattice$n_agree, ncol(m),
                                format(round(1e6 / ev$lattice$factor_estimate), big.mark = ",", scientific = FALSE)))
    }
    return(at_result("NOT PERMITTED", "These are CPM or TPM values, not raw counts.",
                     lines, at_consequence(), na, ev))
  }

  # --- non-integer: FPKM-like, or estimated counts -------------------------
  if (ev$frac_noninteger > 0) {
    sub_poisson   <- isTRUE(ev$dispersion_index < 0.8)
    depth_removed <- isTRUE(ev$col_sum_spread < 1.1)
    if (sub_poisson || depth_removed) {
      lines <- character(0)
      if (sub_poisson) {
        lines <- c(lines, sprintf("variance across samples is smaller than the mean (%.3f); counting events cannot behave that way, so the values were divided by depth and/or gene length",
                                  ev$dispersion_index))
      }
      if (depth_removed) {
        lines <- c(lines, sprintf("sample totals differ by only %.2fx, so sequencing depth has already been divided out",
                                  ev$col_sum_spread))
      }
      return(at_result("NOT PERMITTED",
                       "These are FPKM/RPKM or another normalised measure, not raw counts.",
                       lines,
                       paste(at_consequence(),
                             "If they are salmon/kallisto estimates from an unusually even run, import them with tximport instead of as a matrix."),
                       na, ev))
    }
    return(at_result(
      "CAUTION",
      "These look like estimated counts from salmon or kallisto.",
      c(sprintf("%.0f%% of values are not whole numbers", 100 * ev$frac_noninteger),
        sprintf("but their size is count-like (variance/mean %.2f, sample totals differ %.2fx)",
                ev$dispersion_index, ev$col_sum_spread)),
      "Import them with tximport, which passes the right offsets to DESeq2/edgeR. Rounding them by hand discards that information.",
      na, ev))
  }

  # --- already normalised --------------------------------------------------
  # Three ways a normalised matrix shows itself. A fixed threshold on the size-factor
  # spread alone is depth-dependent: rounding noise put the 1M-read airway fixture at
  # 1.0154 and the 22M-read full matrix at 1.0057. The ratio below is scale-free -
  # measured 0.003-0.072 for normalised matrices against 1.0-1.3 for raw counts.
  flat_factors <- is.finite(ev$sf_spread) && ev$sf_spread < sf_spread_max
  crushed      <- isTRUE(ev$sf_vs_colsum < 0.3) && isTRUE(ev$sf_spread < 1.25)
  flat_totals  <- ev$col_sum_spread < 1.01
  if (flat_factors || crushed || flat_totals) {
    why <- if (flat_totals)
      sprintf("every sample totals the same to within %.2f%%, which sequencing does not produce: the values were either scaled to a common total (then they are not counts) or the libraries were downsampled to equal depth (then they are counts, but reads were thrown away)",
              100 * (ev$col_sum_spread - 1))
    else if (crushed)
      sprintf("size factors vary only %.3fx while sample totals vary %.2fx; in raw counts the two move together, so a ratio this small means depth was divided out and composition left behind",
              ev$sf_spread, ev$col_sum_spread)
    else
      sprintf("median-of-ratios size factors vary only %.3fx across samples, where raw libraries differ by more even in a tightly balanced run",
              ev$sf_spread)
    return(at_result(
      "CAUTION",
      "This matrix appears to have been normalised already.",
      why,
      "Normalising twice makes samples look more alike than they are, which slightly inflates significance. Use the original counts if you still have them.",
      na, ev))
  }

  at_result(
    "PERMITTED",
    "Consistent with raw counts.",
    c(if (is.finite(ev$sf_spread))
        sprintf("all values are whole numbers; sample totals differ %.2fx and size factors %.2fx, as expected when depth has not been divided out",
                ev$col_sum_spread, ev$sf_spread)
      else
        sprintf("all values are whole numbers and sample totals differ %.2fx, as expected when depth has not been divided out",
                ev$col_sum_spread),
      if (is.finite(ev$dispersion_index))
        sprintf("variance/mean is %.2f, above the Poisson floor", ev$dispersion_index)),
    NULL, na, ev)
}

at_consequence <- function() {
  "DESeq2, edgeR and limma will run on this without complaint and report far fewer differences than the data contain (3,993 -> 877 genes in our airway test)."
}

# ---- helpers ---------------------------------------------------------------

at_as_matrix <- function(x) {
  m <- NULL
  if (is.matrix(x) && is.numeric(x)) {
    m <- x
  } else if (is.data.frame(x)) {
    num <- vapply(x, is.numeric, logical(1))
    if (any(num)) m <- as.matrix(x[, num, drop = FALSE])
  } else if (!is.null(attr(class(x), "package")) || isS4(x)) {
    m <- tryCatch({
      if (methods::existsMethod("assay", class(x))) as.matrix(SummarizedExperiment::assay(x))
      else if (!is.null(x$counts)) as.matrix(x$counts)  # DGEList
      else NULL
    }, error = function(e) NULL)
  } else if (is.list(x) && !is.null(x$counts)) {
    m <- as.matrix(x$counts)
  }
  if (is.null(m) || !is.numeric(m)) NULL else m
}

at_result <- function(verdict, headline, evidence, consequence, not_assessed, measurements) {
  structure(list(verdict = verdict,
                 headline = headline,
                 evidence = evidence[!vapply(evidence, is.null, logical(1))],
                 consequence = consequence,
                 not_assessed = not_assessed,
                 measurements = measurements),
            class = "attest_check")
}

# Poisson floor: counts have variance >= mean; values divided by depth or gene
# length become sub-Poisson. Mid-expression genes only.
at_dispersion_index <- function(m) {
  if (ncol(m) < 3) return(NA_real_)
  rm_ <- rowMeans(m)
  pos <- rm_[rm_ > 0]
  if (length(pos) < 200) return(NA_real_)
  qs <- stats::quantile(pos, c(0.2, 0.8))
  keep <- rm_ > qs[1] & rm_ < qs[2]
  y <- m[keep, , drop = FALSE]
  if (nrow(y) < 100) return(NA_real_)
  mu <- rowMeans(y)
  n <- ncol(y)
  v <- rowSums((y - mu)^2) / (n - 1)          # vectorised row variance
  stats::median(v / mu)
}

# median-of-ratios size factors (DESeq2's estimator), base R only
at_size_factors <- function(m) {
  keep <- rowSums(m > 0) == ncol(m)
  if (sum(keep) < 10) return(rep(NA_real_, ncol(m)))
  lm_ <- log(m[keep, , drop = FALSE])
  ref <- rowMeans(lm_)
  apply(lm_, 2, function(col) exp(stats::median(col - ref)))
}

# Scaled-and-rounded data land on a lattice: the smallest distinct values are
# multiples of the scale factor. Checked per sample; claimed only if >=2 agree.
at_lattice <- function(m, n_small = 30, tol = 0.15) {
  per_sample <- apply(m, 2, function(col) {
    occ <- sort(unique(col[col > 0]))
    if (length(occ) < 6) return(NA_real_)
    occ <- occ[seq_len(min(n_small, length(occ)))]   # the smallest distinct values
    g <- stats::median(diff(occ))
    if (g > 1.5) g else NA_real_
  })
  ok <- per_sample[is.finite(per_sample)]
  if (length(ok) < 2) return(list(per_sample = per_sample, factor_estimate = NULL, n_agree = length(ok)))
  med <- stats::median(ok)
  agree <- sum(abs(ok - med) <= tol * med)
  if (agree < 2) return(list(per_sample = per_sample, factor_estimate = NULL, n_agree = agree))
  list(per_sample = per_sample, factor_estimate = med, n_agree = agree)
}

#' @export
print.attest_check <- function(x, ...) {
  cat("attest check\n")
  cat("verdict:", x$verdict, "-", x$headline, "\n")
  if (length(x$evidence)) {
    cat("why:\n"); for (e in x$evidence) cat("  -", e, "\n")
  }
  if (!is.null(x$consequence)) cat("what it costs:", x$consequence, "\n")
  if (length(x$not_assessed)) {
    cat("not assessed:\n"); for (u in x$not_assessed) cat("  -", u, "\n")
  }
  invisible(x)
}
