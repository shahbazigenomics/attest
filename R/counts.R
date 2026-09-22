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
#' @param sf_spread_max size-factor spread (max/min) below which the size
#'   factors are reported as flat. Flat size factors are supporting evidence
#'   only: raw libraries of near-equal depth have them too.
#' @param aliasing_max jaggedness of the small-value histogram above which a
#'   sample is taken to have been divided by a factor and rounded. Default 1.5,
#'   between the most jagged raw matrix measured (0.87; 0.44 over 180 simulated
#'   ones) and the least jagged rounded normalised matrix (2.93), across three
#'   species, full and 3,000-gene matrices, 20k to 30M reads per sample.
#' @return object of class "attest_check": verdict, headline, evidence (text),
#'   consequence, not_assessed, measurements.
#' @examples
#' fx <- readRDS(system.file("extdata", "fixtures.rds", package = "attest"))
#' attest_counts(fx$airway)                       # raw counts
#' m <- fx$airway
#' attest_counts(round(t(t(m) / colSums(m)) * 1e6)) # rounded CPM
#' @export
attest_counts <- function(x, tol_1e6 = 0.01, sf_spread_max = 1.01, aliasing_max = 1.5) {

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
  if (any(is.infinite(m))) {
    # a real count matrix cannot contain Inf/-Inf; seen on published
    # "supplementary counts" files that are actually differential-expression
    # results tables (a fold-change column is Inf wherever the denominator
    # group is all-zero) - past this point every ratio built from column
    # sums or size factors would otherwise be free to become NaN and crash
    # the first unguarded `if` that tests it, rather than being reported.
    return(at_result("UNKNOWN", "The matrix contains infinite values (Inf/-Inf), so it is not a count matrix.",
                     character(0), NULL, character(0), list(n_infinite = sum(is.infinite(m)))))
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
  ev$aliasing          <- at_aliasing(m)
  ev$aliasing_max      <- if (all(is.na(ev$aliasing))) NA_real_ else max(ev$aliasing, na.rm = TRUE)

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
  # Size factors alone cannot tell a normalised matrix from raw libraries of
  # near-equal depth: in both, the typical gene is level across samples and the
  # totals differ only through composition. Measured: depth-balanced raw airway
  # with 20-70% globin reads was called "normalised" by the size-factor rule.
  # So the verdict needs positive evidence of rescaling: dividing counts by a
  # factor that is not 1 and rounding leaves the small-value histogram jagged
  # (aliasing), which sequencing never does. Most jagged sample per matrix,
  # three species, full and 3,000-gene matrices, 20k-30M reads per sample: raw
  # at most 0.87, rounded normalised at least 2.93 (validation/depth_sweep.R).
  flat_factors <- is.finite(ev$sf_spread) && ev$sf_spread < sf_spread_max
  crushed      <- isTRUE(ev$sf_vs_colsum < 0.3) && isTRUE(ev$sf_spread < 1.25)
  flat_totals  <- ev$col_sum_spread < 1.01
  rescaled     <- isTRUE(ev$aliasing_max > aliasing_max)
  if (flat_totals) {
    return(at_result(
      "CAUTION",
      "This matrix appears to have been normalised already.",
      sprintf("every sample totals the same to within %.2f%%, which sequencing does not produce: the values were either scaled to a common total (then they are not counts) or the libraries were downsampled to equal depth (then they are counts, but reads were thrown away)",
              100 * (ev$col_sum_spread - 1)),
      "Normalising twice makes samples look more alike than they are, which slightly inflates significance. Use the original counts if you still have them.",
      na, ev))
  }
  unmeasured <- all(is.na(ev$aliasing))
  if (unmeasured && (flat_factors || crushed)) {
    # fallback: no sample has enough small values to show rescaling directly,
    # so flat size factors are all there is - and they have two readings
    return(at_result(
      "CAUTION",
      "This matrix may have been normalised already.",
      sprintf("size factors vary only %.3fx while sample totals vary %.2fx - what normalisation leaves, but also what raw libraries of near-equal depth show when their totals differ through a few genes (globin in whole blood, for one)",
              ev$sf_spread, ev$col_sum_spread),
      "If it was normalised, normalising again makes samples look more alike than they are, which slightly inflates significance. If these are raw counts from balanced libraries, nothing is wrong.",
      c(na, "rescaling check (the direct evidence of normalisation): no sample has enough small values to judge, so the two readings above cannot be told apart"),
      ev))
  }
  if (rescaled) {
    worst <- names(which.max(ev$aliasing))
    return(at_result(
      "CAUTION",
      "This matrix appears to have been normalised already.",
      c(sprintf("the frequencies of small values are jagged in %d of %d samples (worst: %s, %.2f; raw counts measured at most 0.87): each sample was divided by its own factor and rounded, so some whole numbers collect values from two counts and others from none",
                sum(ev$aliasing > aliasing_max, na.rm = TRUE), ncol(m),
                if (is.null(worst)) "one sample" else worst, ev$aliasing_max),
        if (flat_factors || crushed)
          sprintf("and the median-of-ratios size factors are flat (%.3fx across samples), which is what normalisation leaves",
                  ev$sf_spread)),
      "Normalising twice makes samples look more alike than they are, which slightly inflates significance. Use the original counts if you still have them.",
      na, ev))
  }

  at_result(
    "PERMITTED",
    "Consistent with raw counts.",
    c(if (is.finite(ev$sf_spread) && !(flat_factors || crushed))
        sprintf("all values are whole numbers; sample totals differ %.2fx and size factors %.2fx, as expected when depth has not been divided out",
                ev$col_sum_spread, ev$sf_spread)
      else if (is.finite(ev$sf_spread))
        "all values are whole numbers"
      else
        sprintf("all values are whole numbers and sample totals differ %.2fx, as expected when depth has not been divided out",
                ev$col_sum_spread),
      if (is.finite(ev$dispersion_index))
        sprintf("variance/mean is %.2f, above the Poisson floor", ev$dispersion_index),
      if (flat_factors || crushed)
        sprintf("size factors are flat (%.3fx) for sample totals that differ %.2fx - what a normalised matrix looks like, but also libraries of near-equal depth whose totals differ through a few genes (globin in whole blood, for one); the small-value histograms show no rescaling, so this reads as the second",
                ev$sf_spread, ev$col_sum_spread)),
    NULL, na, ev)
}

# How jagged is each sample's histogram of small values? Counts from sequencing
# give a smooth, decreasing frequency curve over 1, 2, 3, ...; dividing by a
# factor other than 1 and rounding sends two counts to some whole numbers and
# none to others. Measured as the median absolute second difference of
# log(frequency + 1) over values 1-40 that at least 10 genes hold. NA when a
# sample has too few small values to judge (very deep, or very few genes).
at_aliasing <- function(m) {
  if (any(abs(m - round(m)) > 1e-8)) return(rep(NA_real_, ncol(m)))
  a <- apply(m, 2, function(v) {
    v <- v[v > 0 & v <= 40]
    if (length(v) < 200) return(NA_real_)
    f <- tabulate(v, 40)
    k <- which(f >= 10)
    if (length(k) < 6) return(NA_real_)
    lf <- log(f[seq(min(k), max(k))] + 1)
    stats::median(abs(diff(lf, differences = 2)))
  })
  if (!is.null(colnames(m))) names(a) <- colnames(m)
  a
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
  } else if (methods::is(x, "SummarizedExperiment") &&
             requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    # DESeqDataSet and RangedSummarizedExperiment inherit from this. is() follows
    # inheritance; existsMethod("assay", class(x)) - the first version - does
    # not, so every real DESeqDataSet came back unreadable.
    m <- tryCatch({
      an <- SummarizedExperiment::assayNames(x)
      as.matrix(SummarizedExperiment::assay(x, if ("counts" %in% an) "counts" else 1L))
    }, error = function(e) NULL)
  } else if (is.list(x) || isS4(x)) {
    # DGEList (an S4 class built on a list) and tximport's plain list
    m <- tryCatch(if (!is.null(x$counts)) as.matrix(x$counts) else NULL,
                  error = function(e) NULL)
  }
  if (is.null(m) || !is.numeric(m)) NULL else m
}

# `kind` says what sort of statement a check makes:
#   "fault" - something may be wrong with the data; these set the overall verdict
#   "scope" - what the data can support; a limitation, not a defect, so it is
#             reported separately, the way a clinical report separates its
#             limitations section from its result
at_result <- function(verdict, headline, evidence, consequence, not_assessed,
                      measurements, kind = c("fault", "scope")) {
  kind <- match.arg(kind)
  structure(list(verdict = verdict,
                 kind = kind,
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
  cat(if (identical(x$kind, "scope")) "scope:" else "verdict:", x$verdict, "-", x$headline, "\n")
  if (length(x$evidence)) {
    cat("why:\n"); for (e in x$evidence) cat("  -", e, "\n")
  }
  if (!is.null(x$consequence)) cat("what it costs:", x$consequence, "\n")
  if (length(x$not_assessed)) {
    cat("not assessed:\n"); for (u in x$not_assessed) cat("  -", u, "\n")
  }
  invisible(x)
}
