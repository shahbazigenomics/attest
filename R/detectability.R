#' For which genes could a change have been seen at all?
#'
#' Adjudicates the claim "this gene did not change". For every gene it reports
#' the smallest fold change the dataset could have detected, given that gene's
#' own counts and dispersion, the group sizes and the multiple-testing burden.
#' A gene whose detectable change is larger than any effect worth reporting
#' cannot support a negative claim, however small its p-value.
#'
#' The question it answers is what the analysis you will run - DESeq2 - could
#' have detected. Per gene, the smallest fold change that change in *either*
#' direction would reach the given power at the level the gene is judged at:
#'   log FC = z * sqrt( (1/mu + phi)/n + (FC/mu + phi)/n ),  z = z_alpha/2 + z_power
#' The decrease is the harder direction (the lower group has fewer reads), so FC
#' appears on both sides and is solved by iteration. A gene averaging too few
#' reads has no finite answer: its count can only fall to zero.
#'
#' * **Noise (phi).** DESeq2's own dispersion estimates, on the design when one
#'   is given, whenever DESeq2 is installed. On small designs DESeq2 estimates
#'   dispersion 27-36% above the truth, and a formula fed the true value
#'   promises power DESeq2 does not deliver. Without DESeq2, attest's moment
#'   estimate is used and the report says what that costs. Either way, a given
#'   design's residuals define the noise: pooling within the compared groups
#'   counts the donor or cell-line differences a paired design removes.
#' * **Level.** Benjamini-Hochberg rejects p <= q*R/n, so the level depends on
#'   how many genes truly change; R is estimated from these counts with a fast
#'   Wald test, Bonferroni as fallback. Bonferroni throughout would be badly
#'   pessimistic: on airway it put the detectable change above 2-fold for every
#'   gene, where DESeq2 finds about 4,000.
#'
#' **Calibrated against DESeq2 itself.** `validation/detectability_calibration_deseq2.R`
#' simulates experiments from DESeq2's fit to airway (paired and unpaired) and to
#' Kang 2018 pseudobulk (3, 4, 6 and 8 donors, paired), gives each gene exactly
#' its claimed detectable change, and runs DESeq2: a claimed 80% was detected
#' 80-85% of the time (single experiments 78-86%). With the moment estimate
#' instead, 65-78%. The first version of this check, calibrated only against a
#' per-gene `glm.nb` test, gave 0.66 against DESeq2 on an unpaired design.
#'
#' @param x counts: matrix, data.frame, SummarizedExperiment, DESeqDataSet or
#'   DGEList.
#' @param group factor or vector of length `ncol(x)` splitting the samples into
#'   the two groups being compared, or the name of a column in `metadata`.
#' @param metadata optional data.frame holding `group`.
#' @param power power to report the threshold at. Default 0.8.
#' @param alpha family-wise level before the multiple-testing correction.
#'   Default 0.05.
#' @param target_fc the effect size the study cares about, used for the
#'   study-level summary. Default 2 (a doubling).
#' @param design,of_interest optional design formula and the variable being
#'   tested. When given with `metadata`, the calculation uses the design's
#'   effective sample size - group sizes divided by the variance inflation that
#'   the other terms put on this effect - so a partly confounded experiment is
#'   not credited with precision it does not have.
#' @param dispersion where the noise estimate comes from: `"auto"` (DESeq2's
#'   own estimates when DESeq2 is installed and the values are integers, else
#'   the moment estimate), `"deseq2"`, or `"moments"`.
#' @details
#' This check reports *scope*, not a fault: it says what the data can support, so
#' it does not lower the overall verdict of [attest()], the way a clinical report
#' keeps its limitations section separate from its result.
#'
#' @return object of class "attest_check"; `$measurements$min_detectable_fc` is
#'   the per-gene vector.
#' @examples
#' fx <- readRDS(system.file("extdata", "fixtures.rds", package = "attest"))
#' dex <- factor(rep(c("untrt", "trt"), times = 4))
#' d <- attest_detectability(fx$airway, group = dex, dispersion = "moments")
#' d
#' \donttest{
#' # the default: DESeq2's own noise estimate, calibrated against DESeq2
#' cell <- factor(rep(c("A", "B", "C", "D"), each = 2))
#' attest_detectability(fx$airway, group = dex, metadata = data.frame(cell, dex),
#'                      design = ~ cell + dex)
#' }
#' summary(d$measurements$min_detectable_fc[d$measurements$assessable])
#' @export
attest_detectability <- function(x, group = NULL, metadata = NULL, power = 0.8,
                                 alpha = 0.05, target_fc = 2, design = NULL,
                                 of_interest = NULL,
                                 dispersion = c("auto", "deseq2", "moments")) {

  m <- at_as_matrix(x)
  if (is.null(m)) {
    return(at_result("UNKNOWN", "The input could not be read as a numeric matrix.",
                     character(0), NULL, character(0), list(input_class = class(x)[1]), kind = "scope"))
  }
  g <- at_resolve_group(group, metadata, m)
  if (is.null(g)) {
    return(at_result("UNKNOWN",
                     "Detectability needs the two groups being compared: pass group = <factor> or group = \"<column>\" with metadata.",
                     character(0), NULL, character(0), list(dim = dim(m)), kind = "scope"))
  }
  lev <- levels(g)
  if (length(lev) != 2) {
    return(at_result("UNKNOWN",
                     sprintf("Detectability compares two groups; '%s' has %d levels.",
                             deparse(substitute(group)), length(lev)),
                     character(0), NULL, character(0), list(levels = lev), kind = "scope"))
  }
  n1 <- sum(g == lev[1]); n2 <- sum(g == lev[2])
  if (min(n1, n2) < 2) {
    return(at_result("UNKNOWN",
                     sprintf("Detectability needs at least 2 samples per group (%s = %d, %s = %d).",
                             lev[1], n1, lev[2], n2),
                     character(0), NULL, character(0), list(n1 = n1, n2 = n2), kind = "scope"))
  }

  sf   <- at_size_factors(m)
  if (!all(is.finite(sf))) sf <- colSums(m) / mean(colSums(m))
  norm <- t(t(m) / sf)
  mu   <- rowMeans(norm)
  expressed <- mu > 0

  # Noise. The question is what the analysis the user will run could detect,
  # so where DESeq2 is installed its own dispersion estimates are used, on the
  # design when one is given: on these small designs DESeq2 estimates
  # dispersion 27-36% above the truth, and a formula fed the true value
  # promises power DESeq2 does not deliver (0.66-0.70 for a claimed 0.80).
  # Without DESeq2 the moment estimate is used - accurate for the true
  # dispersion, and so optimistic about DESeq2 by a stated, measured margin.
  # Either way a given design's residuals, not the two groups, define the
  # noise: pooling within groups counts donor or cell-line differences the
  # design removes (dispersion 0.039 vs 0.020 on airway).
  dispersion <- match.arg(dispersion)
  mm <- at_design_mm(design, metadata)
  usable_mm <- !is.null(mm) && nrow(mm) == ncol(m) && qr(mm)$rank == ncol(mm) && nrow(mm) - ncol(mm) >= 2
  phi <- NULL
  phi_source <- NULL
  if (dispersion %in% c("auto", "deseq2") && all(m == round(m)) &&
      requireNamespace("DESeq2", quietly = TRUE)) {
    phi <- at_dispersion_deseq2(m, if (usable_mm) design else NULL,
                                if (usable_mm) metadata else NULL, g)
    if (!is.null(phi))
      phi_source <- sprintf("DESeq2's own dispersion estimates for %s",
                            if (usable_mm) paste(deparse(design), collapse = "") else "~ group")
  }
  if (is.null(phi)) {
    phi <- if (usable_mm) at_dispersion_design(norm, mm) else at_dispersion_moments(norm, g)
    phi_source <- sprintf("attest's moment estimate %s",
                          if (usable_mm) sprintf("after fitting %s", paste(deparse(design), collapse = ""))
                          else "within the two groups")
  }
  deseq2_used <- startsWith(phi_source, "DESeq2")
  n_eff <- 2 / (1/n1 + 1/n2)                      # harmonic mean of the group sizes

  # the design's confounding costs precision; charge it here, where precision
  # is judged, rather than in the design check's verdict
  vif <- at_design_vif(design, metadata, of_interest)
  if (!is.null(vif) && !is.finite(vif)) {
    return(at_result("UNKNOWN",
                     "The design cannot estimate this effect at all, so there is no precision to report (see the design check).",
                     character(0), NULL, character(0), list(n1 = n1, n2 = n2, vif = vif), kind = "scope"))
  }
  vif_line <- NULL
  if (!is.null(vif) && vif > 1.05) {
    vif_line <- sprintf("the design shares this effect with its other terms (variance inflation %.2f), so the %d + %d samples are treated as %.1f + %.1f",
                        vif, n1, n2, n1 / vif, n2 / vif)
    n_eff <- n_eff / vif
  }
  n_tested <- sum(expressed)

  # The level each gene is actually judged at. Everyone analyses with
  # Benjamini-Hochberg, not Bonferroni, and BH rejects p <= q*R/n, so the
  # effective per-gene level depends on how many genes really change. We
  # estimate R with a fast Wald test on these counts; with no discoveries we
  # fall back to Bonferroni, which is the conservative answer.
  a_bonf <- alpha / max(n_tested, 1)
  pvals  <- at_quick_wald(norm, g, mu, phi, expressed)
  n_rej  <- if (all(is.na(pvals))) 0 else sum(stats::p.adjust(pvals, "BH") < alpha, na.rm = TRUE)
  a <- if (n_rej > 0) alpha * n_rej / n_tested else a_bonf
  z <- stats::qnorm(1 - a/2) + stats::qnorm(power)
  mdfc <- at_mdfc(mu, phi, n_eff, z)              # minimum detectable fold change
  mdfc[!expressed] <- Inf

  # judge the study on genes with enough counts to be worth a claim; almost every
  # dataset is mostly near-zero genes, and letting those decide the verdict would
  # fail every real experiment
  assessable <- expressed & mu >= 10
  if (!any(assessable)) assessable <- expressed
  frac_ok <- mean(mdfc[assessable] <= target_fc)
  ev <- list(dim = dim(m), n1 = n1, n2 = n2, levels = lev, power = power,
             dispersion_source = phi_source, deseq2_dispersion = deseq2_used,
             design_vif = if (is.null(vif)) NA_real_ else vif, n_effective = n_eff,
             alpha = alpha, per_gene_alpha = a, bonferroni_alpha = a_bonf,
             n_rejected_estimate = n_rej, n_tested = n_tested,
             mean_count = mu,
             target_fc = target_fc, frac_detectable_at_target = frac_ok,
             median_min_detectable_fc = stats::median(mdfc[assessable]),
             assessable = assessable,
             min_detectable_fc = mdfc)

  lines <- c(
    sprintf("of the %s genes with a mean count of 10 or more, %.0f%% could have shown a %.1f-fold change at %.0f%% power",
            format(sum(assessable), big.mark = ","), 100 * frac_ok, target_fc, 100 * power),
    sprintf("among those, the median gene needed %.1f-fold and the quietest quarter %.1f-fold or more",
            stats::median(mdfc[assessable]),
            stats::quantile(mdfc[assessable], 0.75, names = FALSE)),
    {
      low <- mdfc[expressed & !assessable]
      if (is.finite(stats::median(low)))
        sprintf("the %s genes below that count are effectively untestable here (median %.0f-fold needed)",
                format(length(low), big.mark = ","), stats::median(low))
      else
        sprintf("the %s genes below that count are effectively untestable here: for most of them no decrease of any size reaches %.0f%% power, because a count this low can only fall to zero",
                format(length(low), big.mark = ","), 100 * power)
    },
    if (n_rej > 0)
      sprintf("group sizes %d and %d; about %s genes really differ, so Benjamini-Hochberg judges each gene at p < %.2g (Bonferroni would be %.2g)",
              n1, n2, format(n_rej, big.mark = ","), a, a_bonf)
    else
      sprintf("group sizes %d and %d; no gene survives multiple testing on these counts, so each gene is judged at the Bonferroni level p < %.2g",
              n1, n2, a),
    vif_line,
    sprintf("noise: %s", phi_source))

  verdict <- if (frac_ok >= 0.8) "PERMITTED" else if (frac_ok >= 0.4) "CAUTION" else "NOT PERMITTED"
  headline <- switch(verdict,
    "PERMITTED"     = sprintf("A %.1f-fold change was detectable for most genes.", target_fc),
    "CAUTION"       = sprintf("A %.1f-fold change was detectable for under %.0f%% of genes.", target_fc, 80),
    "NOT PERMITTED" = sprintf("Most genes could not have shown a %.1f-fold change.", target_fc))

  at_result(verdict, headline, lines,
    paste("A gene absent from your results list is not evidence that it does not respond;",
          "for genes above their threshold above, it is evidence, and for the rest it is not.",
          "attest_detectability(x, group)$measurements$min_detectable_fc gives the per-gene number,",
          "which is what a negative claim about a named gene has to quote."),
    if (deseq2_used)
      "calibrated against DESeq2 itself on paired and unpaired designs built from airway and Kang 2018 (validation/detectability_calibration_deseq2.R): a claimed 80% gave 80-85% (single simulated experiments 78-86%) across 3-8 replicates. edgeR and limma-voom have not been calibrated"
    else
      paste(if (dispersion == "moments") "The moment estimate was chosen" else "DESeq2 is not installed",
            "so the noise is attest's own estimate of the true dispersion. DESeq2 estimates it higher on small designs and detected 65-78% of changes at these sizes, not 80% (validation/detectability_calibration_deseq2.R); with DESeq2 installed the figures hold"),
    ev, kind = "scope")
}

at_resolve_group <- function(group, metadata, m) {
  if (is.null(group)) return(NULL)
  if (length(group) == 1 && is.character(group)) {
    if (is.null(metadata)) return(NULL)
    metadata <- as.data.frame(metadata)
    if (!group %in% names(metadata) || nrow(metadata) != ncol(m)) return(NULL)
    group <- metadata[[group]]
  }
  if (length(group) != ncol(m)) return(NULL)
  droplevels(as.factor(group))
}

# fast per-gene Wald test, used only to estimate how many genes really differ
at_quick_wald <- function(norm, g, mu, phi, expressed) {
  lev <- levels(g)
  a <- norm[, g == lev[1], drop = FALSE]; b <- norm[, g == lev[2], drop = FALSE]
  ma <- rowMeans(a) + 0.5; mb <- rowMeans(b) + 0.5
  lfc <- log(mb / ma)
  se  <- sqrt((1/ncol(a)) * (1/pmax(ma, 1e-8) + phi) + (1/ncol(b)) * (1/pmax(mb, 1e-8) + phi))
  p <- 2 * stats::pnorm(-abs(lfc / se))
  p[!expressed] <- NA_real_
  p
}

# method-of-moments dispersion per gene, pooled within groups
at_dispersion_moments <- function(norm, g) {
  lev <- levels(g)
  within_var <- 0; within_mu <- 0; df <- 0
  for (l in lev) {
    y <- norm[, g == l, drop = FALSE]
    n <- ncol(y)
    if (n < 2) next
    mu_l <- rowMeans(y)
    v_l  <- rowSums((y - mu_l)^2) / (n - 1)
    within_var <- within_var + v_l * (n - 1)
    within_mu  <- within_mu + mu_l * (n - 1)
    df <- df + (n - 1)
  }
  v  <- within_var / df
  mu <- within_mu / df
  phi <- (v - mu) / pmax(mu^2, .Machine$double.eps)     # var = mu + phi*mu^2
  phi[!is.finite(phi)] <- NA_real_
  phi <- pmax(phi, 0)
  # shrink the noisy per-gene estimates toward the trend, as the count models do
  fill <- stats::median(phi[is.finite(phi) & phi > 0], na.rm = TRUE)
  if (!is.finite(fill)) fill <- 0.1
  phi[is.na(phi)] <- fill
  0.5 * phi + 0.5 * fill
}

# the design's model matrix, or NULL when it cannot be built on these samples
at_design_mm <- function(design, metadata) {
  if (is.null(design) || is.null(metadata) || !inherits(design, "formula")) return(NULL)
  md <- as.data.frame(metadata)
  if (!all(all.vars(design) %in% names(md))) return(NULL)
  mm <- tryCatch(stats::model.matrix(design, md), error = function(e) NULL)
  if (is.null(mm) || nrow(mm) != nrow(md)) return(NULL)
  mm
}

# method-of-moments dispersion from the residuals of the full design, shrunk
# toward the median the same way as the within-group estimate
at_dispersion_design <- function(norm, mm) {
  fit <- stats::lm.fit(mm, t(norm))
  v   <- colSums(fit$residuals^2) / (nrow(mm) - ncol(mm))
  mu  <- rowMeans(norm)
  phi <- (v - mu) / pmax(mu^2, .Machine$double.eps)
  phi[!is.finite(phi)] <- NA_real_
  phi <- pmax(phi, 0)
  fill <- stats::median(phi[is.finite(phi) & phi > 0], na.rm = TRUE)
  if (!is.finite(fill)) fill <- 0.1
  phi[is.na(phi)] <- fill
  0.5 * phi + 0.5 * fill
}

# DESeq2's own dispersion estimates, on the design if given, else on ~ group.
# NULL on any failure, so the caller falls back to the moment estimate.
at_dispersion_deseq2 <- function(m, design, metadata, g) {
  if (is.null(design)) { design <- ~ group; metadata <- data.frame(group = g) }
  metadata <- as.data.frame(metadata)
  rownames(metadata) <- colnames(m)
  storage.mode(m) <- "integer"
  keep <- rowSums(m) > 0
  phi <- tryCatch({
    dds <- suppressWarnings(suppressMessages(
      DESeq2::DESeqDataSetFromMatrix(m[keep, , drop = FALSE], metadata, design)))
    dds <- DESeq2::estimateSizeFactors(dds)
    dds <- suppressWarnings(suppressMessages(DESeq2::estimateDispersions(dds, quiet = TRUE)))
    DESeq2::dispersions(dds)
  }, error = function(e) NULL)
  if (is.null(phi)) return(NULL)
  out <- rep(NA_real_, nrow(m)); out[keep] <- phi
  fill <- stats::median(out[is.finite(out)], na.rm = TRUE)
  out[!is.finite(out)] <- if (is.finite(fill)) fill else 0.1
  out
}

# Smallest fold change detectable at the given z = z_(alpha/2) + z_power, in
# either direction. The decrease is the harder one - the lower group's mean is
# mu / FC, and its counting noise is larger - so FC solves
#   log FC = z * sqrt( (1/mu + phi) / n + (FC/mu + phi) / n ).
# Using the average count for both groups instead understated the noise, and
# realised power against DESeq2 fell to 0.69-0.77.
at_mdfc <- function(mu, phi, n_eff, z) {
  mu <- pmax(mu, .Machine$double.eps)
  fc <- exp(z * sqrt((2 / n_eff) * (1 / mu + phi)))
  for (it in 1:50) {
    nxt <- exp(z * sqrt((1 / n_eff) * (1 / mu + phi) + (1 / n_eff) * (fc / mu + phi)))
    nxt[!is.finite(nxt)] <- Inf
    done <- all(abs(nxt - fc) < 1e-6 * fc | !is.finite(fc))
    fc <- nxt
    if (done) break
  }
  fc
}
