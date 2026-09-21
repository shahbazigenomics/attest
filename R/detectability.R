#' For which genes could a change have been seen at all?
#'
#' Adjudicates the claim "this gene did not change". For every gene it reports
#' the smallest fold change the dataset could have detected, given that gene's
#' own counts and dispersion, the group sizes and the multiple-testing burden.
#' A gene whose detectable change is larger than any effect worth reporting
#' cannot support a negative claim, however small its p-value.
#'
#' The estimate is a two-sample negative-binomial Wald calculation:
#'   log(FC) = (z_alpha/2 + z_power) * sqrt(2/n * (1/mu + phi))
#' The level `alpha` is the one a gene is really judged at. Analyses use
#' Benjamini-Hochberg, which rejects p <= q*R/n, so the level depends on how many
#' genes truly change; R is estimated from these counts with a fast Wald test.
#' With no discoveries it falls back to Bonferroni, the conservative answer.
#' Using Bonferroni throughout would be badly pessimistic: on airway it put the
#' detectable change above 2-fold for every gene, in a dataset where DESeq2 finds
#' 3,993 differentially expressed ones.
#'
#' It is an approximation. `validation/detectability_calibration.R` simulates data
#' with exactly the claimed effect and measures how often it is actually detected;
#' realised power came to 0.78-0.84 against a claim of 0.80 at n = 3, 4, 6 and 10.
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
#' This check reports *scope*, not a fault: it says what the data can support, so
#' it does not lower the overall verdict of [attest()], the way a clinical report
#' keeps its limitations section separate from its result.
#'
#' @return object of class "attest_check"; `$measurements$min_detectable_fc` is
#'   the per-gene vector.
#' @export
attest_detectability <- function(x, group = NULL, metadata = NULL, power = 0.8,
                                 alpha = 0.05, target_fc = 2, design = NULL,
                                 of_interest = NULL) {

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

  phi <- at_dispersion_moments(norm, g)
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
  se_log <- sqrt((2 / n_eff) * (1/pmax(mu, .Machine$double.eps) + phi))
  mdfc <- exp(z * se_log)                         # minimum detectable fold change
  mdfc[!expressed] <- Inf

  # judge the study on genes with enough counts to be worth a claim; almost every
  # dataset is mostly near-zero genes, and letting those decide the verdict would
  # fail every real experiment
  assessable <- expressed & mu >= 10
  if (!any(assessable)) assessable <- expressed
  frac_ok <- mean(mdfc[assessable] <= target_fc)
  ev <- list(dim = dim(m), n1 = n1, n2 = n2, levels = lev, power = power,
             design_vif = if (is.null(vif)) NA_real_ else vif, n_effective = n_eff,
             alpha = alpha, per_gene_alpha = a, bonferroni_alpha = a_bonf,
             n_rejected_estimate = n_rej, n_tested = n_tested,
             mean_count = mu,
             target_fc = target_fc, frac_detectable_at_target = frac_ok,
             median_min_detectable_fc = stats::median(mdfc[expressed]),
             assessable = assessable,
             min_detectable_fc = mdfc)

  lines <- c(
    sprintf("of the %s genes with a mean count of 10 or more, %.0f%% could have shown a %.1f-fold change at %.0f%% power",
            format(sum(assessable), big.mark = ","), 100 * frac_ok, target_fc, 100 * power),
    sprintf("among those, the median gene needed %.1f-fold and the quietest quarter %.1f-fold or more",
            stats::median(mdfc[assessable]),
            stats::quantile(mdfc[assessable], 0.75, names = FALSE)),
    sprintf("the %s genes below that count are effectively untestable here (median %.0f-fold needed)",
            format(sum(expressed & !assessable), big.mark = ","),
            stats::median(mdfc[expressed & !assessable])),
    if (n_rej > 0)
      sprintf("group sizes %d and %d; about %s genes really differ, so Benjamini-Hochberg judges each gene at p < %.2g (Bonferroni would be %.2g)",
              n1, n2, format(n_rej, big.mark = ","), a, a_bonf)
    else
      sprintf("group sizes %d and %d; no gene survives multiple testing on these counts, so each gene is judged at the Bonferroni level p < %.2g",
              n1, n2, a),
    vif_line)

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
    "the estimate is a negative-binomial Wald approximation with method-of-moments dispersion; validation/detectability_calibration.R measures how close the claimed 80% power comes to DESeq2 in simulation",
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
