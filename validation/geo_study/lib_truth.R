# What an author file really is, decided from NCBI's own raw counts for the
# same samples - never from attest's rules. Each function returns a label and
# the measurement behind it, and "undetermined" when the comparison cannot
# decide.

# --- value scale ---------------------------------------------------------------
# Two independent measurements against NCBI's raw counts for the same samples.
#
# Level (primary): for genes NCBI counts >= 50 in every matched sample, the
# median of author value / NCBI count, per sample. The same reads counted by
# two pipelines agree to within a factor of ~2 (read vs fragment counting,
# gene models), so raw and estimated counts sit near 1; CPM, TPM and FPKM sit
# near 1e6 / depth - about 0.01-0.05 at usual depths. On GEO files: raw counts
# 1.00 (GSE161013), FPKM 0.007-0.013 (GSE190775). Values within 3-fold of
# NCBI's are on the count scale; all below 0.2 means depth has been removed.
#
# Depth (secondary, and the only way to tell normalised counts - on the count
# scale, but with depth divided out - from estimated counts): the slope of
# log(author/NCBI level) on log(NCBI depth) across samples: ~0 when depth is
# still in the values, ~-1 when divided out. Used only with >= 4 samples and a
# standard error below 0.15; otherwise it does not decide.
#
# Totals are not used: an author file's total includes features NCBI does not
# count, so on GSE161013 (raw counts) the old slope on totals read 0.09.
truth_value_scale <- function(A, Ag, R) {
  # A: author values, matched columns (as read); Ag: the same by NCBI GeneID;
  # R: NCBI raw counts, same samples in the same order
  n <- ncol(A)
  vals <- A[is.finite(A)]
  integer <- all(abs(vals - round(vals)) < 1e-8)
  base <- list(n = n, integer = integer, level = NA_real_, level_min = NA_real_, level_max = NA_real_,
               depth_slope = NA_real_, depth_slope_se = NA_real_, depth_spread = NA_real_, why = NA_character_)
  done <- function(label, ...) { o <- utils::modifyList(base, list(...)); o$label <- label; o }
  if (any(vals < 0) || (!integer && max(vals) < 50)) return(done("log-transformed"))

  g <- intersect(rownames(Ag), rownames(R))
  g <- g[rowSums(R[g, , drop = FALSE] >= 50) == ncol(R)]
  if (length(g) < 200) return(done("undetermined", why = "fewer than 200 genes counted >= 50 by NCBI in every matched sample"))
  r <- apply(Ag[g, , drop = FALSE] / R[g, , drop = FALSE], 2, stats::median)
  depth <- colSums(R)
  sp <- max(depth) / max(min(depth), 1)
  ds <- dse <- NA_real_
  if (n >= 4 && all(r > 0)) {
    f <- summary(stats::lm(log(r) ~ log(depth)))$coefficients
    if (nrow(f) == 2) { ds <- f[2, 1]; dse <- f[2, 2] }
  }
  slope_ok <- is.finite(dse) && dse < 0.15
  base <- utils::modifyList(base, list(level = stats::median(r), level_min = min(r), level_max = max(r),
                                       depth_slope = ds, depth_slope_se = dse, depth_spread = sp))
  if (all(r < 0.2)) return(done("depth removed"))
  if (all(r > 1/3 & r < 3)) {
    # a decisive slope wins over whole numbers: normalised counts are often
    # rounded before upload ("normalized_readcounts")
    if (slope_ok && ds < -0.7) return(done("depth removed"))
    if (integer) return(done("raw counts"))
    if (slope_ok && ds > -0.3) return(done("estimated counts"))
    return(done("count scale, not whole numbers",
                why = "estimated or size-factor normalised counts; too few samples or too-even depths to tell which"))
  }
  done("ambiguous", why = "levels neither all within 3-fold of NCBI's nor all below 0.2")
}

# --- completeness ------------------------------------------------------------------
# A complete matrix keeps the genes nobody detected. Among NCBI genes with zero
# reads in every matched sample - and that the author's ID system can name -
# what share appear in the author's file at all? Pipelines differ, so some of
# NCBI's zeros are counted by the authors; the thresholds are set wide
# (filtered < 0.3, complete > 0.6) and the continuous share is kept.
truth_completeness <- function(author_gene_ids, R, can_name) {
  zero <- rownames(R)[rowSums(R) == 0]
  zero <- intersect(zero, can_name)
  if (length(zero) < 100)
    return(list(label = "undetermined", present_frac = NA_real_, n_zero = length(zero)))
  f <- mean(zero %in% author_gene_ids)
  list(label = if (f < 0.3) "filtered" else if (f > 0.6) "complete" else "ambiguous",
       present_frac = f, n_zero = length(zero))
}

# --- sex -------------------------------------------------------------------------
# GEO's stated sex against XIST / Y-gene expression in NCBI's counts (not the
# author's file), using attest's calling rule.
truth_sex <- function(meta, R, annot) {
  key <- grep("^ch_(sex|gender)$", names(meta), value = TRUE)
  if (!length(key)) return(NULL)
  stated <- tolower(trimws(meta[[key[1]]])); names(stated) <- meta$gsm
  stated <- setNames(ifelse(stated %in% c("f", "female", "woman", "women"), "female",   # ifelse() drops names
                     ifelse(stated %in% c("m", "male", "man", "men"), "male", NA)), names(stated))
  gs <- intersect(names(stated)[!is.na(stated)], colnames(R))
  if (length(gs) < 2) return(NULL)
  sym <- annot$Symbol[match(rownames(R), annot$GeneID)]
  X <- R[!is.na(sym), gs, drop = FALSE]; rownames(X) <- make.unique(sym[!is.na(sym)])
  call <- at_sex_from_expression(X)
  if (is.null(call$call)) return(NULL)
  data.frame(gsm = gs, stated = stated[gs], expressed = call$call,
             log_ratio = round(call$log_ratio, 2),
             mismatch = call$call != "not determined" & call$call != stated[gs],
             stringsAsFactors = FALSE)
}
