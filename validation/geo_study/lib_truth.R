# What an author file really is, decided from NCBI's own raw counts for the
# same samples - never from attest's rules. Each function returns a label and
# the measurement behind it, and "undetermined" when the comparison cannot
# decide.

# --- value scale ---------------------------------------------------------------
# In raw counts each sample's total rises with its sequencing depth, and NCBI's
# totals measure that depth independently. Regressing log(author total) on
# log(NCBI total) across samples gives a slope near 1 when depth is still in
# the values and near 0 when it has been divided out (CPM, TPM, FPKM, size-
# factor normalised, or rarefied). Needs >= 4 matched samples whose NCBI
# depths differ by >= 1.3x; otherwise undetermined.
truth_value_scale <- function(A, R) {
  n <- ncol(A)
  vals <- A[is.finite(A)]
  integer <- all(abs(vals - round(vals)) < 1e-8)
  if (any(vals < 0) || (!integer && max(vals) < 50))
    return(list(label = "log-transformed", slope = NA_real_, depth_spread = NA_real_,
                integer = integer, n = n))
  r_tot <- colSums(R); a_tot <- colSums(A)
  spread <- max(r_tot) / max(min(r_tot), 1)
  if (n < 4 || spread < 1.3 || any(a_tot <= 0))
    return(list(label = "undetermined", slope = NA_real_, depth_spread = spread,
                integer = integer, n = n,
                why = if (n < 4) "fewer than 4 matched samples" else "sequencing depths too even to tell"))
  b <- unname(stats::coef(stats::lm(log(a_tot) ~ log(r_tot)))[2])
  label <- if (b > 0.7) { if (integer) "raw counts" else "estimated counts" }
           else if (b < 0.3) "depth removed"
           else "ambiguous"
  list(label = label, slope = b, depth_spread = spread, integer = integer, n = n)
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
