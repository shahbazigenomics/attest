#' Are the samples who the sample sheet says they are?
#'
#' Infers each sample's sex from expression (XIST against Y-chromosome genes)
#' and compares it with the sample sheet when one is supplied, and looks for
#' samples that are the same library twice.
#'
#' Sex is inferred from expression, never from X heterozygosity: after
#' X-inactivation a clonal female sample expresses one X, so a DNA-style rule
#' would call it male. Human markers only; for other species the check reports
#' that it could not run.
#'
#' @param x counts: matrix, data.frame, SummarizedExperiment, DESeqDataSet or
#'   DGEList.
#' @param metadata optional data.frame of sample annotation, one row per column
#'   of `x`.
#' @param sex_col optional name of the column in `metadata` holding sex. If
#'   `NULL`, a column named sex/gender is used when present.
#' @param dup_floor absolute correlation above which two samples are reported as
#'   the same library regardless of the cohort. Default 0.999. A pair below it is
#'   still flagged when it stands at least 0.02 above every other pair and above
#'   0.99. Measured: airway replicates reach 0.941, pasilla 0.984, fission 0.9955;
#'   a resequenced library reaches 0.997-0.9999 and an exact copy 1.
#' @return object of class "attest_check".
#' @export
attest_identity <- function(x, metadata = NULL, sex_col = NULL, dup_floor = 0.999) {

  m <- at_as_matrix(x)
  if (is.null(m)) {
    return(at_result("UNKNOWN", "The input could not be read as a numeric matrix.",
                     character(0), NULL, character(0), list(input_class = class(x)[1])))
  }
  if (ncol(m) < 2) {
    return(at_result("UNKNOWN", "Sample identity needs at least 2 samples.",
                     character(0), NULL, character(0), list(dim = dim(m))))
  }
  if (is.null(metadata) && !is.null(attr(x, "class")) && isS4(x)) {
    metadata <- tryCatch(as.data.frame(SummarizedExperiment::colData(x)),
                         error = function(e) NULL)
  }

  ev <- list(dim = dim(m))
  na <- character(0)
  lines <- character(0)
  verdict <- "PERMITTED"
  headline <- "Samples are consistent with their labels."

  # --- sex from expression -------------------------------------------------
  sex <- at_sex_from_expression(m)
  ev$sex_called <- sex$call
  ev$xist_cpm   <- sex$xist
  ev$y_cpm      <- sex$y
  if (is.null(sex$call)) {
    na <- c(na, sprintf("sex check: none of the human markers (XIST, RPS4Y1, KDM5D, DDX3Y, UTY, USP9Y, EIF1AY, NLGN4Y) are in the row names, so either this is not human or the genes were filtered out (%d of 8 found)",
                        sex$n_found))
  } else {
    tab <- table(sex$call)
    lines <- c(lines, sprintf("sex inferred from expression: %s",
                              paste(sprintf("%d %s", tab, names(tab)), collapse = ", ")))
    if (any(sex$call == "not determined")) {
      lines <- c(lines, "'not determined' means XIST and the Y genes disagreed - possible in clonal lines and tumours that have lost XIST")
    }
    stated <- at_stated_sex(metadata, sex_col, ncol(m))
    if (is.null(stated)) {
      na <- c(na, "sex-versus-sample-sheet comparison: no sex column supplied (pass metadata and sex_col)")
    } else {
      cmp <- !is.na(stated) & sex$call != "not determined" & stated != sex$call
      ev$sex_mismatches <- sum(cmp)
      if (any(cmp)) {
        verdict <- "NOT PERMITTED"
        headline <- "A sample does not match its label."
        lines <- c(lines, sprintf("%d sample(s) contradict the sample sheet: %s",
                                  sum(cmp),
                                  paste(sprintf("%s labelled %s, expression says %s",
                                                colnames(m)[cmp], stated[cmp], sex$call[cmp]),
                                        collapse = "; ")))
      } else {
        lines <- c(lines, "every sample's inferred sex matches the sample sheet")
      }
    }
  }

  # --- same library twice --------------------------------------------------
  dup <- at_duplicate_pairs(m, dup_floor)
  ev$max_correlation <- dup$max_cor
  ev$cohort_median   <- dup$median_cor
  ev$duplicate_pairs <- dup$pairs
  if (is.null(dup$max_cor)) {
    na <- c(na, "duplicate-sample check: too few expressed genes to correlate samples")
  } else if (nrow(dup$pairs) > 0) {
    verdict <- if (verdict == "NOT PERMITTED") verdict else "CAUTION"
    if (headline == "Samples are consistent with their labels.")
      headline <- "Two samples look like the same library."
    lines <- c(lines, sprintf("%s (cohort median %.4f)",
                              paste(sprintf("%s and %s correlate at %.4f%s", dup$pairs$a, dup$pairs$b,
                                            dup$pairs$cor,
                                            ifelse(dup$pairs$identical, ", and are identical column for column", "")),
                                    collapse = "; "),
                              dup$median_cor))
  } else {
    lines <- c(lines, sprintf("no pair stands out as the same library (highest %.4f, cohort median %.4f)",
                              dup$max_cor, dup$median_cor))
  }

  na <- c(na, "DNA-RNA match: this check compares samples with each other and with the sample sheet, not with genotypes; for that, run somalier or NGSCheckMate on the BAMs")

  at_result(verdict, headline, lines,
            if (verdict == "NOT PERMITTED")
              "A swapped or mislabelled sample assigns the wrong group to real data, so every comparison that uses it is wrong. Resolve it before analysis."
            else if (verdict == "CAUTION")
              "Two copies of one library act as false replication: the comparison looks better powered than it is."
            else NULL,
            na, ev)
}

# human sex markers, by Ensembl gene id and by symbol
at_sex_markers <- function() list(
  xist = c("ENSG00000229807", "XIST"),
  y    = c("ENSG00000129824", "RPS4Y1", "ENSG00000012817", "KDM5D",
           "ENSG00000067048", "DDX3Y",  "ENSG00000183878", "UTY",
           "ENSG00000114374", "USP9Y",  "ENSG00000198692", "EIF1AY",
           "ENSG00000165246", "NLGN4Y")
)

at_sex_from_expression <- function(m, margin = 2) {
  rn <- toupper(sub("\\.\\d+$", "", rownames(m)))       # drop Ensembl version suffix
  mk <- at_sex_markers()
  ix <- which(rn %in% toupper(mk$xist))
  iy <- which(rn %in% toupper(mk$y))
  n_found <- length(ix) + length(iy)
  if (length(ix) == 0 || length(iy) == 0)
    return(list(call = NULL, n_found = n_found, xist = NA_real_, y = NA_real_))

  cpm  <- t(t(m) / colSums(m)) * 1e6
  xist <- colSums(cpm[ix, , drop = FALSE])
  yexp <- apply(cpm[iy, , drop = FALSE], 2, stats::median)
  lr   <- log2((xist + 0.5) / (yexp + 0.5))            # >0 female-like, <0 male-like
  call <- ifelse(lr >  margin, "female",
          ifelse(lr < -margin, "male", "not determined"))
  list(call = call, n_found = n_found, xist = xist, y = yexp, log_ratio = lr)
}

at_stated_sex <- function(metadata, sex_col, n) {
  if (is.null(metadata)) return(NULL)
  metadata <- as.data.frame(metadata)
  if (nrow(metadata) != n) return(NULL)
  if (is.null(sex_col)) {
    hit <- grep("^(sex|gender)$", names(metadata), ignore.case = TRUE)
    if (!length(hit)) return(NULL)
    sex_col <- names(metadata)[hit[1]]
  }
  if (!sex_col %in% names(metadata)) return(NULL)
  v <- tolower(as.character(metadata[[sex_col]]))
  out <- ifelse(v %in% c("f", "female", "woman"), "female",
         ifelse(v %in% c("m", "male", "man"), "male", NA_character_))
  if (all(is.na(out))) NULL else out
}

# A duplicate is either near-perfect in absolute terms, or far above everything
# else in this cohort. Neither test alone works: yeast replicates reach 0.9955
# (so a fixed 0.99 floor cries wolf), while a duplicate there is only ~1.5 MAD
# above the median (so a MAD rule misses it). Measured on airway/pasilla/fission.
at_duplicate_pairs <- function(m, floor_ = 0.999, gap = 0.02, floor_rel = 0.99,
                               max_genes = 5000) {
  keep <- rowSums(m) > 0
  if (sum(keep) < 200 || ncol(m) < 2)
    return(list(max_cor = NULL, median_cor = NA_real_, pairs = data.frame()))
  y <- m[keep, , drop = FALSE]
  if (nrow(y) > max_genes) y <- y[seq(1, nrow(y), length.out = max_genes), , drop = FALSE]
  l <- log1p(t(t(y) / colSums(y)) * 1e6)
  cc <- suppressWarnings(stats::cor(l, method = "spearman"))
  up <- upper.tri(cc)
  vals <- cc[up]
  idx <- which(up, arr.ind = TRUE)
  nm <- if (is.null(colnames(m))) paste0("sample", seq_len(ncol(m))) else colnames(m)

  identical_cols <- apply(idx, 1, function(k) identical(m[, k[1]], m[, k[2]]))
  ord <- order(vals, decreasing = TRUE)
  flag <- logical(length(vals))
  for (i in ord) {
    rest <- max(vals[-i])                      # the cohort without this pair
    flag[i] <- identical_cols[i] ||
      vals[i] >= floor_ ||
      (vals[i] >= floor_rel && vals[i] - rest >= gap)
  }
  list(max_cor    = max(vals),
       median_cor = stats::median(vals),
       pairs      = data.frame(a = nm[idx[flag, 1]], b = nm[idx[flag, 2]],
                               cor = vals[flag], identical = identical_cols[flag],
                               stringsAsFactors = FALSE))
}
