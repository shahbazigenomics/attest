# What the object already carries
#
# A DESeqDataSet holds the counts, the sample sheet and the design formula; a
# DGEList holds the counts and a sample table with the grouping. Asking the user
# to type them again is friction, and - worse - when they do not, the checks
# that need them are skipped. This returns whatever the object can supply, and
# never raises: an object from a package that is not installed just yields NULL.
at_from_object <- function(x) {
  out <- list(metadata = NULL, design = NULL, group = NULL, source = NULL)
  if (is.null(x)) return(out)

  # --- DGEList (and anything else with a $samples table) -------------------
  if (!isS4(x) && is.list(x) && is.data.frame(x$samples)) {
    if (nrow(x$samples)) out$metadata <- x$samples
    g <- x$samples$group
    if (!is.null(g) && nlevels(droplevels(as.factor(g))) >= 2) out$group <- g
    out$source <- "DGEList"
    return(out)
  }

  if (!isS4(x)) return(out)

  # --- SummarizedExperiment / DESeqDataSet ---------------------------------
  cd <- tryCatch(as.data.frame(SummarizedExperiment::colData(x)),
                 error = function(e) NULL)
  if (!is.null(cd) && nrow(cd)) out$metadata <- cd

  d <- tryCatch(x@design, error = function(e) NULL)      # DESeqDataSet slot
  if (inherits(d, "formula") && length(all.vars(d))) {
    out$design <- d
    out$source <- "DESeqDataSet"
  } else if (!is.null(out$metadata)) {
    out$source <- "SummarizedExperiment"
  }
  out
}

# the variable whose effect is being tested, and the grouping it implies
at_of_interest <- function(design, metadata, given = NULL) {
  if (!is.null(given)) return(given)
  if (!inherits(design, "formula")) return(NULL)
  tv <- all.vars(design)
  if (!length(tv)) return(NULL)
  tv[length(tv)]
}

at_two_level_group <- function(of_interest, metadata) {
  if (is.null(of_interest) || is.null(metadata)) return(NULL)
  md <- as.data.frame(metadata)
  if (!of_interest %in% names(md)) return(NULL)
  v <- droplevels(as.factor(md[[of_interest]]))
  if (nlevels(v) == 2) v else NULL
}
