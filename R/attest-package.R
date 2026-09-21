#' attest: check whether RNA-seq count data can support the analysis run on it
#'
#' `attest` audits a count matrix, its sample sheet and its design before
#' differential expression. It does not analyse expression. Each check
#' adjudicates one claim and returns a verdict - `PERMITTED`, `CAUTION`,
#' `NOT PERMITTED` or `UNKNOWN` - with the evidence behind it and the checks the
#' inputs could not answer. `UNKNOWN` always means "the inputs cannot answer
#' this", never a guess, and nothing throws: malformed input returns a typed
#' result with a reason.
#'
#' Start with [attest()] on a matrix, `SummarizedExperiment`, `DESeqDataSet` or
#' `DGEList`, or [attest_file()] on a count file. The individual checks are
#' [attest_counts()], [attest_identifiers()], [attest_completeness()],
#' [attest_identity()], [attest_design()], [attest_detectability()] and
#' [attest_strandedness()].
#'
#' The package ships `extdata/fixtures.rds`: 3,000-gene subsets of three public
#' datasets - airway (human, 8 samples), pasilla (fly, 7) and fission (yeast,
#' 36) - as raw count matrices, used by the examples, the tests and the vignette.
#' `extdata/strand/` holds STAR `ReadsPerGene.out.tab` files and featureCounts
#' output for libraries of known protocol, simulated from a synthetic genome and
#' counted with the real tools; they back [attest_strandedness()].
#'
#' @keywords internal
"_PACKAGE"
