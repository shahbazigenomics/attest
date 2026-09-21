# attest 0.1.0

First release.

* `attest()` runs every check that its inputs allow and returns one report. It
  takes a matrix, `data.frame`, `SummarizedExperiment`, `DESeqDataSet`,
  `DGEList` or `tximport` list, and reads the sample sheet, design and grouping
  from the object when it carries them. Checks with no inputs are listed as
  "not run".
* `attest_file()` reads a count file - featureCounts output, a GEO supplementary
  table, CSV/TSV, gzip - and reports how it was read.
* Checks: value scale (`attest_counts()`), identifiers (`attest_identifiers()`),
  completeness (`attest_completeness()`), sample identity (`attest_identity()`),
  design (`attest_design()`) and detectability (`attest_detectability()`).
