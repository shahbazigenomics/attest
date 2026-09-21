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
* "Already normalised" is decided from the rounding artefact that dividing by a
  size factor leaves in each sample's small-value histogram, not from flat size
  factors alone, which raw libraries of equal depth also have (whole blood with
  variable globin was a false alarm under the old rule). The new measure also
  catches normalised matrices down to 30k reads per sample, where the old one
  failed below 1M.
