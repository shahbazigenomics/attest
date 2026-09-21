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
* `attest_file()` drops empty rows and reports rows of numbers with no gene
  identifier - the totals row at the foot of a spreadsheet - as NOT PERMITTED,
  since every tool downstream would analyse it as a gene. A `Position` column
  is treated as annotation. `attest_identity()` declines a matrix with missing
  values instead of failing; it crashed on the first such file from GEO.
* `attest_strandedness()`: the library protocol from STAR's `ReadsPerGene.out.tab`,
  and which of its columns a count matrix was built from; a featureCounts
  `.summary` raises a collapse of assigned reads and otherwise says it cannot
  confirm the strand setting. `attest(x, strandedness = files)` adds it to the
  report, and `attest_file()` picks up the summary beside featureCounts output.
* Detectability answers what DESeq2 could have detected on this design, and is
  calibrated against DESeq2 itself: a claimed 80% power was realised 80-85% of
  the time on paired and unpaired designs from airway and Kang 2018. It uses
  DESeq2's own dispersion estimates when DESeq2 is installed, measures noise
  after fitting the design, and takes the harder (decreasing) direction.
