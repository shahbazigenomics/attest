# attest on published GEO data

Does attest find real problems in real published RNA-seq data - problems the
standard DESeq2 / edgeR / limma workflows run straight past - and how often?

Every failure attest was validated on before this was constructed. This study
takes a **random sample of human GEO series**, downloads the count matrices
the authors supplied, runs attest on them exactly as a user would, and checks
each verdict against an **independent truth: NCBI's own raw counts** for the
same samples, which NCBI computes itself from the original reads.

## Running it

From the attest project root, on a machine with internet access and R with
`readxl` and `DESeq2` installed:

```
Rscript validation/geo_study/run_study.R check        # 1. every endpoint, on one real series (1-2 min)
Rscript validation/geo_study/run_study.R pilot        # 2. 20 series, then a summary (~15-30 min)
Rscript validation/geo_study/run_study.R all          # 3. until 200 usable series (hours)
Rscript validation/geo_study/run_study.R consequence  # 4. DESeq2 on flagged files
Rscript validation/geo_study/run_study.R summary      # 5. results/summary.md
```

From RStudio: `setwd("~/Documents/attest"); stage <- "check"; source("validation/geo_study/run_study.R")`.

Everything is resumable: each series is saved as it finishes and skipped on a
re-run. Settings are in `config.R`. Downloads go to `cache/` (author files
nothing flagged are deleted after auditing); results go to `results/`.

**Run `check` first and read its output.** Three things it tests were written
from documentation and could not be tried from where the pipeline was built:
the E-utilities filter `"rnaseq counts"[Filter]` (check prints NCBI's own
reading of the search and the count without the filter), the links on GEO's
download page, and the NCBI annotation file. If `check` fails, stop and
look at its log before running anything else.

## What it measures

**Sample.** The documented E-utilities filter for series with NCBI-generated
counts is tried first; `check` and `sample` confirm from the hit counts that
it actually narrows the search (on the first real run a parsing bug - reading
a per-term count as the total - made it look ignored; the mock now uses
NCBI's real reply layout). If the filter does not cut the search, the frame is
all human series of type "expression profiling by high throughput
sequencing", and whether NCBI has counts is read from each series' own
download page ("no NCBI counts" is then an exclusion with a reported rate).
`check` also confirms that >= 95% of the frame's entries are series, not
samples or platforms. 2000 series are drawn at random with a fixed seed and
processed in that order until 200 have an auditable author matrix. The frame
(search, basis, count, seed, date) is saved. Every series that drops out is
counted with its reason: no NCBI counts, no matrix-like supplementary file,
single-cell only, file too large (> 50 MB), unreadable, samples not matched.
A download page or file NCBI would not serve (bot check, error) is not an
exclusion: it is not saved and is tried again on the next run.

**Author files.** Up to three per series: supplementary `.txt/.tsv/.csv(.gz)`
or `.xlsx` that are not archives, tracks, single-cell or binary. Each is noted
as named like counts (`count`, `raw`, `reads`, `htseq`, ...), like normalised
values (`tpm`, `fpkm`, `cpm`, `norm`, `log`, ...), or neither. Spreadsheets are
converted to text first, as a user would, and that is recorded.

**attest.** `attest_file()` on each file, with no sample sheet - what a user
who has just downloaded the file would get.

**Truth, independent of attest.** Author samples are matched to GSMs (by
accession, then by sample title, then by correlation with a stated margin);
author genes are mapped to NCBI GeneIDs through Ensembl IDs, symbols or
Entrez IDs.
- *Value scale*: slope of log(author sample total) on log(NCBI sample total).
  NCBI's totals measure sequencing depth. Raw counts keep it (slope near 1);
  CPM, TPM, FPKM and normalised values have divided it out (slope near 0).
  Labels: raw counts (> 0.7, whole numbers), estimated counts (> 0.7, not
  whole), depth removed (< 0.3), log-transformed, ambiguous, or
  **undetermined** when fewer than 4 samples match or NCBI depths differ by
  less than 1.3x - never guessed.
- *Completeness*: among genes with zero NCBI reads in every matched sample
  (and nameable in the author's ID system), the share present in the author's
  file. Complete matrices keep them (> 0.6); filtered ones drop them (< 0.3).
- *Sex*: GEO's stated sex against XIST / Y-gene expression in NCBI's counts.

**Consequence.** For files NCBI shows to be non-raw, where the GEO annotation
has a two-level condition with >= 3 samples each: DESeq2 on the author's file
and on NCBI's raw counts, same samples, `~ cond`, padj < 0.05.

## Output

`results/summary.md`: sample flow and exclusions; how often each problem
occurs (files and series, Wilson 95% intervals); attest against the truth
(confusion table, sensitivity, specificity); files named as counts that are
not; completeness against truth; sex mismatches; DESeq2 consequences.
CSVs beside it hold every row.

## Tested before first use

`mock/run_mock.R` builds a mock of GEO and NCBI from real counts (airway; Kang
2018 pseudobulk) with 12 series, each planting a known situation - raw counts,
CPM named `raw_counts`, TPM beside counts, a filtered matrix, htseq summary
rows with unmatched column names, an Excel file with dates for gene names, a
flipped sex label, no matrix, single-cell only, no NCBI counts, DESeq2-
normalised values, log values, and equal depths - and runs every stage.
It also plants NCBI answering a download with a bot-check web page (as
happened on the first real run): the page is never cached as data, and the
annotation falls back to NCBI Gene's `Homo_sapiens.gene_info.gz`; NCBI
ignoring the counts filter (the frame falls back); and a bot check on a
series' download page (retried, not counted as "no NCBI counts").
All 22 checks come out as planted.

## Limits to state with the results

- NCBI's pipeline is not the authors'. The truth rules compare totals and gene
  presence, which survive a change of pipeline; they do not assume the two
  agree gene by gene.
- Completeness and value-scale truth are undetermined for some files; those
  are reported, not dropped.
- The automatic condition choice for the consequence analysis is a heuristic;
  each case is listed so it can be checked by hand.
- Human series only, and only those NCBI has processed.
- If GEO's annotation file cannot be fetched, NCBI Gene's current table is
  used instead; GeneIDs retired since NCBI's counts were made go unmapped.
