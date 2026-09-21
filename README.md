# attest

**Can this data support the analysis you are about to run on it?**

`attest` audits an RNA-seq count matrix and its design *before* differential expression.
It does not analyse expression. It adjudicates one claim at a time, reports the evidence
behind each verdict, and says which claims your inputs could not answer.

```r
library(attest)

attest(dds)                 # a DESeqDataSet already carries the counts, the
                            # sample sheet and the design - nothing else needed
attest_file("counts.txt")   # or the file you were sent, before it is anything else
attest(counts, metadata = colData, design = ~ cell + dex)
```

```
attest - can this data support the analysis?
overall: PERMITTED

[value scale] PERMITTED - Consistent with raw counts.
  - all values are whole numbers; sample totals differ 1.99x and size factors 2.12x, as
    expected when depth has not been divided out
  - variance/mean is 2.29, above the Poisson floor

[identifiers] PERMITTED - The row names identify genes, one row each.
  - all 3,000 gene rows carry Ensembl identifiers

[completeness] PERMITTED - Complete: the matrix still contains genes with no reads at all.
  - 1,446 of 3,000 genes (48.2%) are zero in every sample; filtering removes exactly
    those, so they are the sign that nothing was removed
  not assessed: gene-selection check (protein-coding only, a panel, or another subset):
    needs the annotation size, so pass n_expected = <genes in your GTF>

[identity] PERMITTED - Samples are consistent with their labels.
  - sex inferred from expression: 2 female, 6 male
  - no pair stands out as the same library (highest 0.9415, cohort median 0.9295)
  not assessed: DNA-RNA match: this check compares samples with each other and with the
    sample sheet, not with genotypes; for that, run somalier or NGSCheckMate on the BAMs

[design] PERMITTED - The design can estimate the effect of 'dex'.
  - 'dex' is balanced against the other terms (variance inflation 1.00x, no effective
    loss of samples)
  - group sizes: untrt = 4, trt = 4

-- what this data can support (does not change the verdict above) --

[detectability] CAUTION - A 2.0-fold change was detectable for under 80% of genes.
  - of the 716 genes with a mean count of 10 or more, 60% could have shown a 2.0-fold
    change at 80% power
  - among those, the median gene needed 1.9-fold and the quietest quarter 2.4-fold or more
  - the 838 genes below that count are effectively untestable here (median 65-fold needed)
  - group sizes 4 and 4; about 49 genes really differ, so Benjamini-Hochberg judges each
    gene at p < 0.0016 (Bonferroni would be 3.2e-05)
  cost: A gene absent from your results list is not evidence that it does not respond;
    for genes above their threshold above, it is evidence, and for the rest it is not.
```

(That run is on the 3,000-gene airway fixture shipped with the package, so the gene counts
are smaller than a real annotation would give. Every check is also callable on its own:
`attest_counts()`, `attest_identifiers()`, `attest_completeness()`, `attest_identity()`,
`attest_design()`, `attest_detectability()`. Some `not assessed` lines are trimmed above.)

A check with nothing to work on is **named, not dropped**. `attest(counts)` on a bare matrix
ends with:

```
-- not run, so the verdict above does not cover it --
  - design adequacy (is the effect estimable, and what is it worth after confounding?):
    no sample sheet - pass metadata = <data.frame> and design = ~ <terms>, or hand
    attest() a DESeqDataSet
  - detectability (for which genes could a change have been seen at all?): the two
    groups being compared are not known - pass group = <factor>
```

## Getting your counts in

`attest()` takes a matrix, a `data.frame`, a `SummarizedExperiment`, a `DESeqDataSet`, a
`DGEList` or a `tximport` list, and reads the sample sheet, the design formula and the
grouping straight out of the object when it has them. Anything you pass explicitly wins.

`attest_file()` is for counts that are still a file - a featureCounts table, a GEO
supplementary file, a CSV from a collaborator. It reports how it read the file as the
report's first check: the separator, `#` header lines, gzip, which column held the gene
identifiers, and any duplicated identifiers.

One of those matters more than it looks. featureCounts writes `Chr`, `Start`, `End`,
`Strand` and `Length` between the identifier and the samples. `Start`, `End` and `Length`
are numeric, so a matrix built by hand from that file has three extra "libraries" in it -
and **every check in this package still returns PERMITTED**, because nothing in the numbers
gives it away once the column names are gone. It is catchable at the file and nowhere else,
which is the argument for reading the file here rather than before.

## Why it exists

Count models assume counts, and count-based designs assume the effect is separable from
everything else in the sample sheet. The standard tools check a little of the first and
almost none of the second. We handed the same airway matrix, in four forms, to each tool
and recorded every message it produced (`validation/compare_tools.R`):

| Input | DESeq2 | edgeR | limma-voom | attest |
|---|---|---|---|---|
| raw counts | silent | silent | silent | PERMITTED |
| **rounded CPM** | "converting counts to integer mode" | **silent** | **silent** | NOT PERMITTED |
| **rounded TPM** | "converting counts to integer mode" | **silent** | **silent** | NOT PERMITTED |
| FPKM (non-integer) | ERROR: not integers | **silent** | **silent** | NOT PERMITTED |

Rounding defeats the only check the tools make, and the consequence is not noise, it is
lost findings. On the airway dataset (`~ cell + dex`, padj < 0.05):

| Input | DE genes | Shared with the correct run |
|---|---|---|
| raw counts | 3,993 | — |
| DESeq2-normalised, rounded | 4,011 | ~all (+18 extra) |
| rounded CPM | 877 | 876 |
| rounded TPM | 516 | 513 |

A rounded-CPM matrix returns a *subset* of the true gene list: the analysis looks
successful and quietly loses 78% of the signal.

The same exercise on designs (`validation/compare_design.R`): DESeq2 and edgeR refuse a
*completely* confounded design, limma-voom fits it and returns NA coefficients with a
warning. Partial confounding, and two samples per group, are accepted in silence by all
three.

This is not a criticism of DESeq2, edgeR or limma. They document that they expect raw
counts and a full-rank design; they are not built to police provenance. `attest` checks
the preconditions they reasonably assume.

## The checks

| Check | The claim it adjudicates | What the tools do about it |
|---|---|---|
| `attest_file()` | "this file was read as the count matrix it is" | nothing; annotation columns read as samples are invisible afterwards |
| `attest_counts()` | "this matrix holds raw counts" | DESeq2 rejects non-integers only; rounding defeats it; edgeR and limma silent |
| `attest_identifiers()` | "every row is a gene, named once, from one annotation" | nothing; htseq-count's `__no_feature` row is analysed as a gene like any other |
| `attest_completeness()` | "this matrix holds every gene the pipeline quantified" | nothing; a gene absent from the matrix is indistinguishable from a gene never expressed |
| `attest_identity()` | "these samples are who the sample sheet says" | nothing at the count-matrix level |
| `attest_design()` | "this design can estimate the effect, and nothing in the sample sheet stands in its way" | full-rank refusal only (DESeq2, edgeR), and only for terms in the formula; nothing on columns left out of it |
| `attest_detectability()` | "this gene did not change" | nothing retrospective; PROPER, ssizeRNA and RNASeqPower are prospective planning tools |

### How each one decides

**Value scale.** Six signals, in order of how conclusive they are:

| Signal | Means |
|---|---|
| every column sums to ~1e6 | CPM or TPM |
| variance below the mean across samples | impossible for counts (Poisson floor); the values were divided by depth and/or gene length |
| depth already divided out of non-integer values | FPKM/RPKM or similar |
| median-of-ratios size factors much flatter than the column sums | already normalised |
| negatives, or non-integers with a maximum under 30 | log scale (vst, rlog, log-CPM) |
| regular gaps between small values | scaled and rounded; the gap size recovers the original library size |

**Identifiers.** Three faults, none of them visible in the numbers. htseq-count appends
`__no_feature`, `__ambiguous`, `__too_low_aQual`, `__not_aligned` and
`__alignment_not_unique`; STAR's `ReadsPerGene.out.tab` opens with four `N_*` rows. Left
in, they are analysed as genes, and library sizes, size factors and the multiple-testing
denominator are all computed partly from them. How much that costs depends on the run, so
the check does not quote a figure: it measures the share those rows hold in *your*
libraries and reports it. Then Excel: `SEPT1` becomes
`1-Sep`, `MARCH1` becomes `1-Mar`, which is why HGNC renamed both families in 2020
(`SEPTIN`, `MARCHF`). Then identifiers from two sources - some Ensembl IDs versioned and
some bare, or Ensembl IDs mixed with symbols - where a join keeps one group and drops the
other without saying so.

Two rules here were narrowed by the fixtures before they ever ran on real user data. A
trailing `.N` is not evidence of an Ensembl version: every fission yeast systematic name
ends that way (`SPAC212.09c`, `SPAC977.03`), so the rule requires an Ensembl stem. And an
all-digit identifier is not an Excel serial number, because Entrez gene IDs are all digits;
only the textual date forms are reported. ERCC and SIRV spike-ins are counted separately,
so they do not look like a second annotation.

**Completeness.** A complete matrix keeps genes that are zero in every sample; filtering
removes exactly those, and the floor it leaves on the row totals reveals the threshold
that was used. The *fraction* of all-zero genes is useless as a test — it is 47.4% on
airway, 15.3% on pasilla, 4.0% on fission — so only their presence is used.

**Identity.** Sex is inferred from expression (XIST against seven Y genes), never from X
heterozygosity: after X-inactivation a clonal female sample expresses one X, so a
DNA-style rule calls it male. A duplicate library is flagged when a pair is identical, or
correlates above 0.999, or sits at least 0.02 above every other pair in the cohort.
Neither test alone works: yeast replicates reach 0.9955, so a fixed 0.99 floor cries wolf,
while a true duplicate there is only ~1.5 MAD above the median, so a cohort-relative rule
alone misses it.

**Design.** Rank first (which coefficients are not estimable, named), then group sizes,
then the part no analysis tool can see: **the columns of the sample sheet that the formula
leaves out**. DESeq2 only ever sees the terms it was given, so a sequencing run that holds
every treated sample, or two libraries per donor analysed as four independent replicates,
passes in silence. Both are caught here without any threshold on an association measure,
because at n = 8 such measures are mostly noise:

- a categorical column is flagged when every one of its values falls inside a single group
  (it is *nested* in the condition), which is what batch confounding and repeated
  samples from one subject both look like - the message says it cannot tell which;
- a numeric column (RIN, concentration, date) is flagged when its ranges in the two groups
  do not overlap at all, and the message states how often an unrelated column would do
  that by chance: `2 / choose(n, n1)`, 1 time in 35 at 4 + 4, 1 in 10 at 3 + 3. Measured
  over 40,000 simulated covariates: 0.0273 and 0.0986.

Identifiers, constants, relabelings of the condition, a row index on a sheet sorted by
condition, and columns computed from the counts (`sizeFactor`, `lib.size`) are skipped and
named as skipped.

Confounding *among* the formula's own terms is different: the effect is still estimated
without bias, only less precisely. It is reported as an effective sample size - 8 samples
at variance inflation 2.5 buy the precision of 3.2 - and charged in detectability, which
recomputes every gene's detectable change at that size. It does not lower the verdict:
precision is scope, not a fault.

**Detectability.** For every gene, the smallest fold change this dataset could have
detected, from a two-sample negative-binomial Wald calculation using that gene's own mean
and dispersion, the group sizes, and the level the gene is *really* judged at. That level
is the Benjamini-Hochberg one, `alpha * R / n`, with R estimated from the counts; using
Bonferroni instead put the detectable change above 2-fold for every gene in a dataset
where DESeq2 finds 3,993 differentially expressed ones.

## Faults and scope

The first four checks look for faults: something is wrong with the data or the design, and
the overall verdict is the worst of them. `detectability` is different — nothing is wrong,
the experiment is simply as big as it is — so it reports **scope**, printed separately and
never lowering the overall verdict. FastQC flags that are normal for RNA-seq, MultiQC's
refusal to give an overall pass, and the limitations section of a clinical report all make
the same separation.

Verdicts are `PERMITTED`, `CAUTION`, `NOT PERMITTED` or `UNKNOWN`. `UNKNOWN` means the
inputs cannot answer the question - never a guess. Claims that could not be checked are
listed under `not assessed`, so a `PERMITTED` verdict never hides a check that never ran.

Nothing throws. Malformed input returns a typed result with a reason. There are no bundled
annotation databases.

`attest_as_list()` returns the same content as a plain list for `jsonlite::toJSON()` or a
pipeline step that acts on the verdict.

## Validation

- **39/39** value-scale verdicts correct across airway (human), pasilla (fly) and fission
  (yeast), 13 transformations each: CPM, TPM, FPKM and their rounded forms, log2(CPM+1),
  size-factor normalised, estimated counts, downsampling, gene filtering.
- **0 false alarms** on 180 simulated raw matrices spanning 9 regimes: typical bulk,
  tightly balanced libraries, n = 4, n = 50, shallow 3' (1M and 0.3M reads), sparse UMI
  pseudobulk, and a 60k-gene annotation.
- **Detectability calibrated by simulation**: data generated with exactly the fold change
  the check calls detectable at 80% power was detected 78-84% of the time, at n = 3, 4, 6
  and 10 per group.
- **122 assertions** in the test suite, run against 3,000-gene fixtures shipped with the
  package, so the tests need no Bioconductor data packages.

Several rules above replaced earlier ones that real data falsified: an FPKM rule that
worked on human and failed on fly and yeast, a rounding-lattice rule that does not survive
20M-read libraries, a size-factor threshold that turned out to be depth-dependent. The
design log records each one: what was assumed, which dataset broke it, what replaced it.

## Known limits

- A matrix normalised from shallow libraries (< 1M reads) carries rounding noise that can
  push it past the size-factor test and be read as raw counts.
- Genuine salmon/kallisto estimates from an unusually even experiment can be reported as
  normalised; the message names that reading and points to `tximport`.
- Single-cell and UMI matrices are simulated only, not yet tested on real data.
- A numeric sample-sheet column is tested against two-group comparisons only, and a
  numeric column that is really a batch code (`1, 1, 2, 2, ...`) is treated as a
  number, not as groups; code batches as text to get the nesting test.
- Detectability's effective sample size charges for confounding but not for the precision
  a paired or blocked design gains, so for such designs it errs conservative.
- Detectability is an approximation to what DESeq2 or edgeR would achieve, not a
  reimplementation of either; it is calibrated against simulation, not against their output.
- Sex inference and the DNA-RNA question: the markers are human, and no genotype
  comparison is attempted. For that, run somalier or NGSCheckMate on the BAMs.

## Install

```r
# install.packages("remotes")
remotes::install_github("shahbazigenomics/attest")
```

## Related

[`admissible`](https://github.com/shahbazigenomics/admissible) applies the same idea to
exome data: it does not interpret variants, it audits whether the data can support an
interpretation.

MIT licensed.
