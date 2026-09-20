# attest

**Can this data support the analysis you are about to run on it?**

`attest` audits an RNA-seq count matrix *before* differential expression. It does not
analyse expression. It adjudicates one claim at a time, reports the evidence behind the
verdict, and says which checks your inputs could not answer.

```r
library(attest)
attest_counts(counts)
```

```
attest - is this a raw count matrix?
verdict: NOT PERMITTED - These are CPM or TPM values, not raw counts.
why:
  - every sample totals about 1,000,000 (all slightly below, so the values were scaled
    down then rounded), which happens when each column is divided by its own library size
  - 72% of genes are zero in every sample, so rounding has already destroyed the
    low-expression genes
  - small values are spaced about 1.94 apart in 7 of 8 samples, implying original
    libraries of roughly 515,000 reads
what it costs: DESeq2, edgeR and limma will run on this without complaint and report far
fewer differences than the data contain (3,993 -> 877 genes in our airway test).
```

## Why it exists

Count models assume counts. The standard tools check only that the values are whole
numbers, and **rounding defeats that check**. We handed the same airway matrix, in four
forms, to each tool and recorded every message it produced:

| Input | DESeq2 | edgeR | limma-voom | attest |
|---|---|---|---|---|
| raw counts | silent | silent | silent | PERMITTED |
| **rounded CPM** | "converting counts to integer mode" | **silent** | **silent** | NOT PERMITTED |
| **rounded TPM** | "converting counts to integer mode" | **silent** | **silent** | NOT PERMITTED |
| FPKM (non-integer) | ERROR: not integers | **silent** | **silent** | NOT PERMITTED |

The consequence is not noise, it is lost findings. On the airway dataset
(`~ cell + dex`, padj < 0.05):

| Input | DE genes | Shared with the correct run |
|---|---|---|
| raw counts | 3,993 | — |
| DESeq2-normalised, rounded | 4,011 | ~all (+18 extra) |
| rounded CPM | 877 | 876 |
| rounded TPM | 516 | 513 |

A rounded-CPM matrix returns a *subset* of the true gene list: the analysis looks
successful and quietly loses 78% of the signal.

This is not a criticism of DESeq2, edgeR or limma. They document that they expect raw
counts and are not built to police provenance. `attest` checks the precondition they
reasonably assume.

## What it decides, and how

| Signal | Means |
|---|---|
| every column sums to ~1e6 | CPM or TPM |
| variance below the mean across samples | impossible for counts (Poisson floor); the values were divided by depth and/or gene length |
| depth already divided out of non-integer values | FPKM/RPKM or similar |
| median-of-ratios size factors within 1.01x | already normalised |
| negatives, or non-integers with a maximum under 30 | log scale (vst, rlog, log-CPM) |
| regular gaps between small values | scaled and rounded; the gap size recovers the original library size |

Verdicts are `PERMITTED`, `CAUTION`, `NOT PERMITTED` or `UNKNOWN`. `UNKNOWN` means the
inputs cannot answer the question - never a guess. Checks that could not run are listed
under `not assessed`, so a `PERMITTED` verdict never hides a check that never happened.

Nothing throws. Malformed input returns a typed result with a reason.

## Validation

- **39/39** verdicts correct across airway (human), pasilla (fly) and fission (yeast),
  13 transformations each: CPM, TPM, FPKM and their rounded forms, log2(CPM+1),
  size-factor normalised, estimated counts, downsampling, gene filtering.
- **0 false alarms** on 180 simulated raw matrices spanning 9 regimes: typical bulk,
  tightly balanced libraries, n = 4, n = 50, shallow 3' (1M and 0.3M reads), sparse UMI
  pseudobulk, and a 60k-gene annotation.
- **44 assertions** in the test suite, run against 3,000-gene fixtures shipped with the
  package, so the tests need no Bioconductor data packages.

Several rules in the list above replaced earlier ones that real data falsified. The
design log records each: what was assumed, which dataset broke it, and what replaced it.

## Known limits

- A matrix normalised from shallow libraries (< 1M reads) has rounding noise above the
  1.01x size-factor threshold and can be missed.
- Genuine salmon/kallisto estimates from an unusually even experiment can be reported as
  normalised; the message names that reading and points to `tximport`.
- Single-cell and UMI matrices are only simulated so far, not tested on real data.

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
