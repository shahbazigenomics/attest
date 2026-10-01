# attest on published GEO data

Generated 2026-09-21 from 20 processed series.

Sampling frame: 27,715 series matching `"rnaseq counts"[Filter] AND "Homo sapiens"[Organism] AND gse[ETYP]`; 2000 drawn at random (seed 20260921).

## Flow

| status | Freq |
|---|---|
| ok | 18 |
| error: subscript out of bounds |  1 |
| no matrix file |  1 |


Series with at least one audited author matrix: **18**. Files audited: **21**.


## What attest found (as a user would run it: `attest_file()` on the downloaded file)

| finding | files | series |
|---|---|---|
| value scale NOT PERMITTED (not raw counts) | 4 / 21 (19.0%, 95% CI 7.7-40.0) | 3 / 18 (16.7%, 95% CI 5.8-39.2) |
| value scale CAUTION | 2 / 21 (9.5%, 95% CI 2.7-28.9) | 2 / 18 (11.1%, 95% CI 3.1-32.8) |
| counting summary rows left in (htseq/STAR) | 1 / 21 (4.8%, 95% CI 0.8-22.7) | 1 / 18 (5.6%, 95% CI 1.0-25.8) |
| gene names turned into dates (Excel) | 0 / 21 (0.0%, 95% CI 0.0-15.5) | 0 / 18 (0.0%, 95% CI 0.0-17.6) |
| duplicated gene identifiers | 0 / 21 (0.0%, 95% CI 0.0-15.5) | 0 / 18 (0.0%, 95% CI 0.0-17.6) |
| annotation columns in the table (featureCounts) | 3 / 21 (14.3%, 95% CI 5.0-34.6) | 1 / 18 (5.6%, 95% CI 1.0-25.8) |
| matrix filtered upstream (completeness CAUTION) | 5 / 21 (23.8%, 95% CI 10.6-45.1) | 4 / 18 (22.2%, 95% CI 9.0-45.2) |


## attest against the NCBI-based truth (value scale)

Truth comes from NCBI's own raw counts for the same samples: the slope of log(author total) on log(NCBI total) across samples - near 1 when sequencing depth is still in the values, near 0 when it has been divided out. It uses none of attest's rules.

Files with a determined truth: 1 of 21 audited (the rest: too few matched samples, too-even depths, or samples not matched).

| truth | PERMITTED |
|---|---|
| raw counts | 1 |


- **Sensitivity** - non-raw files attest did not call raw: -

- **Specificity** - raw-count files attest called raw: 1 / 1 (100.0%, 95% CI 20.7-100.0)

- **Named as counts but not counts** (file name says count/raw/reads; NCBI comparison says depth removed or log): 0 / 1 (0.0%, 95% CI 0.0-79.3) of files named as counts; attest flagged 0 of them


## Completeness against the NCBI-based truth

| truth | CAUTION | PERMITTED |
|---|---|---|
| complete | 0 | 3 |
| filtered | 2 | 1 |


## Sex: GEO annotation against XIST / Y-gene expression in NCBI's counts

- series with sex annotated for >= 2 samples: 3

- samples whose expression contradicts their label: 0 / 82 (0.0%, 95% CI 0.0-4.5)

- series with at least one such sample: 0 / 3 (0.0%, 95% CI 0.0-56.2)


## Consequence: DESeq2 on the author's file vs on NCBI's raw counts

For files NCBI shows to be non-raw, with a two-level condition (>= 3 samples each) in the GEO annotation. Same samples, same design (`~ cond`), padj < 0.05.

(none)


## Exclusions

| file_status | Freq |
|---|---|
| attest error |  2 |
| audited | 21 |

