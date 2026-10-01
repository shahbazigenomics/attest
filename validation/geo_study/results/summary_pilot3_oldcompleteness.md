# attest on published GEO data

Generated 2026-09-21 from 20 processed series.

Sampling frame: 27,715 series matching `"rnaseq counts"[Filter] AND "Homo sapiens"[Organism] AND gse[ETYP]`; 2000 drawn at random (seed 20260921).

## Flow

| status | Freq |
|---|---|
| ok | 12 |
| per-sample files only (RAW.tar) |  7 |
| no matrix file |  1 |


Series with at least one audited author matrix: **12**. Files audited: **17**.


## What attest found (as a user would run it: `attest_file()` on the downloaded file)

| finding | files | series |
|---|---|---|
| value scale NOT PERMITTED (not raw counts) | 5 / 17 (29.4%, 95% CI 13.3-53.1) | 4 / 12 (33.3%, 95% CI 13.8-60.9) |
| value scale CAUTION | 2 / 17 (11.8%, 95% CI 3.3-34.3) | 2 / 12 (16.7%, 95% CI 4.7-44.8) |
| counting summary rows left in (htseq/STAR) | 2 / 17 (11.8%, 95% CI 3.3-34.3) | 2 / 12 (16.7%, 95% CI 4.7-44.8) |
| gene names turned into dates (Excel) | 0 / 17 (0.0%, 95% CI 0.0-18.4) | 0 / 12 (0.0%, 95% CI 0.0-24.3) |
| duplicated gene identifiers | 0 / 17 (0.0%, 95% CI 0.0-18.4) | 0 / 12 (0.0%, 95% CI 0.0-24.3) |
| annotation columns in the table (featureCounts) | 8 / 17 (47.1%, 95% CI 26.2-69.0) | 4 / 12 (33.3%, 95% CI 13.8-60.9) |
| matrix filtered upstream (completeness CAUTION) | 5 / 17 (29.4%, 95% CI 13.3-53.1) | 4 / 12 (33.3%, 95% CI 13.8-60.9) |


## attest against the NCBI-based truth (value scale)

Truth comes from NCBI's own raw counts for the same samples, never from attest's rules: per sample, the median of author value / NCBI count over genes NCBI counts >= 50 everywhere. Raw and estimated counts sit within 3-fold of 1; CPM, TPM and FPKM sit near 1e6 / depth (all below 0.2 = depth removed). For non-integer values on the count scale, the slope of that level on NCBI depth (used only if its SE < 0.15) separates estimated counts (~0) from normalised counts (~-1).

Files with a determined truth: 8 of 17 audited.


Why the rest are undetermined:

| reason | Freq |
|---|---|
| author genes could not be mapped to NCBI genes | 3 |
| levels neither all within 3-fold of NCBI's nor all below 0.2 | 1 |
| samples not matched to GSMs (none) | 5 |

| truth | CAUTION | NOT PERMITTED | PERMITTED |
|---|---|---|---|
| depth removed | 0 | 1 | 0 |
| estimated counts | 1 | 0 | 0 |
| raw counts | 0 | 0 | 6 |


- **Sensitivity** - non-raw files attest did not call raw: 1 / 1 (100.0%, 95% CI 20.7-100.0)

- **Specificity** - raw-count files attest called raw: 6 / 6 (100.0%, 95% CI 61.0-100.0)

- **Named as counts but not counts** (file name says count/raw/reads; NCBI comparison says depth removed or log): 0 / 7 (0.0%, 95% CI 0.0-35.4) of files named as counts; attest flagged 0 of them


## Completeness against the NCBI-based truth

| truth | CAUTION | PERMITTED |
|---|---|---|
| complete | 0 | 5 |
| filtered | 2 | 1 |


## Sex: GEO annotation against XIST / Y-gene expression in NCBI's counts

- series with sex annotated for >= 2 samples: 1

- samples whose expression contradicts their label: 0 / 47 (0.0%, 95% CI 0.0-7.6)

- series with at least one such sample: 0 / 1 (0.0%, 95% CI 0.0-79.3)


## Consequence: DESeq2 on the author's file vs on NCBI's raw counts

For files NCBI shows to be non-raw, with a two-level condition (>= 3 samples each) in the GEO annotation. Same samples, same design (`~ cond`), padj < 0.05.

(none)


## Files attest could not read, or failed on

Kept in cache/ for inspection. A file a user would download and that attest cannot read is itself a result.

(none)


## Exclusions

| file_status | Freq |
|---|---|
| audited | 17 |

