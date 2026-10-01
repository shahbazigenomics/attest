# attest on published GEO data

Generated 2026-09-21 from 20 processed series.

Sampling frame: 27,715 series matching `"rnaseq counts"[Filter] AND "Homo sapiens"[Organism] AND gse[ETYP]`; 2000 drawn at random (seed 20260921).

## Flow

| status | Freq |
|---|---|
| ok | 10 |
| per-sample files only (RAW.tar) |  7 |
| ok, but no file could be audited |  2 |
| no matrix file |  1 |


Series with at least one audited author matrix: **10**. Files audited: **15**.


## What attest found (as a user would run it: `attest_file()` on the downloaded file)

| finding | files | series |
|---|---|---|
| value scale NOT PERMITTED (not raw counts) | 4 / 15 (26.7%, 95% CI 10.9-52.0) | 3 / 10 (30.0%, 95% CI 10.8-60.3) |
| value scale CAUTION | 2 / 15 (13.3%, 95% CI 3.7-37.9) | 2 / 10 (20.0%, 95% CI 5.7-51.0) |
| counting summary rows left in (htseq/STAR) | 2 / 15 (13.3%, 95% CI 3.7-37.9) | 2 / 10 (20.0%, 95% CI 5.7-51.0) |
| gene names turned into dates (Excel) | 0 / 15 (0.0%, 95% CI 0.0-20.4) | 0 / 10 (0.0%, 95% CI 0.0-27.8) |
| duplicated gene identifiers | 0 / 15 (0.0%, 95% CI 0.0-20.4) | 0 / 10 (0.0%, 95% CI 0.0-27.8) |
| annotation columns in the table (featureCounts) | 6 / 15 (40.0%, 95% CI 19.8-64.3) | 2 / 10 (20.0%, 95% CI 5.7-51.0) |
| matrix filtered upstream (completeness CAUTION) | 5 / 15 (33.3%, 95% CI 15.2-58.3) | 4 / 10 (40.0%, 95% CI 16.8-68.7) |


## attest against the NCBI-based truth (value scale)

Truth comes from NCBI's own raw counts for the same samples, never from attest's rules: per sample, the median of author value / NCBI count over genes NCBI counts >= 50 everywhere. Raw and estimated counts sit within 3-fold of 1; CPM, TPM and FPKM sit near 1e6 / depth (all below 0.2 = depth removed). For non-integer values on the count scale, the slope of that level on NCBI depth (used only if its SE < 0.15) separates estimated counts (~0) from normalised counts (~-1).

Files with a determined truth: 4 of 15 audited.


Why the rest are undetermined:

| reason | Freq |
|---|---|
| author genes could not be mapped to NCBI genes | 3 |
| levels neither all within 3-fold of NCBI's nor all below 0.2 | 1 |
| samples not matched to GSMs (none) | 7 |

| truth | NOT PERMITTED | PERMITTED |
|---|---|---|
| depth removed | 1 | 0 |
| raw counts | 0 | 3 |


- **Sensitivity** - non-raw files attest did not call raw: 1 / 1 (100.0%, 95% CI 20.7-100.0)

- **Specificity** - raw-count files attest called raw: 3 / 3 (100.0%, 95% CI 43.8-100.0)

- **Named as counts but not counts** (file name says count/raw/reads; NCBI comparison says depth removed or log): 0 / 3 (0.0%, 95% CI 0.0-56.2) of files named as counts; attest flagged 0 of them


## Completeness against the NCBI-based truth

| truth | PERMITTED |
|---|---|
| complete | 3 |
| filtered | 1 |


## Sex: GEO annotation against XIST / Y-gene expression in NCBI's counts

- series with sex annotated for >= 2 samples: 1

- samples whose expression contradicts their label: 0 / 47 (0.0%, 95% CI 0.0-7.6)

- series with at least one such sample: 0 / 1 (0.0%, 95% CI 0.0-79.3)


## Consequence: DESeq2 on the author's file vs on NCBI's raw counts

For files NCBI shows to be non-raw, with a two-level condition (>= 3 samples each) in the GEO annotation. Same samples, same design (`~ cond`), padj < 0.05.

(none)


## Files attest could not read, or failed on

Kept in cache/ for inspection. A file a user would download and that attest cannot read is itself a result.

| gse | file | status | message |
|---|---|---|---|
| GSE115255 | GSE115255_read_counts.tsv.gz | not read by attest | The file could not be parsed as a ';'-separated table. |
| GSE162669 | GSE162669_tpm_summary_exvivo_AMTAM.txt.gz | not read by attest | Only 0 numeric columns left after removing identifiers and annotation, so there is nothing to compare. |


## Exclusions

| file_status | Freq |
|---|---|
| audited | 15 |
| not read by attest |  2 |

