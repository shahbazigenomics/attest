# attest on published GEO data

Generated 2026-09-22 from 369 processed series.

Sampling frame: 27,715 series matching `"rnaseq counts"[Filter] AND "Homo sapiens"[Organism] AND gse[ETYP]`; 2000 drawn at random (seed 20260921).

## Flow

| status | Freq |
|---|---|
| ok | 221 |
| per-sample files only (RAW.tar) | 120 |
| no matrix file |  15 |
| ok, but no file could be audited |  12 |
| single-cell only |   1 |


Series with at least one audited author matrix: **221**. Files audited: **341**.


## What attest found (as a user would run it: `attest_file()` on the downloaded file)

| finding | files | series |
|---|---|---|
| value scale NOT PERMITTED (not raw counts) | 108 / 341 (31.7%, 95% CI 27.0-36.8) | 80 / 221 (36.2%, 95% CI 30.1-42.7) |
| value scale CAUTION | 80 / 341 (23.5%, 95% CI 19.3-28.2) | 65 / 221 (29.4%, 95% CI 23.8-35.7) |
| counting summary rows left in (htseq/STAR) | 8 / 341 (2.3%, 95% CI 1.2-4.6) | 7 / 221 (3.2%, 95% CI 1.5-6.4) |
| gene names turned into dates (Excel) | 31 / 341 (9.1%, 95% CI 6.5-12.6) | 21 / 221 (9.5%, 95% CI 6.3-14.1) |
| duplicated gene identifiers | 0 / 341 (0.0%, 95% CI 0.0-1.1) | 0 / 221 (0.0%, 95% CI 0.0-1.7) |
| annotation columns in the table (featureCounts) | 113 / 341 (33.1%, 95% CI 28.4-38.3) | 77 / 221 (34.8%, 95% CI 28.9-41.3) |
| matrix filtered upstream (completeness CAUTION) | 123 / 341 (36.1%, 95% CI 31.2-41.3) | 87 / 221 (39.4%, 95% CI 33.2-45.9) |


## attest against the NCBI-based truth (value scale)

Truth comes from NCBI's own raw counts for the same samples, never from attest's rules: per sample, the median of author value / NCBI count over genes NCBI counts >= 50 everywhere. Raw and estimated counts sit within 3-fold of 1; CPM, TPM and FPKM sit near 1e6 / depth (all below 0.2 = depth removed). For non-integer values on the count scale, the slope of that level on NCBI depth (used only if its SE < 0.15) separates estimated counts (~0) from normalised counts (~-1).

Files with a determined truth: 129 of 341 audited.


Why the rest are undetermined:

| reason | Freq |
|---|---|
| author genes could not be mapped to NCBI genes |  67 |
| fewer than 200 genes counted >= 50 by NCBI in every matched sample |   4 |
| levels neither all within 3-fold of NCBI's nor all below 0.2 |   9 |
| samples not matched to GSMs (none) | 130 |
| samples not matched to GSMs (title) |   2 |

| truth | CAUTION | NOT PERMITTED | PERMITTED |
|---|---|---|---|
| count scale, not whole numbers |  7 |  2 |  0 |
| depth removed |  5 | 22 |  1 |
| estimated counts | 15 |  0 |  0 |
| log-transformed |  0 |  9 |  0 |
| raw counts |  3 |  0 | 65 |


- **Sensitivity** - non-raw files attest did not call raw: 36 / 37 (97.3%, 95% CI 86.2-99.5)

- **Specificity** - raw-count files attest called raw: 65 / 68 (95.6%, 95% CI 87.8-98.5)

- **Named as counts but not counts** (file name says count/raw/reads; NCBI comparison says depth removed or log): 2 / 78 (2.6%, 95% CI 0.7-8.9) of files named as counts; attest flagged 1 of them


## Completeness against the NCBI-based truth

| truth | CAUTION | PERMITTED |
|---|---|---|
| complete |  7 | 89 |
| filtered | 30 |  9 |


## Sex: GEO annotation against XIST / Y-gene expression in NCBI's counts

- series with sex annotated for >= 2 samples: 20

- samples whose expression contradicts their label: 45 / 519 (8.7%, 95% CI 6.5-11.4)

- series with at least one such sample: 7 / 20 (35.0%, 95% CI 18.1-56.7)


## Consequence: DESeq2 on the author's file vs on NCBI's raw counts

For files NCBI shows to be non-raw, with a two-level condition (>= 3 samples each) in the GEO annotation. Same samples, same design (`~ cond`), padj < 0.05.

| gse | file | labelled_as | attest_value | condition | n_per_group | de_ncbi_raw | de_author_file | de_both | de_ncbi_on_shared_genes |
|---|---|---|---|---|---|---|---|---|---|
| GSE122524 | GSE122524_gene_matrix.txt.gz | unlabelled | NOT PERMITTED | ch_days.of.differentiation | 3+3 | 2116 |  405 |  237 | 1940 |
| GSE131147 | GSE131147_Tagraxofusp_resistance_RNA_rpkm.txt.gz | normalised | NOT PERMITTED | ch_treatment | 6+6 |  886 |   24 |   24 |  736 |
| GSE134886 | GSE134886_gene_all_fpkm.csv.gz | normalised | NOT PERMITTED | ch_tissue | 3+3 |    7 |    0 |    0 |    7 |
| GSE180883 | GSE180883_HTseq-out.txt.gz | counts | PERMITTED | ch_treatment | 3+3 | 6307 |  980 |  558 | 5346 |
| GSE200339 | GSE200339_Normalized_Filtered5CPM_Gene_Counts.txt.gz | normalised | CAUTION | ch_genotype | 3+6 | 1688 |  288 |  269 | 1335 |
| GSE233826 | GSE233826_RNASeq_processed_HT29.txt.gz | unlabelled | NOT PERMITTED | ch_treatment | 3+3 | 4357 | 3313 | 2920 | 3676 |
| GSE86430 | GSE86430_normalized_counts_DESeq.txt.gz | normalised | CAUTION | ch_ethnicity | 4+12 |   26 |   16 |   12 |   21 |


## Files attest could not read, or failed on

Kept in cache/ for inspection. A file a user would download and that attest cannot read is itself a result.

| gse | file | status | message |
|---|---|---|---|
| GSE116580 | GSE116580_TERC-RACE-pileupdata_KOandP.txt.gz | not read by attest | Only 1 numeric column left after removing identifiers and annotation, so there is nothing to compare. |
| GSE122986 | GSE122986_expressed_gene_reads.txt.gz | not read by attest | The file could not be parsed as a 'tab'-separated table. |
| GSE122986 | GSE122986_expressed_gene_FPKM.txt.gz | not read by attest | The file could not be parsed as a 'tab'-separated table. |
| GSE124519 | GSE124519_All_Differentially_Expressed_Transcripts.xlsx | not read by attest | The file could not be parsed as a ','-separated table. |
| GSE133420 | GSE133420_annotated-circRNA_result.xlsx | not read by attest | Only 0 numeric columns left after removing identifiers and annotation, so there is nothing to compare. |
| GSE137480 | GSE137480_Differentially_Expressed_Genes.xlsx | not read by attest | The file could not be parsed as a ','-separated table. |
| GSE137480 | GSE137480_Expression_Gene.xlsx | not read by attest | The file could not be parsed as a ','-separated table. |
| GSE138636 | GSE138636_expressed_gene_reads.txt.gz | not read by attest | The file could not be parsed as a 'tab'-separated table. |
| GSE138636 | GSE138636_expressed_gene_FPKM.txt.gz | not read by attest | The file could not be parsed as a 'tab'-separated table. |
| GSE138691 | GSE138691_R2_Sznajder_Supplementary_Data_1.xlsx | not read by attest | Only 1 numeric column left after removing identifiers and annotation, so there is nothing to compare. |
| GSE141945 | GSE141945_RNAseq.metadata.csv.gz | not read by attest | Only 0 numeric columns left after removing identifiers and annotation, so there is nothing to compare. |
| GSE147657 | GSE147657_data.csv.gz | not read by attest | Only 1 numeric column left after removing identifiers and annotation, so there is nothing to compare. |
| GSE149200 | GSE149200_CircRNA_Expression_Profiling.xlsx | not read by attest | The file could not be parsed as a ' '-separated table. |
| GSE151000 | GSE151000_H2107_ASCL1_ANP_KD_fpkm_whole_genome.txt.gz | not read by attest | The file could not be parsed as a 'tab'-separated table. |
| GSE151000 | GSE151000_hSCLC_PDX_fpkm_whole_genome.txt.gz | not read by attest | The file could not be parsed as a 'tab'-separated table. |
| GSE155432 | GSE155432_polya.length.tail.tools_revised.csv.gz | not read by attest | The file could not be parsed as a ';'-separated table. |
| GSE159220 | GSE159220_Annotated_circRNA.txt.gz | not read by attest | Only 0 numeric columns left after removing identifiers and annotation, so there is nothing to compare. |
| GSE160252 | GSE160252_intron.spanning.rawcounts.txt.gz | not read by attest | The file could not be parsed as a ' '-separated table. |
| GSE186895 | GSE186895_Sample_description.txt.gz | not read by attest | Only 0 numeric columns left after removing identifiers and annotation, so there is nothing to compare. |
| GSE247120 | GSE247120_CircRNA_Expression_Profiling.xlsx | not read by attest | The file could not be parsed as a ' '-separated table. |
| GSE63966 | GSE63966_ProcessedDataMatrix_Final.txt.gz | not read by attest | Only 0 numeric columns left after removing identifiers and annotation, so there is nothing to compare. |
| GSE65655 | GSE65655_G401-sg-HOXB13_FPKM.xlsx | not read by attest | Only 0 numeric columns left after removing identifiers and annotation, so there is nothing to compare. |
| GSE65655 | GSE65655_G401-sg-control_FPKM.xlsx | not read by attest | Only 0 numeric columns left after removing identifiers and annotation, so there is nothing to compare. |


## Exclusions

| file_status | Freq |
|---|---|
| audited | 341 |
| download failed |   1 |
| not read by attest |  23 |

