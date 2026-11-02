# Provenance: GEO-study validation

Quick-reference record for reviewers. Full methodology is in `README.md`;
this file exists so the key reproducibility facts don't have to be dug out
of prose.

## How the sample was drawn

- **Search query (E-utilities, `db=gds`):**
  `"rnaseq counts"[Filter] AND "Homo sapiens"[Organism] AND gse[ETYP]`
  (`config.R`'s `term`; `check`/`sample` confirm from the hit counts that this
  filter actually narrows the search before it is trusted — see README for
  the fallback frame used if it doesn't).
- **Random seed:** `20260921` (`config.R`'s `seed`).
- **Selection rule:** up to 2000 series are drawn at random with that seed
  and then processed **in that fixed order**, not until a fixed count of
  series is processed, but until a target number have an auditable author
  matrix (`config.R`'s `target_usable = 200`). Every series tried along the
  way is logged with its outcome, including ones that drop out (no NCBI
  counts, no matrix file, single-cell only, file too large, unreadable).
- **What this run actually produced:** the `all` stage stopped at the 200-usable
  target after **369 series processed** (log line `series processed: 369, with an
  auditable author matrix: 200`). 20 series had crashed attest during that run;
  after those crashes were fixed they were re-processed within the same 369 series
  (scoped reruns with `target_usable` inflated so nothing outside that scope was
  touched; `_rerun_3files.R` and `_rerun_6misses.R` are two of them), which raised
  the total to **221** series with at least one auditable author matrix. The
  published figures (108/341 not-raw-counts, etc.) are counted over the 341 files
  from those 221 series.
  **Correction:** an earlier version of this file said 221 overshot 200 "because
  usability is only checked between series". That was wrong: `stage_run()` checks
  the count after every series and stops at exactly 200. The overshoot came from
  the scoped reruns, and the results therefore mix code versions (series
  processed before and after the crash fixes).

## Re-run on frozen code (2026-11-01)

The same sample (`sample.csv` copied unchanged, so GEO's growth since 2026-09-21
cannot change which series are drawn) was re-run from scratch with attest frozen
at commit `ce4e4d7`, into `results_rerun_2026-11/`.

- It stopped at exactly 200 usable series, after 376 series: 53 author files
  (37 series) could not be downloaded (NCBI connections cut off at ~60 s, then
  refused with HTTP 403), so the run went further down the fixed order to reach
  200. Its own totals (111/300 not raw, etc.) are **not** comparable to the
  original: the 37 lost series were mostly raw-count files (59% PERMITTED in the
  original run, against 27% for the rest), so losing them inflates the non-raw
  share.
- What it does establish is file-by-file agreement
  (`results_rerun_2026-11/concordance_vs_2026-09-22.csv`). On the 291 files
  audited in both runs: value scale 290/291 identical (GSE82184: NOT PERMITTED
  -> UNKNOWN), identifiers 290/291 (GSE196688: PERMITTED -> CAUTION),
  completeness, identity and the NCBI-based truth 291/291. The other 50 files of
  the original 341 were not re-checked.
- What it does not establish: performance on new data. Several rules were
  refined on this sample, so agreement between two runs on it is a
  reproducibility check, not a held-out test.

Harness changes made for this re-run (study code only; nothing under `R/`):

- `lib_net.R` `list_dir()`: NCBI's HTTPS directory index served empty listings
  for directories that hold files (GSE180883/suppl: index empty, file
  downloads with status 200). An empty index now falls back to GEO's plain-text
  series record (`acc.cgi?...&form=text`), which lists every supplementary file
  and platform. The HTML series page was not used: it returns a reCAPTCHA to
  scripts. If neither source answers, the series is a transient failure.
- `stages.R`: a listing that could not be established, or an author file that
  could not be downloaded, is now a transient failure (retried on the next run),
  not a saved "no supplementary files" / "download failed". The 2026-11-01
  re-run results were produced before the download change, which is why they
  record 53 "download failed" files.
- The `check` stage fails on today's frame because the first series it tests
  (GSE250293) genuinely has empty `matrix/` and `suppl/` directories; the
  series in the sample were checked by hand and are unaffected.

## Rerun time (from `README.md`)

| stage | what it does | rough time |
|---|---|---|
| `check` | every endpoint, on one real series | 1-2 min |
| `pilot` | 20 series, then a summary | ~15-30 min |
| `all` | until 200 usable series (the full run) | hours |
| `consequence` | DESeq2 on flagged files | not separately timed |
| `summary` | rebuild `results/summary.md` from saved per-series results | seconds |

`all` is resumable: each series is saved as it finishes and skipped on a
re-run, so an interrupted run can continue rather than restart.

## Requirements

Running `check`/`pilot`/`all` sources `lib_parse.R`, which uses `data.table`
and `readxl`; `summary.R`/`stages.R` use `DESeq2` for the consequence stage.
All three (plus `edgeR`) are declared in `DESCRIPTION`'s `Suggests`.
(`data.table` and `readxl` were missing from `Suggests` when this file was
first written and have since been added.)

## What was and wasn't verified when this file was written

- Sourcing `run_study.R` end-to-end (package loads, path resolution) was
  confirmed to work from a clean checkout with no hardcoded absolute paths.
- A full network run was **not** completed in this check: the sandbox used
  to verify this blocks outbound requests to `eutils.ncbi.nlm.nih.gov`, so
  `check` got through its package checks and failed at the first live
  E-utilities call. This is an environment restriction on the checking
  machine, not a defect in the script — but it means the live-network path
  (GEO/NCBI fetches) has not been independently re-run end-to-end since the
  2026-09-22 results were generated, only re-read and reasoned about.
