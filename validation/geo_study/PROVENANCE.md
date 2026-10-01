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
- **What this run actually produced:** reaching the 200-usable target took
  **369 series processed** (`results/summary.md`'s "Generated ... from 369
  processed series"), of which **221** ended up with at least one auditable
  author matrix ("ok"), slightly over the 200 target because usability is
  only checked between series, not file-by-file. The published figures
  (108/341 not-raw-counts, etc.) are counted over the 341 files from those
  221 series.
  **Note:** this is not "processed until 369 reached a definitive outcome" —
  369 is how many series it took to reach the 200-usable stopping point in
  this particular run, not the stopping rule itself. A re-run with the same
  seed should reach the same 200-usable point, but the exact number of
  series it takes to get there is not itself fixed by the config (it depends
  on how many series happen to be usable, in encountered order).

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
