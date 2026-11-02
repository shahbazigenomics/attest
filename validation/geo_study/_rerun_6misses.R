# Scoped rerun: get agreement_frac/n_representable for the 6 completeness-miss
# files (GSE155336 x2, GSE198609, GSE233826 x3) computed through the real
# pipeline (stages.R's actual match_columns()/truth_completeness() call)
# instead of the earlier ad-hoc hand-check script, so every number in the
# paper comes from the same code path. Their author files are already on
# disk (kept = TRUE before this fix), so this needs no new downloads - only
# NCBI raw counts / series matrix / annotation are re-touched, all already
# cached. Same safe-rerun pattern as the 24-series and 3-file reruns before
# it: only these 3 series' cached results are deleted, and target_usable is
# inflated so the usable-count early exit can never strand anything else in
# the same 369-candidate scope.
#
# Run from the attest project root:
#   Rscript validation/geo_study/_rerun_6misses.R

targets <- c("GSE155336", "GSE198609", "GSE233826")

for (f in list.files("R", full.names = TRUE)) source(f)
for (f in c("config.R", "lib_net.R", "lib_parse.R", "lib_truth.R", "stages.R", "summary.R"))
  source(file.path("validation/geo_study", f))
options(timeout = study$timeout)
dir.create(study$results, recursive = TRUE, showWarnings = FALSE)
dir.create(study$cache, recursive = TRUE, showWarnings = FALSE)

for (g in targets) {
  f <- series_file(g)
  if (file.exists(f)) { file.remove(f); cat("removed cached result for", g, "\n") }
  else cat("no cached result found for", g, "(already missing?)\n")
}

study$target_usable <- 100000   # never let the usable-count check trigger an early exit
stage_run(limit = 369)          # same 369-candidate scope as the completed 'all' run
stage_summary()

cat("\ndone - check results/files.csv and results/summary.md for the 3 target series.\n")
