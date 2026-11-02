# Scoped rerun: recover completeness truth (and agreement_frac) for the 3
# completeness-miss files whose author downloads were deleted before the
# completeness-agreement_frac fix (kept never considered completeness truth
# until now, so these were dropped before anyone could re-inspect them).
#
# Reprocesses only these 3 series, within the exact same 369-candidate scope
# as the completed 'all' run (target_usable inflated so the usable-count
# early-exit can never trigger before reaching them) - the same safe-rerun
# pattern used for the earlier 24-series crash-fix rerun, so nothing else in
# that scope is silently dropped or re-touched.
#
# Run from the attest project root:
#   Rscript validation/geo_study/_rerun_3files.R

targets <- c("GSE188464", "GSE243702", "GSE75011")

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
