# attest on published GEO data - run from the attest project root:
#
#   Rscript validation/geo_study/run_study.R check        # 1. test every endpoint (1 minute)
#   Rscript validation/geo_study/run_study.R pilot        # 2. 20 series, then a summary
#   Rscript validation/geo_study/run_study.R all          # 3. the full sample (hours; resumable)
#   Rscript validation/geo_study/run_study.R consequence  # 4. DESeq2 on flagged files
#   Rscript validation/geo_study/run_study.R summary      # 5. tables -> results/summary.md
#
# Or from RStudio: setwd("~/Documents/attest"); stage <- "check"; source("validation/geo_study/run_study.R")
# Interrupt any time; re-running continues where it stopped.

args  <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) args[1] else if (exists("stage")) stage else "check"

for (f in list.files("R", full.names = TRUE)) source(f)              # attest itself
for (f in c("config.R", "lib_net.R", "lib_parse.R", "lib_truth.R", "stages.R", "summary.R"))
  source(file.path("validation/geo_study", f))
if (exists("study_override")) study[names(study_override)] <- study_override
options(timeout = study$timeout)
dir.create(study$results, recursive = TRUE, showWarnings = FALSE)
dir.create(study$cache, recursive = TRUE, showWarnings = FALSE)

switch(stage,
  check       = stage_check(),
  sample      = stage_sample(),
  pilot       = { stage_run(limit = study$pilot_n); stage_summary() },
  all         = { stage_run(); stage_summary() },
  consequence = stage_consequence(),
  summary     = stage_summary(),
  stop("unknown stage: ", stage))
