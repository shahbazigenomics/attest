# Settings for the GEO study. Edit here, not in the stage scripts.
# All paths are relative to the attest project root.

study <- list(
  root     = "validation/geo_study",
  cache    = "validation/geo_study/cache",       # downloads; gitignored
  results  = "validation/geo_study/results",     # per-series results and summary

  # which series: human, with NCBI-generated raw counts, series entries only
  term = '"rnaseq counts"[Filter] AND "Homo sapiens"[Organism] AND "gse"[Entry Type]',

  seed          = 20260921,   # the random sample is fixed by this
  n_candidates  = 400,        # series drawn at random, processed in this order
  target_usable = 200,        # stop once this many series have an auditable author matrix
  pilot_n       = 20,

  max_file_mb        = 50,    # author files larger than this are skipped (and counted)
  max_files_per_gse  = 3,
  keep_author_files  = FALSE, # delete author downloads after auditing to save disk

  # NCBI asks E-utilities users for an email address and allows 3 requests/s
  # without an API key (10 with one). Both optional.
  email   = "",
  api_key = "",
  sleep   = 0.4,
  timeout = 600
)
