# Everything that touches the network goes through geo_text(), geo_download()
# and geo_size(), so the whole pipeline can be run against a local mock
# (net$mode = "mock", net$map = named list url -> local file) for testing.

net <- new.env()
net$mode <- "live"
net$map  <- list()

with_retry <- function(f, tries = 3, wait = 2) {
  for (i in seq_len(tries)) {
    out <- tryCatch(f(), error = function(e) e)
    if (!inherits(out, "error")) return(out)
    if (i < tries) Sys.sleep(wait * i)
  }
  stop(conditionMessage(out))
}

mock_lookup <- function(url) {
  f <- net$map[[url]]
  if (!is.null(f)) return(f)
  for (k in grep("^re:", names(net$map), value = TRUE))     # pattern keys, for query URLs
    if (grepl(sub("^re:", "", k), url)) return(net$map[[k]])
  NULL
}

geo_text <- function(url) {
  if (net$mode == "mock") {
    f <- mock_lookup(url)
    if (is.null(f)) stop("mock: no page for ", url)
    return(readLines(f, warn = FALSE))
  }
  Sys.sleep(study$sleep)
  with_retry(function() { con <- url(url); on.exit(close(con)); readLines(con, warn = FALSE) })
}

geo_download <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && file.size(dest) > 0) return(dest)          # resumable
  if (net$mode == "mock") {
    f <- mock_lookup(url)
    if (is.null(f)) stop("mock: no file for ", url)
    file.copy(f, dest, overwrite = TRUE); return(dest)
  }
  Sys.sleep(study$sleep)
  tmp <- paste0(dest, ".part")
  with_retry(function() utils::download.file(url, tmp, mode = "wb", quiet = TRUE))
  file.rename(tmp, dest)
  dest
}

geo_size <- function(url) {                                            # bytes, NA if unknown
  if (net$mode == "mock") { f <- mock_lookup(url); return(if (is.null(f)) NA_real_ else file.size(f)) }
  Sys.sleep(study$sleep)
  h <- tryCatch(curlGetHeaders(url), error = function(e) character(0))
  cl <- grep("^content-length:", h, ignore.case = TRUE, value = TRUE)
  if (!length(cl)) return(NA_real_)
  as.numeric(sub("^[^:]+:\\s*", "", utils::tail(cl, 1)))
}

eutils_url <- function(tool, ...) {
  q <- list(...)
  if (nzchar(study$email))   q$email   <- study$email
  if (nzchar(study$api_key)) q$api_key <- study$api_key
  q$tool <- "attest_study"
  paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", tool, ".fcgi?",
         paste(names(q), vapply(q, function(v) utils::URLencode(as.character(v), reserved = TRUE), ""),
               sep = "=", collapse = "&"))
}

# All GDS UIDs for a search. Series UIDs are 2 followed by the zero-padded GSE
# number, so the accession comes straight from the UID.
esearch_gse <- function(term, page = 10000) {
  first <- geo_text(eutils_url("esearch", db = "gds", term = term, retmax = 0))
  n <- as.integer(sub(".*<Count>([0-9]+)</Count>.*", "\\1", paste(first, collapse = "")))
  if (is.na(n) || n == 0) return(character(0))
  ids <- character(0)
  for (start in seq(0, n - 1, by = page)) {
    x <- paste(geo_text(eutils_url("esearch", db = "gds", term = term,
                                   retmax = page, retstart = start)), collapse = "")
    ids <- c(ids, regmatches(x, gregexpr("(?<=<Id>)[0-9]+(?=</Id>)", x, perl = TRUE))[[1]])
  }
  ids <- ids[startsWith(ids, "2")]
  unique(paste0("GSE", as.numeric(ids) - 2e8))
}

gse_stub <- function(gse) {                        # GSE123456 -> GSE123nnn, GSE123 -> GSEnnn
  num <- sub("^GSE", "", gse)
  paste0("GSE", if (nchar(num) > 3) substr(num, 1, nchar(num) - 3) else "", "nnn")
}
ftp_series <- function(gse, sub) sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/%s/", gse_stub(gse), gse, sub)

# file names in an FTP directory listing (HTML index)
list_dir <- function(url) {
  x <- tryCatch(paste(geo_text(url), collapse = "\n"), error = function(e) "")
  h <- regmatches(x, gregexpr('href="[^"]+"', x))[[1]]
  h <- gsub('^href="|"$', "", h)
  h <- h[!grepl("^[?/]|^\\.\\.|^https?:", h) & !endsWith(h, "/")]
  unique(utils::URLdecode(h))
}

# NCBI-generated counts for a series: the raw counts file and the annotation,
# found on the series' download page, with GEO2R's naming as the fallback
ncbi_counts_urls <- function(gse) {
  base <- "https://www.ncbi.nlm.nih.gov"
  page <- tryCatch(paste(geo_text(sprintf("%s/geo/download/?acc=%s", base, gse)), collapse = "\n"),
                   error = function(e) "")
  links <- gsub("&amp;", "&", regmatches(page, gregexpr('href="[^"]*rnaseq_counts[^"]*"', page))[[1]])
  links <- gsub('^href="|"$', "", links)
  full <- ifelse(startsWith(links, "/"), paste0(base, links), links)
  raw   <- full[grepl("raw_counts", full) & grepl("\\.tsv\\.gz", full)]
  annot <- full[grepl("annot\\.tsv\\.gz", full)]
  list(raw   = if (length(raw)) raw[1] else sprintf(
         "%s/geo/download/?format=file&type=rnaseq_counts&acc=%s&file=%s_raw_counts_GRCh38.p13_NCBI.tsv.gz", base, gse, gse),
       annot = if (length(annot)) annot[1] else sprintf(
         "%s/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz", base),
       from_page = length(raw) > 0)
}
