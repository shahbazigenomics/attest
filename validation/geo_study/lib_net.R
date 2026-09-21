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

# A downloaded file is only kept if it is what was asked for: NCBI sometimes
# answers with a web page (a bot check, an error page) and a 200 status.
# Such a page is never cached - it would be re-read as data on every re-run.
looks_like_page <- function(path) {
  b <- readBin(path, "raw", 512)
  if (!length(b)) return(TRUE)
  if (grepl("\\.gz$", path)) return(!(length(b) >= 2 && b[1] == as.raw(0x1f) && b[2] == as.raw(0x8b)))
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) return(grepl("^\\s*<", rawToChar(b[b != as.raw(0)]), useBytes = TRUE))
  grepl("^\\s*<(!doctype|html|\\?xml|head)", rawToChar(b[b != as.raw(0)]), ignore.case = TRUE, useBytes = TRUE)
}

geo_download <- function(url, dest) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && file.size(dest) > 0) {                       # resumable
    if (!looks_like_page(dest)) return(dest)
    unlink(dest)                                                        # a cached web page: fetch again
  }
  if (net$mode == "mock") {
    f <- mock_lookup(url)
    if (is.null(f)) stop("mock: no file for ", url)
    file.copy(f, dest, overwrite = TRUE)
    if (looks_like_page(dest)) { unlink(dest); stop("NCBI returned a web page instead of the file: ", url) }
    return(dest)
  }
  tmp <- paste0(dest, ".part")
  for (i in 1:3) {
    Sys.sleep(study$sleep * i^2)
    with_retry(function() utils::download.file(url, tmp, mode = "wb", quiet = TRUE))
    if (!looks_like_page(tmp)) { file.rename(tmp, dest); return(dest) }
    unlink(tmp)
    Sys.sleep(5 * i)                                                    # back off before asking again
  }
  stop("NCBI returned a web page (bot check or error) instead of the file: ", url)
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

# NCBI's reading of a search: the count, how it translated the term, and any
# warnings (e.g. a phrase or filter it did not recognise and silently dropped)
esearch_info <- function(term, retmax = 5, retstart = 0) {
  x <- paste(geo_text(eutils_url("esearch", db = "gds", term = term, retmax = retmax, retstart = retstart)), collapse = "")
  tag <- function(t) { m <- regmatches(x, regexpr(sprintf("<%s>.*?</%s>", t, t), x, perl = TRUE))
                       if (length(m)) gsub("<[^>]+>", " ", m) else "" }
  list(count = as.integer(sub(".*<Count>([0-9]+)</Count>.*", "\\1", x)),
       translation = trimws(tag("QueryTranslation")),
       warnings = trimws(gsub("\\s+", " ", paste(tag("WarningList"), tag("ErrorList")))),
       ids = regmatches(x, gregexpr("(?<=<Id>)[0-9]+(?=</Id>)", x, perl = TRUE))[[1]])
}

# Share of a search's UIDs that are series (UIDs 2xxxxxxxx), from 3 x 200 UIDs
# at the start, middle and end of the result list
series_share <- function(term, n) {
  at <- unique(pmax(0, c(0, floor(n / 2), n - 200)))
  ids <- unlist(lapply(at, function(s) esearch_info(term, retmax = 200, retstart = s)$ids))
  if (!length(ids)) return(NA_real_)
  mean(startsWith(ids, "2") & nchar(ids) == 9)
}

# Which search defines the sampling frame: the counts filter if NCBI applies
# it, otherwise term_frame. Decided from counts, never assumed.
resolve_frame <- function() {
  f <- esearch_info(study$term)
  base_term <- trimws(sub('"rnaseq counts"\\[Filter\\]\\s*AND', "", study$term))
  b <- esearch_info(base_term)
  filter_works <- !is.na(f$count) && !is.na(b$count) && f$count > 0 && f$count < 0.9 * b$count
  term <- if (filter_works) study$term else study$term_frame
  n <- if (filter_works) f$count else esearch_info(term)$count
  list(term = term, basis = if (filter_works) "NCBI counts filter" else "all human expression-by-sequencing series; NCBI counts checked per series",
       n = n, filtered = f, unfiltered = b)
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
       from_page = length(raw) > 0,
       # the page itself came back (not a bot check or an error), so a missing link means no counts
       page_ok = nzchar(page) && grepl(gse, page, fixed = TRUE) && !grepl("recaptcha", page, ignore.case = TRUE))
}
