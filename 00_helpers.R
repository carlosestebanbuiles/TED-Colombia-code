# 00_helpers.R : shared helper functions (text cleaning, dates, formatting)
suppressPackageStartupMessages(library(tidyverse))
# ---- helpers (identical to the original Rmd) -----------------------------------------
na_strings <- c("", "na", "n/a", "nan", "sin dato", "no aplica", "sin datos", "nd")
txt <- function(x) { x <- str_squish(as.character(x)); x[str_to_lower(x) %in% na_strings] <- NA_character_; x }
yn  <- function(x) { s <- str_to_lower(txt(x))
  case_when(s %in% c("yes","y","si","sí","1","true") ~ "Yes", s %in% c("no","n","0","false","o") ~ "No", TRUE ~ NA_character_) }
num <- function(x, lo = -Inf, hi = Inf) {
  s <- str_to_lower(txt(x)); s <- str_replace_all(s, ",", "."); s[which(s == "cinco")] <- "5"
  v <- suppressWarnings(as.numeric(s)); need <- is.na(v) & !is.na(s)
  v[need] <- suppressWarnings(as.numeric(str_extract(s[need], "-?\\d+(\\.\\d+)?")))
  v[!is.na(v) & (v < lo | v > hi)] <- NA_real_; v }
safe_date <- function(y, m, d) as.Date(sprintf("%04d-%02d-%02d", y, m, d), format = "%Y-%m-%d")
parse_date_any <- function(x) {
  s <- txt(x); out <- rep(as.Date(NA), length(s)); n <- suppressWarnings(as.numeric(s))
  i <- !is.na(n) & n >= 1900 & n <= 2030 & n == round(n); out[i] <- safe_date(as.integer(n[i]), 7L, 1L)
  i <- !is.na(n) & n > 2030 & n < 60000; out[i] <- as.Date(floor(n[i]), origin = "1899-12-30")
  m <- str_match(s, "^0?0/0?0/(\\d{4})$"); i <- !is.na(m[, 1]); out[i] <- safe_date(as.integer(m[i, 2]), 7L, 1L)
  m <- str_match(s, "^0?0/(0?[1-9]|1[0-2])/(\\d{4})$"); i <- !is.na(m[, 1]); out[i] <- safe_date(as.integer(m[i, 3]), as.integer(m[i, 2]), 15L)
  m <- str_match(s, "^(\\d{1,2})/(\\d{1,2})/(\\d{4})$"); i <- !is.na(m[, 1]) & is.na(out)
  a <- as.integer(m[i, 2]); b <- as.integer(m[i, 3]); y <- as.integer(m[i, 4])
  out[i] <- safe_date(y, ifelse(b > 12, a, b), ifelse(b > 12, b, a))
  m <- str_match(s, "^(\\d{4})-(\\d{1,2})-(\\d{1,2})"); i <- !is.na(m[, 1]) & is.na(out)
  out[i] <- safe_date(as.integer(m[i, 2]), as.integer(m[i, 3]), as.integer(m[i, 4])); out }
months_between <- function(a, b) as.numeric(b - a) / 30.44
any_yes <- function(...) { m <- cbind(...)
  ifelse(rowSums(m == "Yes", na.rm = TRUE) > 0, "Yes", ifelse(rowSums(!is.na(m)) > 0, "No", NA_character_)) }
fmt_p <- function(p) ifelse(is.na(p), "—", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
med_iqr <- function(x, d = 1) { x <- x[!is.na(x)]; if (!length(x)) return("—")
  q <- quantile(x, c(.5, .25, .75), names = FALSE); q[abs(q) < 1e-9] <- 0; sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), q[1], q[2], q[3]) }
n_pct <- function(k, n, d = 1) sprintf(paste0("%d/%d (%.", d, "f%%)"), k, n, 100 * k / n)
wilson <- function(k, n) { ci <- prop.test(k, n, correct = FALSE)$conf.int; sprintf("%.1f–%.1f", 100 * ci[1], 100 * ci[2]) }
clean_names_local <- function(x) {   # same behaviour as janitor::clean_names for these headers
  nm <- names(x); nm <- str_replace_all(nm, "([a-z0-9])([A-Z])", "\\1_\\2")
  nm <- str_replace_all(nm, "([A-Z])([A-Z][a-z])", "\\1_\\2"); nm <- str_to_lower(nm)
  nm <- str_replace_all(nm, "[^a-z0-9]+", "_"); nm <- str_remove_all(nm, "^_|_$"); nm <- make.unique(nm, sep = "_"); names(x) <- nm; x }


# ---- study-specific settings (override with environment variables when running on other data) ----
# Centres in which the registry entry date is the baseline encounter (baseline-to-last-visit interval computable)
BASELINE_CENTRES <- trimws(strsplit(Sys.getenv("TED_BASELINE_CENTRES", "Centre A"), ",")[[1]])
# Centres excluded jointly in the sensitivity analysis of incomplete documentation
EXCLUDE_SET <- trimws(strsplit(Sys.getenv("TED_EXCLUDE_CENTRES", "Centre B,Centre D"), ",")[[1]])
EXCLUDE_LABEL <- if (identical(EXCLUDE_SET, c("Centre B", "Centre D"))) "Centres B and D" else paste(EXCLUDE_SET, collapse = " + ")
