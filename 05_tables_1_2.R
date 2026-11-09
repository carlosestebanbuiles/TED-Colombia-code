# 05_tables_1_2.R : Table 1 (baseline characteristics by activity) and Table 2 (ophthalmic findings)
# Dependency-light re-implementation of the gtsummary tables of the original analysis.
suppressPackageStartupMessages(library(tidyverse))
source("00_helpers.R"); d <- readRDS("outputs/ted_analytic_rev.rds")
# Tables 1 and 2 of the original submission were produced with gtsummary, which uses type-2 quantiles;
# the text and the other tables use R's default (type 7). Both conventions are kept so that every
# published number is reproduced exactly.
med_iqr <- function(x, d = 1) { x <- x[!is.na(x)]; if (!length(x)) return("—")
  q <- quantile(x, c(.5, .25, .75), names = FALSE, type = 2); q[abs(q) < 1e-9] <- 0
  sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), q[1], q[2], q[3]) }

cont_row <- function(v, lab, dig = 1) {
  x <- d[[v]]; g <- d$activity
  p <- tryCatch(wilcox.test(x ~ g, exact = FALSE)$p.value, error = function(e) NA)
  tibble(Characteristic = lab, Overall = med_iqr(x, dig), Inactive = med_iqr(x[g %in% "Inactive"], dig),
         Active = med_iqr(x[g %in% "Active"], dig), p = fmt_p(p))
}
bin_row <- function(v, lab, lvl = "Yes") {
  x <- as.character(d[[v]]); g <- d$activity; y <- x == lvl
  t <- table(g, y); p <- tryCatch(if (any(suppressWarnings(chisq.test(t))$expected < 5)) fisher.test(t)$p.value else chisq.test(t, correct = FALSE)$p.value, error = function(e) NA)
  f <- function(k) n_pct(sum(y[k], na.rm = TRUE), sum(!is.na(y[k])))
  tibble(Characteristic = lab, Overall = f(TRUE), Inactive = f(g %in% "Inactive"), Active = f(g %in% "Active"), p = fmt_p(p))
}
cat_rows <- function(v, lab) {
  x <- as.character(d[[v]]); g <- d$activity; lv <- if (is.factor(d[[v]])) levels(d[[v]]) else sort(unique(na.omit(x)))
  t <- table(g, x); p <- tryCatch(if (any(suppressWarnings(chisq.test(t))$expected < 5)) fisher.test(t, simulate.p.value = TRUE, B = 1e5)$p.value else chisq.test(t)$p.value, error = function(e) NA)
  f <- function(k, l) n_pct(sum(x[k] == l, na.rm = TRUE), sum(!is.na(x[k])))
  bind_rows(tibble(Characteristic = lab, Overall = "", Inactive = "", Active = "", p = fmt_p(p)),
            map_dfr(lv, ~ tibble(Characteristic = paste0("  ", .x), Overall = f(TRUE, .x), Inactive = f(d$activity %in% "Inactive", .x),
                                 Active = f(d$activity %in% "Active", .x), p = "")))
}
set.seed(2026)
table1 <- bind_rows(cont_row("age_go_onset", "Age at TED onset, years", 0), bin_row("sex", "Female sex", "Female"), cont_row("bmi", "BMI, kg/m2"),
  cat_rows("smoking", "Smoking status"), bin_row("diabetes", "Diabetes mellitus"), bin_row("other_autoimmune", "Other autoimmune disease"),
  bin_row("fh_thyroid_disease", "Family history of thyroid disease"), bin_row("visible_goiter", "Visible goitre"),
  bin_row("relapse_after_treatment", "Relapse after thyroid treatment"), bin_row("atd_commenced", "Antithyroid drugs ever"),
  bin_row("radioiodine_treatment", "Radioiodine before evaluation"), bin_row("thyroidectomy", "Thyroidectomy before evaluation"),
  cat_rows("go_timing", "TED onset relative to thyroid diagnosis"), cont_row("ft4", "Free T4 at diagnosis, ng/dL", 2),
  cont_row("tsh", "TSH at diagnosis, mU/L", 3), bin_row("trab_status", "TRAb positive", "Positive"), cont_row("trab", "TRAb, IU/L"))

one <- function(v, lab, lvl = "Yes") { x <- as.character(d[[v]]); tibble(Characteristic = lab, Overall = n_pct(sum(x == lvl, na.rm = TRUE), sum(!is.na(x)))) }
onec <- function(v, lab, dig = 0) tibble(Characteristic = lab, Overall = med_iqr(d[[v]], dig))
table2 <- bind_rows(onec("cas", "CAS (7-point), recorded"), one("activity", "Active disease", "Active"),
  one("painful_oppressive_feeling", "Spontaneous retrobulbar pain"), one("gaze_evoked_pain", "Gaze-evoked pain"), one("active_eyelid_swelling", "Eyelid swelling"),
  one("eyelid_erythema", "Eyelid erythema"), one("conjunctival_redness", "Conjunctival redness"), one("chemosis", "Chemosis"), one("caruncle_swelling", "Caruncle swelling"),
  one("grittiness", "Grittiness"), one("photophobia", "Photophobia"), one("excessive_watering", "Excessive tearing"), one("double_vision", "Diplopia (any)"),
  map_dfr(0:3, ~ tibble(Characteristic = paste("Gorman score", .x), Overall = n_pct(sum(d$gorman == .x, na.rm = TRUE), sum(!is.na(d$gorman))))),
  one("proptosis_clin", "Proptosis on examination"), onec("exo_max", "Exophthalmometry, worse eye, mm", 1), one("exo_ge23", "Exophthalmometry >=23 mm"),
  one("exo_asym", "Inter-eye asymmetry >=2 mm"), onec("ulr_max", "Upper eyelid retraction, worse eye, mm"), one("lid_retr_ge2", "Upper eyelid retraction >=2 mm"),
  one("lagophthalmos", "Lagophthalmos"), one("motility_restriction", "Restricted ocular motility"), one("strabismus", "Manifest strabismus"),
  map_dfr(sort(unique(na.omit(d$cornea))), ~ tibble(Characteristic = .x, Overall = n_pct(sum(d$cornea == .x, na.rm = TRUE), sum(!is.na(d$cornea))))),
  onec("va_worse", "BCVA, worse eye (decimal)", 2), one("color_vision", "Colour vision abnormality"), one("rapd", "RAPD"), one("don", "Dysthyroid optic neuropathy"),
  onec("iop_max", "IOP, higher eye, mmHg"), one("iop_gt21", "IOP >21 mmHg"),
  map_dfr(levels(d$sev4), ~ tibble(Characteristic = .x, Overall = n_pct(sum(d$sev4 == .x), nrow(d)))))
writexl::write_xlsx(list(Table1 = table1, Table2 = table2), "outputs/tables_1_2.xlsx")
print(table1, n = Inf); print(table2, n = Inf)
