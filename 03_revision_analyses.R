# =====================================================================================
# 03_revision_analyses.R : analyses added in response to the reviewers (OED-26-0109)
# Every object stored in the list R is used in the revised manuscript, tables,
# supplementary material or the response letter. Run after 01_derive.R.
# =====================================================================================
suppressPackageStartupMessages({library(tidyverse); library(lme4); library(geepack)})
source("00_helpers.R"); source("00_firth.R")
d <- readRDS("outputs/ted_analytic.rds"); R <- list(); options(width = 200)
cat_h <- function(x) cat("\n==== ", x, " ====\n")

# ---------------------------------------------------------------------------------------
# A. Activity measurement: recorded CAS vs clinician label (Rev 1 #6; Rev 2 Q2; Rev 3 #9)
# ---------------------------------------------------------------------------------------
cat_h("A1 source of baseline activity classification")
R$act_source_bl <- d |> filter(!is.na(activity)) |> count(activity, activity_source) |> pivot_wider(names_from = activity_source, values_from = n, values_fill = 0)
print(R$act_source_bl)
R$act_source_paired <- d |> filter(paired_activity) |> count(activity_source, activity_source_fu)
print(R$act_source_paired)
R$act_source_paired_centre <- d |> filter(paired_activity) |> count(centre, activity_source) |> pivot_wider(names_from = activity_source, values_from = n, values_fill = 0)
print(R$act_source_paired_centre)
R$label_only_centre <- d |> filter(activity_source == "Clinician label") |> count(centre, activity)
print(R$label_only_centre)

cat_h("A2 CAS values >7 truncated")
R$cas_gt7 <- d |> filter(cas_raw > 7 | cas_fu_raw > 7) |> select(id, centre, cas_raw, cas_fu_raw, activity, activity_fu); print(R$cas_gt7)

cat_h("A3 CAS items documented")
R$cas_items_n <- d |> count(cas_items_n); print(R$cas_items_n)
R$cas_items_6 <- d |> filter(cas_items_n == 6) |> select(id, centre, cas, cas_items, activity); print(R$cas_items_6)

# weighted kappa (own implementation; identical to irr::kappa2 with weight = "squared")
kappa_w <- function(a, b, weight = c("unweighted", "squared")) {
  weight <- match.arg(weight); lv <- sort(union(a, b)); k <- length(lv)
  O <- table(factor(a, lv), factor(b, lv)) / length(a); E <- outer(rowSums(O), colSums(O))
  W <- if (weight == "squared") 1 - (outer(seq_len(k), seq_len(k), "-"))^2 / (k - 1)^2 else diag(k)
  (sum(W * O) - sum(W * E)) / (1 - sum(W * E))
}
ag  <- d |> filter(!is.na(cas), !is.na(cas_items))
ag7 <- d |> filter(!is.na(cas), !is.na(cas_items7))
ag2 <- d |> filter(!is.na(cas), !is.na(clinically_active_go))
R$agreement <- tibble(
  comparison = c("Recorded vs item-derived CAS (≥6 items), quadratic-weighted κ", "Recorded vs item-derived CAS ≥3 (≥6 items), κ",
                 "Recorded vs item-derived CAS (all 7 items), quadratic-weighted κ", "Recorded vs item-derived CAS ≥3 (all 7 items), κ",
                 "Recorded CAS ≥3 vs clinician active label, κ"),
  n = c(nrow(ag), nrow(ag), nrow(ag7), nrow(ag7), nrow(ag2)),
  kappa = round(c(kappa_w(ag$cas, ag$cas_items, "squared"), kappa_w(ag$cas >= 3, ag$cas_items >= 3),
                  kappa_w(ag7$cas, ag7$cas_items7, "squared"), kappa_w(ag7$cas >= 3, ag7$cas_items7 >= 3),
                  kappa_w(ag2$cas >= 3, ag2$clinically_active_go == "Yes")), 2),
  exact_pct = round(c(mean(ag$cas == ag$cas_items), NA, mean(ag7$cas == ag7$cas_items7), NA, NA) * 100, 1),
  within1_pct = round(c(mean(abs(ag$cas - ag$cas_items) <= 1), NA, mean(abs(ag7$cas - ag7$cas_items7) <= 1), NA, NA) * 100, 1),
  mean_diff = round(c(mean(ag$cas - ag$cas_items), NA, mean(ag7$cas - ag7$cas_items7), NA, NA), 3))
print(R$agreement)
# discordance at the threshold between CAS >=3 and clinician label
R$cas_vs_label_tab <- table(cas_ge3 = ag2$cas >= 3, label_active = ag2$clinically_active_go == "Yes"); print(R$cas_vs_label_tab)

cat_h("A4 sensitivity analyses restricted to directly recorded CAS")
R$cas_only_baseline <- d |> filter(!is.na(activity_cas)) |> summarise(n = n(), active = sum(activity_cas == "Active"), pct = n_pct(active, n), ci = wilson(active, n)); print(R$cas_only_baseline)
pcas <- d |> filter(!is.na(activity_cas), !is.na(activity_cas_fu))
R$cas_only_paired <- pcas |> summarise(n = n(), active_bl = sum(activity_cas == "Active"), active_fu = sum(activity_cas_fu == "Active"),
                                       inact_n = sum(activity_cas == "Active"), inactivated = sum(inactivated_cas == 1, na.rm = TRUE),
                                       inact_pct = n_pct(inactivated, inact_n), inact_ci = wilson(inactivated, inact_n),
                                       reactivated = sum(activity_cas == "Inactive" & activity_cas_fu == "Active"),
                                       p_mcnemar = fmt_p(binom.test(sum(activity_cas == "Active" & activity_cas_fu == "Inactive"),
                                                                    sum(activity_cas != activity_cas_fu))$p.value))
print(R$cas_only_paired)
dd_cas <- drop_na(d, y_active_cas, age10, male, current_smoker, rai)
R$firth_cas_only <- tidy_firth(firth_fit(y_active_cas ~ age10 + male + current_smoker + rai, dd_cas), "Recorded CAS only (no clinician labels)")
print(R$firth_cas_only |> mutate(across(c(OR, lo, hi), ~ round(.x, 2)), p = fmt_p(p)))
# item-derived CAS with all 7 items as outcome (n with 7 items & no recorded CAS)
d <- d |> mutate(activity_items7 = factor(if_else(cas_items7 >= 3, "Active", "Inactive"), levels = c("Inactive", "Active")),
                 y_active_items7 = as.numeric(activity_items7 == "Active"))
dd_it <- drop_na(d, y_active_items7, age10, male, current_smoker, rai)
R$firth_items7 <- tidy_firth(firth_fit(y_active_items7 ~ age10 + male + current_smoker + rai, dd_it), "Item-derived CAS (all 7 items) as outcome")
print(R$firth_items7 |> mutate(across(c(OR, lo, hi), ~ round(.x, 2)), p = fmt_p(p)))
R$items7_vs_primary <- table(primary = d$activity, items7 = d$activity_items7, useNA = "ifany"); print(R$items7_vs_primary)

# ---------------------------------------------------------------------------------------
# B. Severity denominators (Rev 2 Q3)
# ---------------------------------------------------------------------------------------
cat_h("B severity denominators")
R$don_missing <- d |> filter(is.na(don)) |> select(id, centre, orbitopathy_evidence, cas, activity, sev4); print(R$don_missing)
R$don_cb_overlap <- table(DON = d$don, corneal_breakdown = d$corneal_breakdown, useNA = "ifany"); print(R$don_cb_overlap)
R$cornea_assessed <- sum(!is.na(d$cornea)); cat("Corneal assessment available:", R$cornea_assessed, "\n")
R$sight_threat_components <- d |> filter(sev4 == "Sight-threatening (DON/corneal breakdown)") |> count(don, corneal_breakdown); print(R$sight_threat_components)
R$sev_by_activity <- table(d$sev4, d$activity, useNA = "ifany"); print(R$sev_by_activity)

# ---------------------------------------------------------------------------------------
# C. Follow-up interval, attrition, DON (Rev 1 #4, #10, #11; Rev 2 Q4; Rev 3 #2, #8, #13)
# ---------------------------------------------------------------------------------------
cat_h("C1 baseline-to-last-visit interval (Centre A only)")
R$interval_centreA <- d |> filter(centre %in% BASELINE_CENTRES) |> summarise(n_centre = n(), n_with_interval = sum(!is.na(months_bl_to_fu)),
  median_iqr = med_iqr(months_bl_to_fu), min = round(min(months_bl_to_fu, na.rm = TRUE), 1), max = round(max(months_bl_to_fu, na.rm = TRUE), 1),
  lt6 = sum(months_bl_to_fu < 6, na.rm = TRUE), m6_12 = sum(months_bl_to_fu >= 6 & months_bl_to_fu < 12, na.rm = TRUE),
  m12_24 = sum(months_bl_to_fu >= 12 & months_bl_to_fu < 24, na.rm = TRUE), ge24 = sum(months_bl_to_fu >= 24, na.rm = TRUE))
print(R$interval_centreA)
R$interval_centreA_paired <- d |> filter(centre %in% BASELINE_CENTRES, paired_activity) |> summarise(n = n(), n_int = sum(!is.na(months_bl_to_fu)), median_iqr = med_iqr(months_bl_to_fu)); print(R$interval_centreA_paired)
R$interval_centreA_active <- d |> filter(centre %in% BASELINE_CENTRES, !is.na(inactivated)) |>
  summarise(n = n(), n_int = sum(!is.na(months_bl_to_fu)), median_iqr = med_iqr(months_bl_to_fu), inactivated = sum(inactivated)); print(R$interval_centreA_active)
# Dates actually spanned by the registry (for the Methods: location and period)
R$date_ranges <- tibble(
  item = c("Registry entry (Centre A, prospective at baseline encounter)", "Registry entry / abstraction (Centres B–E)",
           "Last recorded visit (all)", "Reported ocular symptom onset (all)"),
  min = c(min(d$date_of_inclusion[d$centre %in% BASELINE_CENTRES], na.rm = TRUE), min(d$date_of_inclusion[!d$centre %in% BASELINE_CENTRES], na.rm = TRUE),
          min(d$date_of_followup, na.rm = TRUE), min(d$eye_symptom_onset, na.rm = TRUE)),
  max = c(max(d$date_of_inclusion[d$centre %in% BASELINE_CENTRES], na.rm = TRUE), max(d$date_of_inclusion[!d$centre %in% BASELINE_CENTRES], na.rm = TRUE),
          max(d$date_of_followup, na.rm = TRUE), max(d$eye_symptom_onset, na.rm = TRUE)))
print(R$date_ranges)
R$fu_year_centre <- table(centre = d$centre, year = year(d$date_of_followup), useNA = "ifany"); print(R$fu_year_centre)

cat_h("C2 disease duration and inactivation by duration strata")
act <- d |> filter(!is.na(inactivated))
R$duration_active <- act |> summarise(n = n(), n_dur = sum(!is.na(months_go_to_fu)), med = med_iqr(months_go_to_fu)); print(R$duration_active)
R$inact_by_duration <- act |> filter(!is.na(months_go_to_fu)) |>
  mutate(stratum = cut(months_go_to_fu, c(-Inf, 12, 36, Inf), labels = c("<12 months", "12–36 months", ">36 months"))) |>
  group_by(stratum) |> summarise(n = n(), inactivated = sum(inactivated), pct = n_pct(inactivated, n), ci = wilson(inactivated, n), .groups = "drop")
print(R$inact_by_duration)
R$inact_by_duration_all <- d |> filter(paired_activity, !is.na(months_go_to_fu)) |>
  mutate(stratum = cut(months_go_to_fu, c(-Inf, 12, 36, Inf), labels = c("<12 months", "12–36 months", ">36 months"))) |>
  group_by(stratum) |> summarise(n = n(), active_bl = sum(activity == "Active"), active_fu = sum(activity_fu == "Active"), .groups = "drop")
print(R$inact_by_duration_all)
# Centre A: inactivation by true interval
R$inact_by_interval_A <- act |> filter(!is.na(months_bl_to_fu)) |>
  mutate(stratum = cut(months_bl_to_fu, c(-Inf, 6, 12, Inf), labels = c("<6 months", "6–12 months", ">12 months"))) |>
  group_by(stratum) |> summarise(n = n(), inactivated = sum(inactivated), .groups = "drop"); print(R$inact_by_interval_A)
# Exploratory Firth: inactivation ~ log(disease duration)
R$firth_inact_duration <- tidy_firth(firth_fit(inactivated ~ log(months_go_to_fu), drop_na(act, months_go_to_fu)), "Inactivation ~ log(months since symptom onset)")
print(R$firth_inact_duration |> mutate(across(c(OR, lo, hi), ~ round(.x, 2)), p = fmt_p(p)))

cat_h("C3 paired (n=109) vs not paired (n=38)")
cmp_tbl <- function(data, grp) {
  g <- data[[grp]]
  cont <- function(x, lab, dig = 0) { p <- tryCatch(wilcox.test(x ~ g, exact = FALSE)$p.value, error = function(e) NA)
    tibble(Characteristic = lab, !!levels(g)[1] := med_iqr(x[g == levels(g)[1]], dig), !!levels(g)[2] := med_iqr(x[g == levels(g)[2]], dig), p = fmt_p(p)) }
  bin <- function(x, lab, lvl = "Yes") { xx <- x == lvl; t <- table(g, xx)
    p <- tryCatch(if (any(t < 5)) fisher.test(t)$p.value else chisq.test(t)$p.value, error = function(e) NA)
    k1 <- sum(xx[g == levels(g)[1]], na.rm = TRUE); n1 <- sum(!is.na(xx[g == levels(g)[1]])); k2 <- sum(xx[g == levels(g)[2]], na.rm = TRUE); n2 <- sum(!is.na(xx[g == levels(g)[2]]))
    tibble(Characteristic = lab, !!levels(g)[1] := n_pct(k1, n1), !!levels(g)[2] := n_pct(k2, n2), p = fmt_p(p)) }
  cat_ <- function(x, lab) { t <- table(g, x); p <- tryCatch(if (any(t < 5)) fisher.test(t, simulate.p.value = TRUE, B = 1e5)$p.value else chisq.test(t)$p.value, error = function(e) NA)
    bind_rows(tibble(Characteristic = lab, !!levels(g)[1] := "", !!levels(g)[2] := "", p = fmt_p(p)),
              map_dfr(colnames(t), function(l) tibble(Characteristic = paste0("  ", l), !!levels(g)[1] := n_pct(t[1, l], sum(t[1, ])), !!levels(g)[2] := n_pct(t[2, l], sum(t[2, ])), p = ""))) }
  bind_rows(cont(data$age_go_onset, "Age at TED onset, years"), bin(data$sex, "Female sex", "Female"), cat_(data$centre, "Centre"),
            cat_(as.character(data$smoking), "Smoking status"), bin(data$radioiodine_treatment, "Prior radioiodine"),
            bin(as.character(data$activity), "Active disease at baseline", "Active"), cont(data$cas, "Recorded CAS", 0),
            bin(data$activity_source, "Activity from clinician label", "Clinician label"),
            cat_(as.character(data$sev4), "Severity category"), bin(data$don, "DON at baseline"), bin(data$gc_iv, "IV methylprednisolone (administered or planned)"),
            bin(data$gc_any, "Any systemic glucocorticoid"), bin(data$second_line, "Second-line therapy"), bin(data$orbital_surgery, "Surgery (done or planned)"),
            bin(if_else(is.na(data$exo_max), "No", "Yes"), "Exophthalmometry recorded"), bin(if_else(is.na(data$trab_status), "No", "Yes"), "TRAb recorded"))
}
d <- d |> mutate(paired_grp = factor(if_else(paired_activity, "Paired follow-up (n=109)", "No paired follow-up (n=38)"), levels = c("Paired follow-up (n=109)", "No paired follow-up (n=38)")))
R$paired_vs_not <- cmp_tbl(d, "paired_grp"); print(R$paired_vs_not, n = Inf)

cat_h("C4 DON 2x2 and new DON cases")
R$don_2x2 <- d |> filter(!is.na(don), !is.na(optic_neuropathy_followup)) |> count(don, optic_neuropathy_followup) |>
  mutate(status = case_when(don == "Yes" & optic_neuropathy_followup == "Yes" ~ "Persistent", don == "Yes" ~ "Resolved/not recorded at last visit",
                            optic_neuropathy_followup == "Yes" ~ "New-onset", TRUE ~ "Absent at both")); print(R$don_2x2)
R$don_bl_no_fu <- d |> filter(don == "Yes", is.na(optic_neuropathy_followup)) |> nrow(); cat("Baseline DON without follow-up DON record:", R$don_bl_no_fu, "\n")
R$new_don_cases <- d |> filter(new_don == "Yes") |>
  transmute(case = paste0("Case ", row_number()), centre, sex, age = age_go_onset, baseline_CAS = cas, baseline_activity = activity, baseline_severity = sev4,
            baseline_BCVA_worse = va_worse, lastvisit_BCVA_worse = va_worse_fu, colour_vision_fu = color_vision_followup, rapd_fu = rapd_followup,
            lastvisit_CAS = cas_fu, lastvisit_activity = activity_fu, clinician_rated_activity = chg_activity, clinician_rated_severity = chg_severity,
            prior_systemic_GC = gc_prior, IVMP_planned_baseline = ivmp_planned, GC_during_followup = details_go_treatment_since_last,
            months_onset_to_last_visit = round(months_go_to_fu, 1), decompression, radiotherapy, mycophenolate, tocilizumab)
print(as.data.frame(R$new_don_cases))
R$don_baseline_treat <- d |> filter(don == "Yes") |> summarise(n = n(), ivmp = sum(gc_iv == "Yes"), decompression = sum(decompression == "Yes"), fu_don = paste(optic_neuropathy_followup, collapse = ",")); print(R$don_baseline_treat)
# visual function in DON-related items at baseline (for diagnostic parameters statement)
R$don_bl_features <- d |> filter(don == "Yes") |> summarise(n = n(), va_recorded = sum(!is.na(va_worse)), va_lt05 = sum(va_worse < 0.5, na.rm = TRUE),
  colour_recorded = sum(!is.na(color_vision)), colour_abn = sum(color_vision == "Yes", na.rm = TRUE), rapd_recorded = sum(!is.na(rapd)), rapd_yes = sum(rapd == "Yes", na.rm = TRUE)); print(R$don_bl_features)

cat_h("C5 reactivation / clinician-rated worsening")
R$reactivation <- d |> filter(paired_activity, activity == "Inactive") |> summarise(n = n(), reactivated = sum(activity_fu == "Active"), ci = wilson(reactivated, n),
  clinician_worse_activity = sum(chg_activity == "Worse", na.rm = TRUE), rated_n = sum(!is.na(chg_activity)), fu_treatment_yes = sum(active_go_treatment_since_last == "Yes", na.rm = TRUE)); print(R$reactivation)

# ---------------------------------------------------------------------------------------
# D. Treatment exposure windows (Rev 1 #8, #9; Rev 2 Q5; Rev 3 #10)
# ---------------------------------------------------------------------------------------
cat_h("D1 exposure windows")
win <- function(x, lab) tibble(Treatment = lab, n_yes = sum(x == "Yes", na.rm = TRUE), n_known = sum(!is.na(x)), `n/N (%)` = n_pct(n_yes, n_known))
R$exposure_windows <- bind_rows(
  win(d$gc_prior, "Systemic glucocorticoids: recorded in history before/at referral (route not specified)"),
  win(d$gc_planned, "Systemic glucocorticoids: planned at referral visit (any route)"),
  win(d$ivmp_planned, "Intravenous methylprednisolone: planned at referral visit"),
  win(if_else(d$gc_planned == "Yes" & d$gc_route == "Oral", "Yes", if_else(is.na(d$gc_planned), NA, "No")), "Oral glucocorticoids: planned at referral visit"),
  win(d$gc_followup, "Systemic glucocorticoids: administration documented in follow-up notes"),
  win(d$ivmp_followup, "Intravenous methylprednisolone: administration documented in follow-up notes"),
  win(d$gc_any, "Systemic glucocorticoids: composite (recorded or planned, any window)"),
  win(d$gc_iv, "Intravenous methylprednisolone: composite (recorded or planned, any window)"),
  win(if_else(d$gc_iv == "Yes" & (d$ivmp_followup %in% "Yes" | d$active_go_treatment_since_last %in% "Yes"), "Yes", "No"), "  of which active TED treatment documented after referral"),
  win(if_else(d$gc_iv == "Yes" & d$active_go_treatment_since_last %in% "No" & !(d$ivmp_followup %in% "Yes"), "Yes", "No"), "  of which no active TED treatment documented after referral"),
  win(d$topical_eye_treatment, "Topical ocular treatment (history)"),
  win(d$selenium_initial, "Selenium (initial medication field)"), win(d$selenium, "Selenium (composite incl. notes)"),
  win(d$orbital_irradiation, "Orbital radiotherapy: history before referral"), win(d$irradiation_plan, "Orbital radiotherapy: planned at referral"),
  win(if_else(str_detect(d$det_fu, "radioterap"), "Yes", if_else(is.na(d$active_go_treatment_since_last), NA, "No")), "Orbital radiotherapy: documented in follow-up notes"),
  win(d$radiotherapy, "Orbital radiotherapy: composite"),
  win(d$surgery_go, "Surgery for TED: history before referral"), win(d$surgery_plan, "Surgery for TED: planned at referral"),
  win(if_else(str_detect(d$det_fu, "descompres|tarsorraf|estrabismo|blefarotom|oculopl|evisc"), "Yes", if_else(is.na(d$active_go_treatment_since_last), NA, "No")), "Surgery for TED: documented in follow-up notes"),
  win(d$orbital_surgery, "Surgery for TED: composite"),
  win(d$mycophenolate, "Mycophenolate (notes, planned or given)"), win(d$tocilizumab, "Tocilizumab (notes, planned or given)"))
print(R$exposure_windows, n = Inf)
R$ivmp_sources <- d |> filter(gc_iv == "Yes") |> count(gc_prior, ivmp_followup) |> mutate(gc_prior = paste("history of systemic GC:", gc_prior), ivmp_followup = paste("IVMP in follow-up notes:", ivmp_followup)); print(R$ivmp_sources)

cat_h("D2 IVMP dose documentation")
R$ivmp_dose <- d |> filter(gc_iv == "Yes") |> summarise(n = n(), dose_recorded = sum(!is.na(gc_dose)), pulse_dose_valid = sum(!is.na(ivmp_pulse_dose)),
  med_pulse = med_iqr(ivmp_pulse_dose, 0), duration_recorded = sum(!is.na(gc_duration_months)), med_duration = med_iqr(gc_duration_months, 1)); print(R$ivmp_dose)
R$ivmp_pulse_dist <- d |> filter(!is.na(ivmp_pulse_dose)) |> count(ivmp_pulse_dose); print(R$ivmp_pulse_dist)
R$ivmp_dose_other <- d |> filter(gc_iv == "Yes", !is.na(gc_dose), is.na(ivmp_pulse_dose)) |> select(id, centre, gc_dose, gc_duration_months); print(R$ivmp_dose_other)
R$cumulative_dose_text <- d |> filter(str_detect(det_fu, "g[mr]? de metilpred|dosis total|\\d[\\.,]\\d ?gm|metilprednisolona +\\d[\\.,]\\d")) |> select(id, centre, details_go_treatment_since_last); print(as.data.frame(R$cumulative_dose_text))

cat_h("D3 glucocorticoids in patients inactive at baseline")
R$gc_inactive <- d |> filter(activity == "Inactive", gc_any == "Yes") |>
  mutate(window = case_when(gc_planned == "Yes" & gc_prior == "Yes" ~ "Planned at referral, with prior course in history",
                            gc_planned == "Yes" ~ "Planned at referral, no prior course", TRUE ~ "History only (treated before referral)")) |>
  count(window, gc_iv) |> arrange(window); print(R$gc_inactive)
R$ivmp_inactive_sev <- d |> filter(activity == "Inactive", gc_iv == "Yes") |> count(sev4, don, diplopia_sig); print(R$ivmp_inactive_sev)
R$ivmp_inactive_centre <- d |> filter(activity == "Inactive", gc_iv == "Yes") |> count(centre); print(R$ivmp_inactive_centre)

cat_h("D4 surgery sub-types")
R$surgery_types <- d |> filter(orbital_surgery == "Yes") |> summarise(n = n(), decompression = sum(surg_decompression), strabismus = sum(surg_strabismus),
  eyelid = sum(surg_eyelid), unspecified = sum(surg_other_unspec), history_yes = sum(surgery_go == "Yes"), plan_yes = sum(surgery_plan == "Yes", na.rm = TRUE)); print(R$surgery_types)
R$surgery_text <- d |> filter(orbital_surgery == "Yes") |> select(id, centre, surgery_go, surgery_plan, other_treatment_plan, details_go_treatment_since_last); print(as.data.frame(R$surgery_text))

cat_h("D5 inactivation by IVMP documentation status (exploratory)")
R$inact_by_ivmp_doc <- act |> mutate(grp = case_when(gc_iv == "No" ~ "No IVMP recorded", ivmp_followup %in% "Yes" | active_go_treatment_since_last %in% "Yes" ~ "IVMP; active treatment documented after referral", active_go_treatment_since_last %in% "No" ~ "IVMP; no active treatment documented after referral", TRUE ~ "IVMP; not recorded")) |>
  group_by(grp) |> summarise(n = n(), inactivated = sum(inactivated), pct = n_pct(inactivated, n), modsev = sum(modsev), .groups = "drop"); print(R$inact_by_ivmp_doc)
R$inact_by_ivmp_sev <- act |> group_by(gc_iv, modsev) |> summarise(n = n(), inactivated = sum(inactivated), pct = n_pct(inactivated, n), .groups = "drop"); print(R$inact_by_ivmp_sev)
R$inact_ivmp_fisher <- fisher.test(table(act$gc_iv, act$inactivated))$p.value

# ---------------------------------------------------------------------------------------
# E. Between-centre heterogeneity and modelling (Rev 1 #7, #13; Rev 2 Q6; Rev 3 #11)
# ---------------------------------------------------------------------------------------
cat_h("E1 models for active disease at baseline")
dd <- drop_na(d, y_active, age10, male, current_smoker, rai) |> arrange(centre)
f_prim <- firth_fit(y_active ~ age10 + male + current_smoker + rai, dd)
f_cent <- firth_fit(y_active ~ age10 + male + current_smoker + rai + centre, dd)
# crude (univariable) Firth
crude <- map_dfr(c("age10", "male", "current_smoker", "rai"), function(v) tidy_firth(firth_fit(as.formula(paste("y_active ~", v)), dd), "Crude (univariable, same 108 patients)"))
# mixed-effects logistic regression (random intercept for centre)
mm <- glmer(y_active ~ age10 + male + current_smoker + rai + (1 | centre), data = dd, family = binomial, control = glmerControl(optimizer = "bobyqa"))
mm_ci <- tryCatch(confint(mm, parm = "beta_", method = "profile", quiet = TRUE), error = function(e) confint(mm, parm = "beta_", method = "Wald"))
R$mixed_model <- tibble(model = "Mixed-effects logistic (random intercept for centre)", term = names(fixef(mm)), OR = exp(fixef(mm)),
                        lo = exp(mm_ci[, 1]), hi = exp(mm_ci[, 2]), p = coef(summary(mm))[, "Pr(>|z|)"], n = nobs(mm), events = sum(dd$y_active)) |> filter(term != "(Intercept)")
R$mixed_model_var <- as.numeric(VarCorr(mm)$centre[1]); cat("Random-intercept variance:", R$mixed_model_var, "\n")
# GEE with exchangeable correlation, clustered by centre (5 clusters: interpret with caution)
gee <- geeglm(y_active ~ age10 + male + current_smoker + rai, id = centre, data = dd, family = binomial, corstr = "exchangeable")
gs <- coef(summary(gee))
R$gee_model <- tibble(model = "GEE, exchangeable correlation, clustered by centre", term = rownames(gs), OR = exp(gs$Estimate),
                      lo = exp(gs$Estimate - qnorm(0.975) * gs$Std.err), hi = exp(gs$Estimate + qnorm(0.975) * gs$Std.err), p = gs$`Pr(>|W|)`, n = nrow(dd), events = sum(dd$y_active)) |> filter(term != "(Intercept)")
# leave-one-centre-out
R$loco <- map_dfr(sort(unique(dd$centre)), function(cc) { s <- filter(dd, centre != cc)
  tidy_firth(firth_fit(y_active ~ age10 + male + current_smoker + rai, s), paste("Excluding", cc)) |> mutate(excluded = cc) })
R$excl_BD <- tidy_firth(firth_fit(y_active ~ age10 + male + current_smoker + rai, filter(dd, !centre %in% EXCLUDE_SET)), paste("Excluding", EXCLUDE_LABEL))
d_mi <- d |> mutate(smk_missing = as.numeric(is.na(current_smoker)), cs_f = coalesce(current_smoker, 0), rai_missing = as.numeric(is.na(rai)), rai_f = coalesce(rai, 0))
f_mi <- firth_fit(y_active ~ age10 + male + cs_f + smk_missing + rai_f + rai_missing, drop_na(d_mi, y_active, age10))
f_ever <- firth_fit(y_active ~ age10 + male + ever_smoker + rai, drop_na(d, y_active, age10, male, ever_smoker, rai))
f_sens <- firth_fit(y_active_sens ~ age10 + male + current_smoker + rai, drop_na(d, y_active_sens, age10, male, current_smoker, rai))
R$models_orig_sens <- bind_rows(tidy_firth(f_mi, "Missing-indicator model"), tidy_firth(f_sens, "Outcome incl. item-derived CAS when recorded CAS missing"), tidy_firth(f_ever, "Ever versus never smoking"))
R$models_all <- bind_rows(tidy_firth(f_prim, "Primary Firth (complete case)"), crude, tidy_firth(f_cent, "Firth adjusted for centre (fixed effects)") |> filter(!str_detect(term, "centre")),
                          R$mixed_model, R$gee_model, R$models_orig_sens, R$firth_cas_only, R$firth_items7, R$loco |> select(-excluded), R$excl_BD)
print(R$models_all |> mutate(across(c(OR, lo, hi), ~ round(.x, 2)), p = fmt_p(p)), n = Inf)
# centre fixed effects themselves
R$centre_effects <- tidy_firth(f_cent, "centre") |> filter(str_detect(term, "centre")); print(R$centre_effects |> mutate(across(c(OR, lo, hi), ~ round(.x, 2)), p = fmt_p(p)))

cat_h("E2 included (n=108) vs excluded (n=31) from the primary model")
d <- d |> mutate(in_model = !is.na(y_active) & !is.na(age10) & !is.na(male) & !is.na(current_smoker) & !is.na(rai),
                 model_grp = factor(if_else(in_model, "Included in model (n=108)", "Excluded (n=31)"), levels = c("Included in model (n=108)", "Excluded (n=31)")))
R$excluded_reasons <- d |> filter(!is.na(activity), !in_model) |> summarise(n = n(), missing_age = sum(is.na(age10)), missing_smoking = sum(is.na(current_smoker)), missing_rai = sum(is.na(rai)),
  missing_smoking_only = sum(is.na(current_smoker) & !is.na(rai) & !is.na(age10)), missing_rai_only = sum(!is.na(current_smoker) & is.na(rai) & !is.na(age10)), both = sum(is.na(current_smoker) & is.na(rai))); print(R$excluded_reasons)
R$excluded_centre <- d |> filter(!is.na(activity)) |> count(centre, in_model) |> pivot_wider(names_from = in_model, values_from = n, values_fill = 0); print(R$excluded_centre)
R$included_vs_excluded <- cmp_tbl(filter(d, !is.na(activity)), "model_grp"); print(R$included_vs_excluded, n = Inf)

cat_h("E3 missingness by centre and by subgroup")
miss_vars <- c(age_go_onset = "Age at TED onset", smoking = "Smoking status", radioiodine_treatment = "Prior radioiodine", cas = "Recorded CAS", activity = "Activity classification",
               exo_max = "Exophthalmometry", ulr_max = "Upper eyelid retraction", gorman = "Gorman score", cornea = "Corneal assessment", va_worse = "Visual acuity",
               don = "DON record", trab_status = "TRAb", date_of_followup = "Follow-up date", activity_fu = "Follow-up activity", cas_fu = "Follow-up CAS")
R$missing_by_centre <- map_dfr(names(miss_vars), function(v) {
  x <- d[[v]]; tibble(Variable = miss_vars[[v]], Overall = n_pct(sum(is.na(x)), nrow(d)),
                      !!!set_names(map(sort(unique(d$centre)), ~ n_pct(sum(is.na(x[d$centre == .x])), sum(d$centre == .x))), sort(unique(d$centre)))) })
print(R$missing_by_centre, n = Inf)
R$missing_by_activity <- map_dfr(c("exo_max", "ulr_max", "gorman", "cornea", "trab_status", "activity_fu"), function(v) { x <- d[[v]]
  tibble(Variable = miss_vars[[v]], Inactive = n_pct(sum(is.na(x) & d$activity == "Inactive", na.rm = TRUE), sum(d$activity == "Inactive", na.rm = TRUE)),
         Active = n_pct(sum(is.na(x) & d$activity == "Active", na.rm = TRUE), sum(d$activity == "Active", na.rm = TRUE)),
         p = fmt_p(fisher.test(table(is.na(x)[!is.na(d$activity)], d$activity[!is.na(d$activity)]))$p.value)) }); print(R$missing_by_activity)

cat_h("E4 expanded centre table")
cen <- function(f, lab, dig = 1) tibble(Characteristic = lab, Overall = f(d), !!!set_names(map(sort(unique(d$centre)), ~ f(filter(d, centre == .x))), sort(unique(d$centre))))
bin_c <- function(var, lvl = "Yes") function(s) { x <- s[[var]]; n_pct(sum(x == lvl, na.rm = TRUE), sum(!is.na(x))) }
R$centre_table <- bind_rows(
  tibble(Characteristic = "Patients, n", Overall = as.character(nrow(d)), !!!set_names(map(sort(unique(d$centre)), ~ as.character(sum(d$centre == .x))), sort(unique(d$centre)))),
  cen(function(s) med_iqr(s$age_go_onset, 0), "Age at TED onset, years"), cen(bin_c("sex", "Female"), "Female"),
  cen(bin_c("smoking", "Current"), "Current smoking"), cen(bin_c("radioiodine_treatment"), "Prior radioiodine"),
  cen(bin_c("activity", "Active"), "Active disease at baseline"), cen(bin_c("activity_source", "Clinician label"), "Activity from clinician label (CAS missing)"),
  cen(function(s) med_iqr(s$cas, 0), "Recorded CAS"),
  cen(bin_c("sev4", "Mild"), "Mild (fully assessed)"), cen(bin_c("sev4", "Mild features, incomplete assessment"), "Mild features, incomplete assessment"),
  cen(bin_c("sev4", "Moderate-to-severe"), "Moderate-to-severe"), cen(bin_c("sev4", "Sight-threatening (DON/corneal breakdown)"), "Sight-threatening"),
  cen(bin_c("don"), "DON at baseline"), cen(function(s) n_pct(sum(!is.na(s$exo_max)), nrow(s)), "Exophthalmometry recorded"),
  cen(function(s) n_pct(sum(!is.na(s$ulr_max)), nrow(s)), "Upper eyelid retraction recorded"), cen(function(s) n_pct(sum(!is.na(s$cornea)), nrow(s)), "Corneal assessment recorded"),
  cen(function(s) n_pct(sum(!is.na(s$trab_status)), nrow(s)), "TRAb recorded"),
  cen(bin_c("gc_prior"), "Systemic GC in history before referral"), cen(bin_c("ivmp_planned"), "IVMP planned at referral"), cen(bin_c("gc_iv"), "IVMP composite (recorded or planned)"),
  cen(bin_c("second_line"), "Second-line therapy"), cen(bin_c("orbital_surgery"), "Surgery (done or planned)"),
  cen(function(s) n_pct(sum(s$has_followup), nrow(s)), "Any follow-up information"), cen(function(s) n_pct(sum(s$paired_activity), nrow(s)), "Paired activity assessment"),
  cen(function(s) n_pct(sum(!is.na(s$cas) & !is.na(s$cas_fu)), nrow(s)), "Paired recorded CAS"),
  cen(function(s) med_iqr(s$months_go_to_fu, 1), "Months from symptom onset to last visit"),
  cen(function(s) n_pct(sum(s$inactivated == 1, na.rm = TRUE), sum(!is.na(s$inactivated))), "Inactivation among active with follow-up"))
print(R$centre_table, n = Inf)

# ---------------------------------------------------------------------------------------
# F. Statistical reporting details (Rev 2 Q7; Rev 3 #4)
# ---------------------------------------------------------------------------------------
cat_h("F paired CAS: Hodges-Lehmann and effect sizes")
pc <- d |> filter(!is.na(cas), !is.na(cas_fu))
wt <- wilcox.test(pc$cas_fu, pc$cas, paired = TRUE, exact = FALSE, correct = TRUE, conf.int = TRUE)
tt <- t.test(pc$cas_fu - pc$cas)
R$cas_paired <- tibble(n = nrow(pc), baseline = med_iqr(pc$cas, 0), last = med_iqr(pc$cas_fu, 0), mean_change = round(mean(pc$cas_fu - pc$cas), 2),
                       mean_change_ci = sprintf("%.2f to %.2f", tt$conf.int[1], tt$conf.int[2]), median_change = median(pc$cas_fu - pc$cas),
                       hodges_lehmann = round(unname(wt$estimate), 2), hl_ci = sprintf("%.2f to %.2f", wt$conf.int[1], wt$conf.int[2]), p = fmt_p(wt$p.value),
                       improved_ge2 = sum(pc$cas_fu - pc$cas <= -2), worsened_ge2 = sum(pc$cas_fu - pc$cas >= 2), unchanged = sum(pc$cas_fu == pc$cas))
print(as.data.frame(R$cas_paired))
R$cas_paired_cas_only_active <- pc |> filter(cas >= 3) |> summarise(n = n(), inactive_fu = sum(cas_fu < 3), pct = n_pct(inactive_fu, n), ci = wilson(inactive_fu, n), med_bl = med_iqr(cas, 0), med_fu = med_iqr(cas_fu, 0)); print(R$cas_paired_cas_only_active)

# paired binary outcomes with absolute change and 95% CI (Wald for paired difference)
mcn <- function(a, b, lvl = "Yes", label = a) { s <- d |> select(x = all_of(a), y = all_of(b)) |> filter(!is.na(x), !is.na(y)) |> mutate(x = x == lvl, y = y == lvl)
  n <- nrow(s); b01 <- sum(s$x & !s$y); b10 <- sum(!s$x & s$y); diff <- (sum(s$y) - sum(s$x)) / n; se <- sqrt((b01 + b10) / n^2 - (b01 - b10)^2 / n^3)
  tibble(Outcome = label, `Paired n` = n, Baseline = sprintf("%d (%.1f%%)", sum(s$x), 100 * mean(s$x)), `Last visit` = sprintf("%d (%.1f%%)", sum(s$y), 100 * mean(s$y)),
         `Absolute change, % (95% CI)` = sprintf("%.1f (%.1f to %.1f)", 100 * diff, 100 * (diff - 1.96 * se), 100 * (diff + 1.96 * se)),
         Improved = b01, Worsened = b10, p = fmt_p(if (b01 + b10 > 0) binom.test(b01, b01 + b10)$p.value else NA)) }
R$table4a <- bind_rows(mcn("activity", "activity_fu", "Active", "Active disease (CAS ≥3 or clinician label)"), mcn("activity_cas", "activity_cas_fu", "Active", "Active disease, recorded CAS ≥3 only"),
  mcn("painful_oppressive_feeling", "painful_oppressive_feeling_followup", label = "Spontaneous retrobulbar pain"), mcn("gaze_evoked_pain", "gaze_evoked_pain_followup", label = "Gaze-evoked pain"),
  mcn("grittiness", "grittiness_followup", label = "Grittiness"), mcn("photophobia", "photophobia_followup", label = "Photophobia"), mcn("excessive_watering", "excessive_watering_followup", label = "Excessive tearing"),
  mcn("double_vision", "double_vision_followup", label = "Diplopia (any)"), mcn("gorman_ge2", "gorman_ge2_fu", label = "Inconstant/constant diplopia (Gorman ≥2)"),
  mcn("color_vision", "color_vision_followup", label = "Colour vision abnormality"), mcn("rapd", "rapd_followup", label = "RAPD"), mcn("don", "optic_neuropathy_followup", label = "Dysthyroid optic neuropathy"))
print(R$table4a, n = Inf, width = Inf)
paired_cont <- function(a, b, label, dig = 2) { s <- d |> select(x = all_of(a), y = all_of(b)) |> drop_na()
  w <- wilcox.test(s$y, s$x, paired = TRUE, exact = FALSE, correct = TRUE, conf.int = TRUE)
  tibble(Variable = label, n = nrow(s), Baseline = med_iqr(s$x + 0, dig), `Last visit` = med_iqr(s$y + 0, dig), `Hodges–Lehmann change (95% CI)` = sprintf(paste0("%.", dig, "f (%.", dig, "f to %.", dig, "f)"), w$estimate, w$conf.int[1], w$conf.int[2]), p = fmt_p(w$p.value)) }
R$table4b <- bind_rows(paired_cont("cas", "cas_fu", "CAS (7-point)", 1), paired_cont("va_worse", "va_worse_fu", "BCVA, worse eye (decimal)"), paired_cont("logmar_worse", "logmar_worse_fu", "logMAR, worse eye"),
                       paired_cont("iop_max", "iop_max_fu", "IOP, higher eye, mmHg", 1), paired_cont("ft4", "ft4_fu", "Free T4, ng/dL"), paired_cont("tsh", "tsh_fu", "TSH, mU/L", 3))
print(R$table4b, width = Inf)

cat_h("F TRAb")
R$trab <- d |> filter(!is.na(activity)) |> group_by(activity) |> summarise(n = n(), tested = sum(!is.na(trab_status)), positive = sum(trab_status == "Positive", na.rm = TRUE)); print(R$trab)
R$trab_fisher_p <- fisher.test(table(d$trab_status, d$activity))$p.value; cat("Fisher p:", R$trab_fisher_p, "\n")

# Centre A (prospective entry): time-to-inactivation feasibility
cat_h("Centre A time-to-event feasibility")
R$centreA_tte <- act |> filter(centre %in% BASELINE_CENTRES) |> summarise(n = n(), with_interval = sum(!is.na(months_bl_to_fu)), inactivated = sum(inactivated)); print(R$centreA_tte)

# ---------------------------------------------------------------------------------------
# G. Longitudinal outcome excluding single centres (Rev 2 Q6) and abstract numbers
# ---------------------------------------------------------------------------------------
cat_h("G inactivation excluding centres")
R$inact_excl_centre <- map_dfr(c("None", sort(unique(d$centre)), EXCLUDE_LABEL), function(cc) {
  s <- if (cc == "None") act else if (cc == EXCLUDE_LABEL) filter(act, !centre %in% EXCLUDE_SET) else filter(act, centre != cc)
  p <- d |> filter(paired_activity); p <- if (cc == "None") p else if (cc == EXCLUDE_LABEL) filter(p, !centre %in% EXCLUDE_SET) else filter(p, centre != cc)
  tibble(Excluded = cc, `Paired n` = nrow(p), `Active at baseline` = sprintf("%d (%.1f%%)", sum(p$activity == "Active"), 100 * mean(p$activity == "Active")),
         `Active at last visit` = sprintf("%d (%.1f%%)", sum(p$activity_fu == "Active"), 100 * mean(p$activity_fu == "Active")),
         `Inactivation` = n_pct(sum(s$inactivated), nrow(s)), `95% CI` = wilson(sum(s$inactivated), nrow(s))) })
print(R$inact_excl_centre)
R$abstract_numbers <- list(n = nrow(d), women = sum(d$sex == "Female"), active = "65/139", sight_threat = "12/147", don_bl = "10/145", cb = "3/71",
  ivmp_planned = "59/137", ivmp_tx_after_referral = 23, ivmp_no_tx_after_referral = 21, ivmp_tx_not_recorded = 15, paired = 109, active_bl = 51, active_fu = 14, inact = "37/51", new_don = "3/107", don_fu = "8/116",
  interval_A = R$interval_centreA$median_iqr, duration = med_iqr(d$months_go_to_fu))
print(R$abstract_numbers)
# ---------------------------------------------------------------------------------------
# H. Second-round additions: protocol deviations, IVMP documentation, sex-specific Hertel
# ---------------------------------------------------------------------------------------
cat_h("H1 protocol deviations")
R$age_lt18 <- d |> filter(age_go_onset < 18 | (year(date_of_inclusion) - yob) < 18) |> select(id, centre, age_go_onset)
R$fu_before_2020 <- d |> filter(year(date_of_followup) < 2020) |> count(centre)
print(R$age_lt18); print(R$fu_before_2020)
act_fu4 <- with(d, case_when(!is.na(cas_fu) & cas_fu >= 4 ~ "Active", !is.na(cas_fu) ~ "Inactive",
                             clinically_active_go_followup == "Yes" ~ "Active", clinically_active_go_followup == "No" ~ "Inactive"))
pp <- d$paired_activity
R$fu_threshold4 <- tibble(paired = sum(pp), active_fu_cas3 = sum(d$activity_fu[pp] == "Active"), active_fu_cas4 = sum(act_fu4[pp] == "Active"),
  active_bl = sum(d$activity[pp] == "Active"), inact_cas3 = sum(d$inactivated == 1, na.rm = TRUE),
  inact_cas4 = sum(d$activity[pp] == "Active" & act_fu4[pp] == "Inactive"))
print(R$fu_threshold4)
aA <- act |> filter(!is.na(months_bl_to_fu))
R$centreA_ge6 <- tibble(n_ge6 = sum(aA$months_bl_to_fu >= 6), inact_ge6 = sum(aA$inactivated[aA$months_bl_to_fu >= 6]))
print(R$centreA_ge6)

cat_h("H2 IVMP documentation after referral")
d <- d |> mutate(ivmp_grp = case_when(gc_iv == "No" ~ "No IVMP recorded or planned",
  ivmp_followup %in% "Yes" | active_go_treatment_since_last %in% "Yes" ~ "IVMP planned; active TED treatment documented after referral",
  active_go_treatment_since_last %in% "No" ~ "IVMP planned; no active TED treatment documented after referral",
  TRUE ~ "IVMP planned; treatment after referral not recorded"))
R$ivmp_grp_all <- d |> count(ivmp_grp)
R$ivmp_grp_prior <- d |> filter(gc_iv == "Yes") |> count(ivmp_grp, gc_prior)
R$ivmp_grp_inact <- d |> filter(!is.na(inactivated)) |> group_by(ivmp_grp) |> summarise(n = n(), inactivated = sum(inactivated), pct = n_pct(inactivated, n), modsev = sum(modsev), .groups = "drop")
R$tx_since_last <- d |> count(active_go_treatment_since_last)
R$ivmp_fu_notes <- sum(d$ivmp_followup %in% "Yes")
R$ivmp_planned_txfield <- d |> filter(gc_iv == "Yes") |> count(active_go_treatment_since_last)
print(R$ivmp_grp_all); print(R$ivmp_grp_inact); print(R$ivmp_planned_txfield)

cat_h("H3 exophthalmometry: sex and population-specific thresholds")
sevx <- function(lab, thr_f, thr_m) { thr <- ifelse(d$sex == "Male", thr_m, thr_f)
  s <- case_when(d$don == "Yes" | d$corneal_breakdown == "Yes" ~ "ST",
                 d$lid_retr_ge2 == "Yes" | (!is.na(d$exo_max) & d$exo_max >= thr) | d$diplopia_sig == "Yes" ~ "MS", TRUE ~ "Mild")
  tibble(Threshold = lab, `Women, mm` = thr_f, `Men, mm` = thr_m, `Moderate-to-severe` = n_pct(sum(s == "MS"), nrow(d)),
         `Sight-threatening` = n_pct(sum(s == "ST"), nrow(d)), `Exophthalmometry above threshold` = n_pct(sum(!is.na(d$exo_max) & d$exo_max >= thr), sum(!is.na(d$exo_max)))) }
R$exo_sex_thresholds <- bind_rows(sevx("Operational (both sexes)", 23, 23),
  sevx("Mexican adults: ULN + 3 mm (approx.)", 22, 23), sevx("ATA/ETA White ULN + 3 mm", 22, 24),
  sevx("Hispanic (New York) ULN + 3 mm (approx.)", 24, 25))
R$exo_by_sex <- d |> filter(!is.na(exo_max)) |> group_by(sex) |> summarise(n = n(), median_iqr = med_iqr(exo_max, 1), ge23 = sum(exo_max >= 23), .groups = "drop")
R$ms_exo_only <- sum(d$sev4 == "Moderate-to-severe" & !(d$lid_retr_ge2 %in% "Yes") & !(d$diplopia_sig %in% "Yes"))
print(R$exo_sex_thresholds); print(R$exo_by_sex); cat("MS solely by exophthalmometry:", R$ms_exo_only, "\n")

saveRDS(R, "outputs/revision_results.rds")
saveRDS(R, "outputs/revision_results.rds"); saveRDS(d, "outputs/ted_analytic_rev.rds")
# flat export of all tables to xlsx for the response letter
tbls <- R[map_lgl(R, is.data.frame)]
writexl::write_xlsx(map(tbls, ~ mutate(as_tibble(.x), across(where(is.factor), as.character))), "outputs/revision_results.xlsx")
cat("\nSaved outputs/revision_results.rds and .xlsx (", length(tbls), "tables )\n")
