# =====================================================================================
# TED Colombia — revision analysis (Therapeutic Advances in Ophthalmology, OED-26-0109)
# 01_derive.R : import BD.xlsx and derive the analytic data set.
#
# This file reproduces, line by line, the derivation rules of the original
# TED_Colombia_analysis.Rmd (section 2) so that every number in the revised
# manuscript can be traced to BD.xlsx. New variables added for the revision are
# marked with "# [REV]".
#
# Dependencies: tidyverse, readxl (janitor is NOT required: a local clean_names() is used)
# =====================================================================================
suppressPackageStartupMessages({library(tidyverse); library(readxl)})

bd_path <- Sys.getenv("TED_BD", unset = "BD.xlsx")

source("00_helpers.R")

# ---- import ----------------------------------------------------------------------------
raw <- read_excel(bd_path, skip = 1, col_types = "text") |> clean_names_local()
# janitor splits "ThYroidectomY" -> "th_yroidectom_y"; local version yields same tokens:
names(raw)[names(raw) %in% c("th_yroidectom_y")] <- "thyroidectomy"
names(raw)[names(raw) %in% c("th_yroid_dermopath_y")] <- "thyroid_dermopathy"


# Centre pseudonyms. The mapping from institution names to "Centre A"–"Centre E" is kept in a
# private file (config/centre_map.csv, columns: centre_code,centre) that is not part of the public
# deposit. If the file is absent (e.g. synthetic data, which already use pseudonyms), codes are used as is.
centre_map_file <- Sys.getenv("TED_CENTRE_MAP", unset = "config/centre_map.csv")
if (file.exists(centre_map_file)) {
  cm <- read.csv(centre_map_file, stringsAsFactors = FALSE, encoding = "UTF-8")
  centre_map <- setNames(cm$centre, cm$centre_code)
} else {
  centre_map <- NULL
}

yn_vars <- c("relapse_after_treatment","atd_commenced","atd_current","atd_stopped","radioiodine_treatment",
  "thyroidectomy","visible_goiter","thyroid_dermopathy","diabetes","addisons","pernicious_anaemia","vitiligo",
  "rheumatoid_arthritis","other_autoimmune","fh_thyroid_disease","fh_diabetes","fh_other_autoimmune",
  "metimazol_initial","propiltiuracilo_initial","lt4_initial","selenium_initial","topical_eye_treatment",
  "systemic_steroids","orbital_irradiation","surgery_go","painful_oppressive_feeling","gaze_evoked_pain",
  "excessive_watering","photophobia","grittiness","double_vision","rapd","color_vision",
  "active_eyelid_swelling","eyelid_erythema","conjunctival_redness","chemosis","caruncle_swelling",
  "plical_swelling","redness_lat_rect_insertion","superior_limbic_keratoconjunctivitis","lagophthalmos",
  "lateral_flare","proptosis_od","proptosis_os","motility_abnormal_head_posture","esotropia","exotropia",
  "hypotropia","hypertropia","binocular_single_vision_without_prism","optic_neuropathy_evidence",
  "ocular_comorbidities","orbitopathy_evidence","clinically_active_go","immunosuppression_plan",
  "irradiation_plan","surgery_plan","steroids_planned","other_immunosuppression_planned",
  "thyroid_treatment_changes","new_radioiodine","new_thyroidectomy","painful_oppressive_feeling_followup",
  "gaze_evoked_pain_followup","excessive_watering_followup","photophobia_followup","grittiness_followup",
  "double_vision_followup","rapd_followup","color_vision_followup","optic_neuropathy_followup",
  "ocular_comorbidities_followup","clinically_active_go_followup","active_go_treatment_since_last")
stopifnot(all(yn_vars %in% names(raw)))
date_vars <- c("date_of_inclusion","onset_thyroid_symptoms","date_of_diagnosis","eye_symptom_onset",
               "radioiodine_treatment_date_1","thyroidectomy_date","date_of_followup")

ft4_conv <- function(x) { v <- num(x); ifelse(!is.na(v) & v > 10, v / 12.87, v) }
trab_status <- function(x) { s <- str_to_lower(txt(x)); v <- num(x)
  case_when(str_detect(s, "posit|poait") ~ "Positive", str_detect(s, "negat") ~ "Negative",
            !is.na(v) & v > 1.75 ~ "Positive", !is.na(v) ~ "Negative", TRUE ~ NA_character_) }
parse_exo <- function(x) { s <- txt(x); n <- suppressWarnings(as.numeric(s)); ser <- !is.na(n) & n > 30000
  s[ser] <- format(as.Date(n[ser], origin = "1899-12-30"), "%m/%d")
  s <- str_replace_all(s, ",", ".") |> str_replace_all(" y ", "/"); parts <- str_split_fixed(s, "/", 2)
  od <- num(parts[, 1], 10, 35); os <- ifelse(parts[, 2] == "" & !is.na(s), od, num(parts[, 2], 10, 35))
  tibble(exo_od = od, exo_os = os) }
motility_code <- function(x) { s <- str_to_lower(txt(x))
  case_when(s %in% c("no","o") ~ 0, s %in% c("si","sí","yes") ~ 1, !is.na(num(s)) ~ as.numeric(num(s) != 0), TRUE ~ NA_real_) }

d <- raw |>
  mutate(
    id = sprintf("P%03d", row_number()),
    centre = if (is.null(centre_map)) txt(centre_code) else unname(centre_map[txt(centre_code)]),
    across(all_of(yn_vars), yn),
    across(all_of(date_vars), parse_date_any),
    date_of_inclusion = if_else(year(date_of_inclusion) < 2000, as.Date(NA), date_of_inclusion),
    yob = num(year_of_birth, 1900, 2015), sex = txt(sex), male = as.numeric(sex == "Male"),
    age_go_onset = year(eye_symptom_onset) - yob, age_inclusion = year(date_of_inclusion) - yob,
    weight = num(body_weight, 30, 200), height = num(height, 120, 210), bmi = weight / (height / 100)^2,
    race = case_when(str_detect(str_to_lower(race), "mestiz") ~ "Mestizo (Hispanic)", !is.na(txt(race)) ~ "Other", TRUE ~ NA_character_),
    ft4 = ft4_conv(f_t4), ft4_fu = ft4_conv(f_t4_followup), tsh = num(tsh), tsh_fu = num(tsh_followup),
    trab = num(tr_ab), trab_status = trab_status(tr_ab), trab_fu = num(tr_ab_followup),
    tpo_status = case_when(is.na(txt(tpo_ab)) ~ NA_character_, str_detect(str_to_lower(tpo_ab), "negat") ~ "Negative", TRUE ~ "Positive"),
    smoking = recode(txt(smoker), "Never smoked" = "Never", "Ex-smoker" = "Former", "Current smoker" = "Current"),
    smoking = factor(smoking, levels = c("Never","Former","Current")),
    current_smoker = as.numeric(smoking == "Current"),
    ever_smoker = as.numeric(smoking %in% c("Former","Current")) + ifelse(is.na(smoking), NA, 0),
    months_thyroid_symptoms_to_dx = months_between(onset_thyroid_symptoms, date_of_diagnosis),
    months_thyroid_symptoms_to_dx = if_else(months_thyroid_symptoms_to_dx < 0, NA_real_, months_thyroid_symptoms_to_dx),
    months_dx_to_go = months_between(date_of_diagnosis, eye_symptom_onset),
    go_timing = case_when(months_dx_to_go < -6 ~ "GO before thyroid diagnosis (>6 mo)", months_dx_to_go <= 6 ~ "Concurrent (±6 mo)",
                          months_dx_to_go > 6 ~ "GO after thyroid diagnosis (>6 mo)"),
    go_timing = factor(go_timing, levels = c("GO before thyroid diagnosis (>6 mo)","Concurrent (±6 mo)","GO after thyroid diagnosis (>6 mo)")),
    go_after_rai = case_when(radioiodine_treatment != "Yes" | is.na(radioiodine_treatment) ~ NA_character_,
                             is.na(radioiodine_treatment_date_1) | is.na(eye_symptom_onset) ~ NA_character_,
                             eye_symptom_onset > radioiodine_treatment_date_1 ~ "Yes", TRUE ~ "No"),
    months_go_to_fu = months_between(eye_symptom_onset, date_of_followup),
    months_go_to_fu = if_else(months_go_to_fu < 0, NA_real_, months_go_to_fu),
    cas_raw = num(cas_score, 0, 10), cas_fu_raw = num(cas_score_followup, 0, 10),         # [REV] keep raw for truncation count
    cas = pmin(cas_raw, 7), cas_fu = pmin(cas_fu_raw, 7),
    car_pli = case_when(caruncle_swelling == "Yes" | plical_swelling == "Yes" ~ 1,
                        is.na(caruncle_swelling) & is.na(plical_swelling) ~ NA_real_, TRUE ~ 0))

cas_items_mat <- d |> transmute(across(c(painful_oppressive_feeling, gaze_evoked_pain, eyelid_erythema, active_eyelid_swelling,
                                         conjunctival_redness, chemosis), ~ as.numeric(.x == "Yes")), car_pli)
d$cas_items_n     <- rowSums(!is.na(cas_items_mat))                                                   # [REV]
d$cas_items       <- ifelse(d$cas_items_n >= 6, rowSums(cas_items_mat, na.rm = TRUE), NA_real_)
d$cas_items7      <- ifelse(d$cas_items_n == 7, rowSums(cas_items_mat), NA_real_)                      # [REV] all 7 items documented

act_from <- function(cas, clin) case_when(!is.na(cas) & cas >= 3 ~ "Active", !is.na(cas) ~ "Inactive",
                                          clin == "Yes" ~ "Active", clin == "No" ~ "Inactive", TRUE ~ NA_character_)
d <- d |>
  mutate(
    activity = act_from(cas, clinically_active_go), activity_fu = act_from(cas_fu, clinically_active_go_followup),
    activity_source    = case_when(!is.na(cas) ~ "Recorded CAS", !is.na(clinically_active_go) ~ "Clinician label", TRUE ~ NA_character_),       # [REV]
    activity_source_fu = case_when(!is.na(cas_fu) ~ "Recorded CAS", !is.na(clinically_active_go_followup) ~ "Clinician label", TRUE ~ NA_character_), # [REV]
    activity_cas    = case_when(cas >= 3 ~ "Active", !is.na(cas) ~ "Inactive"),            # [REV] recorded CAS only
    activity_cas_fu = case_when(cas_fu >= 3 ~ "Active", !is.na(cas_fu) ~ "Inactive"),      # [REV]
    activity_sens = coalesce(activity, if_else(cas_items >= 3, "Active", "Inactive")),
    across(c(activity, activity_fu, activity_sens, activity_cas, activity_cas_fu), ~ factor(.x, levels = c("Inactive","Active")))) |>
  bind_cols(parse_exo(raw$exophthalmometer)) |>
  mutate(
    exo_max = pmax(exo_od, exo_os, na.rm = TRUE),
    exo_ge23 = if_else(is.na(exo_max), NA_character_, if_else(exo_max >= 23, "Yes", "No")),
    exo_asym = if_else(is.na(exo_od) | is.na(exo_os), NA_character_, if_else(abs(exo_od - exo_os) >= 2, "Yes", "No")),
    proptosis_clin = any_yes(proptosis_od, proptosis_os),
    across(c(upper_lid_retraction_od, upper_lid_retraction_os, lower_lid_retraction_od, lower_lid_retraction_os), ~ num(.x, 0, 15)),
    ulr_max = pmax(upper_lid_retraction_od, upper_lid_retraction_os, na.rm = TRUE),
    lid_retr_ge2 = case_when(ulr_max >= 2 ~ "Yes", !is.na(ulr_max) | !is.na(lower_lid_retraction_od) ~ "No"),
    across(c(adduction_od, adduction_os, abduction_od, abduction_os, elevation_90_od, elevation_90_os,
             depression_270_od, depression_270_os), motility_code, .names = "mot_{.col}")) |>
  mutate(
    motility_restriction = { m <- as.matrix(pick(starts_with("mot_")))
      ifelse(rowSums(!is.na(m)) == 0, NA_character_, ifelse(rowSums(m == 1, na.rm = TRUE) > 0, "Yes", "No")) },
    strabismus = any_yes(esotropia, exotropia, hypotropia, hypertropia),
    gorman = num(gorman_score, 0, 3),
    gorman_fu = num(if_else(str_to_lower(txt(gorman_score_followup)) %in% "o", "0", gorman_score_followup), 0, 3),
    diplopia_sig = case_when(!is.na(gorman) ~ if_else(gorman >= 2, "Yes", "No"), double_vision == "Yes" ~ "Yes", double_vision == "No" ~ "No"),
    across(c(best_visual_acuity_od, best_visual_acuity_os, best_visual_acuity_followup_od, best_visual_acuity_followup_os), ~ num(.x, 0, 2)),
    va_worse = pmin(best_visual_acuity_od, best_visual_acuity_os, na.rm = TRUE),
    va_worse_fu = pmin(best_visual_acuity_followup_od, best_visual_acuity_followup_os, na.rm = TRUE),
    logmar_worse = -log10(pmax(va_worse, 0.001)), logmar_worse_fu = -log10(pmax(va_worse_fu, 0.001)),
    across(c(intraocular_pressure_od, intraocular_pressure_os, intraocular_pressure_followup_od, intraocular_pressure_followup_os), ~ num(.x, 5, 60)),
    iop_max = pmax(intraocular_pressure_od, intraocular_pressure_os, na.rm = TRUE),
    iop_max_fu = pmax(intraocular_pressure_followup_od, intraocular_pressure_followup_os, na.rm = TRUE),
    iop_gt21 = if_else(is.na(iop_max), NA_character_, if_else(iop_max > 21, "Yes", "No")),
    cornea_txt = str_to_lower(txt(cornea_findings)),
    cornea_txt2 = str_remove_all(cornea_txt, "sin queratitis|sin tinci[oó]n(es)?"),
    cornea = case_when(is.na(cornea_txt) ~ NA_character_,
                       str_detect(cornea_txt, "úlcera|ulcera|defecto epitelial") ~ "Corneal breakdown (ulcer/epithelial defect)",
                       str_detect(cornea_txt2, "queratit|punteado|punteada|queratopat|tinci[oó]n") ~ "Keratopathy (punctate staining/keratitis)",
                       TRUE ~ "No GO-related corneal findings"),
    corneal_breakdown = if_else(is.na(cornea), NA_character_, if_else(str_starts(cornea, "Corneal breakdown"), "Yes", "No")),
    don = optic_neuropathy_evidence,
    severity = case_when(don == "Yes" | corneal_breakdown == "Yes" ~ "Sight-threatening (DON/corneal breakdown)",
                         lid_retr_ge2 == "Yes" | exo_ge23 == "Yes" | diplopia_sig == "Yes" ~ "Moderate-to-severe", TRUE ~ "Mild"),
    sev4 = if_else(severity == "Mild" & is.na(lid_retr_ge2) & is.na(exo_ge23), "Mild features, incomplete assessment", severity),
    sev4 = factor(sev4, levels = c("Mild","Mild features, incomplete assessment","Moderate-to-severe","Sight-threatening (DON/corneal breakdown)")),
    modsev = sev4 %in% c("Moderate-to-severe","Sight-threatening (DON/corneal breakdown)"))

# ---- treatments (original composite + [REV] exposure windows) ----------------------------
d <- d |>
  mutate(
    det_fu   = str_to_lower(coalesce(txt(details_go_treatment_since_last), "")),   # free text: treatment DURING follow-up
    det_plan = str_to_lower(paste(coalesce(txt(other_treatment_plan), ""), coalesce(txt(raw$other_immunosuppression_planned), ""), sep = " | ")),
    det = paste(det_fu, det_plan, sep = " | "),
    gc_route = txt(steroids_if_yes),
    # original composite "administered or planned"
    gc_any = if_else(systemic_steroids %in% "Yes" | steroids_planned %in% "Yes" | str_detect(det, "metilpred|esteroide sist|bolos"), "Yes", "No"),
    gc_iv = if_else(gc_any == "Yes" & (gc_route %in% "Systemic" | str_detect(det, "metilpred|bolos|iv ")), "Yes", "No"),
    gc_oral_only = if_else(gc_any == "Yes" & gc_iv == "No", "Yes", "No"),
    gc_dose = num(steroids_initial_dose),
    ivmp_pulse_dose = if_else(gc_iv == "Yes" & gc_dose %in% c(250, 500, 750, 1000), gc_dose, NA_real_),
    gc_duration_months = num(steroids_planned_duration, 0.1, 24),
    mycophenolate = if_else(str_detect(det, "micofenol"), "Yes", "No"),
    tocilizumab = if_else(str_detect(det, "tocil") & !str_detect(det, "no fue necesario"), "Yes", "No"),
    radiotherapy = if_else(orbital_irradiation %in% "Yes" | irradiation_plan %in% "Yes" | str_detect(det, "radioterap"), "Yes", "No"),
    second_line = if_else(mycophenolate == "Yes" | tocilizumab == "Yes" | radiotherapy == "Yes", "Yes", "No"),
    selenium = if_else(selenium_initial %in% "Yes" | str_detect(det, "selen"), "Yes", "No"),
    orbital_surgery = if_else(surgery_go %in% "Yes" | surgery_plan %in% "Yes" | str_detect(det, "descompres|tarsorraf|estrabismo|blefarotom|oculopl|evisc"), "Yes", "No"),
    decompression = if_else(str_detect(det, "descompres"), "Yes", "No"),
    # [REV] exposure windows for systemic glucocorticoids / IVMP
    #   (1) before/at referral: the "Historia de orbitopatía" field Systemic_Steroids (history of systemic steroids)
    #   (2) planned at the baseline visit: Steroids_Planned (+ route in Steroids_If_Yes)
    #   (3) administered during follow-up: Active_GO_Treatment_Since_Last + free-text details
    gc_prior   = systemic_steroids,
    gc_planned = steroids_planned,
    ivmp_planned = if_else(steroids_planned %in% "Yes" & gc_route %in% "Systemic", "Yes", if_else(is.na(steroids_planned), NA_character_, "No")),
    gc_fu_text   = str_detect(det_fu, "metilpred|esteroide|bolos|corticoide|prednis|dexametasona|deflazacort"),
    ivmp_fu_text = str_detect(det_fu, "metilpred|bolos|iv "),
    gc_followup  = case_when(active_go_treatment_since_last %in% "Yes" & gc_fu_text ~ "Yes",
                             !is.na(active_go_treatment_since_last) ~ "No", gc_fu_text ~ "Yes", TRUE ~ NA_character_),
    ivmp_followup = case_when(active_go_treatment_since_last %in% "Yes" & ivmp_fu_text ~ "Yes",
                              !is.na(active_go_treatment_since_last) ~ "No", ivmp_fu_text ~ "Yes", TRUE ~ NA_character_),
    # [REV] surgery sub-types from free text (plan + follow-up)
    surg_decompression = str_detect(det, "descompres"),
    surg_strabismus    = str_detect(det, "estrabismo|muscul|recto"),
    surg_eyelid        = str_detect(det, "tarsorraf|blefaro|párpado|parpado|retracci|mullerect|oculopl"),
    surg_other_unspec  = orbital_surgery == "Yes" & !(surg_decompression | surg_strabismus | surg_eyelid))

# ---- follow-up outcomes -------------------------------------------------------------------
worst_of <- function(a, b) case_when(is.na(a) & is.na(b) ~ NA_character_, a %in% "Worse" | b %in% "Worse" ~ "Worse",
                                     a %in% "Better" | b %in% "Better" ~ "Better", TRUE ~ "Same")
d <- d |>
  mutate(
    chg_activity = txt(changes_in_signs_activity), chg_severity = txt(changes_in_signs_severity),
    chg_proptosis = worst_of(txt(proptosis_followup_od), txt(proptosis_followup_os)), chg_motility = txt(motility_findings_followup),
    has_followup = !is.na(date_of_followup) | !is.na(cas_fu) | !is.na(clinically_active_go_followup) | !is.na(chg_activity),
    paired_activity = !is.na(activity) & !is.na(activity_fu),                                          # [REV]
    delta_cas = cas_fu - cas,
    inactivated = case_when(activity == "Active" & !is.na(activity_fu) ~ as.numeric(activity_fu == "Inactive")),
    inactivated_cas = case_when(activity_cas == "Active" & !is.na(activity_cas_fu) ~ as.numeric(activity_cas_fu == "Inactive")),  # [REV]
    new_don = case_when(don == "No" & !is.na(optic_neuropathy_followup) ~ optic_neuropathy_followup),
    unfavourable = case_when(is.na(chg_activity) & is.na(chg_severity) & is.na(activity_fu) ~ NA_character_,
      chg_activity %in% "Worse" | chg_severity %in% "Worse" | new_don %in% "Yes" | (activity %in% "Inactive" & activity_fu %in% "Active") ~ "Yes", TRUE ~ "No"),
    age10 = age_go_onset / 10,
    rai = case_when(radioiodine_treatment == "Yes" ~ 1, radioiodine_treatment == "No" ~ 0),
    y_active = as.numeric(activity == "Active"), y_active_sens = as.numeric(activity_sens == "Active"),
    y_active_cas = as.numeric(activity_cas == "Active"),                                              # [REV]
    ivgc = as.numeric(gc_iv == "Yes"),
    gorman_ge2 = if_else(is.na(gorman), NA_character_, if_else(gorman >= 2, "Yes", "No")),
    gorman_ge2_fu = if_else(is.na(gorman_fu), NA_character_, if_else(gorman_fu >= 2, "Yes", "No")),
    # [REV] true baseline-to-last-visit interval: only in Centre A does Date_of_Inclusion correspond to the baseline encounter
    months_bl_to_fu = if_else(centre %in% BASELINE_CENTRES, months_between(date_of_inclusion, date_of_followup), NA_real_),
    months_bl_to_fu = if_else(months_bl_to_fu < 0, NA_real_, months_bl_to_fu))

dir.create("outputs", showWarnings = FALSE)
saveRDS(d, "outputs/ted_analytic.rds")
cat("Derived data set:", nrow(d), "patients;", ncol(d), "variables\n")
