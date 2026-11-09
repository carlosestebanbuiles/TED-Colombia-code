# 04_figures.R : publication figures for the revised manuscript (PNG 300 dpi + TIFF 600 dpi)
suppressPackageStartupMessages({library(tidyverse); library(patchwork); library(scales)})
source("00_helpers.R"); d <- readRDS("outputs/ted_analytic_rev.rds"); R <- readRDS("outputs/revision_results.rds")
fig_dir <- "outputs/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

pal_activity <- c(Inactive = "#2166AC", Active = "#B2182B"); pal_time <- c(Baseline = "#B2182B", `Last visit` = "#2166AC")
pal_treat <- c(`IV methylprednisolone` = "#D55E00", `Oral GC only` = "#E69F00", `Second-line therapy` = "#009E73", Surgery = "#0072B2")
pal_model <- c(`Primary Firth (complete case)` = "#111111", `Firth adjusted for centre (fixed effects)` = "#8C8C8C", `Mixed-effects logistic (random intercept for centre)` = "#D55E00")
theme_ted <- theme_classic(base_size = 11, base_family = "sans") +
  theme(axis.text = element_text(colour = "grey25"), axis.title = element_text(colour = "grey10"),
        axis.line = element_line(colour = "grey45", linewidth = 0.4), axis.ticks = element_line(colour = "grey45", linewidth = 0.4),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3), legend.position = "bottom", legend.title = element_blank(),
        legend.key.size = unit(4, "mm"), plot.title = element_text(face = "bold", size = 11, hjust = 0), plot.title.position = "plot")

# ---- Figure 1. Study flow (analysis sets, with source of activity classification) -------------
n_lab <- sum(d$activity_source == "Clinician label", na.rm = TRUE); n_cas <- sum(d$activity_source == "Recorded CAS", na.rm = TRUE)
flow <- tibble(y = c(6, 5, 4, 3, 2, 1), lab = c(
  sprintf("Patients in the multicentre TED registry\n(5 centres, 6 sites): n = %d\nNumber of records screened not documented", nrow(d)),
  sprintf("Baseline activity classifiable: n = %d\nrecorded CAS n = %d (active %d); clinician label only n = %d (active %d)\nnot classifiable n = %d",
          sum(!is.na(d$activity)), n_cas, sum(d$activity_source == "Recorded CAS" & d$activity == "Active", na.rm = TRUE), n_lab,
          sum(d$activity_source == "Clinician label" & d$activity == "Active", na.rm = TRUE), sum(is.na(d$activity))),
  sprintf("Any follow-up information: n = %d\n(no follow-up information: n = %d)", sum(d$has_followup), sum(!d$has_followup)),
  sprintf("Paired activity status (baseline and last visit): n = %d\nboth from recorded CAS n = %d; clinician label at one or both visits n = %d\nPaired recorded CAS: n = %d",
          sum(d$paired_activity), sum(d$paired_activity & d$activity_source == "Recorded CAS" & d$activity_source_fu == "Recorded CAS"),
          sum(d$paired_activity & (d$activity_source == "Clinician label" | d$activity_source_fu == "Clinician label")), sum(!is.na(d$cas) & !is.na(d$cas_fu))),
  sprintf("Active at baseline with follow-up activity: n = %d\n(primary follow-up outcome: inactivation, %d/%d)", sum(!is.na(d$inactivated)), sum(d$inactivated == 1, na.rm = TRUE), sum(!is.na(d$inactivated))),
  sprintf("Baseline-to-last-visit interval computable: n = %d (Centre A only)\nOther centres: last visit dated, baseline encounter date not recorded", sum(!is.na(d$months_bl_to_fu)))))
fig1 <- ggplot(flow, aes(x = 0, y = y)) +
  geom_segment(data = flow[-6, ], aes(x = 0, xend = 0, y = y - 0.3, yend = y - 0.7), arrow = arrow(length = unit(2, "mm")), colour = "grey40") +
  geom_label(aes(label = lab), size = 3, label.padding = unit(3, "mm"), fill = "white") + xlim(-1, 1) + ylim(0.5, 6.5) + theme_void()

# ---- Figure 2. Treatment by activity/severity (A) and exposure windows for glucocorticoids (B) ----
fig2a <- d |> filter(!is.na(activity)) |>
  mutate(grp = if_else(modsev, "Moderate-to-severe/\nsight-threatening", "Mild or\nincompletely assessed"), grp = paste(activity, grp, sep = "\n")) |>
  group_by(grp) |> summarise(n = n(), `IV methylprednisolone` = mean(gc_iv == "Yes"), `Oral GC only` = mean(gc_oral_only == "Yes"),
                             `Second-line therapy` = mean(second_line == "Yes"), Surgery = mean(orbital_surgery == "Yes"), .groups = "drop") |>
  pivot_longer(-c(grp, n)) |> mutate(grp = paste0(grp, "\n(n = ", n, ")"), name = factor(name, levels = names(pal_treat))) |>
  ggplot(aes(grp, value, fill = name)) + geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = round(100 * value), group = name), position = position_dodge(0.8), vjust = -0.4, size = 2.8, colour = "grey20") +
  scale_y_continuous(labels = percent, limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) + scale_fill_manual(values = pal_treat) +
  labs(x = NULL, y = "Patients (%)", title = "A. Recorded or planned treatment by baseline activity and severity") +
  theme_ted + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))
win <- tibble(Window = factor(c("History before/at referral\n(any systemic GC)", "Planned at referral visit\n(IV methylprednisolone)", "Planned at referral visit\n(oral GC)",
                                "Documented in follow-up notes\n(IV methylprednisolone)", "Documented in follow-up notes\n(any systemic GC)"),
                              levels = rev(c("History before/at referral\n(any systemic GC)", "Planned at referral visit\n(IV methylprednisolone)", "Planned at referral visit\n(oral GC)",
                                             "Documented in follow-up notes\n(IV methylprednisolone)", "Documented in follow-up notes\n(any systemic GC)"))),
              k = R$exposure_windows$n_yes[c(1, 3, 4, 6, 5)], N = R$exposure_windows$n_known[c(1, 3, 4, 6, 5)]) |>
  mutate(p = k / N, lab = sprintf("%d/%d (%.0f%%)", k, N, 100 * p))   # rows of Table 3A computed in 03_revision_analyses.R
fig2b <- ggplot(win, aes(p, Window)) + geom_col(fill = "#8C8C8C", width = 0.65) + geom_text(aes(label = lab), hjust = -0.1, size = 3, colour = "grey20") +
  scale_x_continuous(labels = percent, limits = c(0, 0.75), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Patients with the window documented (%)", y = NULL, title = "B. Systemic glucocorticoids by documentation window") + theme_ted
fig2 <- fig2a / fig2b + plot_layout(heights = c(1.1, 0.9))

# ---- Figure 3. Activity and symptoms, baseline vs last visit ------------------------------------
fig3a <- d |> filter(paired_activity) |> select(id, Baseline = activity, `Last visit` = activity_fu) |>
  pivot_longer(-id, names_to = "time", values_to = "status") |>
  ggplot(aes(time, fill = status)) + geom_bar(position = "fill", width = 0.6, colour = "white") +
  scale_y_continuous(labels = percent) + scale_fill_manual(values = pal_activity) +
  labs(x = NULL, y = "Patients (paired, %)", title = sprintf("A. Activity status (n = %d pairs)", sum(d$paired_activity))) + theme_ted
sym <- R$table4a |> filter(!str_detect(Outcome, "Active disease|Gorman|RAPD|Colour|Dysthyroid")) |>
  mutate(b = as.numeric(str_extract(Baseline, "(?<=\\()[0-9.]+")), f = as.numeric(str_extract(`Last visit`, "(?<=\\()[0-9.]+")), Outcome = fct_reorder(Outcome, b))
fig3b <- ggplot(sym) + geom_segment(aes(x = f, xend = b, y = Outcome, yend = Outcome), colour = "grey60") +
  geom_point(aes(b, Outcome, colour = "Baseline"), size = 3) + geom_point(aes(f, Outcome, colour = "Last visit"), size = 3) +
  scale_colour_manual(values = pal_time) + labs(x = "Prevalence among paired patients (%)", y = NULL, title = "B. Symptoms and signs") + theme_ted
fig3 <- fig3a + fig3b + plot_layout(widths = c(1, 1.6))

# ---- Figure 4. Forest plot: primary, centre fixed effects, mixed model ------------------------
fp <- R$models_all |> filter(model %in% names(pal_model)) |>
  mutate(term = recode(term, age10 = "Age (per 10 years)", male = "Male sex", current_smoker = "Current smoking", rai = "Prior radioiodine"),
         model = factor(model, levels = names(pal_model)))
fig4 <- ggplot(fp, aes(OR, fct_rev(factor(term, levels = c("Age (per 10 years)", "Male sex", "Current smoking", "Prior radioiodine"))), colour = model)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") + geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6)) +
  scale_x_log10() + scale_colour_manual(values = pal_model) + guides(colour = guide_legend(nrow = 3)) +
  labs(x = "Odds ratio for active disease at baseline (log scale)", y = NULL) + theme_ted

# ---- Supplementary Figure S1. Leave-one-centre-out and other sensitivity models (age) --------
sfp <- R$models_all |> filter(term == "age10", !str_detect(model, "GEE")) |>
  mutate(model = factor(model, levels = rev(unique(model))))
figS1 <- ggplot(sfp, aes(OR, model)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = "#111111") + geom_text(aes(label = sprintf("n = %d", n)), x = log10(4.2), size = 2.8, colour = "grey30") +
  scale_x_log10(limits = c(0.6, 5)) + labs(x = "Odds ratio for active disease per 10-year increase in age at TED onset (log scale)", y = NULL) + theme_ted

figs <- list(list(fig1, "Figure1_flow", 7.5, 8), list(fig2, "Figure2_management", 9, 8.5), list(fig3, "Figure3_clinical_course", 10, 5),
             list(fig4, "Figure4_forest", 7.5, 4.5), list(figS1, "FigureS1_sensitivity_age", 9, 5))
for (f in figs) { base <- file.path(fig_dir, f[[2]])
  ggsave(paste0(base, ".png"), f[[1]], width = f[[3]], height = f[[4]], dpi = 300, bg = "white")
  ggsave(paste0(base, ".tiff"), f[[1]], width = f[[3]], height = f[[4]], dpi = 600, device = "tiff", compression = "lzw", bg = "white") }
print(tibble(file = list.files(fig_dir)) |> mutate(kb = round(file.size(file.path(fig_dir, file)) / 1024)))
