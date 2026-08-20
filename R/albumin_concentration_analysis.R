options(stringsAsFactors = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")

library(ggplot2)

dir.create(file.path("output", "concentration"), recursive = TRUE, showWarnings = FALSE)

parse_dt <- function(x) {
  x <- as.character(x)
  x[x == "" | is.na(x)] <- NA_character_
  out <- as.POSIXct(rep(NA_character_, length(x)), tz = "UTC")
  formats <- c("%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M", "%d/%m/%Y %H:%M:%S", "%Y-%m-%d")
  for (fmt in formats) {
    idx <- is.na(out) & !is.na(x)
    if (!any(idx)) break
    parsed <- as.POSIXct(strptime(x[idx], format = fmt, tz = "UTC"))
    idx_pos <- which(idx)
    out[idx_pos[!is.na(parsed)]] <- parsed[!is.na(parsed)]
  }
  out
}

read_clean <- function(path) {
  read.csv(path, check.names = FALSE, na.strings = c("", "NA"))
}

df <- read_clean(file.path("data", "df_final_v2_lt3_raw.csv"))
events <- read_clean(file.path("data", "albumin_events_48h_lt3.csv"))

df$icu_intime_dt <- parse_dt(df$icu_intime)
df$first_albumin_dt <- parse_dt(df$first_albumin_time)
df$first_albumin_h <- as.numeric(difftime(df$first_albumin_dt, df$icu_intime_dt, units = "hours"))

landmark <- df[
  df$prior_albumin_use == 0 &
    df$los_icu * 24 >= 6 &
    (is.na(df$first_albumin_h) | df$first_albumin_h >= 6),
]
landmark$exposed_6_48 <- !is.na(landmark$first_albumin_h) &
  landmark$first_albumin_h >= 6 &
  landmark$first_albumin_h < 48

exposed_ids <- landmark$stay_id[landmark$exposed_6_48]
events_window <- events[
  events$stay_id %in% exposed_ids &
    events$hours_from_icu_intime >= 6 &
    events$hours_from_icu_intime < 48,
]

events_window$concentration <- factor(
  events_window$concentration,
  levels = c("4%", "5%", "20%", "25%")
)
events_window$tonicity <- ifelse(
  events_window$concentration %in% c("20%", "25%"),
  "Hyperoncotic (20%/25%)",
  "Iso-/low-oncotic (4%/5%)"
)

event_summary <- aggregate(
  cbind(n_events = rep(1, nrow(events_window)), total_dose_g = events_window$dose_g, total_volume_ml = events_window$amount_ml) ~ concentration + tonicity,
  data = events_window,
  FUN = sum,
  na.rm = TRUE
)
event_summary$n_patients <- vapply(seq_len(nrow(event_summary)), function(i) {
  z <- events_window[
    events_window$concentration == event_summary$concentration[i] &
      events_window$tonicity == event_summary$tonicity[i],
  ]
  length(unique(z$stay_id))
}, integer(1))
event_summary$event_pct <- round(event_summary$n_events / sum(event_summary$n_events) * 100, 1)
event_summary$patient_pct <- round(event_summary$n_patients / length(unique(events_window$stay_id)) * 100, 1)
event_summary$dose_pct <- round(event_summary$total_dose_g / sum(event_summary$total_dose_g, na.rm = TRUE) * 100, 1)
event_summary$volume_pct <- round(event_summary$total_volume_ml / sum(event_summary$total_volume_ml, na.rm = TRUE) * 100, 1)
event_summary <- event_summary[order(event_summary$concentration), ]

first_events <- events_window[order(events_window$stay_id, events_window$hours_from_icu_intime), ]
first_events <- first_events[!duplicated(first_events$stay_id), ]
first_conc_summary <- as.data.frame(table(first_events$concentration, useNA = "ifany"))
names(first_conc_summary) <- c("first_concentration", "n_patients")
first_conc_summary$pct_patients <- round(first_conc_summary$n_patients / sum(first_conc_summary$n_patients) * 100, 1)

patient_conc <- split(events_window, events_window$stay_id)
patient_profile <- do.call(rbind, lapply(names(patient_conc), function(id) {
  z <- patient_conc[[id]]
  has_iso_low <- any(z$concentration %in% c("4%", "5%"), na.rm = TRUE)
  has_hyper <- any(z$concentration %in% c("20%", "25%"), na.rm = TRUE)
  profile <- if (has_iso_low && has_hyper) {
    "Mixed"
  } else if (has_hyper) {
    "Hyperoncotic only"
  } else {
    "Iso-/low-oncotic only"
  }
  data.frame(
    stay_id = as.integer(id),
    profile = profile,
    n_events = nrow(z),
    total_dose_g = sum(z$dose_g, na.rm = TRUE),
    total_volume_ml = sum(z$amount_ml, na.rm = TRUE)
  )
}))
profile_summary <- as.data.frame(table(patient_profile$profile))
names(profile_summary) <- c("profile", "n_patients")
profile_summary$pct_patients <- round(profile_summary$n_patients / sum(profile_summary$n_patients) * 100, 1)

patient_conc_dose <- aggregate(dose_g ~ stay_id + concentration + tonicity, data = events_window, FUN = sum, na.rm = TRUE)

write.csv(event_summary, file.path("output", "concentration", "albumin_concentration_event_summary.csv"), row.names = FALSE)
write.csv(first_conc_summary, file.path("output", "concentration", "albumin_first_concentration_summary.csv"), row.names = FALSE)
write.csv(profile_summary, file.path("output", "concentration", "albumin_patient_profile_summary.csv"), row.names = FALSE)
write.csv(patient_conc_dose, file.path("output", "concentration", "albumin_patient_concentration_dose.csv"), row.names = FALSE)

count_label <- function(n, pct) paste0(n, " (", pct, "%)")

p1_data <- first_conc_summary
p1_data$label <- count_label(p1_data$n_patients, p1_data$pct_patients)
p1 <- ggplot(p1_data, aes(x = first_concentration, y = n_patients, fill = first_concentration)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = label), vjust = -0.35, size = 3.6) +
  scale_fill_manual(values = c("4%" = "#7EA6E0", "5%" = "#2C7FB8", "20%" = "#D99152", "25%" = "#C0392B"), guide = "none", drop = FALSE) +
  labs(
    x = "First albumin concentration in 6-<48 h",
    y = "Patients",
    title = "First Albumin Concentration Among Treated Landmark Patients",
    subtitle = paste0("N = ", length(exposed_ids), " treated patients")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path("output", "concentration", "fig_albumin_first_concentration.png"), p1, width = 7, height = 5, dpi = 300, device = grDevices::png)
ggsave(file.path("output", "concentration", "fig_albumin_first_concentration.pdf"), p1, width = 7, height = 5)

p2_data <- profile_summary
p2_data$profile <- factor(p2_data$profile, levels = c("Iso-/low-oncotic only", "Hyperoncotic only", "Mixed"))
p2_data$label <- count_label(p2_data$n_patients, p2_data$pct_patients)
p2 <- ggplot(p2_data, aes(x = profile, y = n_patients, fill = profile)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(aes(label = label), vjust = -0.35, size = 3.6) +
  scale_fill_manual(values = c("Iso-/low-oncotic only" = "#2C7FB8", "Hyperoncotic only" = "#C0392B", "Mixed" = "#6A4C93"), guide = "none") +
  labs(
    x = NULL,
    y = "Patients",
    title = "Patient-Level Albumin Concentration Pattern",
    subtitle = "All albumin events during 6-<48 h"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path("output", "concentration", "fig_albumin_patient_profile.png"), p2, width = 7, height = 5, dpi = 300, device = grDevices::png)
ggsave(file.path("output", "concentration", "fig_albumin_patient_profile.pdf"), p2, width = 7, height = 5)

p3 <- ggplot(patient_conc_dose, aes(x = concentration, y = dose_g, fill = concentration)) +
  geom_boxplot(width = 0.55, outlier.alpha = 0.35) +
  geom_jitter(width = 0.12, alpha = 0.45, size = 1.8) +
  scale_fill_manual(values = c("4%" = "#7EA6E0", "5%" = "#2C7FB8", "20%" = "#D99152", "25%" = "#C0392B"), guide = "none", drop = FALSE) +
  labs(
    x = "Albumin concentration",
    y = "Patient-level dose in 6-<48 h, g",
    title = "Dose Distribution by Albumin Concentration",
    subtitle = "A patient can contribute to more than one concentration if mixed products were used"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path("output", "concentration", "fig_albumin_dose_by_concentration.png"), p3, width = 7, height = 5, dpi = 300, device = grDevices::png)
ggsave(file.path("output", "concentration", "fig_albumin_dose_by_concentration.pdf"), p3, width = 7, height = 5)

summary_lines <- c(
  "# Albumin Concentration Analysis",
  "",
  paste0("Landmark exposed patients (6-<48 h): ", length(exposed_ids)),
  paste0("Albumin events in 6-<48 h: ", nrow(events_window)),
  "",
  "## Event-Level Summary",
  paste(capture.output(print(event_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## First Concentration Summary",
  paste(capture.output(print(first_conc_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Patient-Level Concentration Pattern",
  paste(capture.output(print(profile_summary, row.names = FALSE)), collapse = "\n")
)
writeLines(summary_lines, file.path("output", "concentration", "albumin_concentration_summary.md"))

cat("Done. Outputs written to output/concentration\n")
