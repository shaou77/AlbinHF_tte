options(stringsAsFactors = FALSE)

cat("\n============================================================\n")
cat("  v2 supplementary output parity\n")
cat("============================================================\n\n")

write_table_pair <- function(x, filename) {
  write.csv(x, file.path(dirs$tables, filename), row.names = FALSE)
  v2_name <- sub("\\.csv$", "_v2.csv", filename)
  if (!identical(v2_name, filename)) {
    write.csv(x, file.path(dirs$tables, v2_name), row.names = FALSE)
  }
}

calc_ess <- function(w) {
  w <- w[is.finite(w)]
  if (length(w) == 0 || sum(w^2) == 0) return(NA_real_)
  (sum(w)^2) / sum(w^2)
}

format_ci_table <- function(x) {
  rd <- if ("rd" %in% names(x)) x$rd else x$RD
  rd_low <- if ("rd_low" %in% names(x)) x$rd_low else x$RD_lo
  rd_high <- if ("rd_high" %in% names(x)) x$rd_high else x$RD_hi
  hr <- if ("hr" %in% names(x)) x$hr else x$HR
  hr_low <- if ("hr_low" %in% names(x)) x$hr_low else x$HR_lo
  hr_high <- if ("hr_high" %in% names(x)) x$hr_high else x$HR_hi
  data.frame(
    Analysis = x$analysis,
    Endpoint = paste0(x$endpoint, "d"),
    `RD (95% CI)` = sprintf("%.1f%% (%.1f%% to %.1f%%)", 100 * rd, 100 * rd_low, 100 * rd_high),
    `HR (95% CI)` = sprintf("%.3f (%.3f to %.3f)", hr, hr_low, hr_high),
    check.names = FALSE
  )
}

result_to_old_names <- function(x, label) {
  data.frame(
    analysis = label,
    endpoint = x$endpoint,
    RD = x$rd,
    RD_lo = x$rd_low,
    RD_hi = x$rd_high,
    HR = x$hr,
    HR_lo = x$hr_low,
    HR_hi = x$hr_high,
    check.names = FALSE
  )
}

calc_evalue <- function(hr) {
  if (!is.finite(hr) || hr <= 0) return(NA_real_)
  if (hr < 1) hr <- 1 / hr
  hr + sqrt(hr * (hr - 1))
}

calc_evalue_ci <- function(hr_lo, hr_hi) {
  if (!is.finite(hr_lo) || !is.finite(hr_hi)) return(NA_real_)
  if (hr_lo <= 1 && hr_hi >= 1) return(1)
  calc_evalue(if (abs(hr_lo - 1) < abs(hr_hi - 1)) hr_lo else hr_hi)
}

make_clones_custom <- function(df, horizon_days, a_complier_col = "exposed_6_48",
                               b_deviation_col = "exposed_6_48",
                               window_start_h = 6, window_end_h = 48) {
  out <- make_outcome(df, horizon_days)
  base <- df
  base$out_time_days <- out$time
  base$out_status <- out$status
  base$first_albumin_landmark_h <- base$first_albumin_h - window_start_h
  grace_h <- window_end_h - window_start_h

  a_ok <- as.logical(base[[a_complier_col]])
  a_ok[is.na(a_ok)] <- FALSE
  b_dev <- as.logical(base[[b_deviation_col]])
  b_dev[is.na(b_dev)] <- FALSE

  a <- base
  a$strategy <- "A"
  a$art_censored <- ifelse(a_ok, 0L, 1L)
  a$censor_h <- ifelse(a$art_censored == 1L, grace_h, NA_real_)

  b <- base
  b$strategy <- "B"
  b$art_censored <- ifelse(b_dev, 1L, 0L)
  b$censor_h <- ifelse(b$art_censored == 1L, b$first_albumin_landmark_h, NA_real_)

  cloned <- rbind(a, b)
  cloned$strategy <- factor(cloned$strategy, levels = c("B", "A"))
  cloned$clone_id <- seq_len(nrow(cloned))

  event_h <- cloned$out_time_days * 24
  has_art <- cloned$art_censored == 1L & !is.na(cloned$censor_h)
  event_before_censor <- has_art & cloned$out_status == 1L & event_h <= cloned$censor_h
  censored_before_event <- has_art & !event_before_censor

  cloned$follow_time_days <- cloned$out_time_days
  cloned$event <- cloned$out_status
  cloned$follow_time_days[censored_before_event] <- pmax(cloned$censor_h[censored_before_event] / 24, 1e-4)
  cloned$event[censored_before_event] <- 0L
  cloned$pp_exit_h <- pmin(grace_h, ifelse(cloned$out_status == 1L, event_h, Inf), ifelse(has_art, cloned$censor_h, Inf), na.rm = TRUE)
  cloned$pp_exit_h[!is.finite(cloned$pp_exit_h)] <- grace_h
  cloned$grace_h <- grace_h
  cloned
}

estimate_endpoint_custom <- function(df, covariates, horizon_days, label,
                                     time_df = 3, trim = c(0.01, 0.99),
                                     a_complier_col = "exposed_6_48",
                                     b_deviation_col = "exposed_6_48") {
  cloned <- make_clones_custom(df, horizon_days, a_complier_col, b_deviation_col)
  weighted <- fit_ipcw(cloned, covariates, time_df = time_df, trim = trim)

  cox <- survival::coxph(
    survival::Surv(follow_time_days, event) ~ strategy + cluster(stay_id),
    data = weighted,
    weights = sw_trunc,
    robust = TRUE
  )
  ci <- suppressMessages(confint(cox))
  log_hr <- coef(cox)[["strategyA"]]
  se_log_hr <- sqrt(diag(vcov(cox)))[["strategyA"]]

  sf <- survival::survfit(
    survival::Surv(follow_time_days, event) ~ strategy,
    data = weighted,
    weights = sw_trunc
  )
  ss <- summary(sf, times = horizon_days, extend = TRUE)
  surv_by <- setNames(ss$surv, gsub("strategy=", "", ss$strata))
  se_by <- setNames(ss$std.err, gsub("strategy=", "", ss$strata))
  risk_a <- 1 - surv_by[["A"]]
  risk_b <- 1 - surv_by[["B"]]
  rd <- risk_a - risk_b
  rd_se <- sqrt(se_by[["A"]]^2 + se_by[["B"]]^2)

  data.frame(
    analysis = label,
    endpoint = horizon_days,
    n_patients = length(unique(df$stay_id)),
    exposed_6_48 = sum(df$exposed_6_48, na.rm = TRUE),
    risk_a = as.numeric(risk_a),
    risk_b = as.numeric(risk_b),
    rd = as.numeric(rd),
    rd_se = as.numeric(rd_se),
    log_hr = as.numeric(log_hr),
    se_log_hr = as.numeric(se_log_hr),
    hr = as.numeric(exp(log_hr)),
    hr_low = as.numeric(exp(ci["strategyA", 1])),
    hr_high = as.numeric(exp(ci["strategyA", 2])),
    ph_p = tryCatch(as.numeric(survival::cox.zph(cox)$table["strategy", "p"]), error = function(e) NA_real_)
  )
}

run_custom_pooled <- function(label, data_list, covariates, endpoints = c(28, 90, 365),
                              time_df = 3, trim = c(0.01, 0.99),
                              a_complier_col = "exposed_6_48",
                              b_deviation_col = "exposed_6_48") {
  rows <- list()
  k <- 1L
  for (i in seq_along(data_list)) {
    for (ep in endpoints) {
      res <- tryCatch(
        estimate_endpoint_custom(data_list[[i]], covariates, ep, label, time_df, trim, a_complier_col, b_deviation_col),
        error = function(e) {
          cat(sprintf("  %s imputation %d endpoint %s failed: %s\n", label, i, ep, e$message))
          NULL
        }
      )
      if (!is.null(res)) {
        res$imputation <- i
        rows[[k]] <- res
        k <- k + 1L
      }
    }
  }
  if (!length(rows)) return(NULL)
  pool_results(do.call(rbind, rows))
}

primary_runs <- all_main_runs[vapply(all_main_runs, function(x) x$result$endpoint[1] == 28, logical(1))]
if (length(primary_runs) > 0) {
  ess_rows <- lapply(seq_along(primary_runs), function(i) {
    z <- primary_runs[[i]]$weighted
    data.frame(
      imputation = i,
      n_truncated = sum(abs(z$sw - z$sw_trunc) > 1e-12, na.rm = TRUE),
      pct_truncated = round(100 * mean(abs(z$sw - z$sw_trunc) > 1e-12, na.rm = TRUE), 2),
      lo_bound = min(z$sw_trunc, na.rm = TRUE),
      hi_bound = max(z$sw_trunc, na.rm = TRUE),
      ess_A = round(calc_ess(z$sw_trunc[z$strategy == "A"]), 1),
      ess_B = round(calc_ess(z$sw_trunc[z$strategy == "B"]), 1),
      n_A = sum(z$strategy == "A"),
      n_B = sum(z$strategy == "B")
    )
  })
  write_table_pair(do.call(rbind, ess_rows), "table_s_weight_ess.csv")

  z <- primary_runs[[1]]$weighted
  strategy_deviation <- data.frame(
    Strategy = c("A (initiate albumin 6-<48h)", "B (avoid albumin 6-<48h)"),
    `Total clones` = c(sum(z$strategy == "A"), sum(z$strategy == "B")),
    `Compliers (n)` = c(sum(z$strategy == "A" & z$art_censored == 0), sum(z$strategy == "B" & z$art_censored == 0)),
    `Deviators (n)` = c(sum(z$strategy == "A" & z$art_censored == 1), sum(z$strategy == "B" & z$art_censored == 1)),
    `Deviation %` = round(c(
      mean(z$art_censored[z$strategy == "A"] == 1) * 100,
      mean(z$art_censored[z$strategy == "B"] == 1) * 100
    ), 1),
    `Censoring time` = c(
      "42h after landmark (fixed)",
      sprintf(
        "Median %.1f (IQR %.1f-%.1f) h after landmark",
        median(z$censor_h[z$strategy == "B" & z$art_censored == 1], na.rm = TRUE),
        quantile(z$censor_h[z$strategy == "B" & z$art_censored == 1], 0.25, na.rm = TRUE),
        quantile(z$censor_h[z$strategy == "B" & z$art_censored == 1], 0.75, na.rm = TRUE)
      )
    ),
    check.names = FALSE
  )
  write_table_pair(strategy_deviation, "table_s_strategy_deviation.csv")
}

evalue_table <- data.frame(
  endpoint = main_results$endpoint,
  HR = main_results$hr,
  HR_95CI = sprintf("%.3f to %.3f", main_results$hr_low, main_results$hr_high),
  Evalue_point = round(vapply(main_results$hr, calc_evalue, numeric(1)), 2),
  Evalue_CI = round(mapply(calc_evalue_ci, main_results$hr_low, main_results$hr_high), 2)
)
write_table_pair(evalue_table, "table_s_evalue.csv")

bootstrap_table_path <- file.path(dirs$tables, "table_s_bootstrap_ci.csv")
write_bootstrap_placeholder <- TRUE
if (file.exists(bootstrap_table_path)) {
  existing_boot <- tryCatch(read.csv(bootstrap_table_path, check.names = FALSE), error = function(e) NULL)
  write_bootstrap_placeholder <- is.null(existing_boot) ||
    !"boot_RD_lo" %in% names(existing_boot) ||
    all(is.na(existing_boot$boot_RD_lo))
}
if (write_bootstrap_placeholder) {
  bootstrap_note <- data.frame(
    endpoint = main_results$endpoint,
    rubin_RD = main_results$rd,
    rubin_RD_lo = main_results$rd_low,
    rubin_RD_hi = main_results$rd_high,
    rubin_HR = main_results$hr,
    rubin_HR_lo = main_results$hr_low,
    rubin_HR_hi = main_results$hr_high,
    boot_RD_lo = NA_real_,
    boot_RD_hi = NA_real_,
    boot_RD_median = NA_real_,
    boot_HR_lo = NA_real_,
    boot_HR_hi = NA_real_,
    boot_HR_median = NA_real_,
    note = "Run bootstrap_v2.R to populate fixed-imputation bootstrap intervals."
  )
  write_table_pair(bootstrap_note, "table_s_bootstrap_ci.csv")
  saveRDS(data.frame(note = bootstrap_note$note[1]), file.path(dirs$data, "boot_results.rds"))
}

outlier_log <- data.frame(
  variable = c("all covariates"),
  action = c("No additional v2 winsorization or outlier recoding was applied in this script."),
  note = c("v2 uses the extracted landmark dataset and records this audit row for parity with the original output.")
)
write_table_pair(outlier_log, "outlier_handling_log.csv")

if (!exists("miss_tbl")) {
  miss_tbl <- data.frame(
    variable = imp_covs,
    label = vapply(imp_covs, label_var, character(1)),
    n_missing = vapply(imp_covs, function(v) sum(is.na(landmark[[v]])), integer(1)),
    pct_missing = vapply(imp_covs, function(v) mean(is.na(landmark[[v]])) * 100, numeric(1))
  )
  miss_tbl <- miss_tbl[order(-miss_tbl$pct_missing), ]
}
missing_summary <- data.frame(
  variable = miss_tbl$variable,
  label = miss_tbl$label,
  n_missing = miss_tbl$n_missing,
  pct_missing = round(miss_tbl$pct_missing, 1),
  imputation_role = ifelse(miss_tbl$n_missing > 0, "included in imputation model", "complete in landmark cohort")
)
write_table_pair(missing_summary, "table_s_missing_summary.csv")

table1_smd <- data.frame(
  variable = bal$variable,
  label = vapply(bal$variable, label_var, character(1)),
  smd_unweighted = bal$smd_unweighted,
  smd_weighted = bal$smd_weighted
)
write_table_pair(table1_smd, "table1_smd_column.csv")

events_path <- file.path("data", "albumin_events_48h_lt3.csv")
landmark$albumin_dose_6_48h <- 0
if (file.exists(events_path)) {
  ev <- read.csv(events_path, check.names = FALSE)
  ev <- ev[ev$hours_from_icu_intime >= 6 & ev$hours_from_icu_intime < 48, , drop = FALSE]
  if (nrow(ev) > 0) {
    dose_by <- aggregate(dose_g ~ stay_id, data = ev, FUN = sum, na.rm = TRUE)
    dose_map <- setNames(dose_by$dose_g, dose_by$stay_id)
    landmark$albumin_dose_6_48h <- ifelse(
      as.character(landmark$stay_id) %in% names(dose_map),
      dose_map[as.character(landmark$stay_id)],
      0
    )
  }
}
landmark$dose_ge_25g_6_48h <- landmark$albumin_dose_6_48h >= 25
landmark$dose_ge_50g_6_48h <- landmark$albumin_dose_6_48h >= 50
for (i in seq_along(imputed_datasets)) {
  imputed_datasets[[i]]$albumin_dose_6_48h <- landmark$albumin_dose_6_48h
  imputed_datasets[[i]]$dose_ge_25g_6_48h <- landmark$dose_ge_25g_6_48h
  imputed_datasets[[i]]$dose_ge_50g_6_48h <- landmark$dose_ge_50g_6_48h
}

main_display <- result_to_old_names(main_results, "Main analysis")
early_display <- result_to_old_names(sens_results, "Exclude early treatment covariates")
icu48_list <- lapply(imputed_datasets, function(x) x[x$los_icu * 24 >= 48, , drop = FALSE])

greater_equal <- intToUtf8(0x2265)
extra_sens <- list(
  run_custom_pooled(paste("ICU", greater_equal, "48h"), icu48_list, main_covs),
  run_custom_pooled(paste("Dose", greater_equal, "25g vs 0g"), imputed_datasets, main_covs, a_complier_col = "dose_ge_25g_6_48h"),
  run_custom_pooled(paste("Dose", greater_equal, "50g vs 0g"), imputed_datasets, main_covs, a_complier_col = "dose_ge_50g_6_48h"),
  run_custom_pooled("Spline df = 2", imputed_datasets, main_covs, time_df = 2),
  run_custom_pooled("Spline df = 5", imputed_datasets, main_covs, time_df = 5),
  run_custom_pooled("Trim 0.5/99.5%", imputed_datasets, main_covs, trim = c(0.005, 0.995)),
  run_custom_pooled("Trim 2.5/97.5%", imputed_datasets, main_covs, trim = c(0.025, 0.975))
)
extra_sens <- Filter(Negate(is.null), extra_sens)
all_sensitivity <- do.call(rbind, c(list(main_display, early_display), lapply(extra_sens, function(x) {
  result_to_old_names(x, x$analysis[1])
})))
write_table_pair(all_sensitivity, "table_s_sensitivity_all_raw.csv")
write_table_pair(format_ci_table(all_sensitivity), "table_s_sensitivity_all.csv")
saveRDS(all_sensitivity, file.path(dirs$data, "all_sensitivity.rds"))
saveRDS(evalue_table, file.path(dirs$data, "evalue_table.rds"))

subgroup_defs <- list(
  list(name = "Albumin < 2.5 g/dL", rows = which(imputed_datasets[[1]]$albumin_closest < 2.5)),
  list(name = paste("Albumin", greater_equal, "2.5 g/dL"), rows = which(imputed_datasets[[1]]$albumin_closest >= 2.5)),
  list(name = "Vasopressor use", rows = which(imputed_datasets[[1]]$vaso_ne_dose_max_6h > 0)),
  list(name = "No vasopressor", rows = which(imputed_datasets[[1]]$vaso_ne_dose_max_6h == 0))
)

subgroup_rows <- lapply(subgroup_defs, function(sg) {
  data_list <- lapply(imputed_datasets, function(x) x[sg$rows, , drop = FALSE])
  n_total <- nrow(data_list[[1]])
  n_exposed <- sum(data_list[[1]]$exposed_6_48, na.rm = TRUE)
  n_unexposed <- n_total - n_exposed
  feasible <- n_exposed >= 10 && n_unexposed >= 10
  if (!feasible) {
    return(data.frame(
      subgroup = sg$name, n = n_total, n_exposed = n_exposed, feasible = FALSE,
      RD = NA_real_, RD_lo = NA_real_, RD_hi = NA_real_,
      HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_,
      note = sprintf("Insufficient: exposed=%d, unexposed=%d", n_exposed, n_unexposed)
    ))
  }
  pooled <- run_custom_pooled(sg$name, data_list, main_covs, endpoints = 28)
  if (is.null(pooled) || nrow(pooled) == 0) {
    return(data.frame(
      subgroup = sg$name, n = n_total, n_exposed = n_exposed, feasible = FALSE,
      RD = NA_real_, RD_lo = NA_real_, RD_hi = NA_real_,
      HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_, note = "Pipeline failed"
    ))
  }
  data.frame(
    subgroup = sg$name, n = n_total, n_exposed = n_exposed, feasible = TRUE,
    RD = pooled$rd[1], RD_lo = pooled$rd_low[1], RD_hi = pooled$rd_high[1],
    HR = pooled$hr[1], HR_lo = pooled$hr_low[1], HR_hi = pooled$hr_high[1],
    note = "Completed"
  )
})
subgroup_table <- do.call(rbind, subgroup_rows)
subgroup_display <- data.frame(
  Subgroup = subgroup_table$subgroup,
  N = subgroup_table$n,
  `N exposed` = subgroup_table$n_exposed,
  Feasible = ifelse(subgroup_table$feasible, "Yes", "No"),
  `RD (95% CI)` = ifelse(
    subgroup_table$feasible,
    sprintf("%.1f%% (%.1f%% to %.1f%%)", 100 * subgroup_table$RD, 100 * subgroup_table$RD_lo, 100 * subgroup_table$RD_hi),
    subgroup_table$note
  ),
  `HR (95% CI)` = ifelse(
    subgroup_table$feasible,
    sprintf("%.3f (%.3f to %.3f)", subgroup_table$HR, subgroup_table$HR_lo, subgroup_table$HR_hi),
    "-"
  ),
  check.names = FALSE
)
write_table_pair(subgroup_table, "table_s_subgroup_results_raw.csv")
write_table_pair(subgroup_display, "table_s_subgroup_results.csv")
write_table_pair(data.frame(
  subgroup = subgroup_table$subgroup,
  n_total = subgroup_table$n,
  n_exposed = subgroup_table$n_exposed,
  pct_exposed = round(100 * subgroup_table$n_exposed / subgroup_table$n, 1),
  note = subgroup_table$note
), "subgroup_feasibility_notes.csv")

summary_main <- data.frame(
  category = "Main", analysis = "Main analysis",
  RD = main_results$rd[main_results$endpoint == 28],
  RD_lo = main_results$rd_low[main_results$endpoint == 28],
  RD_hi = main_results$rd_high[main_results$endpoint == 28],
  HR = main_results$hr[main_results$endpoint == 28],
  HR_lo = main_results$hr_low[main_results$endpoint == 28],
  HR_hi = main_results$hr_high[main_results$endpoint == 28]
)
summary_sens <- all_sensitivity[all_sensitivity$endpoint == 28 & all_sensitivity$analysis != "Main analysis", ]
summary_sens$category <- "Sensitivity"
summary_sens <- summary_sens[, c("category", "analysis", "RD", "RD_lo", "RD_hi", "HR", "HR_lo", "HR_hi")]
summary_sub <- subgroup_table[subgroup_table$feasible, c("subgroup", "RD", "RD_lo", "RD_hi", "HR", "HR_lo", "HR_hi")]
names(summary_sub)[1] <- "analysis"
summary_sub$category <- "Subgroup"
summary_sub <- summary_sub[, c("category", "analysis", "RD", "RD_lo", "RD_hi", "HR", "HR_lo", "HR_hi")]
full_summary <- rbind(summary_main, summary_sens, summary_sub)
write_table_pair(full_summary, "table_s_comprehensive_summary_raw.csv")
write_table_pair(data.frame(
  Category = full_summary$category,
  Analysis = full_summary$analysis,
  `RD (95% CI)` = sprintf("%.1f%% (%.1f%% to %.1f%%)", 100 * full_summary$RD, 100 * full_summary$RD_lo, 100 * full_summary$RD_hi),
  `HR (95% CI)` = sprintf("%.3f (%.3f to %.3f)", full_summary$HR, full_summary$HR_lo, full_summary$HR_hi),
  check.names = FALSE
), "table_s_comprehensive_summary.csv")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  forest28 <- all_sensitivity[all_sensitivity$endpoint == 28, ]
  forest28$analysis <- factor(forest28$analysis, levels = rev(unique(forest28$analysis)))
  p_sens <- ggplot(forest28, aes(x = RD * 100, y = analysis, xmin = RD_lo * 100, xmax = RD_hi * 100)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(height = 0.25, linewidth = 0.6, colour = "#2C3E50") +
    geom_point(size = 3, colour = "#C0392B") +
    labs(x = "28-Day Risk Difference (%)", y = NULL, title = "Sensitivity Analyses: 28-Day Risk Difference") +
    theme_minimal(base_size = 11)
  ggsave(file.path(dirs$fig_suppl, "fig_s_sensitivity_forest.pdf"), p_sens, width = 9, height = max(4, 0.55 * nrow(forest28)))

  sub_forest <- subgroup_table[subgroup_table$feasible, ]
  if (nrow(sub_forest) > 0) {
    forest_sub <- rbind(
      data.frame(subgroup = "Overall (Main analysis)", RD = summary_main$RD, RD_lo = summary_main$RD_lo, RD_hi = summary_main$RD_hi, is_main = TRUE),
      data.frame(subgroup = sub_forest$subgroup, RD = sub_forest$RD, RD_lo = sub_forest$RD_lo, RD_hi = sub_forest$RD_hi, is_main = FALSE)
    )
    forest_sub$subgroup <- factor(forest_sub$subgroup, levels = rev(forest_sub$subgroup))
    p_sub <- ggplot(forest_sub, aes(x = RD * 100, y = subgroup, xmin = RD_lo * 100, xmax = RD_hi * 100)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_errorbarh(height = 0.25, linewidth = 0.6, colour = "#2C3E50") +
      geom_point(aes(shape = is_main, size = is_main), colour = "#C0392B") +
      scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16), guide = "none") +
      scale_size_manual(values = c("TRUE" = 4, "FALSE" = 3), guide = "none") +
      labs(x = "28-Day Risk Difference (%)", y = NULL, title = "Subgroup Analyses: 28-Day Risk Difference") +
      theme_minimal(base_size = 11)
    ggsave(file.path(dirs$fig_main, "fig_subgroup_forest.pdf"), p_sub, width = 9, height = max(3.5, 0.7 * nrow(forest_sub)))
    ggsave(file.path(dirs$fig_main, "fig_subgroup_forest.png"), p_sub, width = 9, height = max(3.5, 0.7 * nrow(forest_sub)), dpi = 300, device = grDevices::png)
  }

  full_forest <- full_summary
  full_forest$analysis <- factor(full_forest$analysis, levels = rev(full_forest$analysis))
  p_full <- ggplot(full_forest, aes(x = RD * 100, y = analysis, xmin = RD_lo * 100, xmax = RD_hi * 100, colour = category)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(height = 0.25, linewidth = 0.6) +
    geom_point(size = 3) +
    scale_colour_manual(values = c(Main = "#C0392B", Sensitivity = "#2C3E50", Subgroup = "#287D3C")) +
    labs(x = "28-Day Risk Difference (%)", y = NULL, title = "Comprehensive Summary: 28-Day Risk Difference", colour = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
  ggsave(file.path(dirs$fig_suppl, "fig_s_comprehensive_forest.pdf"), p_full, width = 10, height = max(5, 0.55 * nrow(full_forest)), limitsize = FALSE)
}

saveRDS(landmark, file.path(dirs$data, "df_analytic.rds"))
saveRDS(counts, file.path(dirs$data, "flowchart_data.rds"))
saveRDS(list(window_start_h = 6, window_end_h = 48, n_impute = length(imputed_datasets), endpoints_days = endpoints, weight_trim = c(0.01, 0.99), time_spline_df = 3), file.path(dirs$data, "params.rds"))
saveRDS(main_covs, file.path(dirs$data, "model_covariates.rds"))
saveRDS(main_covs, file.path(dirs$data, "dag_covariates.rds"))
saveRDS(var_labels, file.path(dirs$data, "var_labels.rds"))
saveRDS(main_results, file.path(dirs$data, "pooled_table.rds"))
saveRDS(main_results_by_imp, file.path(dirs$data, "all_effects.rds"))
saveRDS(all_main_runs, file.path(dirs$data, "imp_results_step12.rds"))
saveRDS(lapply(primary_runs, `[[`, "weighted"), file.path(dirs$data, "weighted_cloned_list.rds"))
saveRDS(lapply(primary_runs, function(x) make_clones(imputed_datasets[[x$result$imputation[1]]], 28)), file.path(dirs$data, "cloned_list.rds"))
if (exists("mids_obj")) saveRDS(mids_obj, file.path(dirs$data, "mids_obj.rds"))
save.image(file.path(dirs$data, "full_workspace.RData"))

expected_outputs <- c(
  "tables/outlier_handling_log.csv",
  "tables/table_s_bootstrap_ci.csv",
  "tables/table_s_comprehensive_summary.csv",
  "tables/table_s_covariate_balance_v2.csv",
  "tables/table_s_evalue.csv",
  "tables/table_s_missing_summary.csv",
  "tables/table_s_positivity_v2.csv",
  "tables/table_s_sensitivity_all.csv",
  "tables/table_s_strategy_deviation.csv",
  "tables/table_s_subgroup_results.csv",
  "tables/table_s_weight_diagnostics_v2.csv",
  "tables/table_s_weight_ess.csv",
  "tables/table1_baseline_characteristics_v2.csv",
  "tables/table1_imputed_supplement_v2.csv",
  "tables/table1_smd_column.csv",
  "tables/table2_pooled_results_v2.csv",
  "figures/main/fig_km_28d.pdf",
  "figures/main/fig_km_28d.png",
  "figures/main/fig_subgroup_forest.pdf",
  "figures/supplement/fig_s_comprehensive_forest.pdf",
  "figures/supplement/fig_s_sensitivity_forest.pdf"
)
audit <- data.frame(
  expected_output = expected_outputs,
  exists = file.exists(file.path("output", expected_outputs))
)
write.csv(audit, file.path("output", "v2_original_output_parity_audit.csv"), row.names = FALSE)

cat("Supplementary parity outputs complete.\n\n")
