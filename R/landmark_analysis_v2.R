options(stringsAsFactors = FALSE)

required <- c("survival", "splines", "sandwich", "lmtest")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}

library(survival)
library(splines)

dir.create("output", showWarnings = FALSE, recursive = TRUE)
dirs <- list(
  data = file.path("output", "data"),
  tables = file.path("output", "tables"),
  fig_main = file.path("output", "figures", "main"),
  fig_suppl = file.path("output", "figures", "supplement")
)
for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

parse_dt <- function(x) {
  x <- as.character(x)
  x[x == "" | is.na(x)] <- NA_character_
  out <- as.POSIXct(rep(NA_character_, length(x)), tz = "UTC")
  formats <- c(
    "%Y-%m-%d %H:%M:%S",
    "%Y/%m/%d %H:%M",
    "%d/%m/%Y %H:%M:%S",
    "%Y-%m-%d"
  )
  for (fmt in formats) {
    idx <- is.na(out) & !is.na(x)
    if (!any(idx)) break
    parsed <- as.POSIXct(strptime(x[idx], format = fmt, tz = "UTC"))
    idx_pos <- which(idx)
    out[idx_pos[!is.na(parsed)]] <- parsed[!is.na(parsed)]
  }
  out
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

impute_median_mode <- function(df, covariates) {
  for (v in covariates) {
    if (!v %in% names(df)) next
    if (is.numeric(df[[v]]) || is.integer(df[[v]])) {
      med <- suppressWarnings(stats::median(df[[v]], na.rm = TRUE))
      if (!is.finite(med)) med <- 0
      df[[v]][is.na(df[[v]])] <- med
    } else {
      m <- mode_value(df[[v]])
      if (is.na(m)) m <- "Missing"
      df[[v]][is.na(df[[v]]) | df[[v]] == ""] <- m
      df[[v]] <- factor(df[[v]])
    }
  }
  df
}

compress_categories <- function(x, min_n = 20) {
  tab <- table(x, useNA = "no")
  rare <- names(tab)[tab < min_n]
  x[x %in% rare] <- "Other"
  factor(x)
}

weighted_mean <- function(x, w) sum(w * x, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)
weighted_var <- function(x, w) {
  mu <- weighted_mean(x, w)
  sum(w * (x - mu)^2, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)
}

smd_numeric <- function(x, g, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  x0 <- x[g == "B"]; x1 <- x[g == "A"]
  w0 <- w[g == "B"]; w1 <- w[g == "A"]
  m0 <- weighted_mean(x0, w0); m1 <- weighted_mean(x1, w1)
  v0 <- weighted_var(x0, w0); v1 <- weighted_var(x1, w1)
  den <- sqrt((v0 + v1) / 2)
  if (!is.finite(den) || den == 0) return(0)
  (m1 - m0) / den
}

smd_factor_max <- function(x, g, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  x <- factor(x)
  vals <- levels(x)
  if (length(vals) <= 1) return(0)
  max(abs(vapply(vals, function(v) smd_numeric(as.numeric(x == v), g, w), numeric(1))), na.rm = TRUE)
}

make_outcome <- function(df, horizon_days) {
  dod <- parse_dt(df$dod)
  landmark <- df$icu_intime_dt + 6 * 3600
  days <- ceiling(as.numeric(difftime(dod, landmark, units = "days")))
  event <- ifelse(!is.na(days) & days >= 0 & days <= horizon_days, 1L, 0L)
  time <- ifelse(event == 1L, days, horizon_days)
  data.frame(time = pmax(as.numeric(time), 1e-4), status = event)
}

make_clones <- function(df, horizon_days, window_start_h = 6, window_end_h = 48) {
  out <- make_outcome(df, horizon_days)
  base <- df
  base$out_time_days <- out$time
  base$out_status <- out$status
  base$first_albumin_landmark_h <- base$first_albumin_h - window_start_h
  base$exposed_window <- !is.na(base$first_albumin_h) &
    base$first_albumin_h >= window_start_h &
    base$first_albumin_h < window_end_h
  grace_h <- window_end_h - window_start_h

  a <- base
  a$strategy <- "A"
  a$art_censored <- ifelse(a$exposed_window, 0L, 1L)
  a$censor_h <- ifelse(a$art_censored == 1L, grace_h, NA_real_)

  b <- base
  b$strategy <- "B"
  b$art_censored <- ifelse(b$exposed_window, 1L, 0L)
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

make_person_period <- function(cloned) {
  pieces <- vector("list", nrow(cloned))
  for (i in seq_len(nrow(cloned))) {
    n_int <- max(1L, min(as.integer(ceiling(cloned$pp_exit_h[i])), as.integer(cloned$grace_h[i])))
    z <- cloned[rep(i, n_int), , drop = FALSE]
    z$interval <- seq_len(n_int) - 1L
    z$interval_start <- z$interval
    z$interval_end <- z$interval + 1
    z$censor_event <- 0L
    if (cloned$art_censored[i] == 1L && !is.na(cloned$censor_h[i])) {
      if (!(cloned$out_status[i] == 1L && cloned$out_time_days[i] * 24 <= cloned$censor_h[i])) {
        censor_interval <- min(max(floor(cloned$censor_h[i]), 0), n_int - 1L)
        z$censor_event[z$interval == censor_interval] <- 1L
      }
    }
    pieces[[i]] <- z
  }
  do.call(rbind, pieces)
}

fit_ipcw <- function(cloned, covariates, time_df = 3, trim = c(0.01, 0.99)) {
  keep_cols <- unique(c(
    "clone_id", "stay_id", "strategy", "art_censored", "censor_h",
    "out_status", "out_time_days", "pp_exit_h", "grace_h",
    "follow_time_days", "event", covariates
  ))
  lean_cloned <- cloned[, keep_cols[keep_cols %in% names(cloned)], drop = FALSE]
  pp <- make_person_period(lean_cloned)
  covariates <- covariates[covariates %in% names(pp)]
  rhs_time <- paste0("ns(interval, df = ", time_df, ")")
  f_num <- as.formula(paste("censor_event ~ strategy +", rhs_time, "+ strategy:", rhs_time))
  f_den <- as.formula(paste("censor_event ~ strategy +", rhs_time, "+ strategy:", rhs_time, "+", paste(covariates, collapse = " + ")))

  num <- glm(f_num, data = pp, family = binomial())
  den <- glm(f_den, data = pp, family = binomial(), control = list(maxit = 100))

  p_num <- pmin(pmax(1 - predict(num, type = "response"), 1e-6), 1)
  p_den <- pmin(pmax(1 - predict(den, type = "response"), 1e-6), 1)
  pp$sw_step <- p_num / p_den
  pp$sw <- ave(pp$sw_step, pp$clone_id, FUN = cumprod)
  final <- pp[!duplicated(pp$clone_id, fromLast = TRUE), c("clone_id", "sw")]
  q <- stats::quantile(final$sw, probs = trim, na.rm = TRUE)
  final$sw_trunc <- pmin(pmax(final$sw, q[1]), q[2])
  merged <- merge(cloned, final, by = "clone_id", all.x = TRUE, sort = FALSE)
  attr(merged, "weight_limits") <- q
  attr(merged, "person_period_n") <- nrow(pp)
  merged
}

estimate_endpoint <- function(df, covariates, horizon_days, label) {
  cat("Running", label, "endpoint", horizon_days, "days...\n")
  cloned <- make_clones(df, horizon_days)
  weighted <- fit_ipcw(cloned, covariates)

  cox <- survival::coxph(
    survival::Surv(follow_time_days, event) ~ strategy + cluster(stay_id),
    data = weighted,
    weights = sw_trunc,
    robust = TRUE
  )
  ci <- suppressMessages(confint(cox))
  hr <- exp(coef(cox)[["strategyA"]])
  log_hr <- coef(cox)[["strategyA"]]
  se_log_hr <- sqrt(diag(vcov(cox)))[["strategyA"]]
  hr_low <- exp(ci["strategyA", 1])
  hr_high <- exp(ci["strategyA", 2])

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

  ph <- tryCatch(survival::cox.zph(cox)$table["strategy", "p"], error = function(e) NA_real_)

  weight_limits <- attr(weighted, "weight_limits")
  weight_summary <- data.frame(
    analysis = label,
    endpoint = horizon_days,
    n_clones = nrow(weighted),
    person_period_rows = attr(weighted, "person_period_n"),
    weight_min = min(weighted$sw, na.rm = TRUE),
    weight_p01 = unname(weight_limits[1]),
    weight_median = median(weighted$sw, na.rm = TRUE),
    weight_p99 = unname(weight_limits[2]),
    weight_max = max(weighted$sw, na.rm = TRUE),
    weight_trunc_min = min(weighted$sw_trunc, na.rm = TRUE),
    weight_trunc_max = max(weighted$sw_trunc, na.rm = TRUE)
  )

  result <- data.frame(
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
    hr = as.numeric(hr),
    hr_low = as.numeric(hr_low),
    hr_high = as.numeric(hr_high),
    ph_p = as.numeric(ph)
  )

  list(result = result, weights = weight_summary, weighted = weighted)
}

pool_scalar <- function(q, se) {
  ok <- is.finite(q) & is.finite(se)
  q <- q[ok]
  se <- se[ok]
  m <- length(q)
  if (m == 0) return(c(qbar = NA, se = NA, low = NA, high = NA))
  qbar <- mean(q)
  ubar <- mean(se^2)
  b <- if (m > 1) stats::var(q) else 0
  total <- ubar + (1 + 1 / m) * b
  pooled_se <- sqrt(total)
  if (!is.finite(pooled_se)) pooled_se <- NA_real_
  crit <- 1.96
  c(qbar = qbar, se = pooled_se, low = qbar - crit * pooled_se, high = qbar + crit * pooled_se)
}

pool_results <- function(results) {
  keys <- unique(results[, c("analysis", "endpoint")])
  rows <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    sub <- results[results$analysis == keys$analysis[i] & results$endpoint == keys$endpoint[i], ]
    rd_p <- pool_scalar(sub$rd, sub$rd_se)
    hr_p <- pool_scalar(sub$log_hr, sub$se_log_hr)
    rows[[i]] <- data.frame(
      analysis = keys$analysis[i],
      endpoint = keys$endpoint[i],
      n_patients = sub$n_patients[1],
      exposed_6_48 = sub$exposed_6_48[1],
      risk_a = mean(sub$risk_a, na.rm = TRUE),
      risk_b = mean(sub$risk_b, na.rm = TRUE),
      rd = rd_p[["qbar"]],
      rd_low = rd_p[["low"]],
      rd_high = rd_p[["high"]],
      hr = exp(hr_p[["qbar"]]),
      hr_low = exp(hr_p[["low"]]),
      hr_high = exp(hr_p[["high"]]),
      ph_p_min = min(sub$ph_p, na.rm = TRUE),
      ph_p_max = max(sub$ph_p, na.rm = TRUE),
      m = nrow(sub)
    )
  }
  do.call(rbind, rows)
}

balance_table <- function(weighted, covariates) {
  rows <- lapply(covariates[covariates %in% names(weighted)], function(v) {
    x <- weighted[[v]]
    unweighted <- if (is.numeric(x) || is.integer(x)) {
      smd_numeric(as.numeric(x), weighted$strategy)
    } else {
      smd_factor_max(x, weighted$strategy)
    }
    weighted_smd <- if (is.numeric(x) || is.integer(x)) {
      smd_numeric(as.numeric(x), weighted$strategy, weighted$sw_trunc)
    } else {
      smd_factor_max(x, weighted$strategy, weighted$sw_trunc)
    }
    data.frame(variable = v, smd_unweighted = unweighted, smd_weighted = weighted_smd)
  })
  do.call(rbind, rows)
}

df <- read.csv("data/df_final_v2_lt3_raw.csv", check.names = FALSE)
df$icu_intime_dt <- parse_dt(df$icu_intime)
df$first_albumin_dt <- parse_dt(df$first_albumin_time)
df$first_albumin_h <- as.numeric(difftime(df$first_albumin_dt, df$icu_intime_dt, units = "hours"))

df$race <- compress_categories(df$race, min_n = 30)
df$icu_type <- compress_categories(df$icu_type, min_n = 30)
df$gender <- factor(df$gender)

landmark <- df[
  df$prior_albumin_use == 0 &
    df$los_icu * 24 >= 6 &
    (is.na(df$first_albumin_h) | df$first_albumin_h >= 6),
]
landmark$exposed_6_48 <- !is.na(landmark$first_albumin_h) &
  landmark$first_albumin_h >= 6 &
  landmark$first_albumin_h < 48

events54_path <- "data/albumin_events_54h_lt3.csv"
first54 <- data.frame(stay_id = integer(), first_albumin_h_54 = numeric())
if (file.exists(events54_path)) {
  events54 <- read.csv(events54_path, check.names = FALSE)
  if (nrow(events54) > 0) {
    first54 <- aggregate(hours_from_icu_intime ~ stay_id, data = events54, FUN = min)
    names(first54)[2] <- "first_albumin_h_54"
  }
}
landmark <- merge(landmark, first54, by = "stay_id", all.x = TRUE, sort = FALSE)
landmark$exposed_48_54 <- !is.na(landmark$first_albumin_h_54) &
  landmark$first_albumin_h_54 >= 48 &
  landmark$first_albumin_h_54 < 54

for (ep in c(28, 90, 365)) {
  outcome_ep <- make_outcome(landmark, ep)
  landmark[[paste0("mortality_", ep, "d_aux")]] <- outcome_ep$status
  landmark[[paste0("follow_time_", ep, "d_aux")]] <- outcome_ep$time
}
landmark$exposed_6_48_aux <- as.integer(landmark$exposed_6_48)

counts <- data.frame(
  metric = c(
    "source_rows",
    "no_prior_albumin",
    "icu_stay_ge_6h_no_prior",
    "albumin_0_to_lt6h_among_no_prior_icu_ge6h",
    "landmark_cohort",
    "albumin_6_to_lt48h",
    "albumin_48_to_lt54h",
    "albumin_6_to_lt54h"
  ),
  value = c(
    nrow(df),
    sum(df$prior_albumin_use == 0),
    sum(df$prior_albumin_use == 0 & df$los_icu * 24 >= 6),
    sum(df$prior_albumin_use == 0 & df$los_icu * 24 >= 6 & !is.na(df$first_albumin_h) & df$first_albumin_h >= 0 & df$first_albumin_h < 6),
    nrow(landmark),
    sum(landmark$exposed_6_48),
    sum(landmark$exposed_48_54, na.rm = TRUE),
    sum(landmark$exposed_6_48 | landmark$exposed_48_54)
  )
)
write.csv(counts, "output/landmark_counts.csv", row.names = FALSE)

main_covs <- c(
  "gender", "admission_age", "race", "weight", "charlson", "icu_type", "gcs",
  "albumin_closest", "creatinine_max", "bun_baseline", "glucose_baseline",
  "sodium_baseline", "potassium_baseline", "chloride_baseline",
  "bicarbonate_baseline", "aniongap_baseline", "calcium_baseline",
  "wbc_baseline", "hemoglobin_baseline", "platelet_baseline",
  "inr_baseline", "pt_baseline", "aptt_baseline",
  "alt_baseline", "ast_baseline", "alp_baseline", "bilirubin_total_enz_baseline",
  "furosemide_eq_dose_6h", "diuretics_use_baseline",
  "vaso_use_at_admission", "vaso_ne_dose_max_6h",
  "mech_vent_baseline", "rrt_baseline", "urineoutput_6h_rate"
)
main_covs <- main_covs[main_covs %in% names(landmark)]

sensitivity_covs <- setdiff(main_covs, c(
  "furosemide_eq_dose_6h", "diuretics_use_baseline",
  "vaso_use_at_admission", "vaso_ne_dose_max_6h",
  "mech_vent_baseline", "rrt_baseline", "urineoutput_6h_rate"
))

imp_covs <- unique(c(main_covs, sensitivity_covs))
imputation_aux_vars <- c(
  "exposed_6_48_aux",
  "mortality_28d_aux", "follow_time_28d_aux",
  "mortality_90d_aux", "follow_time_90d_aux",
  "mortality_365d_aux", "follow_time_365d_aux"
)
imputation_aux_vars <- imputation_aux_vars[imputation_aux_vars %in% names(landmark)]
use_mice <- requireNamespace("mice", quietly = TRUE)
if (use_mice) {
  cat("Running MICE multiple imputation (m=5, maxit=20, exposure/outcome auxiliary predictors)...\n")
  imp_input <- landmark[, unique(c(imp_covs, imputation_aux_vars)), drop = FALSE]
  for (v in names(imp_input)) {
    if (is.character(imp_input[[v]])) imp_input[[v]] <- factor(imp_input[[v]])
  }
  imp_method <- mice::make.method(imp_input)
  imp_pred <- mice::make.predictorMatrix(imp_input)
  imp_method[imputation_aux_vars] <- ""
  imp_pred[imputation_aux_vars, ] <- 0
  set.seed(20260614)
  mids_obj <- mice::mice(
    imp_input,
    m = 5,
    maxit = 20,
    method = imp_method,
    predictorMatrix = imp_pred,
    printFlag = FALSE
  )
  saveRDS(mids_obj, "output/mids_landmark_v2.rds")
  imputed_datasets <- lapply(seq_len(mids_obj$m), function(i) {
    tmp <- landmark
    completed <- mice::complete(mids_obj, i)
    tmp[, imp_covs] <- completed[, imp_covs]
    tmp
  })
  imputation_label <- "mice_m5_maxit20_aux_exposure_outcome"
} else {
  cat("mice not available; using deterministic median/mode imputation.\n")
  imputed_datasets <- list(impute_median_mode(landmark, imp_covs))
  imputation_label <- "deterministic_median_mode"
}

endpoints <- c(28, 90, 365)
all_main_runs <- list()
all_sens_runs <- list()
for (i in seq_along(imputed_datasets)) {
  cat("Imputation", i, "of", length(imputed_datasets), "\n")
  main_runs_i <- lapply(endpoints, function(ep) estimate_endpoint(imputed_datasets[[i]], main_covs, ep, "main_6h_landmark_6_48h"))
  sens_runs_i <- lapply(endpoints, function(ep) estimate_endpoint(imputed_datasets[[i]], sensitivity_covs, ep, "sensitivity_no_early_treatment_covariates"))
  for (j in seq_along(main_runs_i)) main_runs_i[[j]]$result$imputation <- i
  for (j in seq_along(sens_runs_i)) sens_runs_i[[j]]$result$imputation <- i
  all_main_runs <- c(all_main_runs, main_runs_i)
  all_sens_runs <- c(all_sens_runs, sens_runs_i)
}

main_results_by_imp <- do.call(rbind, lapply(all_main_runs, `[[`, "result"))
sens_results_by_imp <- do.call(rbind, lapply(all_sens_runs, `[[`, "result"))
main_results <- pool_results(main_results_by_imp)
sens_results <- pool_results(sens_results_by_imp)
write.csv(main_results_by_imp, "output/main_results_by_imputation.csv", row.names = FALSE)
write.csv(sens_results_by_imp, "output/sensitivity_no_early_treatment_results_by_imputation.csv", row.names = FALSE)
write.csv(main_results, "output/main_results.csv", row.names = FALSE)
write.csv(sens_results, "output/sensitivity_no_early_treatment_results.csv", row.names = FALSE)

weights <- rbind(
  do.call(rbind, lapply(all_main_runs, `[[`, "weights")),
  do.call(rbind, lapply(all_sens_runs, `[[`, "weights"))
)
write.csv(weights, "output/weight_summary.csv", row.names = FALSE)

bal <- balance_table(all_main_runs[[1]]$weighted, main_covs)
original_group <- ifelse(imputed_datasets[[1]]$exposed_6_48, "A", "B")
bal$smd_unweighted <- vapply(bal$variable, function(v) {
  x <- imputed_datasets[[1]][[v]]
  if (is.numeric(x) || is.integer(x)) {
    smd_numeric(as.numeric(x), original_group)
  } else {
    smd_factor_max(x, original_group)
  }
}, numeric(1))
write.csv(bal, "output/balance_main.csv", row.names = FALSE)

ph_tests <- rbind(
  main_results_by_imp[, c("analysis", "endpoint", "imputation", "ph_p")],
  sens_results_by_imp[, c("analysis", "endpoint", "imputation", "ph_p")]
)
write.csv(ph_tests, "output/ph_tests.csv", row.names = FALSE)

multiply_symbol <- intToUtf8(0x00D7)
superscript_three <- intToUtf8(0x00B3)
micro_symbol <- intToUtf8(0x00B5)

var_labels <- c(
  stay_id = "Stay ID",
  gender = "Sex",
  admission_age = "Age, years",
  race = "Race/Ethnicity",
  weight = "Weight, kg",
  charlson = "Charlson Comorbidity Index",
  icu_type = "ICU Type",
  gcs = "GCS Score",
  albumin_closest = "Baseline Albumin, g/dL",
  albumin_min = "Minimum Albumin, g/dL",
  creatinine_max = "Creatinine, mg/dL",
  bun_baseline = "BUN, mg/dL",
  glucose_baseline = "Glucose, mg/dL",
  sodium_baseline = "Sodium, mEq/L",
  potassium_baseline = "Potassium, mEq/L",
  chloride_baseline = "Chloride, mEq/L",
  bicarbonate_baseline = "Bicarbonate, mEq/L",
  aniongap_baseline = "Anion Gap, mEq/L",
  calcium_baseline = "Calcium, mg/dL",
  wbc_baseline = paste0("WBC, ", multiply_symbol, "10", superscript_three, "/", micro_symbol, "L"),
  hemoglobin_baseline = "Hemoglobin, g/dL",
  platelet_baseline = paste0("Platelets, ", multiply_symbol, "10", superscript_three, "/", micro_symbol, "L"),
  inr_baseline = "INR",
  pt_baseline = "PT, s",
  aptt_baseline = "aPTT, s",
  alt_baseline = "ALT, U/L",
  ast_baseline = "AST, U/L",
  alp_baseline = "ALP, U/L",
  bilirubin_total_enz_baseline = "Total Bilirubin, mg/dL",
  furosemide_eq_dose_6h = "Furosemide Equiv. Dose 0-6 h, mg",
  diuretics_use_baseline = "Diuretic Use 0-6 h",
  vaso_use_at_admission = "Vasopressor Use 0-6 h",
  vaso_ne_dose_max_6h = paste0("NE Equiv. Dose 0-6 h, ", micro_symbol, "g/kg/min"),
  mech_vent_baseline = "Mechanical Ventilation 0-6 h",
  rrt_baseline = "RRT 0-6 h",
  urineoutput_6h_rate = "Urine Output 0-6 h, mL/h",
  los_icu = "ICU LOS, days",
  los_hospital = "Hospital LOS, days"
)

label_var <- function(v) {
  out <- unname(var_labels[v])
  ifelse(is.na(out), v, out)
}

format_results_table <- function(x) {
  data.frame(
    Endpoint = paste0(x$endpoint, "-day"),
    `Risk A` = sprintf("%.1f%%", 100 * x$risk_a),
    `Risk B` = sprintf("%.1f%%", 100 * x$risk_b),
    `RD (95% CI)` = sprintf(
      "%.1f%% (%.1f%% to %.1f%%)",
      100 * x$rd, 100 * x$rd_low, 100 * x$rd_high
    ),
    `HR (95% CI)` = sprintf("%.3f (%.3f to %.3f)", x$hr, x$hr_low, x$hr_high),
    check.names = FALSE
  )
}

factor_for_display <- function(x) {
  if (is.factor(x)) return(x)
  if (is.character(x)) return(factor(x))
  ux <- sort(unique(x[!is.na(x)]))
  if (length(ux) <= 2 && all(ux %in% c(0, 1))) {
    return(factor(x, levels = c(0, 1), labels = c("No", "Yes")))
  }
  x
}

make_table1 <- function(dat, path) {
  if (!requireNamespace("tableone", quietly = TRUE)) return(invisible(NULL))
  table1_cont <- intersect(c(
    "admission_age", "weight", "charlson", "gcs",
    "albumin_closest", "albumin_min",
    "creatinine_max", "bun_baseline", "glucose_baseline",
    "sodium_baseline", "potassium_baseline", "chloride_baseline",
    "bicarbonate_baseline", "aniongap_baseline", "calcium_baseline",
    "hemoglobin_baseline", "platelet_baseline", "wbc_baseline",
    "inr_baseline", "pt_baseline", "aptt_baseline",
    "alt_baseline", "ast_baseline", "alp_baseline",
    "bilirubin_total_enz_baseline",
    "furosemide_eq_dose_6h", "vaso_ne_dose_max_6h",
    "urineoutput_6h_rate", "los_icu", "los_hospital"
  ), names(dat))
  table1_cat <- intersect(c(
    "gender", "race", "icu_type",
    "diuretics_use_baseline", "vaso_use_at_admission",
    "mech_vent_baseline", "rrt_baseline"
  ), names(dat))
  dat$exposed_6_48_f <- factor(
    as.integer(dat$exposed_6_48),
    levels = c(0, 1),
    labels = c("No Albumin", "Albumin")
  )
  for (v in table1_cat) dat[[v]] <- factor_for_display(dat[[v]])
  tab <- tableone::CreateTableOne(
    vars = c(table1_cat, table1_cont),
    strata = "exposed_6_48_f",
    data = dat,
    factorVars = table1_cat,
    test = TRUE,
    smd = TRUE
  )
  tab_print <- print(
    tab,
    nonnormal = table1_cont,
    exact = table1_cat,
    smd = TRUE,
    showAllLevels = TRUE,
    printToggle = FALSE,
    noSpaces = TRUE
  )
  rn <- rownames(tab_print)
  for (v in names(var_labels)) rn <- gsub(v, var_labels[[v]], rn, fixed = TRUE)
  rownames(tab_print) <- rn
  write.csv(tab_print, path, fileEncoding = "UTF-8")
}

make_positivity_table <- function(dat, path) {
  greater_equal <- intToUtf8(0x2265)
  cut_albumin <- function(x) cut(x, c(-Inf, 2.0, 2.5, Inf), labels = c("<2.0", "2.0-2.5", paste0(greater_equal, "2.5")), right = FALSE)
  cut_charlson <- function(x) cut(x, c(-Inf, 3, 6, Inf), labels = c("0-3", "4-6", paste0(greater_equal, "7")))
  cut_gcs <- function(x) cut(x, c(-Inf, 8, 12, Inf), labels = c("3-8", "9-12", "13-15"))
  cut_vaso <- function(x) factor(ifelse(x == 0, "No vasopressor", ifelse(x <= 0.1, "<=0.1", ifelse(x <= 0.3, "0.1-0.3", ">0.3"))))
  tmp <- dat
  tmp$alb_strata <- cut_albumin(tmp$albumin_closest)
  tmp$charlson_strata <- cut_charlson(tmp$charlson)
  tmp$gcs_strata <- cut_gcs(tmp$gcs)
  tmp$vaso_strata <- cut_vaso(tmp$vaso_ne_dose_max_6h)
  vars <- intersect(c("alb_strata", "charlson_strata", "gcs_strata", "icu_type", "vaso_strata", "mech_vent_baseline", "rrt_baseline"), names(tmp))
  rows <- lapply(vars, function(v) {
    pieces <- split(tmp, tmp[[v]], drop = TRUE)
    do.call(rbind, lapply(names(pieces), function(level) {
      z <- pieces[[level]]
      data.frame(
        variable = v,
        stratum = level,
        n = nrow(z),
        n_exposed = sum(z$exposed_6_48, na.rm = TRUE),
        n_unexposed = sum(!z$exposed_6_48, na.rm = TRUE),
        pct_exposed = round(mean(z$exposed_6_48, na.rm = TRUE) * 100, 1)
      )
    }))
  })
  write.csv(do.call(rbind, rows), path, row.names = FALSE)
}

make_surv_curve <- function(weighted, endpoint) {
  sf <- survival::survfit(
    survival::Surv(follow_time_days, event) ~ strategy,
    data = weighted,
    weights = sw_trunc
  )
  ss <- summary(sf, times = seq(0, endpoint, by = 0.5), extend = TRUE)
  data.frame(
    time = ss$time,
    cif = 1 - ss$surv,
    strategy = gsub("strategy=", "", ss$strata)
  )
}

make_km_plot <- function(curve, weighted, pooled_row, endpoint) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  library(ggplot2)
  risk_times <- if (endpoint <= 28) seq(0, endpoint, 7) else if (endpoint <= 90) seq(0, endpoint, 15) else seq(0, endpoint, 60)
  nar <- do.call(rbind, lapply(c("A", "B"), function(s) {
    z <- weighted[weighted$strategy == s, ]
    data.frame(strategy = s, time = risk_times, n_risk = vapply(risk_times, function(t) sum(z$follow_time_days >= t), integer(1)))
  }))
  colors <- c(A = "#C0392B", B = "#2C7FB8")
  labels <- c(A = "Strategy A (albumin 6-<48 h)", B = "Strategy B (no albumin 6-<48 h)")
  p_main <- ggplot(curve, aes(time, cif * 100, color = strategy, linetype = strategy)) +
    geom_step(linewidth = 0.9) +
    scale_color_manual(values = colors, labels = labels) +
    scale_linetype_manual(values = c(A = "solid", B = "dashed"), labels = labels) +
    scale_x_continuous(limits = c(0, endpoint), breaks = risk_times) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    annotate(
      "text",
      x = endpoint * 0.50,
      y = max(curve$cif * 100, na.rm = TRUE) * 0.95,
      hjust = 0,
      vjust = 1,
      size = 3.2,
      label = sprintf(
        "RD: %.1f%% (%.1f%% to %.1f%%)\nHR: %.3f (%.3f to %.3f)",
        100 * pooled_row$rd, 100 * pooled_row$rd_low, 100 * pooled_row$rd_high,
        pooled_row$hr, pooled_row$hr_low, pooled_row$hr_high
      )
    ) +
    labs(x = "Days from 6-hour landmark", y = "Cumulative mortality", color = NULL, linetype = NULL,
         title = sprintf("Weighted Kaplan-Meier Estimate of %d-Day Mortality", endpoint)) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  nar_wide <- reshape(nar, idvar = "time", timevar = "strategy", direction = "wide")
  p_nar <- ggplot(nar_wide, aes(time)) +
    geom_text(aes(y = 1.4, label = n_risk.A), color = colors[["A"]], size = 2.8) +
    geom_text(aes(y = 0.5, label = n_risk.B), color = colors[["B"]], size = 2.8) +
    annotate("text", x = -endpoint * 0.04, y = 1.4, label = "A", hjust = 1, color = colors[["A"]]) +
    annotate("text", x = -endpoint * 0.04, y = 0.5, label = "B", hjust = 1, color = colors[["B"]]) +
    scale_x_continuous(limits = c(-endpoint * 0.08, endpoint), breaks = risk_times) +
    scale_y_continuous(limits = c(0, 2)) +
    labs(title = "Number at risk") +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 9))
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_main / p_nar + patchwork::plot_layout(heights = c(4, 1))
  } else {
    p_main
  }
}

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  miss_tbl <- data.frame(
    variable = imp_covs,
    label = vapply(imp_covs, label_var, character(1)),
    n_missing = vapply(imp_covs, function(v) sum(is.na(landmark[[v]])), integer(1)),
    pct_missing = vapply(imp_covs, function(v) mean(is.na(landmark[[v]])) * 100, numeric(1))
  )
  miss_tbl <- miss_tbl[order(-miss_tbl$pct_missing), ]
  write.csv(miss_tbl, file.path(dirs$tables, "table_s_missing_rates_v2.csv"), row.names = FALSE)
  p_miss <- ggplot(miss_tbl, aes(x = reorder(label, pct_missing), y = pct_missing, fill = pct_missing > 20)) +
    geom_col() +
    geom_hline(yintercept = 20, linetype = "dashed", color = "#B91C1C") +
    coord_flip() +
    scale_fill_manual(values = c("FALSE" = "#2C7FB8", "TRUE" = "#C0392B"), guide = "none") +
    labs(x = NULL, y = "Missing (%)", title = "Missing Rates of Baseline Covariates") +
    theme_minimal(base_size = 11)
  ggsave(file.path(dirs$fig_suppl, "fig_s_missing_rates.pdf"), p_miss, width = 8, height = max(6, 0.25 * nrow(miss_tbl)), limitsize = FALSE)

  if (requireNamespace("naniar", quietly = TRUE)) {
    p_pattern <- naniar::vis_miss(landmark[, imp_covs, drop = FALSE], sort_miss = TRUE, cluster = TRUE) +
      labs(title = "Missing Data Pattern") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
    ggsave(file.path(dirs$fig_suppl, "fig_s_missing_pattern.pdf"), p_pattern, width = 12, height = 6, limitsize = FALSE)
  }

  if (exists("mids_obj")) {
    imputed_var_names <- names(mids_obj$method[mids_obj$method != ""])
    key_diag_vars <- head(intersect(c("creatinine_max", "albumin_closest", "gcs", "hemoglobin_baseline", "platelet_baseline", "bun_baseline"), imputed_var_names), 6)
    if (length(key_diag_vars) > 0) {
      grDevices::pdf(file.path(dirs$fig_suppl, "fig_s_mice_convergence.pdf"), width = 7, height = 4)
      for (vv in key_diag_vars) try(print(plot(mids_obj, y = vv)), silent = TRUE)
      grDevices::dev.off()
    }

    # Generate the actual publication-quality density diagnostic. The previous
    # logic reused key_diag_vars selected for convergence plots; those variables
    # had only 2-4 missing observations, so no density figure was produced and a
    # later legacy script silently filled the slot with a missing-rate bar chart.
    figure_s4_script <- "E:/workstation/tte/Figure_S4_mice_density.R"
    if (file.exists(figure_s4_script)) {
      source(figure_s4_script, local = FALSE)
    } else {
      warning("Figure S4 source not found: ", figure_s4_script)
    }
  }

  ps_dat <- imputed_datasets[[1]]
  ps_formula <- as.formula(paste("as.integer(exposed_6_48) ~", paste(main_covs, collapse = " + ")))
  ps_fit <- tryCatch(glm(ps_formula, data = ps_dat, family = binomial(), control = list(maxit = 100)), error = function(e) NULL)
  if (!is.null(ps_fit)) {
    ps_plot_data <- data.frame(
      ps = predict(ps_fit, type = "response"),
      group = factor(ps_dat$exposed_6_48, levels = c(FALSE, TRUE), labels = c("No Albumin", "Albumin"))
    )
    p_ps <- ggplot(ps_plot_data, aes(ps, fill = group, color = group)) +
      geom_density(alpha = 0.30, linewidth = 0.6) +
      scale_fill_manual(values = c("No Albumin" = "#2C7FB8", "Albumin" = "#C0392B")) +
      scale_color_manual(values = c("No Albumin" = "#2C7FB8", "Albumin" = "#C0392B")) +
      labs(x = "Predicted probability of albumin initiation", y = "Density", title = "Propensity Score Distribution", fill = NULL, color = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")
    ggsave(file.path(dirs$fig_suppl, "fig_s_propensity_score.pdf"), p_ps, width = 7, height = 5)
  }

  censor_dat <- all_main_runs[[1]]$weighted
  censor_b <- censor_dat[censor_dat$strategy == "B" & censor_dat$art_censored == 1L, ]
  if (nrow(censor_b) > 0) {
    p_censor <- ggplot(censor_b, aes(censor_h + 6)) +
      geom_histogram(aes(y = after_stat(density)), bins = 36, fill = "#2C7FB8", color = "white", alpha = 0.75) +
      geom_density(color = "#C0392B", linewidth = 0.8) +
      labs(x = "Time of artificial censoring (hours from ICU admission)", y = "Density", title = "Strategy B: Artificial Censoring Times") +
      theme_minimal(base_size = 11)
    ggsave(file.path(dirs$fig_suppl, "fig_s_censoring_times_stratB.pdf"), p_censor, width = 7, height = 4.5)
  }

  wt1 <- all_main_runs[[1]]$weighted
  wt_plot <- rbind(
    data.frame(stage = "Before truncation", strategy = wt1$strategy, sw = wt1$sw),
    data.frame(stage = "After truncation", strategy = wt1$strategy, sw = wt1$sw_trunc)
  )
  p_w_hist <- ggplot(wt_plot, aes(sw, fill = strategy)) +
    geom_histogram(bins = 50, alpha = 0.65, position = "identity") +
    facet_wrap(~stage, scales = "free_x") +
    scale_fill_manual(values = c(A = "#C0392B", B = "#2C7FB8")) +
    labs(x = "Stabilized weight", y = "Count", title = "IPCW Weight Distribution", fill = "Strategy") +
    theme_minimal(base_size = 11)
  ggsave(file.path(dirs$fig_suppl, "fig_s_weight_distribution.pdf"), p_w_hist, width = 10, height = 5)
  p_w_box <- ggplot(wt_plot, aes(strategy, sw, fill = strategy)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.25) +
    facet_wrap(~stage) +
    scale_fill_manual(values = c(A = "#C0392B", B = "#2C7FB8"), guide = "none") +
    labs(x = "Strategy", y = "Stabilized weight", title = "Weight Distribution by Strategy") +
    theme_minimal(base_size = 11)
  ggsave(file.path(dirs$fig_suppl, "fig_s_weight_boxplot.pdf"), p_w_box, width = 8, height = 5)

  bal_plot <- bal
  bal_plot$label <- vapply(bal_plot$variable, label_var, character(1))
  bal_long <- rbind(
    data.frame(label = bal_plot$label, smd = abs(bal_plot$smd_unweighted), stage = "Unweighted"),
    data.frame(label = bal_plot$label, smd = abs(bal_plot$smd_weighted), stage = "IPCW-weighted")
  )
  bal_long$label <- factor(bal_long$label, levels = rev(bal_plot$label))
  p_love <- ggplot(bal_long, aes(smd, label, color = stage, shape = stage)) +
    geom_point(size = 2.3, alpha = 0.85) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Unweighted" = "#C0392B", "IPCW-weighted" = "#2C7FB8")) +
    labs(x = "|Standardized Mean Difference|", y = NULL, title = "Covariate Balance Before and After IPCW", color = NULL, shape = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
  ggsave(file.path(dirs$fig_suppl, "fig_s_love_plot.pdf"), p_love, width = 8, height = max(6, 0.3 * nrow(bal_plot)), limitsize = FALSE)

  pp1 <- make_person_period(wt1)
  key_balance_vars <- intersect(c("albumin_closest", "creatinine_max", "charlson", "vaso_ne_dose_max_6h", "hemoglobin_baseline", "gcs"), main_covs)
  key_balance_vars <- key_balance_vars[sapply(key_balance_vars, function(v) is.numeric(pp1[[v]]))]
  if (length(key_balance_vars) > 0) {
    pp1$time_bin <- floor(pp1$interval / 4) * 4 + 6
    split_pp <- split(pp1, list(pp1$time_bin, pp1$strategy), drop = TRUE)
    means <- do.call(rbind, lapply(names(split_pp), function(nm) {
      z <- split_pp[[nm]]
      keys <- strsplit(nm, "\\.")[[1]]
      data.frame(
        time_bin = as.numeric(keys[1]),
        strategy = keys[2],
        variable = key_balance_vars,
        mean_unweighted = vapply(key_balance_vars, function(v) mean(z[[v]], na.rm = TRUE), numeric(1)),
        mean_weighted = vapply(key_balance_vars, function(v) weighted_mean(z[[v]], z$sw_trunc), numeric(1))
      )
    }))
    means_long <- rbind(
      data.frame(time_bin = means$time_bin, strategy = means$strategy, variable = means$variable, mean_val = means$mean_unweighted, weighting = "Unweighted"),
      data.frame(time_bin = means$time_bin, strategy = means$strategy, variable = means$variable, mean_val = means$mean_weighted, weighting = "IPCW-weighted")
    )
    means_long$variable_label <- vapply(means_long$variable, label_var, character(1))
    p_balance_time <- ggplot(means_long, aes(time_bin, mean_val, color = strategy, linetype = strategy)) +
      geom_line(linewidth = 0.7) +
      facet_grid(variable_label ~ weighting, scales = "free_y") +
      scale_color_manual(values = c(A = "#C0392B", B = "#2C7FB8")) +
      labs(x = "Hours from ICU admission", y = "Mean", title = "Covariate Means Over Time by Strategy", color = "Strategy", linetype = "Strategy") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom", strip.text.y = element_text(size = 7, angle = 0, hjust = 0))
    ggsave(file.path(dirs$fig_suppl, "fig_s_balance_over_time.pdf"), p_balance_time, width = 10, height = max(8, 2 * length(key_balance_vars)), limitsize = FALSE)
  }

  km_curve_data <- list()
  for (ep in endpoints) {
    run_index <- which(vapply(all_main_runs, function(x) x$result$endpoint[1] == ep && x$result$imputation[1] == 1, logical(1)))[1]
    curve_list <- lapply(seq_along(imputed_datasets), function(i) {
      idx <- which(vapply(all_main_runs, function(x) x$result$endpoint[1] == ep && x$result$imputation[1] == i, logical(1)))[1]
      make_surv_curve(all_main_runs[[idx]]$weighted, ep)
    })
    grid <- seq(0, ep, by = 0.5)
    pooled_curve <- do.call(rbind, lapply(c("A", "B"), function(s) {
      mat <- sapply(curve_list, function(cu) {
        zz <- cu[cu$strategy == s, ]
        approx(zz$time, zz$cif, xout = grid, method = "constant", f = 0, rule = 2)$y
      })
      data.frame(time = grid, cif = rowMeans(mat, na.rm = TRUE), strategy = s)
    }))
    pooled_row <- main_results[main_results$endpoint == ep, ]
    p_km <- make_km_plot(pooled_curve, all_main_runs[[run_index]]$weighted, pooled_row, ep)
    outpath <- if (ep == 28) file.path(dirs$fig_main, "fig_km_28d.pdf") else file.path(dirs$fig_suppl, sprintf("fig_s_km_%dd.pdf", ep))
    ggsave(outpath, p_km, width = 8, height = 6.5)
    ggsave(sub("\\.pdf$", ".png", outpath), p_km, width = 8, height = 6.5, dpi = 300, device = grDevices::png)
    km_curve_data[[as.character(ep)]] <- pooled_curve
  }
  saveRDS(km_curve_data, file.path(dirs$data, "km_curve_data.rds"))

  figure_s9_script <- "E:/workstation/tte/Figure_S9_schoenfeld_28d_all_imputations.R"
  if (file.exists(figure_s9_script)) {
    source(figure_s9_script, local = FALSE)
  } else {
    warning("Figure S9 script not found: ", figure_s9_script)
  }

  figure_s10_script <- "E:/workstation/tte/Figure_S10_km_90d.R"
  if (file.exists(figure_s10_script)) {
    source(figure_s10_script, local = FALSE)
  } else {
    warning("Figure S10 script not found: ", figure_s10_script)
  }

  figure_s11_script <- "E:/workstation/tte/Figure_S11_km_365d.R"
  if (file.exists(figure_s11_script)) {
    source(figure_s11_script, local = FALSE)
  } else {
    warning("Figure S11 script not found: ", figure_s11_script)
  }
}

write.csv(counts, file.path(dirs$tables, "table_flow_v2.csv"), row.names = FALSE)
write.csv(main_results, file.path(dirs$tables, "table2_pooled_results_raw_v2.csv"), row.names = FALSE)
write.csv(format_results_table(main_results), file.path(dirs$tables, "table2_pooled_results_v2.csv"), row.names = FALSE)
write.csv(sens_results, file.path(dirs$tables, "table_s_sensitivity_no_early_treatment_raw_v2.csv"), row.names = FALSE)
write.csv(format_results_table(sens_results), file.path(dirs$tables, "table_s_sensitivity_no_early_treatment_v2.csv"), row.names = FALSE)
write.csv(weights, file.path(dirs$tables, "table_s_weight_diagnostics_v2.csv"), row.names = FALSE)
write.csv(bal, file.path(dirs$tables, "table_s_covariate_balance_v2.csv"), row.names = FALSE)
write.csv(ph_tests, file.path(dirs$tables, "table_s_ph_tests_v2.csv"), row.names = FALSE)
write.csv(data.frame(
  window = c("6-48 h", "48-54 h", "6-54 h"),
  exposed = c(
    sum(landmark$exposed_6_48, na.rm = TRUE),
    sum(landmark$exposed_48_54, na.rm = TRUE),
    sum(landmark$exposed_6_48 | landmark$exposed_48_54, na.rm = TRUE)
  )
), file.path(dirs$tables, "table_s_window_48_54_check_v2.csv"), row.names = FALSE)
make_table1(landmark, file.path(dirs$tables, "table1_baseline_characteristics_v2.csv"))
make_table1(imputed_datasets[[1]], file.path(dirs$tables, "table1_imputed_supplement_v2.csv"))
make_positivity_table(imputed_datasets[[1]], file.path(dirs$tables, "table_s_positivity_v2.csv"))

if (file.exists("supplementary_outputs_v2.R")) {
  source("supplementary_outputs_v2.R")
}
if (file.exists("legacy_display_outputs_v2.R")) {
  source("legacy_display_outputs_v2.R")
}

coverage_lines <- c(
  "# v2 Output Coverage",
  "",
  "## Status",
  "The v2 output now mirrors the main classes of v1 outputs: flow counts, baseline Table 1, pooled effect Table 2, sensitivity results, missing-data diagnostics, MICE diagnostics, positivity check, IPCW diagnostics, covariate balance, PH checks, and KM figures.",
  "",
  "## Design changes from v1",
  "- Time zero is the 6-hour landmark rather than ICU admission.",
  "- Exposure is albumin initiation at 6-<48 hours rather than 0-48 hours.",
  "- Patients with albumin before 6 hours are excluded from the landmark cohort.",
  "- All result text, table captions, and figure legends should use the 6-hour landmark wording.",
  "",
  "## Core files for manuscript update",
  "- Main Table 1: output/tables/table1_baseline_characteristics_v2.csv",
  "- Main Table 2: output/tables/table2_pooled_results_v2.csv",
  "- Main Figure 2 or KM figure: output/figures/main/fig_km_28d.pdf",
  "- Supplementary sensitivity table: output/tables/table_s_sensitivity_no_early_treatment_v2.csv",
  "- Supplementary diagnostics: output/tables/table_s_*.csv and output/figures/supplement/*.pdf"
)
writeLines(coverage_lines, file.path("output", "v2_output_coverage.md"))

summary_lines <- c(
  "# Landmark Analysis v2 Summary",
  "",
  "## Counts",
  paste(capture.output(print(counts, row.names = FALSE)), collapse = "\n"),
  "",
  "## Main Results",
  paste(capture.output(print(main_results, row.names = FALSE)), collapse = "\n"),
  "",
  "## Sensitivity Results",
  paste(capture.output(print(sens_results, row.names = FALSE)), collapse = "\n"),
  "",
  "## Notes",
  paste0("- Imputation mode: ", imputation_label, "."),
  "- RD confidence intervals use an approximate Rubin-style combination of weighted Kaplan-Meier standard errors.",
  "- HR confidence intervals pool log-HR estimates using Rubin-style rules."
)
writeLines(summary_lines, "output/analysis_summary.md")

cat("Done. Outputs written to E:/workstation/tte/v2/output\n")
