suppressPackageStartupMessages(library(survival))

project_dir <- "E:/workstation/tte/v2"
weighted_path <- file.path(project_dir, "output", "data", "weighted_cloned_list.rds")
out_table <- file.path(project_dir, "output", "tables", "table_s_interaction_tests_v2.csv")
out_by_imp <- file.path(project_dir, "output", "data", "interaction_tests_by_imputation_v2.csv")
out_validation <- file.path(project_dir, "output", "tables", "table_s_interaction_validation_v2.csv")

weighted_list <- readRDS(weighted_path)
if (length(weighted_list) < 2L) {
  stop("At least two imputed weighted datasets are required.")
}

pool_scalar <- function(q, se) {
  ok <- is.finite(q) & is.finite(se) & se >= 0
  q <- q[ok]
  se <- se[ok]
  m <- length(q)
  if (m == 0L) {
    return(c(
      estimate = NA_real_, se = NA_real_, lower = NA_real_,
      upper = NA_real_, z = NA_real_, p = NA_real_, m = 0
    ))
  }
  qbar <- mean(q)
  ubar <- mean(se^2)
  b <- if (m > 1L) stats::var(q) else 0
  total <- ubar + (1 + 1 / m) * b
  pooled_se <- sqrt(total)
  z <- qbar / pooled_se
  p <- stats::pchisq(z^2, df = 1, lower.tail = FALSE)
  c(
    estimate = qbar,
    se = pooled_se,
    lower = qbar - 1.96 * pooled_se,
    upper = qbar + 1.96 * pooled_se,
    z = z,
    p = p,
    m = m
  )
}

weighted_km_rd <- function(data, horizon = 28) {
  fit <- survival::survfit(
    survival::Surv(follow_time_days, event) ~ strategy,
    data = data,
    weights = sw_trunc
  )
  sm <- summary(fit, times = horizon, extend = TRUE)
  risk <- setNames(1 - sm$surv, sub("strategy=", "", sm$strata, fixed = TRUE))
  risk_se <- setNames(sm$std.err, sub("strategy=", "", sm$strata, fixed = TRUE))
  if (!all(c("A", "B") %in% names(risk))) {
    stop("Both strategies are required to estimate a risk difference.")
  }
  c(
    risk_a = unname(risk[["A"]]),
    risk_b = unname(risk[["B"]]),
    rd = unname(risk[["A"]] - risk[["B"]]),
    rd_se = unname(sqrt(risk_se[["A"]]^2 + risk_se[["B"]]^2))
  )
}

subgroup_definitions <- list(
  vasopressor = list(
    label = "Vasopressor use during first 6 ICU hours",
    reference = "No vasopressor",
    comparison = "Vasopressor use",
    derive = function(data) as.integer(data$vaso_ne_dose_max_6h > 0)
  ),
  albumin = list(
    label = "Baseline albumin",
    reference = ">=2.5 g/dL",
    comparison = "<2.5 g/dL",
    derive = function(data) as.integer(data$albumin_closest < 2.5)
  )
)

by_imp_rows <- list()
validation_rows <- list()
k <- 1L
v <- 1L

for (subgroup_name in names(subgroup_definitions)) {
  definition <- subgroup_definitions[[subgroup_name]]

  for (imputation in seq_along(weighted_list)) {
    data <- weighted_list[[imputation]]
    data$trt <- as.integer(data$strategy == "A")
    data$subgroup_indicator <- definition$derive(data)

    if (anyNA(data$subgroup_indicator)) {
      stop(sprintf(
        "Missing subgroup membership for %s in imputation %d.",
        subgroup_name, imputation
      ))
    }
    if (!identical(sort(unique(data$subgroup_indicator)), c(0L, 1L))) {
      stop(sprintf(
        "Both subgroup levels are required for %s in imputation %d.",
        subgroup_name, imputation
      ))
    }

    cox <- survival::coxph(
      survival::Surv(follow_time_days, event) ~
        trt + trt:subgroup_indicator + strata(subgroup_indicator) + cluster(stay_id),
      data = data,
      weights = sw_trunc,
      robust = TRUE
    )
    coefficient_names <- names(stats::coef(cox))
    interaction_name <- coefficient_names[
      grepl("trt:subgroup_indicator", coefficient_names, fixed = TRUE)
    ]
    if (length(interaction_name) != 1L) {
      stop(sprintf(
        "Could not identify the interaction coefficient for %s, imputation %d.",
        subgroup_name, imputation
      ))
    }

    beta <- stats::coef(cox)
    variance <- stats::vcov(cox)
    interaction_log_hr <- unname(beta[[interaction_name]])
    interaction_se <- unname(sqrt(diag(variance))[[interaction_name]])
    log_hr_reference <- unname(beta[["trt"]])
    log_hr_comparison <- log_hr_reference + interaction_log_hr

    data_reference <- data[data$subgroup_indicator == 0L, , drop = FALSE]
    data_comparison <- data[data$subgroup_indicator == 1L, , drop = FALSE]
    rd_reference <- weighted_km_rd(data_reference)
    rd_comparison <- weighted_km_rd(data_comparison)
    rd_interaction <- rd_comparison[["rd"]] - rd_reference[["rd"]]
    rd_interaction_se <- sqrt(
      rd_comparison[["rd_se"]]^2 + rd_reference[["rd_se"]]^2
    )

    patient_data <- data[!duplicated(data$stay_id), , drop = FALSE]
    counts <- table(
      factor(patient_data$subgroup_indicator, levels = c(0, 1)),
      factor(patient_data$exposed_6_48, levels = c(FALSE, TRUE))
    )

    by_imp_rows[[k]] <- data.frame(
      subgroup = definition$label,
      reference = definition$reference,
      comparison = definition$comparison,
      imputation = imputation,
      n_reference = sum(counts[1, ]),
      exposed_reference = counts[1, 2],
      n_comparison = sum(counts[2, ]),
      exposed_comparison = counts[2, 2],
      log_hr_reference = log_hr_reference,
      hr_reference = exp(log_hr_reference),
      log_hr_comparison = log_hr_comparison,
      hr_comparison = exp(log_hr_comparison),
      interaction_log_hr = interaction_log_hr,
      interaction_se_log_hr = interaction_se,
      interaction_ratio_hr = exp(interaction_log_hr),
      interaction_p_hr = stats::pchisq(
        (interaction_log_hr / interaction_se)^2,
        df = 1,
        lower.tail = FALSE
      ),
      rd_reference = rd_reference[["rd"]],
      rd_reference_se = rd_reference[["rd_se"]],
      rd_comparison = rd_comparison[["rd"]],
      rd_comparison_se = rd_comparison[["rd_se"]],
      interaction_rd = rd_interaction,
      interaction_se_rd = rd_interaction_se,
      interaction_p_rd = stats::pchisq(
        (rd_interaction / rd_interaction_se)^2,
        df = 1,
        lower.tail = FALSE
      ),
      check.names = FALSE
    )
    k <- k + 1L
  }

  subgroup_by_imp <- do.call(rbind, by_imp_rows)
  subgroup_by_imp <- subgroup_by_imp[
    subgroup_by_imp$subgroup == definition$label,
    ,
    drop = FALSE
  ]

  pooled_hr_interaction <- pool_scalar(
    subgroup_by_imp$interaction_log_hr,
    subgroup_by_imp$interaction_se_log_hr
  )
  pooled_rd_interaction <- pool_scalar(
    subgroup_by_imp$interaction_rd,
    subgroup_by_imp$interaction_se_rd
  )
  pooled_rd_reference <- mean(subgroup_by_imp$rd_reference)
  pooled_rd_comparison <- mean(subgroup_by_imp$rd_comparison)

  validation_rows[[v]] <- data.frame(
    subgroup = definition$label,
    level = definition$reference,
    n = subgroup_by_imp$n_reference[1],
    n_exposed = subgroup_by_imp$exposed_reference[1],
    pooled_RD = pooled_rd_reference,
    pooled_HR = exp(mean(subgroup_by_imp$log_hr_reference)),
    check.names = FALSE
  )
  v <- v + 1L
  validation_rows[[v]] <- data.frame(
    subgroup = definition$label,
    level = definition$comparison,
    n = subgroup_by_imp$n_comparison[1],
    n_exposed = subgroup_by_imp$exposed_comparison[1],
    pooled_RD = pooled_rd_comparison,
    pooled_HR = exp(mean(subgroup_by_imp$log_hr_comparison)),
    check.names = FALSE
  )
  v <- v + 1L

  validation_rows[[v]] <- data.frame(
    subgroup = definition$label,
    level = "Interaction HR ratio",
    n = NA_integer_,
    n_exposed = NA_integer_,
    pooled_RD = NA_real_,
    pooled_HR = exp(pooled_hr_interaction[["estimate"]]),
    check.names = FALSE
  )
  v <- v + 1L

  validation_rows[[v]] <- data.frame(
    subgroup = definition$label,
    level = "Interaction RD difference",
    n = NA_integer_,
    n_exposed = NA_integer_,
    pooled_RD = pooled_rd_interaction[["estimate"]],
    pooled_HR = NA_real_,
    check.names = FALSE
  )
  v <- v + 1L
}

by_imp <- do.call(rbind, by_imp_rows)
write.csv(by_imp, out_by_imp, row.names = FALSE)

summary_rows <- list()
s <- 1L
for (subgroup_label in unique(by_imp$subgroup)) {
  data <- by_imp[by_imp$subgroup == subgroup_label, , drop = FALSE]
  pooled_hr <- pool_scalar(data$interaction_log_hr, data$interaction_se_log_hr)
  pooled_rd <- pool_scalar(data$interaction_rd, data$interaction_se_rd)

  summary_rows[[s]] <- data.frame(
    subgroup = subgroup_label,
    scale = "Hazard ratio",
    contrast = sprintf("%s vs %s", data$comparison[1], data$reference[1]),
    interaction_estimate = exp(pooled_hr[["estimate"]]),
    lower_95 = exp(pooled_hr[["lower"]]),
    upper_95 = exp(pooled_hr[["upper"]]),
    p_interaction = pooled_hr[["p"]],
    m = pooled_hr[["m"]],
    interpretation = "Ratio of subgroup-specific hazard ratios",
    check.names = FALSE
  )
  s <- s + 1L

  summary_rows[[s]] <- data.frame(
    subgroup = subgroup_label,
    scale = "Risk difference",
    contrast = sprintf("%s minus %s", data$comparison[1], data$reference[1]),
    interaction_estimate = 100 * pooled_rd[["estimate"]],
    lower_95 = 100 * pooled_rd[["lower"]],
    upper_95 = 100 * pooled_rd[["upper"]],
    p_interaction = pooled_rd[["p"]],
    m = pooled_rd[["m"]],
    interpretation = "Difference in 28-day risk differences, percentage points",
    check.names = FALSE
  )
  s <- s + 1L
}

summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, out_table, row.names = FALSE)

validation <- do.call(rbind, validation_rows)
write.csv(validation, out_validation, row.names = FALSE)

cat("\nFormal interaction tests pooled across imputations\n")
print(summary_table, row.names = FALSE)
cat("\nSubgroup estimate validation\n")
print(validation, row.names = FALSE)
cat("\nOutputs:\n")
cat(out_table, "\n")
cat(out_by_imp, "\n")
cat(out_validation, "\n")
