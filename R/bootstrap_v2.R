options(stringsAsFactors = FALSE)

required <- c("survival")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}
library(survival)

dirs <- list(
  data = file.path("output", "data"),
  tables = file.path("output", "tables")
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

n_boot <- as.integer(Sys.getenv("V2_BOOT_N", "200"))
if (!is.finite(n_boot) || n_boot < 20) n_boot <- 200L
seed <- as.integer(Sys.getenv("V2_BOOT_SEED", "20260619"))
set.seed(seed)

main_runs_path <- file.path(dirs$data, "imp_results_step12.rds")
main_results_path <- file.path("output", "main_results.csv")
if (!file.exists(main_runs_path)) stop("Missing ", main_runs_path, "; run landmark_analysis_v2.R first.")
if (!file.exists(main_results_path)) stop("Missing ", main_results_path, "; run landmark_analysis_v2.R first.")

main_runs <- readRDS(main_runs_path)
main_results <- read.csv(main_results_path, check.names = FALSE)
endpoints <- sort(unique(vapply(main_runs, function(x) x$result$endpoint[1], numeric(1))))
imputations <- sort(unique(vapply(main_runs, function(x) x$result$imputation[1], numeric(1))))

pool_scalar <- function(q, se) {
  ok <- is.finite(q) & is.finite(se)
  q <- q[ok]
  se <- se[ok]
  m <- length(q)
  if (m == 0) return(c(qbar = NA_real_, se = NA_real_, low = NA_real_, high = NA_real_))
  qbar <- mean(q)
  ubar <- mean(se^2)
  b <- if (m > 1) stats::var(q) else 0
  total <- ubar + (1 + 1 / m) * b
  pooled_se <- sqrt(total)
  c(qbar = qbar, se = pooled_se, low = qbar - 1.96 * pooled_se, high = qbar + 1.96 * pooled_se)
}

find_run <- function(endpoint, imputation) {
  idx <- which(vapply(
    main_runs,
    function(x) x$result$endpoint[1] == endpoint && x$result$imputation[1] == imputation,
    logical(1)
  ))
  if (length(idx) != 1) stop("Cannot locate run for endpoint=", endpoint, ", imputation=", imputation)
  main_runs[[idx]]
}

bootstrap_one_fit <- function(weighted, endpoint, boot_id) {
  ids <- unique(weighted$stay_id)
  sampled <- sample(ids, length(ids), replace = TRUE)
  boot_pieces <- vector("list", length(sampled))
  for (i in seq_along(sampled)) {
    z <- weighted[weighted$stay_id == sampled[i], , drop = FALSE]
    z$boot_cluster <- paste0(sampled[i], "_", i)
    boot_pieces[[i]] <- z
  }
  boot_dat <- do.call(rbind, boot_pieces)
  boot_dat$strategy <- factor(boot_dat$strategy, levels = c("B", "A"))

  cox <- survival::coxph(
    survival::Surv(follow_time_days, event) ~ strategy + cluster(boot_cluster),
    data = boot_dat,
    weights = sw_trunc,
    robust = TRUE
  )
  log_hr <- unname(coef(cox)[["strategyA"]])
  se_log_hr <- unname(sqrt(diag(vcov(cox)))[["strategyA"]])

  sf <- survival::survfit(
    survival::Surv(follow_time_days, event) ~ strategy,
    data = boot_dat,
    weights = sw_trunc
  )
  ss <- summary(sf, times = endpoint, extend = TRUE)
  surv_by <- setNames(ss$surv, gsub("strategy=", "", ss$strata))
  se_by <- setNames(ss$std.err, gsub("strategy=", "", ss$strata))
  risk_a <- 1 - surv_by[["A"]]
  risk_b <- 1 - surv_by[["B"]]
  rd <- risk_a - risk_b
  rd_se <- sqrt(se_by[["A"]]^2 + se_by[["B"]]^2)

  data.frame(
    boot_id = boot_id,
    endpoint = endpoint,
    rd = as.numeric(rd),
    rd_se = as.numeric(rd_se),
    log_hr = as.numeric(log_hr),
    se_log_hr = as.numeric(se_log_hr)
  )
}

boot_rows <- vector("list", n_boot * length(endpoints))
k <- 1L
cat(sprintf("Running fixed-imputation bootstrap: %d replicates x %d endpoints x %d imputations\n",
            n_boot, length(endpoints), length(imputations)))

for (b in seq_len(n_boot)) {
  if (b %% 10 == 0 || b == 1) cat("  bootstrap", b, "of", n_boot, "\n")
  for (ep in endpoints) {
    imp_rows <- lapply(imputations, function(m) {
      run <- find_run(ep, m)
      tryCatch(
        {
          z <- bootstrap_one_fit(run$weighted, ep, b)
          z$imputation <- m
          z
        },
        error = function(e) {
          data.frame(
            boot_id = b,
            endpoint = ep,
            rd = NA_real_,
            rd_se = NA_real_,
            log_hr = NA_real_,
            se_log_hr = NA_real_,
            imputation = m,
            error = e$message
          )
        }
      )
    })
    imp_df <- do.call(rbind, imp_rows)
    rd_pool <- pool_scalar(imp_df$rd, imp_df$rd_se)
    hr_pool <- pool_scalar(imp_df$log_hr, imp_df$se_log_hr)
    boot_rows[[k]] <- data.frame(
      boot_id = b,
      endpoint = ep,
      rd = rd_pool[["qbar"]],
      log_hr = hr_pool[["qbar"]],
      hr = exp(hr_pool[["qbar"]]),
      n_imputations = sum(is.finite(imp_df$rd)),
      failed_imputations = sum(!is.finite(imp_df$rd))
    )
    k <- k + 1L
  }
}

boot_df <- do.call(rbind, boot_rows)
saveRDS(boot_df, file.path(dirs$data, "boot_results.rds"))
write.csv(boot_df, file.path(dirs$data, "boot_results_v2.csv"), row.names = FALSE)

ci_rows <- lapply(endpoints, function(ep) {
  z <- boot_df[boot_df$endpoint == ep & is.finite(boot_df$rd) & is.finite(boot_df$log_hr), , drop = FALSE]
  if (nrow(z) < 20) {
    return(data.frame(
      endpoint = ep,
      boot_RD_lo = NA_real_,
      boot_RD_hi = NA_real_,
      boot_RD_median = NA_real_,
      boot_HR_lo = NA_real_,
      boot_HR_hi = NA_real_,
      boot_HR_median = NA_real_,
      n_boot_success = nrow(z)
    ))
  }
  rd_ci <- stats::quantile(z$rd, c(0.025, 0.975), na.rm = TRUE)
  loghr_ci <- stats::quantile(z$log_hr, c(0.025, 0.975), na.rm = TRUE)
  data.frame(
    endpoint = ep,
    boot_RD_lo = unname(rd_ci[1]),
    boot_RD_hi = unname(rd_ci[2]),
    boot_RD_median = stats::median(z$rd, na.rm = TRUE),
    boot_HR_lo = exp(unname(loghr_ci[1])),
    boot_HR_hi = exp(unname(loghr_ci[2])),
    boot_HR_median = exp(stats::median(z$log_hr, na.rm = TRUE)),
    n_boot_success = nrow(z)
  )
})
boot_ci <- do.call(rbind, ci_rows)

comparison <- merge(
  main_results[, c("endpoint", "rd", "rd_low", "rd_high", "hr", "hr_low", "hr_high")],
  boot_ci,
  by = "endpoint",
  all.x = TRUE,
  sort = FALSE
)
names(comparison)[names(comparison) == "rd"] <- "rubin_RD"
names(comparison)[names(comparison) == "rd_low"] <- "rubin_RD_lo"
names(comparison)[names(comparison) == "rd_high"] <- "rubin_RD_hi"
names(comparison)[names(comparison) == "hr"] <- "rubin_HR"
names(comparison)[names(comparison) == "hr_low"] <- "rubin_HR_lo"
names(comparison)[names(comparison) == "hr_high"] <- "rubin_HR_hi"
comparison$method <- "Fixed imputation; patient bootstrap applied to final IPCW-weighted cloned datasets"
comparison$n_boot_requested <- n_boot
comparison$seed <- seed

write.csv(comparison, file.path(dirs$tables, "table_s_bootstrap_ci.csv"), row.names = FALSE)
write.csv(comparison, file.path(dirs$tables, "table_s_bootstrap_ci_v2.csv"), row.names = FALSE)
cat("Bootstrap outputs written to output/tables/table_s_bootstrap_ci.csv and output/data/boot_results.rds\n")
