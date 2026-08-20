options(stringsAsFactors = FALSE)

if (!requireNamespace("survival", quietly = TRUE)) stop("survival is required")
library(survival)

dirs <- list(
  data = file.path("output", "data"),
  tables = file.path("output", "tables"),
  fig_suppl = file.path("output", "figures", "supplement")
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

main_runs_path <- file.path(dirs$data, "imp_results_step12.rds")
if (!file.exists(main_runs_path)) stop("Missing ", main_runs_path, "; run landmark_analysis_v2.R first.")
main_runs <- readRDS(main_runs_path)

ph_rows <- vector("list", length(main_runs))
pdf_path <- file.path(dirs$fig_suppl, "fig_s_schoenfeld_all.pdf")
grDevices::pdf(pdf_path, width = 7, height = 5)
on.exit(grDevices::dev.off(), add = TRUE)

for (i in seq_along(main_runs)) {
  run <- main_runs[[i]]
  endpoint <- run$result$endpoint[1]
  imputation <- run$result$imputation[1]
  weighted <- run$weighted
  weighted$strategy <- factor(weighted$strategy, levels = c("B", "A"))
  title <- sprintf("Schoenfeld Residuals: Strategy, %s-day, imputation %s", endpoint, imputation)
  ph <- tryCatch(
    {
      cox <- survival::coxph(
        survival::Surv(follow_time_days, event) ~ strategy + cluster(stay_id),
        data = weighted,
        weights = sw_trunc,
        robust = TRUE
      )
      survival::cox.zph(cox)
    },
    error = function(e) e
  )
  if (inherits(ph, "error")) {
    plot.new()
    text(0.5, 0.5, paste(title, "\nFailed:", ph$message))
    ph_rows[[i]] <- data.frame(
      analysis = "main_6h_landmark_6_48h",
      endpoint = endpoint,
      imputation = imputation,
      ph_p = NA_real_,
      status = ph$message
    )
  } else {
    plot(ph, main = title)
    abline(h = 0, lty = 2, col = "grey50")
    p_val <- tryCatch(as.numeric(ph$table["strategy", "p"]), error = function(e) NA_real_)
    ph_rows[[i]] <- data.frame(
      analysis = "main_6h_landmark_6_48h",
      endpoint = endpoint,
      imputation = imputation,
      ph_p = p_val,
      status = "ok"
    )
  }
}

ph_all <- do.call(rbind, ph_rows)
write.csv(ph_all, file.path(dirs$tables, "table_s_ph_tests_main_all_v2.csv"), row.names = FALSE)
cat("Schoenfeld outputs written to ", pdf_path, "\n", sep = "")
