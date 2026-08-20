options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

root <- "E:/workstation/tte/v2"
output_root <- file.path(root, "output")
figure_dir <- file.path(output_root, "figures", "main", "final")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

km_all <- readRDS(file.path(output_root, "data", "km_curve_data.rds"))
cloned <- readRDS(file.path(output_root, "data", "weighted_cloned_list.rds"))
pooled <- read.csv(
  file.path(output_root, "tables", "table2_pooled_results_raw_v2.csv"),
  check.names = FALSE
)

stopifnot(length(cloned) == 5L, "28" %in% names(km_all))
result_28 <- pooled[pooled$endpoint == 28, , drop = FALSE]
stopifnot(nrow(result_28) == 1L, result_28$n_patients == 1308L, result_28$exposed_6_48 == 80L)

km <- km_all[["28"]]
km$strategy <- factor(km$strategy, levels = c("A", "B"))
stopifnot(all(c("time", "cif", "strategy") %in% names(km)))

risk_times <- c(0, 7, 14, 21, 28)
first_imputation <- cloned[[1]]
risk <- do.call(
  rbind,
  lapply(c("A", "B"), function(strategy_i) {
    d <- first_imputation[first_imputation$strategy == strategy_i, , drop = FALSE]
    data.frame(
      strategy = strategy_i,
      time = risk_times,
      n_risk = vapply(risk_times, function(t) sum(d$follow_time_days >= t), integer(1))
    )
  })
)

expected_risk <- data.frame(
  strategy = rep(c("A", "B"), each = 5),
  time = rep(risk_times, 2),
  n_risk = c(1308, 69, 60, 57, 55, 1308, 1042, 955, 900, 853)
)
if (!isTRUE(all.equal(risk, expected_risk, check.attributes = FALSE))) {
  print(risk)
  stop("Number-at-risk counts no longer match the validated final figure.")
}

strategy_labels <- c(
  A = "Strategy A (albumin 6\u2013<48 h)",
  B = "Strategy B (no albumin 6\u2013<48 h)"
)
strategy_colours <- c(A = "#C93A2B", B = "#2C7FB8")
strategy_linetypes <- c(A = "solid", B = "22")

unicode_minus <- function(x) gsub("-", "\u2212", x, fixed = TRUE)
rd_label <- sprintf(
  "RD: %.1f%% (%s to %.1f%%)",
  100 * result_28$rd,
  unicode_minus(sprintf("%.1f%%", 100 * result_28$rd_low)),
  100 * result_28$rd_high
)
hr_label <- sprintf(
  "HR: %.3f (%.3f to %.3f)",
  result_28$hr,
  result_28$hr_low,
  result_28$hr_high
)

base_family <- "Arial"

p_curve <- ggplot(km, aes(x = time, y = 100 * cif, colour = strategy, linetype = strategy)) +
  geom_step(linewidth = 1.05, direction = "hv") +
  annotate(
    "label", x = 15.6, y = 8.6,
    label = paste(rd_label, hr_label, sep = "\n"),
    hjust = 0, vjust = 0, size = 3.6, family = base_family,
    linewidth = 0.25, label.padding = unit(0.22, "lines"),
    colour = "#111111", fill = "white"
  ) +
  scale_colour_manual(values = strategy_colours, labels = strategy_labels, name = NULL) +
  scale_linetype_manual(values = strategy_linetypes, labels = strategy_labels, name = NULL) +
  scale_x_continuous(breaks = risk_times, expand = expansion(mult = c(0, 0.005))) +
  scale_y_continuous(
    limits = c(0, 36), breaks = c(0, 10, 20, 30),
    expand = expansion(mult = c(0, 0.015))
  ) +
  coord_cartesian(xlim = c(0, 28.4)) +
  labs(x = NULL, y = "Cumulative mortality (%)") +
  guides(
    colour = guide_legend(override.aes = list(linewidth = 1.2, linetype = c("solid", "22"))),
    linetype = "none"
  ) +
  theme_minimal(base_size = 11, base_family = base_family) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#E1E4E8", linewidth = 0.45),
    axis.line = element_line(colour = "#111111", linewidth = 0.55),
    axis.ticks = element_line(colour = "#111111", linewidth = 0.45),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 12, margin = margin(r = 9)),
    axis.text.y = element_text(size = 10.5, colour = "#111111"),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.94), colour = NA),
    legend.key.width = unit(1.25, "cm"),
    legend.text = element_text(size = 10.5, colour = "#111111"),
    plot.margin = margin(t = 7, r = 8, b = 0, l = 8)
  )

risk$strategy <- factor(risk$strategy, levels = c("A", "B"))
risk$y <- ifelse(risk$strategy == "A", 1.35, 0.48)

p_risk <- ggplot(risk, aes(x = time, y = y, label = comma(n_risk), colour = strategy)) +
  geom_text(size = 3.7, family = base_family, show.legend = FALSE) +
  annotate("text", x = 14, y = 2.28, label = "Number at risk", fontface = "bold", size = 3.9, family = base_family) +
  annotate("text", x = -1.55, y = 1.35, label = "A", fontface = "bold", size = 3.7, family = base_family, colour = strategy_colours[["A"]]) +
  annotate("text", x = -1.55, y = 0.48, label = "B", fontface = "bold", size = 3.7, family = base_family, colour = strategy_colours[["B"]]) +
  scale_colour_manual(values = strategy_colours) +
  scale_x_continuous(breaks = risk_times, expand = expansion(mult = c(0, 0.005))) +
  scale_y_continuous(limits = c(0.1, 2.45), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 28.4), clip = "off") +
  labs(x = "Days from 6-hour landmark", y = NULL) +
  theme_classic(base_size = 11, base_family = base_family) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line(colour = "#111111", linewidth = 0.55),
    axis.ticks.x = element_line(colour = "#111111", linewidth = 0.45),
    axis.text.x = element_text(size = 10.5, colour = "#111111"),
    axis.title.x = element_text(size = 12, margin = margin(t = 7)),
    plot.margin = margin(t = 0, r = 8, b = 7, l = 8)
  )

final_plot <- p_curve / p_risk + plot_layout(heights = c(3.45, 1.15))

png_path <- file.path(figure_dir, "Figure2.png")
pdf_path <- file.path(figure_dir, "Figure2.pdf")
tiff_path <- file.path(figure_dir, "Figure2.tiff")

ggsave(png_path, final_plot, width = 7.0, height = 5.05, units = "in", dpi = 600, bg = "white", device = grDevices::png)
ggsave(tiff_path, final_plot, width = 7.0, height = 5.05, units = "in", dpi = 600, bg = "white", device = grDevices::tiff, compression = "lzw")
ggsave(pdf_path, final_plot, width = 7.0, height = 5.05, units = "in", device = cairo_pdf, bg = "white")

write.csv(risk, file.path(figure_dir, "Figure2_number_at_risk.csv"), row.names = FALSE)

cat("Created:\n", png_path, "\n", pdf_path, "\n", tiff_path, "\n", sep = "")
cat(rd_label, "\n", hr_label, "\n", sep = "")
