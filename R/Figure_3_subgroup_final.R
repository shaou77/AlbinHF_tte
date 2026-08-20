options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

root <- "E:/workstation/tte/v2"
output_root <- file.path(root, "output")
figure_dir <- file.path(output_root, "figures", "main", "final")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

pooled <- read.csv(
  file.path(output_root, "tables", "table2_pooled_results_raw_v2.csv"),
  check.names = FALSE
)
sub <- read.csv(
  file.path(output_root, "tables", "table_s_subgroup_results_raw_v2.csv"),
  check.names = FALSE
)

overall <- pooled[pooled$endpoint == 28, , drop = FALSE]
stopifnot(nrow(overall) == 1L, nrow(sub) == 4L, all(sub$feasible))

rows <- data.frame(
  type = c("estimate", "section", "estimate", "estimate", "section", "estimate", "estimate"),
  label = c(
    "Overall (main analysis)",
    "Baseline serum albumin",
    "Albumin < 2.5 g/dL",
    "Albumin \u2265 2.5 g/dL",
    "Vasopressor use within first 6 h",
    "Vasopressor use",
    "No vasopressor"
  ),
  y = c(6, 5, 4, 3, 2, 1, 0),
  N = c(overall$n_patients, NA, sub$n[1], sub$n[2], NA, sub$n[3], sub$n[4]),
  n = c(overall$exposed_6_48, NA, sub$n_exposed[1], sub$n_exposed[2], NA, sub$n_exposed[3], sub$n_exposed[4]),
  RD = c(overall$rd, NA, sub$RD[1], sub$RD[2], NA, sub$RD[3], sub$RD[4]),
  RD_lo = c(overall$rd_low, NA, sub$RD_lo[1], sub$RD_lo[2], NA, sub$RD_lo[3], sub$RD_lo[4]),
  RD_hi = c(overall$rd_high, NA, sub$RD_hi[1], sub$RD_hi[2], NA, sub$RD_hi[3], sub$RD_hi[4]),
  HR = c(overall$hr, NA, sub$HR[1], sub$HR[2], NA, sub$HR[3], sub$HR[4]),
  HR_lo = c(overall$hr_low, NA, sub$HR_lo[1], sub$HR_lo[2], NA, sub$HR_lo[3], sub$HR_lo[4]),
  HR_hi = c(overall$hr_high, NA, sub$HR_hi[1], sub$HR_hi[2], NA, sub$HR_hi[3], sub$HR_hi[4])
)

est <- rows[rows$type == "estimate", , drop = FALSE]
est$is_overall <- est$label == "Overall (main analysis)"
stopifnot(
  isTRUE(all.equal(est$N, c(1308, 395, 913, 313, 995), check.attributes = FALSE)),
  isTRUE(all.equal(est$n, c(80, 40, 40, 43, 37), check.attributes = FALSE)),
  abs(est$RD[1] - 0.0221339811103647) < 1e-12,
  abs(est$HR[1] - 1.036819) < 1e-6
)

unicode_minus <- function(x) gsub("-", "\u2212", x, fixed = TRUE)
fmt_rd <- function(est, low, high) {
  unicode_minus(sprintf("%.1f%% (%.1f%% to %.1f%%)", 100 * est, 100 * low, 100 * high))
}
fmt_hr <- function(est, low, high) sprintf("%.3f (%.3f to %.3f)", est, low, high)

est$rd_text <- mapply(fmt_rd, est$RD, est$RD_lo, est$RD_hi)
est$hr_text <- mapply(fmt_hr, est$HR, est$HR_lo, est$HR_hi)
est$n_text <- format(est$N, big.mark = ",", scientific = FALSE)
est$exposed_text <- format(est$n, big.mark = ",", scientific = FALSE)

base_family <- "Arial"
navy <- "#243B53"
red <- "#C93A2B"
blue <- "#2C7FB8"
row_guides <- c(6, 4, 3, 1, 0)
section_guides <- c(4.55, 1.55)

base_void <- theme_void(base_family = base_family) +
  theme(plot.margin = margin(4, 2, 15, 2))

p_left <- ggplot() +
  geom_hline(yintercept = row_guides, colour = "#EDF0F2", linewidth = 0.4) +
  geom_hline(yintercept = section_guides, colour = "#E5E7EB", linewidth = 0.5) +
  annotate("text", x = 0.00, y = 7.0, label = "Subgroup", hjust = 0, fontface = "bold", size = 3.8, family = base_family) +
  annotate("text", x = 1.22, y = 7.0, label = "N", hjust = 0.5, fontface = "bold", size = 3.8, family = base_family) +
  annotate("text", x = 1.52, y = 7.30, label = "Strategy A", hjust = 0.5, fontface = "bold", size = 3.8, family = base_family) +
  annotate("text", x = 1.52, y = 6.95, label = "n", hjust = 0.5, fontface = "bold", size = 3.6, family = base_family) +
  geom_text(
    data = rows[rows$type == "section", ],
    aes(x = 0, y = y, label = label),
    hjust = 0, fontface = "bold", size = 3.75, family = base_family
  ) +
  geom_text(
    data = est[!est$is_overall, ],
    aes(x = 0.02, y = y, label = label),
    hjust = 0, size = 3.30, family = base_family
  ) +
  geom_text(
    data = est[est$is_overall, ],
    aes(x = 0.02, y = y, label = label),
    hjust = 0, fontface = "bold", size = 3.30, family = base_family
  ) +
  geom_text(data = est, aes(x = 1.22, y = y, label = n_text), size = 3.35, family = base_family) +
  geom_text(data = est, aes(x = 1.52, y = y, label = exposed_text), size = 3.35, family = base_family) +
  coord_cartesian(xlim = c(0, 1.66), ylim = c(-0.6, 7.55), clip = "off") +
  base_void

p_forest <- ggplot(est, aes(y = y)) +
  geom_hline(yintercept = row_guides, colour = "#EDF0F2", linewidth = 0.4) +
  geom_vline(xintercept = seq(-20, 20, 10), colour = "#E1E4E8", linewidth = 0.42) +
  geom_vline(xintercept = 0, colour = "#AAB2BD", linewidth = 0.6, linetype = "22") +
  geom_segment(aes(x = 100 * RD_lo, xend = 100 * RD_hi, yend = y), colour = navy, linewidth = 0.75) +
  geom_segment(aes(x = 100 * RD_lo, xend = 100 * RD_lo, y = y - 0.11, yend = y + 0.11), colour = navy, linewidth = 0.75) +
  geom_segment(aes(x = 100 * RD_hi, xend = 100 * RD_hi, y = y - 0.11, yend = y + 0.11), colour = navy, linewidth = 0.75) +
  geom_point(data = est[!est$is_overall, ], aes(x = 100 * RD), shape = 15, size = 2.55, colour = blue) +
  geom_point(data = est[est$is_overall, ], aes(x = 100 * RD), shape = 18, size = 3.8, colour = red) +
  scale_x_continuous(limits = c(-25, 30), breaks = seq(-20, 20, 10), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.6, 7.55), expand = c(0, 0)) +
  labs(x = "28-day risk difference (%)", y = NULL) +
  theme_classic(base_size = 10.5, base_family = base_family) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(colour = "#111111", linewidth = 0.55),
    axis.ticks.x = element_line(colour = "#111111", linewidth = 0.45),
    axis.text.x = element_text(size = 9.5, colour = "#111111"),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 6)),
    plot.margin = margin(4, 6, 15, 6)
  )

p_rd <- ggplot() +
  geom_hline(yintercept = row_guides, colour = "#EDF0F2", linewidth = 0.4) +
  annotate("text", x = 0, y = 7.0, label = "RD (95% CI)", hjust = 0, fontface = "bold", size = 3.8, family = base_family) +
  geom_text(data = est, aes(x = 0, y = y, label = rd_text), hjust = 0, size = 3.35, family = base_family) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.6, 7.55), clip = "off") +
  base_void

p_hr <- ggplot() +
  geom_hline(yintercept = row_guides, colour = "#EDF0F2", linewidth = 0.4) +
  annotate("text", x = 0, y = 7.0, label = "HR (95% CI)", hjust = 0, fontface = "bold", size = 3.8, family = base_family) +
  geom_text(data = est, aes(x = 0, y = y, label = hr_text), hjust = 0, size = 3.35, family = base_family) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.6, 7.55), clip = "off") +
  base_void

final_plot <- p_left + p_forest + p_rd + p_hr +
  plot_layout(widths = c(3.45, 2.85, 2.35, 2.15))

png_path <- file.path(figure_dir, "Figure3.png")
pdf_path <- file.path(figure_dir, "Figure3.pdf")
tiff_path <- file.path(figure_dir, "Figure3.tiff")

ggsave(png_path, final_plot, width = 8.4, height = 3.65, units = "in", dpi = 600, bg = "white", device = grDevices::png)
ggsave(tiff_path, final_plot, width = 8.4, height = 3.65, units = "in", dpi = 600, bg = "white", device = grDevices::tiff, compression = "lzw")
ggsave(pdf_path, final_plot, width = 8.4, height = 3.65, units = "in", device = cairo_pdf, bg = "white")

write.csv(est, file.path(figure_dir, "Figure3_display_values.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("Created:\n", png_path, "\n", pdf_path, "\n", tiff_path, "\n", sep = "")
print(est[, c("label", "N", "n", "rd_text", "hr_text")], row.names = FALSE)
