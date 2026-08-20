options(stringsAsFactors = FALSE)

if (!exists("dirs")) {
  dirs <- list(
    data = file.path("output", "data"),
    tables = file.path("output", "tables"),
    fig_main = file.path("output", "figures", "main"),
    fig_suppl = file.path("output", "figures", "supplement")
  )
}
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

csv_to_html_table <- function(dat) {
  header <- paste0("<tr>", paste0("<th>", html_escape(names(dat)), "</th>", collapse = ""), "</tr>")
  rows <- apply(dat, 1, function(row) {
    paste0("<tr>", paste0("<td>", html_escape(row), "</td>", collapse = ""), "</tr>")
  })
  paste0("<table>\n", header, "\n", paste(rows, collapse = "\n"), "\n</table>")
}

write_docx_table <- function(path, title, dat) {
  if (!requireNamespace("officer", quietly = TRUE)) return(FALSE)
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, title, style = "heading 1")
  doc <- officer::body_add_par(
    doc,
    "Generated from the revised v2 landmark analysis outputs.",
    style = "Normal"
  )
  doc <- officer::body_add_table(doc, value = dat, style = "table_template")
  print(doc, target = path)
  TRUE
}

table1_path <- file.path(dirs$tables, "table1_baseline_characteristics_v2.csv")
if (file.exists(table1_path)) {
  table1 <- read.csv(table1_path, check.names = FALSE)
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<title>v2 Table 1</title>",
    "<style>body{font-family:Arial,sans-serif;margin:24px;}table{border-collapse:collapse;font-size:12px;}th,td{border:1px solid #ccc;padding:4px 6px;}th{background:#f2f2f2;}</style>",
    "</head><body><h1>Table 1. Baseline Characteristics, v2 Landmark Cohort</h1>",
    csv_to_html_table(table1),
    "</body></html>"
  )
  writeLines(html, file.path(dirs$tables, "table1_gtsummary.html"), useBytes = TRUE)
  write_docx_table(file.path(dirs$tables, "table1_gtsummary.docx"), "Table 1. Baseline Characteristics", table1)
  write_docx_table(file.path(dirs$tables, "table1_gtsummary_with_smd.docx"), "Table 1. Baseline Characteristics with SMD", table1)
}

missing_path <- file.path(dirs$tables, "table_s_missing_summary.csv")
if (file.exists(missing_path)) {
  missing_summary <- read.csv(missing_path, check.names = FALSE)
  write_docx_table(file.path(dirs$tables, "table_s_missing_summary.docx"), "Supplementary Missingness Summary", missing_summary)
}

density_pdf <- file.path(dirs$fig_suppl, "fig_s_mice_density.pdf")
if (!file.exists(density_pdf)) {
  stop(
    "True MICE density diagnostic is missing. Run ",
    "E:/workstation/tte/Figure_S4_mice_density.R; refusing to create a ",
    "missing-rate placeholder in the Figure S4 output slot."
  )
}

copy_if_exists <- function(from, to) {
  if (file.exists(from)) file.copy(from, to, overwrite = TRUE)
}

copy_if_exists(
  file.path(dirs$fig_main, "fig_km_28d.png"),
  file.path(dirs$fig_main, "abstract_km_clean_axis.png")
)
copy_if_exists(
  file.path(dirs$fig_main, "fig_km_28d.png"),
  file.path(dirs$fig_main, "fig_km_28d_transparent.png")
)
copy_if_exists(
  file.path(dirs$fig_main, "fig_subgroup_forest.png"),
  file.path(dirs$fig_main, "abstract_subgroup_forest_clean.png")
)

cat("Legacy display outputs complete.\n")
