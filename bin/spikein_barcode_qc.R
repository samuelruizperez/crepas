#!/usr/bin/env Rscript

# ===============================================================================
# spikein_barcode_qc.R
#
# Created for the grothlab/crepas pipeline by:
#     - Samuel Ruiz-Perez <samper@cancer.dk>
#     https://github.com/grothlab/crepas/
#
# Description:
#     Read one or more spike-in count tables and matched uniq_totals values.
#     Write:
#       1) a long table with sample_id per barcode row
#       2) a summary table with A+B counts per target and per-sample
#          on-target normalization
# ===============================================================================

options(show.error.locations = TRUE)

required.libs <- c("readr", "dplyr", "argparse", "ComplexHeatmap", "circlize")

for (lib in required.libs) {
  suppressPackageStartupMessages({
    if (!require(lib, character.only = TRUE)) {
      stop("Missing required R package: ", lib)
    }
  })
}

# ===============================================================================
# Argument parsing
# ===============================================================================

parser <- argparse::ArgumentParser(
  description = "Generate spike-in QC long and summary tables"
)

parser$add_argument(
  "-c", "--count_tables",
  action = "store",
  type = "character",
  nargs = "+",
  required = TRUE,
  help = "One or more count tables with columns: barcode_target, barcode_id, barcode_sequence, R1_count, R2_count"
)

parser$add_argument(
  "-u", "--uniq_totals",
  action = "store",
  type = "double",
  nargs = "+",
  required = TRUE,
  help = "One or more uniq total values, provided in the same order as count_tables"
)

parser$add_argument(
  "-p", "--prefix",
  action = "store",
  type = "character",
  default = "spikein_barcodes_qc",
  help = "Prefix for output files [default: spikein_barcodes_qc]"
)

parser$add_argument(
  "-o", "--outdir",
  action = "store",
  type = "character",
  default = "./",
  help = "Output directory [default: ./]"
)

opt <- parser$parse_args()
opt_count_tables <- opt$count_tables
opt_uniq_totals <- opt$uniq_totals
opt_prefix <- opt$prefix
opt_outdir <- opt$outdir

if (length(opt_count_tables) != length(opt_uniq_totals)) {
  stop(
    "count_tables and uniq_totals must have the same length. Got ",
    length(opt_count_tables), " count_tables and ", length(opt_uniq_totals), " uniq_totals."
  )
}

# ===============================================================================
# Helpers
# ===============================================================================

required_cols <- c(
  "barcode_target",
  "barcode_id",
  "barcode_sequence",
  "R1_count",
  "R2_count"
)

sample_id_from_path <- function(path) {
  tools::file_path_sans_ext(basename(path))
}

read_count_table <- function(path, sample_id, uniq_total) {
  if (!file.exists(path)) {
    stop("Count table not found: ", path)
  }

  df <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Table '", path, "' is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df <- df %>%
    dplyr::select(dplyr::all_of(required_cols)) %>%
    dplyr::mutate(
      barcode_target = as.character(barcode_target),
      barcode_id = as.character(barcode_id),
      barcode_sequence = as.character(barcode_sequence),
      R1_count = as.numeric(R1_count),
      R2_count = as.numeric(R2_count),
      sample_id = sample_id,
      uniq_total = uniq_total
    )

  bad_rows <- which(is.na(df$R1_count) | is.na(df$R2_count))
  if (length(bad_rows) > 0) {
    stop(
      "Table '", path, "' has non-numeric values in R1_count/R2_count at rows: ",
      paste(head(bad_rows, 10), collapse = ", "),
      if (length(bad_rows) > 10) " ..." else ""
    )
  }

  df
}

# ===============================================================================
# Build long table
# ===============================================================================

message("[", Sys.time(), "] Reading ", length(opt_count_tables), " count table(s)...")

table_list <- lapply(seq_along(opt_count_tables), function(i) {
  read_count_table(
    path = opt_count_tables[[i]],
    sample_id = sample_id_from_path(opt_count_tables[[i]]),
    uniq_total = opt_uniq_totals[[i]]
  )
})

long_df <- dplyr::bind_rows(table_list) %>%
  dplyr::select(sample_id, uniq_total, dplyr::all_of(required_cols))

# ===============================================================================
# Build summary table
# ===============================================================================

summary_df <- long_df %>%
  dplyr::mutate(total_count = R1_count + R2_count) %>%
  dplyr::group_by(sample_id, barcode_target) %>%
  dplyr::summarise(
    barcode_ab_count = sum(total_count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(
    on_target_normalization = barcode_ab_count / sum(barcode_ab_count, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::rename(target = barcode_target) %>%
  dplyr::select(sample_id, target, barcode_ab_count, on_target_normalization)

# ===============================================================================
# Build heatmap
# ===============================================================================

recovery_col_label <- "Total barcode /\nuniq reads"

sample_recovery_df <- long_df %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(
    target = recovery_col_label,
    recovery_ratio = sum(R1_count + R2_count, na.rm = TRUE) / dplyr::first(uniq_total),
    .groups = "drop"
  )

target_levels <- c(sort(unique(summary_df$target)), recovery_col_label)
sample_levels <- summary_df$sample_id %>% unique()

target_heatmap_df <- summary_df %>%
  dplyr::select(sample_id, target, on_target_normalization) %>%
  dplyr::mutate(
    sample_id = sub("\\..*$", "", sample_id)
  )

recovery_heatmap_df <- sample_recovery_df %>%
  dplyr::mutate(
    sample_id = sub("\\..*$", "", sample_id)
  )

sample_levels_clean <- sub("\\..*$", "", sample_levels)
target_only_levels <- sort(unique(summary_df$target))

target_matrix <- with(
  target_heatmap_df,
  tapply(
    on_target_normalization,
    list(
      factor(sample_id, levels = sample_levels_clean),
      factor(target, levels = target_only_levels)
    ),
    identity
  )
)

recovery_by_sample <- setNames(recovery_heatmap_df$recovery_ratio, recovery_heatmap_df$sample_id)
recovery_matrix <- matrix(
  recovery_by_sample[sample_levels_clean],
  ncol = 1,
  dimnames = list(sample_levels_clean, recovery_col_label)
)

target_col_fun <- circlize::colorRamp2(c(0, 1), c("#f8c22dff", "#54a7d3ff"))
max_recovery <- max(recovery_matrix, na.rm = TRUE)
if (!is.finite(max_recovery) || max_recovery <= 0) {
  max_recovery <- 1
}
recovery_col_fun <- circlize::colorRamp2(c(0, max_recovery), c("#fff5f0", "#cb181d"))

target_ht <- ComplexHeatmap::Heatmap(
  target_matrix,
  name = "On-target recovery",
  col = target_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 9),
  rect_gp = grid::gpar(col = "white", lwd = 0.5),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid::grid.text(sprintf("%.3f", target_matrix[i, j]), x, y, gp = grid::gpar(fontsize = 7))
  },
  column_names_rot = 45,
  width = grid::unit(ncol(target_matrix) * 10.0, "mm")
)

recovery_ht <- ComplexHeatmap::Heatmap(
  recovery_matrix,
  name = "Total barcode /\nuniq reads",
  col = recovery_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = grid::gpar(fontsize = 11),
  rect_gp = grid::gpar(col = "white", lwd = 0.5),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid::grid.text(sprintf("%.3f", recovery_matrix[i, j]), x, y, gp = grid::gpar(fontsize = 7))
  },
  width = grid::unit(2.4, "cm")
)

heatmap_obj <- target_ht + recovery_ht

# ===============================================================================
# Write outputs
# ===============================================================================

if (!dir.exists(opt_outdir)) {
  dir.create(opt_outdir, recursive = TRUE)
}

long_file <- file.path(opt_outdir, paste0(opt_prefix, ".long.tsv"))
summary_file <- file.path(opt_outdir, paste0(opt_prefix, ".summary.tsv"))
heatmap_pdf <- file.path(opt_outdir, paste0(opt_prefix, ".heatmap.pdf"))
heatmap_png <- file.path(opt_outdir, paste0(opt_prefix, ".heatmap.png"))

readr::write_tsv(long_df, long_file, col_names = TRUE)
readr::write_tsv(summary_df, summary_file, col_names = TRUE)

grDevices::pdf(heatmap_pdf, width = 18, height = 6)
ComplexHeatmap::draw(heatmap_obj, heatmap_legend_side = "right", merge_legends = FALSE)
grDevices::dev.off()

grDevices::png(heatmap_png, width = 18, height = 6, units = "in", res = 300)
ComplexHeatmap::draw(heatmap_obj, heatmap_legend_side = "right", merge_legends = FALSE)
grDevices::dev.off()

message("[", Sys.time(), "] Wrote long table: ", long_file)
message("[", Sys.time(), "] Long table rows: ", nrow(long_df))
message("[", Sys.time(), "] Wrote summary table: ", summary_file)
message("[", Sys.time(), "] Summary table rows: ", nrow(summary_df))
message("[", Sys.time(), "] Wrote heatmap PDF: ", heatmap_pdf)
message("[", Sys.time(), "] Wrote heatmap PNG: ", heatmap_png)
