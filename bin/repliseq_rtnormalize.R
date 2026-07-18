#!/usr/bin/env Rscript

# ===============================================================================
# repliseq_rtnormalize.R
#
# Written for the grothlab/crepas pipeline by:
#     - Samuel Ruiz-Pérez <samper@cancer.dk>
#     https://github.com/grothlab/crepas/
#
# Description:
#     Calculates and normalizes a replication-timing (E/L ratio) track for one
#     Repli-seq sample from a multi-replicate, multi-phase raw read-count table
#     (chr, start, end, <phase_brep> columns): CPM normalization, early/late
#     replicate combination (paired per biological replicate, pooled across
#     replicates, or auto-chosen via a correlation heuristic), optional
#     quantile normalization between replicates, and loess or rolling-window
#     smoothing.
# ===============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(preprocessCore)
  library(zoo)
})

# ===============================================================================
# Argument parsing
# ===============================================================================

parser <- ArgumentParser()

parser$add_argument("-c", "--counts", action = "store",
                    type = "character",
                    required = TRUE,
                    help = "Raw read-count table (chr, start, end, then one column per replicate), e.g. deepTools multiBamSummary --outRawCounts output; its one header line is always skipped [required]")

parser$add_argument("--phases", action = "store",
                    type = "character",
                    nargs = "+",
                    required = TRUE,
                    help = "Phase ('early' or 'late') of each count column, in the same order as --counts [required]")

parser$add_argument("--breps", action = "store",
                    type = "character",
                    nargs = "+",
                    required = TRUE,
                    help = "Biological replicate identifier of each count column, in the same order as --counts. Used to match early/late replicate pairs for --method paired [required]")

parser$add_argument("-m", "--method", action = "store",
                    default = "auto",
                    choices = c("auto", "paired", "pooled"),
                    help = "How to combine replicates: 'paired' computes RT per matched early/late biological-replicate pair and averages; 'pooled' averages replicates within each phase first, then computes one RT ratio; 'auto' chooses between the two via a correlation heuristic [default: %(default)s]")

parser$add_argument("-n", "--normalization", action = "store",
                    default = "cpm",
                    choices = c("cpm", "qn"),
                    help = "'cpm' only CPM-normalizes each replicate; 'qn' additionally quantile-normalizes between replicates (on top of CPM), within the paired RT columns or within each phase's CPM columns depending on --method [default: %(default)s]")

parser$add_argument("-s", "--smooth", action = "store",
                    default = "loess",
                    choices = c("loess", "roll", "none"),
                    help = "Smoothing strategy applied to the raw RT track [default: %(default)s]")

parser$add_argument("--loess_span", action = "store",
                    default = 0.002,
                    type = "double",
                    help = "Loess smoothing span, as a fraction of points per chromosome [default: %(default)s]")

parser$add_argument("--roll_k", action = "store",
                    default = 5,
                    type = "integer",
                    help = "Rolling-mean window size, in bins [default: %(default)s]")

parser$add_argument("--corr", action = "store",
                    default = "spearman",
                    choices = c("spearman", "pearson"),
                    help = "Correlation method used by the 'auto' replicate-combination heuristic [default: %(default)s]")

parser$add_argument("--pairing_threshold", action = "store",
                    default = 0.02,
                    type = "double",
                    help = "Under --method auto, 'paired' is chosen when (matched-pair avg. correlation - mismatched-pair avg. correlation) is at least this [default: %(default)s]")

parser$add_argument("--pseudocount", action = "store",
                    default = 1,
                    type = "double",
                    help = "Pseudocount added to CPM values before taking log2(early/late) or log10 (for correlation) [default: %(default)s]")

parser$add_argument("--min_points_for_cor", action = "store",
                    default = 50,
                    type = "integer",
                    help = "Minimum number of finite paired observations required to compute a correlation under --method auto; pairs with fewer are treated as NA [default: %(default)s]")

parser$add_argument("-o", "--outdir", action = "store",
                    default = ".",
                    type = "character",
                    help = "Path to output directory [default: %(default)s]")

parser$add_argument("-p", "--prefix", action = "store",
                    default = "",
                    type = "character",
                    help = "Prefix prepended to every output filename, as '<prefix>.<suffix>' [default: none, i.e. just '<suffix>']")

opt <- parser$parse_args()

if (length(opt$phases) != length(opt$breps)) {
  stop("--phases and --breps must have the same length")
}
if (!all(opt$phases %in% c("early", "late"))) {
  stop("--phases values must be 'early' or 'late'")
}

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

out_path <- function(suffix) {
  filename <- if (nzchar(opt$prefix)) paste(opt$prefix, suffix, sep = ".") else suffix
  file.path(opt$outdir, filename)
}

# ===============================================================================
# Helpers (ported from RepliSeq_process.R, generalized to N replicates)
# ===============================================================================

cpm_vec <- function(x) {
  total <- sum(x, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) {
    return(rep(NA_real_, length(x)))
  }
  (x / total) * 1e6
}

safe_cor <- function(x, y, method, min_points) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_points) {
    return(NA_real_)
  }
  cor(x[ok], y[ok], method = method)
}

qn_matrix <- function(mat) {
  mat <- as.matrix(mat)
  if (ncol(mat) < 2) {
    return(mat)
  }
  ok_rows <- apply(mat, 1, function(x) all(is.finite(x)))
  out <- mat
  if (sum(ok_rows) > 1) {
    out[ok_rows, ] <- preprocessCore::normalize.quantiles(mat[ok_rows, , drop = FALSE])
  }
  out
}

smooth_loess_chr <- function(y, x, span) {
  ok <- is.finite(y) & is.finite(x)
  out <- rep(NA_real_, length(y))
  if (sum(ok) < 4) {
    return(out)
  }
  fit <- loess(y[ok] ~ x[ok], span = span, degree = 1, family = "symmetric",
              control = loess.control(surface = "direct"))
  out[ok] <- predict(fit, x[ok])
  out
}

smooth_roll_chr <- function(y, k) {
  zoo::rollapply(y, width = k, FUN = mean, fill = NA_real_, align = "center", na.rm = TRUE)
}

write_bedgraph <- function(dt, value_col, path) {
  out <- dt[is.finite(get(value_col)), .(chr, start, end, value = get(value_col))]
  setorder(out, chr, start)
  fwrite(out, file = path, sep = "\t", col.names = FALSE)
}

# ===============================================================================
# Read the raw counts table
# ===============================================================================

message("[", Sys.time(), "] Reading counts table: ", opt$counts)
col_ids <- paste(opt$phases, opt$breps, sep = "_")
if (anyDuplicated(col_ids)) {
  stop("Duplicate (phase, brep) combinations in --phases/--breps: ", paste(col_ids[duplicated(col_ids)], collapse = ", "))
}

# --counts is deepTools multiBamSummary's --outRawCounts output, which always has one
# header line (e.g. "#'chr'\t'start'\t'end'\t'label1'..."); skip it unconditionally
# rather than trusting its quoting, and use --phases/--breps for column identity instead.
dt <- fread(opt$counts, header = FALSE, sep = "\t", skip = 1, col.names = c("chr", "start", "end", col_ids))
dt[, start := as.integer(start)]
dt[, end := as.integer(end)]
setorder(dt, chr, start)
dt[, mid := (start + end) / 2]

early_ids <- col_ids[opt$phases == "early"]
late_ids <- col_ids[opt$phases == "late"]
early_breps <- opt$breps[opt$phases == "early"]
late_breps <- opt$breps[opt$phases == "late"]

if (length(early_ids) < 1 || length(late_ids) < 1) {
  stop("At least one 'early' and one 'late' replicate are required")
}

# ===============================================================================
# CPM normalization per replicate
# ===============================================================================

for (id in col_ids) {
  dt[, (paste0(id, "_cpm")) := cpm_vec(get(id))]
}

pc <- opt$pseudocount

# ===============================================================================
# Decide how to combine replicates: paired, pooled, or auto (correlation heuristic)
# ===============================================================================

matched_breps <- intersect(early_breps, late_breps)

chosen_method <- opt$method
paired_avg <- NA_real_
crossed_avg <- NA_real_
pairing_delta <- NA_real_

if (opt$method %in% c("auto", "paired")) {
  if (length(matched_breps) == 0) {
    if (opt$method == "paired") {
      warning("--method paired requested, but no biological replicate is shared between early and late phases (early breps: ",
              paste(early_breps, collapse = ", "), "; late breps: ", paste(late_breps, collapse = ", "),
              "); falling back to 'pooled'.")
    }
    chosen_method <- "pooled"
  } else {
    for (id in col_ids) {
      dt[, (paste0(id, "_logcpm")) := log10(get(paste0(id, "_cpm")) + pc)]
    }

    matched_cors <- c()
    crossed_cors <- c()
    for (e in seq_along(early_ids)) {
      for (l in seq_along(late_ids)) {
        r <- safe_cor(
          dt[[paste0(early_ids[e], "_logcpm")]],
          dt[[paste0(late_ids[l], "_logcpm")]],
          opt$corr, opt$min_points_for_cor
        )
        if (early_breps[e] == late_breps[l]) {
          matched_cors <- c(matched_cors, r)
        } else {
          crossed_cors <- c(crossed_cors, r)
        }
      }
    }
    paired_avg <- mean(matched_cors, na.rm = TRUE)
    crossed_avg <- if (length(crossed_cors) > 0) mean(crossed_cors, na.rm = TRUE) else NA_real_
    pairing_delta <- paired_avg - crossed_avg

    if (opt$method == "auto") {
      chosen_method <- if (is.finite(pairing_delta) && pairing_delta >= opt$pairing_threshold) "paired" else "pooled"
    }
  }
}

message("[", Sys.time(), "] Replicate-combination method: ", chosen_method)

# ===============================================================================
# Compute the raw RT ratio
# ===============================================================================

if (chosen_method == "paired") {
  rt_cols <- c()
  for (brep in matched_breps) {
    e_id <- early_ids[early_breps == brep][1]
    l_id <- late_ids[late_breps == brep][1]
    rt_col <- paste0("RT_", brep)
    dt[, (rt_col) := log2((get(paste0(e_id, "_cpm")) + pc) / (get(paste0(l_id, "_cpm")) + pc))]
    rt_cols <- c(rt_cols, rt_col)
  }

  if (opt$normalization == "qn" && length(rt_cols) > 1) {
    dt[, (rt_cols) := as.data.table(qn_matrix(dt[, ..rt_cols]))]
  }

  dt[, RT_raw := rowMeans(dt[, ..rt_cols], na.rm = TRUE)]
} else {
  early_cpm_cols <- paste0(early_ids, "_cpm")
  late_cpm_cols <- paste0(late_ids, "_cpm")
  early_mat <- dt[, ..early_cpm_cols]
  late_mat <- dt[, ..late_cpm_cols]

  if (opt$normalization == "qn") {
    early_mat <- qn_matrix(early_mat)
    late_mat <- qn_matrix(late_mat)
  }

  dt[, Early_agg := rowMeans(early_mat, na.rm = TRUE)]
  dt[, Late_agg := rowMeans(late_mat, na.rm = TRUE)]
  dt[, RT_raw := log2((Early_agg + pc) / (Late_agg + pc))]
}

# ===============================================================================
# Smoothing
# ===============================================================================

dt[, RT_smooth := RT_raw]
if (opt$smooth == "loess") {
  message("[", Sys.time(), "] Loess smoothing (span = ", opt$loess_span, ")...")
  dt[, RT_smooth := smooth_loess_chr(RT_raw, mid, span = opt$loess_span), by = chr]
} else if (opt$smooth == "roll") {
  message("[", Sys.time(), "] Rolling-mean smoothing (k = ", opt$roll_k, " bins)...")
  dt[, RT_smooth := smooth_roll_chr(RT_raw, k = opt$roll_k), by = chr]
}

# ===============================================================================
# Write outputs
# ===============================================================================

write_bedgraph(dt, "RT_raw", out_path("RT.raw.bedGraph"))
write_bedgraph(dt, "RT_smooth", out_path("RT.smooth.bedGraph"))

qc_lines <- c(
  paste0("counts:\t", opt$counts),
  paste0("early_replicates:\t", paste(early_ids, collapse = ", ")),
  paste0("late_replicates:\t", paste(late_ids, collapse = ", ")),
  paste0("requested_method:\t", opt$method),
  paste0("chosen_method:\t", chosen_method),
  paste0("corr_method:\t", opt$corr),
  paste0("matched_replicate_pairs_avg_corr:\t", paired_avg),
  paste0("mismatched_replicate_pairs_avg_corr:\t", crossed_avg),
  paste0("pairing_delta:\t", pairing_delta),
  paste0("pairing_threshold:\t", opt$pairing_threshold),
  paste0("normalization:\t", opt$normalization),
  paste0("pseudocount:\t", opt$pseudocount),
  paste0("smoothing:\t", opt$smooth),
  if (opt$smooth == "loess") paste0("loess_span:\t", opt$loess_span) else "",
  if (opt$smooth == "roll") paste0("roll_k:\t", opt$roll_k) else "",
  paste0("n_bins:\t", nrow(dt))
)
qc_lines <- qc_lines[nzchar(qc_lines)]
writeLines(qc_lines, con = out_path("qc.txt"))

message("[", Sys.time(), "] Done!")
