# ranking_analysis.R
# Layer 2 Global Ranking — Professor Kugiumtzis revised methodology.
#
# METHODOLOGY (per supervisor's instruction, email 19/6/2026):
#   - Uses Layer 1 paired t-test results (layer1_summary.csv) as input.
#   - Assigns TIED ranks per configuration: methods not significantly
#     different from the best (sig_worse = 0, including the reference)
#     share averaged positions (Spearman convention).
#   - Computes mean tied-rank per method across all configurations.
#   - Runs pairwise Wilcoxon signed-rank tests with Holm correction
#     on the per-config tied-rank vectors. No Friedman test.
#   - Produces CD diagram and mean MCC heatmap.
#
# INPUT:  layer1_results/layer1_summary.csv  (from layer1_analysis.m)
#         mcc_long.csv                        (from ranking_analysis.m, for heatmap)
# OUTPUT: layer2_results/
#           global_ranking_table.csv
#           dgp_ranking_table_S1.csv / S2 / S3 / S4 / S5
#           global_cd_diagram.pdf
#           dgp_cd_diagram_S1.pdf / S2 / S3 / S4 / S5
#           wilcoxon_results.txt
#           mcc_heatmap.pdf
#
# Required packages: scmamp, ggplot2, dplyr, reshape2
# -------------------------------------------------------------------------

library(scmamp)
library(ggplot2)
library(dplyr)
library(reshape2)

# ---- 0. PATHS ------------------------------------------------------------

results_dir <- "C:/Users/dimit/OneDrive/ΗΜΜΥ/Διπλωματική/Chapter 2 - methods software/simulations/results/task3/"
layer1_csv  <- file.path(results_dir, "layer1_results", "layer1_summary.csv")
mcc_csv     <- file.path(results_dir, "mcc_long.csv")
output_dir  <- file.path(results_dir, "layer2_results")

if (!dir.exists(output_dir)) dir.create(output_dir)

sink_file <- file.path(output_dir, "wilcoxon_results.txt")
if (file.exists(sink_file)) file.remove(sink_file)

cat("Loading Layer 1 results:", layer1_csv, "\n")
df_l1 <- read.csv(layer1_csv, stringsAsFactors = FALSE)
cat("Rows:", nrow(df_l1), "\n")
cat("Columns:", names(df_l1), "\n\n")

# ---- 1. COMPUTE TIED RANKS PER CONFIGURATION ----------------------------
# For each configuration:
#   - Methods with sig_worse = 0 (not significantly worse than best,
#     INCLUDING the reference itself) form the top cluster.
#   - Methods with sig_worse = 1 are ranked below.
#   - Methods with sig_worse = NaN (missing data) get rank NA.
#   - Ties within each group receive averaged (mid) ranks.
#
# Sort order within config: by mean_MCC descending.
# The tied-rank assignment:
#   1. Sort all methods by mean_MCC descending within config.
#   2. Assign integer positions 1..n_methods.
#   3. Methods not significantly worse from the best (sig_worse=0)
#      share the average of their integer positions.
#   4. Methods with sig_worse=1 keep their individual integer positions
#      (they are already separated from the top cluster).
#   Note: within the sig_worse=1 group, further ties are possible if
#   two methods have identical mean_MCC — these are also averaged.

assign_tied_ranks <- function(df_config) {
  df_sorted <- df_config[order(-df_config$mean_MCC, na.last = TRUE), ]
  n <- nrow(df_sorted)
  
  # Integer positions (NaN methods go to end)
  df_sorted$int_rank <- seq_len(n)
  df_sorted$tied_rank <- NA_real_
  
  # Find the cluster of non-significantly-worse methods (sig_worse = 0 or is_reference = 1)
  top_cluster_idx <- which(df_sorted$sig_worse == 0 | df_sorted$is_reference == 1)
  
  if (length(top_cluster_idx) > 0) {
    shared_rank <- mean(df_sorted$int_rank[top_cluster_idx])
    df_sorted$tied_rank[top_cluster_idx] <- shared_rank
  }
  
  # Remaining methods (sig_worse = 1): assign their integer positions
  # but also average ties within them if mean_MCC is identical
  remaining_idx <- which(df_sorted$sig_worse == 1)
  if (length(remaining_idx) > 0) {
    remaining_mccs <- df_sorted$mean_MCC[remaining_idx]
    unique_mccs <- unique(remaining_mccs)
    for (mcc_val in unique_mccs) {
      tied_idx <- remaining_idx[remaining_mccs == mcc_val]
      df_sorted$tied_rank[tied_idx] <- mean(df_sorted$int_rank[tied_idx])
    }
  }
  
  # NaN methods: leave as NA
  nan_idx <- which(is.na(df_sorted$sig_worse))
  # tied_rank already NA for these
  
  return(df_sorted[, c("method_label", "mean_MCC", "tied_rank", "sig_worse", "is_reference")])
}

# Apply to all configurations
configs      <- unique(df_l1$config)
dgp_families <- unique(df_l1$dgp_family)

# Build tied-rank matrix: rows = configs, cols = methods
# First get all method labels
all_methods <- unique(df_l1$method_label)

tied_rank_list <- list()
for (cfg in configs) {
  df_cfg <- df_l1[df_l1$config == cfg, ]
  ranked <- assign_tied_ranks(df_cfg)
  tied_rank_list[[cfg]] <- setNames(ranked$tied_rank, ranked$method_label)
}

# Convert to matrix
tied_mat_full <- do.call(rbind, lapply(configs, function(cfg) {
  row <- tied_rank_list[[cfg]]
  row[all_methods]  # ensure consistent column order
}))
rownames(tied_mat_full) <- configs
colnames(tied_mat_full) <- all_methods

cat("Tied rank matrix dimensions:", dim(tied_mat_full), "\n")
cat("Configurations:", nrow(tied_mat_full), "\n")
cat("Methods:", ncol(tied_mat_full), "\n\n")

# ---- 1b. HELPER: CD DIAGRAM DRIVEN BY THE ACTUAL WILCOXON+HOLM RESULT ---
# Draws a classic Demsar-style critical-difference diagram, but the
# connecting bars are built directly from p_adj_mat (the pairwise Wilcoxon
# signed-rank test with Holm correction from section 2b below), instead of
# letting plotCD() recompute its own Nemenyi-based grouping internally on
# the raw MCC matrix. Two methods are joined by a bar iff every pair within
# the group has p_adj >= alpha (maximal non-significant cliques in rank
# order, same convention as the original Demsar (2006) CD diagram).

plot_wilcoxon_cd <- function(avg_ranks, p_adj_mat, alpha = 0.05, title = "") {
  ord     <- order(avg_ranks)
  methods <- names(avg_ranks)[ord]
  ranks   <- avg_ranks[ord]
  n       <- length(methods)
  p_adj_mat <- p_adj_mat[methods, methods]
  
  # -- maximal cliques of mutually non-significant methods, in rank order --
  cliques <- list()
  for (i in 1:n) {
    j <- i
    while (j < n) {
      seg <- i:(j + 1)
      ok  <- TRUE
      for (a in seq_along(seg)[-length(seg)]) {
        for (b in (a + 1):length(seg)) {
          if (p_adj_mat[methods[seg[a]], methods[seg[b]]] < alpha) ok <- FALSE
        }
      }
      if (ok) j <- j + 1 else break
    }
    if (j > i) cliques[[length(cliques) + 1]] <- i:j
  }
  keep <- rep(TRUE, length(cliques))
  if (length(cliques) > 0) {
    for (a in seq_along(cliques)) for (b in seq_along(cliques)) {
      if (a != b && all(cliques[[a]] %in% cliques[[b]]) &&
          length(cliques[[a]]) < length(cliques[[b]])) keep[a] <- FALSE
    }
  }
  cliques <- cliques[keep]
  
  # -- layout: axis on top, best-ranked half labelled left, rest right ------
  rng   <- pretty(range(ranks))
  xlim  <- range(rng)
  half  <- ceiling(n / 2)
  left  <- 1:half
  right <- if (half < n) (half + 1):n else integer(0)
  
  # top margin holds two separate rows: the axis (numbers + ticks) and,
  # further out, the title on its own line, so the two never collide
  op <- par(mar = c(1, 9, 5, 9), xpd = NA)
  on.exit(par(op))
  plot(NA, xlim = xlim, ylim = c(0, n + 3), axes = FALSE,
       xlab = "", ylab = "", main = "")
  axis(3, at = rng, line = 0)
  title(main = title, line = 3, cex.main = 1.0, font.main = 2)
  
  y0 <- n + 1.7
  segments(min(rng), y0, max(rng), y0)
  segments(rng, y0 - 0.15, rng, y0 + 0.15)
  
  # stub length scales with the axis span, and the label sits one full
  # space-character's width beyond the stub, both measured in the plot's
  # own coordinate system, so labels clear the line on every diagram
  # regardless of rank span or label length
  stub <- 0.05 * diff(range(rng))
  gap  <- strwidth(" ", cex = 0.9) * 1.5
  
  for (k in seq_along(left)) {
    m <- methods[left[k]]; y <- n - k + 1
    x_end <- min(rng) - stub
    segments(ranks[m], y0, ranks[m], y)
    segments(ranks[m], y, x_end, y)
    text(x_end - gap, y, m, adj = 1, cex = 0.9)
  }
  for (k in seq_along(right)) {
    m <- methods[right[k]]; y <- n - (half + k) + 1
    x_end <- max(rng) + stub
    segments(ranks[m], y0, ranks[m], y)
    segments(ranks[m], y, x_end, y)
    text(x_end + gap, y, m, adj = 0, cex = 0.9)
  }
  
  # -- connecting bars, one row per clique, stacked below the axis ----------
  bar_y <- y0 - 0.5
  for (cl in cliques) {
    r <- ranks[cl]
    segments(min(r), bar_y, max(r), bar_y, lwd = 3)
    bar_y <- bar_y - 0.35
  }
}

# ---- 2. HELPER: RUN PIPELINE ON A TIED-RANK MATRIX ----------------------

run_tied_rank_pipeline <- function(tied_mat, mcc_mat, label, output_prefix) {
  
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("ANALYSIS:", label, "\n")
  cat(rep("=", 60), "\n\n")
  
  n_configs <- nrow(tied_mat)
  n_methods <- ncol(tied_mat)
  
  cat("Configurations:", n_configs, "\n")
  cat("Methods:", n_methods, "\n\n")
  
  # -- 2a. Average tied ranks --------------------------------------------
  avg_tied_ranks <- colMeans(tied_mat, na.rm = TRUE)
  mean_mcc_global <- colMeans(mcc_mat, na.rm = TRUE)
  
  ranks_df <- data.frame(
    Method       = names(avg_tied_ranks),
    Avg_TiedRank = round(avg_tied_ranks, 3),
    Mean_MCC     = round(mean_mcc_global[names(avg_tied_ranks)], 4)
  )
  ranks_df <- ranks_df[order(ranks_df$Avg_TiedRank), ]
  
  cat("Average tied ranks (lower = better):\n")
  print(ranks_df, row.names = FALSE)
  
  rank_csv <- file.path(output_dir, paste0(output_prefix, "_ranking_table.csv"))
  write.csv(ranks_df, rank_csv, row.names = FALSE)
  cat("\nRanking table saved:", rank_csv, "\n")
  
  # -- 2b. Pairwise Wilcoxon + Holm correction ---------------------------
  # Operate on tied-rank vectors (not mean MCC) per supervisor's instruction
  cat("\n--- Pairwise Wilcoxon + Holm correction (on tied ranks) ---\n")
  
  method_names <- colnames(tied_mat)
  p_raw_mat <- matrix(1, nrow = n_methods, ncol = n_methods,
                      dimnames = list(method_names, method_names))
  
  for (i in 1:(n_methods - 1)) {
    for (j in (i + 1):n_methods) {
      xi <- tied_mat[, i]
      xj <- tied_mat[, j]
      ok <- !is.na(xi) & !is.na(xj)
      if (sum(ok) >= 3) {
        wt <- wilcox.test(xi[ok], xj[ok], paired = TRUE, exact = FALSE)
        p_raw_mat[i, j] <- wt$p.value
        p_raw_mat[j, i] <- wt$p.value
      }
    }
  }
  
  p_upper <- p_raw_mat[upper.tri(p_raw_mat)]
  p_adj   <- p.adjust(p_upper, method = "holm")
  
  p_adj_mat <- matrix(1, nrow = n_methods, ncol = n_methods,
                      dimnames = list(method_names, method_names))
  p_adj_mat[upper.tri(p_adj_mat)] <- p_adj
  p_adj_mat[lower.tri(p_adj_mat)] <- t(p_adj_mat)[lower.tri(p_adj_mat)]
  
  n_sig <- sum(p_adj_mat[upper.tri(p_adj_mat)] < 0.05)
  cat("Significant pairs (alpha=0.05):", n_sig, "\n")
  
  # Write to results file
  sink(sink_file, append = TRUE)
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("ANALYSIS:", label, "\n")
  cat(rep("=", 60), "\n\n")
  cat("Average tied ranks:\n")
  print(ranks_df, row.names = FALSE)
  cat("\nPairwise Wilcoxon (Holm-corrected) — significant pairs (alpha=0.05):\n")
  sig_pairs <- which(p_adj_mat < 0.05 & upper.tri(p_adj_mat), arr.ind = TRUE)
  if (nrow(sig_pairs) == 0) {
    cat("None\n")
  } else {
    for (k in seq_len(nrow(sig_pairs))) {
      r <- sig_pairs[k, 1]; co <- sig_pairs[k, 2]
      cat(sprintf("  %s vs %s: p_adj = %.4f\n",
                  rownames(p_adj_mat)[r],
                  colnames(p_adj_mat)[co],
                  p_adj_mat[r, co]))
    }
  }
  sink()
  
  # -- 2c. CD Diagram (driven by the actual Wilcoxon+Holm result) --------
  # Bars now connect methods per p_adj_mat computed in 2b above, instead of
  # plotCD()'s internal Nemenyi recomputation on the raw MCC matrix, so the
  # figure and the reported significant pairs are guaranteed to agree.
  pdf_path <- file.path(output_dir, paste0(output_prefix, "_cd_diagram.pdf"))
  pdf(pdf_path, width = 12, height = 5)
  tryCatch({
    plot_wilcoxon_cd(
      avg_ranks = avg_tied_ranks,
      p_adj_mat = p_adj_mat,
      alpha     = 0.05,
      title     = paste0("Critical Difference Diagram (Wilcoxon+Holm): ", label)
    )
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("CD diagram error:", e$message), cex = 0.8)
  })
  dev.off()
  cat("CD diagram saved:", pdf_path, "\n")
  
  return(list(
    ranks     = ranks_df,
    p_adj_mat = p_adj_mat
  ))
}

# ---- 3. BUILD MCC MATRIX (for CD diagram and heatmap) -------------------

cat("Loading mcc_long.csv for heatmap and CD diagrams...\n")
df_mcc <- read.csv(mcc_csv, stringsAsFactors = FALSE)

build_mean_mcc_matrix <- function(data) {
  agg <- data %>%
    group_by(config, method_label) %>%
    summarise(mean_MCC = mean(MCC, na.rm = TRUE), .groups = "drop")
  mat <- acast(agg, config ~ method_label, value.var = "mean_MCC")
  return(mat)
}

# ---- 4. GLOBAL ANALYSIS: ALL 16 MAIN CONFIGS (S1-S4) --------------------

main_configs <- configs[!grepl("^S5", configs)]

tied_mat_main <- tied_mat_full[main_configs, , drop = FALSE]
# Remove methods that are all NA in main configs
valid_methods_main <- colnames(tied_mat_main)[colSums(!is.na(tied_mat_main)) > 0]
tied_mat_main <- tied_mat_main[, valid_methods_main, drop = FALSE]

mcc_mat_main <- build_mean_mcc_matrix(df_mcc[df_mcc$dgp_family != "S5", ])
# Align columns
common_methods_main <- intersect(valid_methods_main, colnames(mcc_mat_main))
tied_mat_main <- tied_mat_main[, common_methods_main, drop = FALSE]
mcc_mat_main  <- mcc_mat_main[, common_methods_main, drop = FALSE]

results_global <- run_tied_rank_pipeline(
  tied_mat      = tied_mat_main,
  mcc_mat       = mcc_mat_main,
  label         = "Global — All 16 Main Configurations (S1-S4)",
  output_prefix = "global"
)

# ---- 5. PER-DGP ANALYSIS ------------------------------------------------

for (dgp in c("S1", "S2", "S3", "S4")) {
  cfg_dgp <- configs[grepl(paste0("^", dgp), configs)]
  
  cat("\n--- DGP:", dgp, "(", length(cfg_dgp), "configs) ---\n")
  
  if (length(cfg_dgp) < 2) {
    cat("Skipping: need at least 2 configurations.\n")
    next
  }
  
  tied_dgp <- tied_mat_full[cfg_dgp, , drop = FALSE]
  valid_m  <- colnames(tied_dgp)[colSums(!is.na(tied_dgp)) > 0]
  tied_dgp <- tied_dgp[, valid_m, drop = FALSE]
  
  mcc_dgp  <- build_mean_mcc_matrix(df_mcc[df_mcc$dgp_family == dgp, ])
  common_m <- intersect(valid_m, colnames(mcc_dgp))
  tied_dgp <- tied_dgp[, common_m, drop = FALSE]
  mcc_dgp  <- mcc_dgp[, common_m, drop = FALSE]
  
  run_tied_rank_pipeline(
    tied_mat      = tied_dgp,
    mcc_mat       = mcc_dgp,
    label         = paste("DGP Family:", dgp),
    output_prefix = paste0("dgp_", dgp)
  )
}

# ---- 6. S5 ANALYSIS -----------------------------------------------------

s5_configs <- configs[grepl("^S5", configs)]

if (length(s5_configs) >= 2) {
  tied_s5 <- tied_mat_full[s5_configs, , drop = FALSE]
  valid_s5 <- colnames(tied_s5)[colSums(!is.na(tied_s5)) > 0]
  tied_s5  <- tied_s5[, valid_s5, drop = FALSE]
  
  mcc_s5   <- build_mean_mcc_matrix(df_mcc[df_mcc$dgp_family == "S5", ])
  common_s5 <- intersect(valid_s5, colnames(mcc_s5))
  tied_s5  <- tied_s5[, common_s5, drop = FALSE]
  mcc_s5   <- mcc_s5[, common_s5, drop = FALSE]
  
  run_tied_rank_pipeline(
    tied_mat      = tied_s5,
    mcc_mat       = mcc_s5,
    label         = "S5 — Sparse VAR(2), K=50 (200 realizations)",
    output_prefix = "s5"
  )
}

# ---- 7. MEAN MCC HEATMAP (all 18 configurations, including S5) ----------

cat("\nGenerating mean MCC heatmap...\n")

df_heatmap <- df_mcc %>%
  group_by(config, method_label) %>%
  summarise(mean_MCC = mean(MCC, na.rm = TRUE), .groups = "drop")

method_order <- results_global$ranks$Method
df_heatmap$method_label <- factor(df_heatmap$method_label, levels = rev(method_order))
config_order_all <- c(rownames(mcc_mat_main), s5_configs)
df_heatmap$config <- factor(df_heatmap$config, levels = config_order_all)

heatmap_plot <- ggplot(df_heatmap, aes(x = config, y = method_label, fill = mean_MCC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", mean_MCC)), size = 2.2, color = "black") +
  scale_fill_gradient2(
    low      = "#d73027",
    mid      = "#ffffbf",
    high     = "#1a9850",
    midpoint = 0.5,
    limits   = c(-0.1, 1),
    name     = "Mean MCC"
  ) +
  labs(
    title    = "Mean MCC Heatmap — 14 Methods × 18 Configurations",
    subtitle = "Green = high performance, Red = low performance | Methods ordered by global tied-rank",
    x        = "Configuration",
    y        = "Method (ordered by global tied-rank)"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 8),
    plot.title      = element_text(face = "bold", size = 11),
    legend.position = "right"
  )

heatmap_path <- file.path(output_dir, "mcc_heatmap.pdf")
ggsave(heatmap_path, heatmap_plot, width = 18, height = 7)
cat("Heatmap saved:", heatmap_path, "\n")

# ---- 8. FINAL SUMMARY ---------------------------------------------------

cat("\n", rep("=", 60), "\n", sep = "")
cat("ALL ANALYSES COMPLETE\n")
cat("Output files in:", output_dir, "\n")
cat(rep("=", 60), "\n")