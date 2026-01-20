# ============================================================
# TSNE + KMEANS (Statistica-like) + LABELS
# - Read CSV (';' separator)
# - K-means (Hartigan-Wong) with init centers = first K rows
# - t-SNE 2D visualization
# - Plot colored by cluster + participant labels (ggrepel)
# - Save outputs
#
# IMPORTANT:
# - If your CSV already contains Z-values, do NOT standardize again.
# ============================================================

# -----------------------------
# PARAMETERS
# -----------------------------
CSV_PATH <- "Z_values_v2.csv"
SEP      <- ";"
OUT_DIR  <- "kmeans_outputs_R"

SUBJECT_COL <- "Subject"
FEATURE_COLS <- c(
  "corr_Duration_ImpactSpeed",
  "corr_Duration_Amplitude",
  "corr_ImpactSpeed_Amplitude"
)

K_STAT <- 3
ITER_MAX_STAT <- 10
SET_SEED <- 0

# If you sorted in Statistica before clustering, set TRUE
SORT_BY_SUBJECT <- FALSE

# t-SNE parameters
TSNE_DIMS <- 2
TSNE_PERPLEXITY <- 5      # For ~33 samples, 5 is safe. Must be < (n-1)/3
TSNE_MAX_ITER <- 1000

# Plot options
POINT_SIZE <- 3
LABEL_SIZE <- 3.5
SAVE_WIDTH <- 8
SAVE_HEIGHT <- 7
SAVE_DPI <- 300

# -----------------------------
# SETUP
# -----------------------------
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Packages
if (!requireNamespace("Rtsne", quietly = TRUE)) install.packages("Rtsne")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")

library(Rtsne)
library(ggplot2)
library(ggrepel)

# -----------------------------
# LOAD DATA
# -----------------------------
df <- read.csv(CSV_PATH, sep = SEP, header = TRUE, stringsAsFactors = FALSE)

required_cols <- c(SUBJECT_COL, FEATURE_COLS)
missing <- setdiff(required_cols, colnames(df))
if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse = ", ")))

if (SORT_BY_SUBJECT) {
  df <- df[order(df[[SUBJECT_COL]]), ]
}

subjects <- df[[SUBJECT_COL]]

# Feature matrix
X <- as.matrix(df[, FEATURE_COLS])
storage.mode(X) <- "double"

n <- nrow(X)
if (TSNE_PERPLEXITY >= (n - 1) / 3) {
  stop(sprintf("Perplexity too high for n=%d. Choose perplexity < (n-1)/3 = %.2f",
               n, (n - 1) / 3))
}

cat("[INFO] Data loaded. n =", n, "| p =", ncol(X), "\n")
cat("[INFO] First K rows used as init centers for k-means:\n")
print(df[1:K_STAT, c(SUBJECT_COL, FEATURE_COLS)])

# -----------------------------
# K-MEANS (Statistica-like init)
# -----------------------------
init_centers <- X[1:K_STAT, , drop = FALSE]

set.seed(SET_SEED)
km <- kmeans(
  X,
  centers   = init_centers,       # init = first K observations
  iter.max  = ITER_MAX_STAT,      # Statistica iterations = 10
  algorithm = "Hartigan-Wong"
)

clusters <- km$cluster  # 1..K

# Save assignments
df_assign <- data.frame(Subject = subjects, Cluster = clusters)
df_assign <- df_assign[order(df_assign$Cluster, df_assign$Subject), ]
write.table(
  df_assign,
  file = file.path(OUT_DIR, sprintf("assignments_K%d.csv", K_STAT)),
  sep = SEP, row.names = FALSE, col.names = TRUE, quote = FALSE
)

# Save cluster lists
split_subjects <- split(df_assign$Subject, df_assign$Cluster)
df_counts <- data.frame(
  Cluster  = as.integer(names(split_subjects)),
  N        = as.integer(sapply(split_subjects, length)),
  Subjects = I(unname(split_subjects))
)
write.table(
  df_counts,
  file = file.path(OUT_DIR, sprintf("counts_K%d.csv", K_STAT)),
  sep = SEP, row.names = FALSE, col.names = TRUE, quote = FALSE
)

cat("\n===== KMEANS RESULT =====\n")
for (k in 1:K_STAT) {
  subs <- split_subjects[[as.character(k)]]
  cat(sprintf("Cluster %d: n=%d -> %s\n", k, length(subs), paste(subs, collapse = ", ")))
}

# -----------------------------
# t-SNE (2D) for visualization
# -----------------------------
set.seed(SET_SEED)
tsne_res <- Rtsne(
  X,
  dims = TSNE_DIMS,
  perplexity = TSNE_PERPLEXITY,
  verbose = TRUE,
  max_iter = TSNE_MAX_ITER,
  check_duplicates = FALSE
)

df_tsne <- data.frame(
  Subject = subjects,
  Cluster = factor(clusters),
  TSNE1 = tsne_res$Y[, 1],
  TSNE2 = tsne_res$Y[, 2]
)

# Save t-SNE coordinates
write.table(
  df_tsne,
  file = file.path(OUT_DIR, sprintf("tsne_coords_K%d.csv", K_STAT)),
  sep = SEP, row.names = FALSE, col.names = TRUE, quote = FALSE
)

# -----------------------------
# Plot t-SNE with labels (ggrepel)
# -----------------------------
p_tsne <- ggplot(df_tsne, aes(x = TSNE1, y = TSNE2, color = Cluster)) +
  geom_point(size = POINT_SIZE, alpha = 0.9) +
  geom_text_repel(
    aes(label = Subject),
    size = LABEL_SIZE,
    max.overlaps = Inf,      # show all labels
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.3,
    show.legend = FALSE
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "t-SNE visualization of K-means clusters",
    subtitle = paste("K =", K_STAT,
                     "| kmeans = Hartigan–Wong (init = first K obs)",
                     "| perplexity =", TSNE_PERPLEXITY),
    x = "t-SNE dimension 1",
    y = "t-SNE dimension 2",
    color = "Cluster"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_tsne)

# Save plot
ggsave(
  filename = file.path(OUT_DIR, sprintf("tSNE_clusters_K%d_labeled.png", K_STAT)),
  plot = p_tsne,
  width = SAVE_WIDTH,
  height = SAVE_HEIGHT,
  dpi = SAVE_DPI
)

cat("\n[INFO] Saved outputs to:", OUT_DIR, "\n")
cat("[INFO] t-SNE figure saved as:",
    file.path(OUT_DIR, sprintf("tSNE_clusters_K%d_labeled.png", K_STAT)), "\n")


# ============================================================
# Within-Cluster Sum of Squares (WCSS) analysis
# ============================================================

# km$withinss : WSS per cluster (within-cluster sum of squares)
# km$tot.withinss : sum of WSS values (global WCSS)
# km$totss : total sum of squares (around the global mean)
# km$betweenss : between-cluster sum of squares

wcss_per_cluster <- km$withinss
wcss_total <- km$tot.withinss
tss_total <- km$totss
bss_total <- km$betweenss

# Variance explained (R² analogue for k-means)
variance_explained <- bss_total / tss_total

cat("\n===== WITHIN-CLUSTER SUM OF SQUARES (WCSS) ANALYSIS =====\n")
cat(sprintf("Total SS (TSS)        = %.6f\n", tss_total))
cat(sprintf("Between SS (BSS)      = %.6f\n", bss_total))
cat(sprintf("Total Within SS (WCSS)= %.6f\n", wcss_total))
cat(sprintf("Variance explained    = %.3f (BSS/TSS)\n", variance_explained))

# Detailed table per cluster
df_wcss <- data.frame(
  Cluster = 1:K_STAT,
  N = as.integer(table(clusters)),
  WithinSS = as.numeric(wcss_per_cluster)
)

# Add proportions (useful for reporting in papers)
df_wcss$WithinSS_Percent <- 100 * df_wcss$WithinSS / wcss_total

# Export CSV
write.table(
  df_wcss,
  file = file.path(OUT_DIR, sprintf("wcss_analysis_K%d.csv", K_STAT)),
  sep = SEP, row.names = FALSE, col.names = TRUE, quote = FALSE
)

# Figure: WSS per cluster
png(filename = file.path(OUT_DIR, sprintf("wcss_by_cluster_K%d.png", K_STAT)),
    width = 900, height = 600)

barplot(
  height = df_wcss$WithinSS,
  names.arg = paste0("C", df_wcss$Cluster),
  xlab = "Cluster",
  ylab = "Within-cluster sum of squares (WithinSS)",
  main = sprintf(
    "WCSS by cluster (K=%d) | Total WCSS=%.2f | Explained=%.2f%%",
    K_STAT, wcss_total, 100 * variance_explained
  )
)
grid()
dev.off()

cat("\n[INFO] WCSS table saved -> ",
    file.path(OUT_DIR, sprintf("wcss_analysis_K%d.csv", K_STAT)), "\n")
cat("[INFO] WCSS figure saved -> ",
    file.path(OUT_DIR, sprintf("wcss_by_cluster_K%d.png", K_STAT)), "\n")

