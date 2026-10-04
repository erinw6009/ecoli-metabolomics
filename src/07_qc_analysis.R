
# ============================================================
# Task 5 — QC analysis
# Assess analytical stability and reproducibility
# ============================================================

library(ggplot2)

# ------------------------------------------------------------
# 1. Load intensity matrix
# ------------------------------------------------------------

feature_matrix <- read.csv(
  "feature_intensity_matrix.csv",
  row.names = 1,
  check.names = FALSE
)

cat("Feature matrix:", nrow(feature_matrix), "features x",
    ncol(feature_matrix), "samples\n")

# ------------------------------------------------------------
# 2. Identify QC samples
# ------------------------------------------------------------

qc_cols <- grep(
  "QCPool",
  colnames(feature_matrix),
  ignore.case = TRUE
)

qc_matrix <- feature_matrix[, qc_cols, drop = FALSE]

cat("QC samples:", ncol(qc_matrix), "\n")

# ------------------------------------------------------------
# 3. Log transform
# ------------------------------------------------------------
# Log transformation reduces the effect of extremely
# abundant features dominating the analysis.

qc_log <- log10(qc_matrix + 1)

# ------------------------------------------------------------
# 4. Pairwise QC correlations
# ------------------------------------------------------------

qc_cor <- cor(
  qc_log,
  method = "pearson",
  use = "pairwise.complete.obs"
)

write.csv(
  qc_cor,
  "qc_pairwise_correlations.csv"
)

cat("\nPairwise QC correlations:\n")
print(round(qc_cor, 3))

# ------------------------------------------------------------
# 5. Calculate feature-level QC CV
# ------------------------------------------------------------

qc_mean <- rowMeans(
  qc_matrix,
  na.rm = TRUE
)

qc_sd <- apply(
  qc_matrix,
  1,
  sd,
  na.rm = TRUE
)

qc_cv <- qc_sd / qc_mean

qc_summary <- data.frame(
  feature_id = rownames(qc_matrix),
  mean_intensity = qc_mean,
  sd_intensity = qc_sd,
  CV = qc_cv
)

qc_summary <- qc_summary[
  is.finite(qc_summary$CV),
]

write.csv(
  qc_summary,
  "qc_feature_cv.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Summarise QC CV
# ------------------------------------------------------------

cat("\nQC feature CV summary:\n")
print(summary(qc_summary$CV))

# ------------------------------------------------------------
# 7. Plot CV distribution
# ------------------------------------------------------------

p_cv <- ggplot(
  qc_summary,
  aes(x = CV)
) +
  geom_histogram(
    bins = 50
  ) +
  labs(
    title = "Distribution of Feature CVs Across QC Samples",
    x = "Coefficient of Variation",
    y = "Number of Features"
  ) +
  theme_minimal()

ggsave(
  "qc_cv_distribution.png",
  p_cv,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 8. PCA of QC samples
# ------------------------------------------------------------

# Transpose so samples are rows and features are columns.
pca_input <- t(qc_log)

# Keep only features with finite values in every QC sample.
finite_features <- apply(
  pca_input,
  2,
  function(x) all(is.finite(x))
)

pca_input <- pca_input[
  ,
  finite_features,
  drop = FALSE
]

# Remove features with zero variance.
feature_variance <- apply(
  pca_input,
  2,
  var
)

pca_input <- pca_input[
  ,
  is.finite(feature_variance) &
  feature_variance > 0,
  drop = FALSE
]

cat(
  "\nFeatures retained for QC PCA:",
  ncol(pca_input),
  "\n"
)

# PCA
pca <- prcomp(
  pca_input,
  center = TRUE,
  scale. = FALSE
)

pca_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

write.csv(
  pca_df,
  "qc_pca_coordinates.csv",
  row.names = FALSE
)

p_qc_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    label = Sample
  )
) +
  geom_point(size = 3) +
  geom_text(
    vjust = -0.7,
    size = 3
  ) +
  labs(
    title = "PCA of QC Samples",
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)")
  ) +
  theme_minimal()

ggsave(
  "qc_pca.png",
  p_qc_pca,
  width = 8,
  height = 6,
  dpi = 300
)

cat("\nQC ANALYSIS COMPLETE\n")