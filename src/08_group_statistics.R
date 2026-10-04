library(ggplot2)

 
# STEP 5: COMPREHENSIVE STATISTICAL ANALYSIS
 

 
# 1. Load processed intensity matrix
 

df <- read.csv(
  "feature_intensity_matrix.csv",
  row.names = 1,
  check.names = FALSE
)

cat("Feature matrix:", nrow(df), "features x", ncol(df), "samples\n")


 
# 2. Identify biological samples
 

# Groups are encoded in the sample names:
# STDMIX_A1 = Group A
# STDMIX_B1 = Group B
# STDMIX_C1 = Group C
# STDMIX_D1 = Group D

biological_cols <- grep(
  "STDMIX_[ABCD]1_",
  colnames(df),
  value = TRUE
)

bio <- df[, biological_cols, drop = FALSE]

groups <- sub(
  ".*STDMIX_([ABCD])1_.*",
  "\\1",
  biological_cols
)

groups <- factor(groups, levels = c("A", "B", "C", "D"))

cat("\nBiological samples:", ncol(bio), "\n")
print(table(groups))


 
# 3. Log transformation
 

# Add 1 to avoid log2(0)

bio_log <- log2(bio + 1)


 
# 4. Remove problematic features
 

# Remove features containing NA/Inf

valid <- apply(
  bio_log,
  1,
  function(x) all(is.finite(x))
)

bio_log <- bio_log[valid, ]

# Remove features with zero variance

variable <- apply(
  bio_log,
  1,
  function(x) sd(x) > 0
)

bio_log <- bio_log[variable, ]

cat(
  "\nFeatures retained for statistical analysis:",
  nrow(bio_log),
  "\n"
)


 
# 5. PCA
 

pca <- prcomp(
  t(bio_log),
  center = TRUE,
  scale. = TRUE
)

pca_scores <- as.data.frame(pca$x[, 1:2])

pca_scores$Sample <- rownames(pca_scores)
pca_scores$Group <- groups

percent_var <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_plot <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    shape = Group
  )
) +
  geom_point(size = 4) +
  geom_text(
    aes(label = Group),
    vjust = -1,
    size = 4
  ) +
  labs(
    title = "PCA of Biological Samples",
    x = paste0("PC1 (", round(percent_var[1], 1), "%)"),
    y = paste0("PC2 (", round(percent_var[2], 1), "%)")
  ) +
  theme_minimal(base_size = 14)

ggsave(
  "biological_pca_final.png",
  pca_plot,
  width = 9,
  height = 7
)

write.csv(
  pca_scores,
  "biological_pca_scores_final.csv",
  row.names = FALSE
)


 
# 6. One-way ANOVA for every feature
 

cat("\nRunning feature-wise ANOVA...\n")

anova_results <- data.frame(
  feature_id = rownames(bio_log),
  p_value = NA_real_,
  F_statistic = NA_real_,
  mean_A = NA_real_,
  mean_B = NA_real_,
  mean_C = NA_real_,
  mean_D = NA_real_
)

for (i in seq_len(nrow(bio_log))) {

  values <- as.numeric(bio_log[i, ])

  model <- aov(
    values ~ groups
  )

  result <- summary(model)[[1]]

  anova_results$p_value[i] <- result$`Pr(>F)`[1]
  anova_results$F_statistic[i] <- result$`F value`[1]

  anova_results$mean_A[i] <- mean(values[groups == "A"])
  anova_results$mean_B[i] <- mean(values[groups == "B"])
  anova_results$mean_C[i] <- mean(values[groups == "C"])
  anova_results$mean_D[i] <- mean(values[groups == "D"])
}


 
# Multiple-testing correction
 

anova_results$FDR <- p.adjust(
  anova_results$p_value,
  method = "BH"
)


 
# Effect size: eta squared
 

anova_results$eta_squared <- NA_real_

for (i in seq_len(nrow(bio_log))) {

  values <- as.numeric(bio_log[i, ])

  grand_mean <- mean(values)

  ss_between <- sum(
    table(groups) *
      (
        tapply(values, groups, mean) -
        grand_mean
      )^2
  )

  ss_total <- sum(
    (values - grand_mean)^2
  )

  anova_results$eta_squared[i] <-
    ss_between / ss_total
}


 
# Maximum group difference
 

anova_results$max_mean_difference <- apply(
  anova_results[, c(
    "mean_A",
    "mean_B",
    "mean_C",
    "mean_D"
  )],
  1,
  function(x) max(x) - min(x)
)

anova_results$significant_FDR_05 <-
  anova_results$FDR < 0.05


anova_results <- anova_results[
  order(anova_results$FDR),
]

write.csv(
  anova_results,
  "group_anova_results_final.csv",
  row.names = FALSE
)


cat(
  "\nFeatures with FDR < 0.05:",
  sum(anova_results$FDR < 0.05),
  "\n"
)


 
# 7. Top significant features
 

top_features <- head(
  anova_results[
    order(
      anova_results$FDR,
      -anova_results$eta_squared
    ),
  ],
  20
)

write.csv(
  top_features,
  "top_20_differential_features.csv",
  row.names = FALSE
)

print(
  top_features[
    ,
    c(
      "feature_id",
      "FDR",
      "F_statistic",
      "eta_squared",
      "max_mean_difference"
    )
  ]
)


 
# 8. Heatmap of top 25 significant features
 

top25 <- head(
  anova_results$feature_id[
    order(anova_results$FDR)
  ],
  25
)

heat <- bio_log[top25, ]

# Standardise each feature so the heatmap shows
# relative differences between groups

heat_scaled <- t(
  scale(t(heat))
)

heat_df <- as.data.frame(heat_scaled)

heat_df$Feature <- rownames(heat_df)

heat_long <- reshape(
  heat_df,
  varying = colnames(heat_df)[
    colnames(heat_df) != "Feature"
  ],
  v.names = "Expression",
  timevar = "Sample",
  times = colnames(heat_df)[
    colnames(heat_df) != "Feature"
  ],
  direction = "long"
)

rownames(heat_long) <- NULL

heat_long$Group <- groups[
  match(
    heat_long$Sample,
    biological_cols
  )
]

heat_plot <- ggplot(
  heat_long,
  aes(
    x = Sample,
    y = Feature,
    fill = Expression
  )
) +
  geom_tile() +
  labs(
    title = "Top 25 Differential Features",
    x = "Biological samples",
    y = "Feature"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  )

ggsave(
  "top25_feature_heatmap.png",
  heat_plot,
  width = 12,
  height = 9
)


 
# 9. Boxplots for top 4 features
 
 
# Plot the top 4 most significant features
 

 
# Plot the top 4 most significant features
 

top4 <- anova_results$feature_id[
  order(anova_results$FDR)
][1:4]

plot_data <- do.call(
  rbind,
  lapply(top4, function(feature) {

    values <- as.numeric(bio_log[feature, ])

    data.frame(
      Feature = feature,
      Sample = colnames(bio_log),
      Intensity = values,
      Group = groups
    )
  })
)

plot_data$Feature <- factor(
  plot_data$Feature,
  levels = top4
)

p <- ggplot(
  plot_data,
  aes(x = Group, y = Intensity)
) +
  geom_boxplot() +
  geom_jitter(
    width = 0.15,
    alpha = 0.7
  ) +
  facet_wrap(
    ~ Feature,
    scales = "free_y"
  ) +
  labs(
    title = "Top Differentially Abundant Features",
    x = "Biological Group",
    y = "Log-transformed intensity"
  ) +
  theme_minimal()

ggsave(
  "top_differential_features.png",
  p,
  width = 10,
  height = 7,
  dpi = 300
)

cat("\nSaved:\n")
cat("  top_differential_features.png\n")
 
# 10. Summary
 

cat("\n")
cat("========================================\n")
cat("STEP 5 STATISTICAL ANALYSIS COMPLETE\n")
cat("========================================\n")
cat(
  "Biological samples:",
  ncol(bio_log),
  "\n"
)
cat(
  "Features analysed:",
  nrow(bio_log),
  "\n"
)
cat(
  "FDR < 0.05:",
  sum(anova_results$FDR < 0.05),
  "\n"
)
cat(
  "Top feature:",
  anova_results$feature_id[1],
  "\n"
)
cat(
  "Top feature FDR:",
  signif(anova_results$FDR[1], 4),
  "\n"
)
cat(
  "Top feature eta-squared:",
  round(anova_results$eta_squared[1], 3),
  "\n"
)

cat("\nSaved:\n")
cat("  biological_pca_final.png\n")
cat("  biological_pca_scores_final.csv\n")
cat("  group_anova_results_final.csv\n")
cat("  top_20_differential_features.csv\n")
cat("  top25_feature_heatmap.png\n")
cat("  top4_feature_boxplots.png\n")