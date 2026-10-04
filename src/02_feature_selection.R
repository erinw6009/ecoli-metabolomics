 
# Step 2: Feature selection for ipaPy2
 
library(xcms)
library(MSnbase)
# Load XCMS result
xdata_filled <- readRDS("xdata_filled.rds")

# Extract intensity matrix
feature_matrix <- featureValues(
  xdata_filled,
  value = "into"
)

# Identify the 7 QC-pool injections
qc_cols <- grep(
  "QCPool",
  colnames(feature_matrix),
  ignore.case = TRUE
)

cat("Number of QC samples:", length(qc_cols), "\n")

# Extract QC intensities
qc_matrix <- feature_matrix[, qc_cols, drop = FALSE]

# Calculate the proportion of QC samples in which
# each feature was detected with non-zero intensity.
qc_detection_rate <- rowMeans(
  qc_matrix > 0,
  na.rm = TRUE
)

# Calculate coefficient of variation across QC samples.
# CV = standard deviation / mean.
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

# Create feature-selection table
feature_selection <- data.frame(
  feature_id = rownames(feature_matrix),
  qc_detection_rate = qc_detection_rate,
  qc_cv = qc_cv,
  qc_mean_intensity = qc_mean
)

# Remove features that cannot be meaningfully assessed
feature_selection <- feature_selection[
  is.finite(feature_selection$qc_cv) &
  is.finite(feature_selection$qc_mean_intensity),
]

# Require detection in at least 5 of the 7 QC injections
reliable_features <- feature_selection[
  feature_selection$qc_detection_rate >= 5 / 7,
]

cat(
  "Features detected in >=5/7 QC samples:",
  nrow(reliable_features),
  "\n"
)

# Rank by QC reproducibility first, then signal intensity.
# Lower CV = more reproducible measurement.
reliable_features <- reliable_features[
  order(
    reliable_features$qc_cv,
    -reliable_features$qc_mean_intensity
  ),
]

# Select at most 1,000 analytically reliable features
n_select <- min(
  1000,
  nrow(reliable_features)
)

selected_features <- reliable_features[
  seq_len(n_select),
]

# Save selection
write.csv(
  selected_features,
  "selected_features_1000.csv",
  row.names = FALSE
)

cat(
  "Selected features for ipaPy2:",
  nrow(selected_features),
  "\n"
)
 
# Prepare selected features for ipaPy2
 

metadata <- read.csv(
  "feature_metadata.csv",
  stringsAsFactors = FALSE
)

ipapy2_input <- merge(
  selected_features,
  metadata,
  by = "feature_id"
)

ipapy2_input <- ipapy2_input[
  ,
  c(
    "feature_id",
    "mz",
    "rt",
    "qc_detection_rate",
    "qc_cv",
    "qc_mean_intensity"
  )
]

cat(
  "Features available for ipaPy2:",
  nrow(ipapy2_input),
  "\n"
)

cat(
  "Features with m/z:",
  sum(!is.na(ipapy2_input$mz)),
  "\n"
)

cat(
  "Features with retention time:",
  sum(!is.na(ipapy2_input$rt)),
  "\n"
)

write.csv(
  ipapy2_input,
  "ipapy2_input_1000.csv",
  row.names = FALSE
)