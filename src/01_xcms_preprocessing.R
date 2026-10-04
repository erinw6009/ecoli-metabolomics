# ============================================================
# Step 1: XCMS preprocessing of LC-MS mzML data
# E. coli metabolomics PhD application
# ============================================================

# Packages ---------------------------------------------------
library(xcms)
library(MSnbase)
library(BiocParallel)

# Use serial processing to keep memory usage manageable
bp <- SerialParam()

# Paths ------------------------------------------------------
raw_dir <- "raw_data_POS"

# Find all positive-mode mzML files
files <- list.files(
  raw_dir,
  pattern = "\\.mzML$",
  full.names = TRUE
)

files <- sort(files)

cat("Number of mzML files:", length(files), "\n")

# ------------------------------------------------------------
# 1. Load raw MS data
# ------------------------------------------------------------

raw_data <- readMSData(
  files,
  mode = "onDisk"
)

# ------------------------------------------------------------
# 2. Peak picking using CentWave
#
# Parameters used:
# ppm       = 15
# peakwidth = 5-20 seconds
# snthresh  = 10
# prefilter = minimum 3 scans and intensity 100
# mzdiff    = -0.001
# ------------------------------------------------------------

cwp <- CentWaveParam(
  ppm = 15,
  peakwidth = c(5, 20),
  snthresh = 10,
  prefilter = c(3, 100),
  mzCenterFun = "wMean",
  integrate = 1,
  mzdiff = -0.001,
  fitgauss = FALSE,
  noise = 0,
  verboseColumns = FALSE,
  firstBaselineCheck = TRUE
)

xdata <- findChromPeaks(
  raw_data,
  param = cwp,
  BPPARAM = bp
)

cat(
  "Chromatographic peaks detected:",
  nrow(chromPeaks(xdata)),
  "\n"
)

saveRDS(xdata, "xdata_peakpicked.rds")

# ------------------------------------------------------------
# 3. Retention-time alignment using Obiwarp
# ------------------------------------------------------------

obiwarp <- ObiwarpParam(
  binSize = 0.6
)

xdata_aligned <- adjustRtime(
  xdata,
  param = obiwarp
)

saveRDS(xdata_aligned, "xdata_aligned.rds")

# ------------------------------------------------------------
# 4. Peak grouping / correspondence
#
# Peaks with similar m/z and retention time are grouped
# into common features across samples.
# ------------------------------------------------------------

pdp <- PeakDensityParam(
  sampleGroups = rep(1, length(files)),
  minFraction = 0.5,
  minSamples = 1,
  binSize = 0.01
)

xdata_grouped <- groupChromPeaks(
  xdata_aligned,
  param = pdp
)

cat(
  "Features after peak grouping:",
  nrow(featureDefinitions(xdata_grouped)),
  "\n"
)

saveRDS(xdata_grouped, "xdata_grouped.rds")

# ------------------------------------------------------------
# 5. Fill missing chromatographic peaks
# ------------------------------------------------------------

xdata_filled <- fillChromPeaks(
  xdata_grouped
)

saveRDS(xdata_filled, "xdata_filled.rds")

# ------------------------------------------------------------
# 6. Extract feature intensity matrix
#
# "into" = integrated peak intensity
# Rows    = features
# Columns = samples
# ------------------------------------------------------------

feature_matrix <- featureValues(
  xdata_filled,
  value = "into"
)

# ------------------------------------------------------------
# 7. Extract feature metadata
#
# Keep the main feature descriptors required downstream:
# feature ID, median m/z and median retention time.
# ------------------------------------------------------------

feature_info <- featureDefinitions(xdata_filled)

feature_info_simple <- data.frame(
  feature_id = rownames(feature_info),
  mz = feature_info$mzmed,
  rt = feature_info$rtmed
)

# ------------------------------------------------------------
# 8. Export outputs for downstream Python / ipaPy2 analysis
# ------------------------------------------------------------

write.csv(
  feature_matrix,
  "feature_intensity_matrix.csv",
  row.names = TRUE
)

write.csv(
  feature_info_simple,
  "feature_metadata.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Basic validation
# ------------------------------------------------------------

cat(
  "Final intensity matrix dimensions:",
  paste(dim(feature_matrix), collapse = " x "),
  "\n"
)

cat(
  "Final feature metadata dimensions:",
  paste(dim(feature_info_simple), collapse = " x "),
  "\n"
)

cat("Step 1 complete.\n")