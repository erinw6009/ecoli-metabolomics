 
# Step 2: ipaPy2 baseline MS1 annotation
 

import os
import pandas as pd
from ipaPy2 import ipa

 
# Paths
 

SELECTED = "ipapy2_input_1000.csv"
INTENSITIES = "feature_intensity_matrix.csv"
DB_DIR = "ipaPy2/DB"

 
# 1. Load selected features
 

selected = pd.read_csv(SELECTED)

selected["feature_id"] = selected["feature_id"].astype(str)

print("Selected features:", len(selected))

 
# 2. Load the full XCMS intensity matrix
 

intensity = pd.read_csv(
    INTENSITIES,
    index_col=0
)

intensity.index = intensity.index.astype(str)

# Keep the same feature order as selected_features_1000.csv
intensity = intensity.loc[selected["feature_id"]]

print("Intensity matrix:", intensity.shape)

 
# 3. Build ipaPy2 input
#
# ipaPy2 expects:
#
#   column 1 = ids
#   column 2 = mzs
#   column 3 = RTs
#   columns 4 onward = intensity for each sample
#
 

df = pd.DataFrame({
    "ids": selected["feature_id"].values,
    "mzs": selected["mz"].values,
    "RTs": selected["rt"].values
})

# Add the 35 sample intensity columns
for sample in intensity.columns:
    df[sample] = intensity[sample].values

print("ipaPy2 dataframe:", df.shape)
print("Number of intensity columns:", len(df.columns) - 3)

 
# 4. Cluster features
 

print("Clustering features...")

df = ipa.clusterFeatures(
    df,
    Cthr=0.8,
    RTwin=1,
    Intmode="max"
)

print("Feature clustering complete.")
print("Clusters/features returned:", len(df))

 
# 5. Map isotope patterns
 

print("Mapping isotope patterns...")

ipa.map_isotope_patterns(
    df,
    isoDiff=1,
    ppm=100,
    ionisation=1,
    MinIsoRatio=0.5
)
print("Isotope mapping complete.")

 
# 6. Load ipaPy2's default repository database
 

adduct_file = os.path.join(
    DB_DIR,
    "adducts.csv"
)

db_file = os.path.join(
    DB_DIR,
    "IPA_MS1.csv"
)

adducts = pd.read_csv(adduct_file)
db = pd.read_csv(db_file)

print("Adduct database rows:", len(adducts))
print("IPA MS1 database rows:", len(db))

 
# 7. Compute possible adducts
 

print("Computing adducts...")

all_adds = ipa.compute_all_adducts(
    adducts,
    db,
    ionisation=1,
    ncores=1
)

print("Adduct calculation complete.")
print("Possible adduct entries:", len(all_adds))

 
# 8. Baseline MS1 annotation
#
# 15 ppm matches the XCMS mass-trace setting used earlier.
 

print("Running ipaPy2 MS1 annotation...")

annotations = ipa.MS1annotation(
    df,
    all_adds,
    ppm=15,
    ncores=1
)

print("MS1 annotation complete.")
print("Annotated feature entries:", len(annotations))

 
# 9. Save processed feature information
 

df.to_csv(
    "ipapy2_features_baseline.csv",
    index=False
)

# Save raw Python annotation dictionary
import pickle

with open(
    "ipapy2_baseline_annotations.pkl",
    "wb"
) as f:
    pickle.dump(annotations, f)

print("Saved:")
print("  ipapy2_features_baseline.csv")
print("  ipapy2_baseline_annotations.pkl")

print("BASELINE IPA-PY2 STEP COMPLETE")