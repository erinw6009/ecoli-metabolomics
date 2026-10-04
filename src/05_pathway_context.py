 
# Step 4: Biochemical-Connection Context
 
#
# Goal:
# Use ECMDB-supported ipaPy2 annotations and known biochemical
# connections between KEGG compounds to add biological context.
#
# Step 3 established the ECMDB mapping.
# Step 4 uses those KEGG IDs to find connections between
# supported metabolites.
 

import pandas as pd


 
# 1. Load Step 3 ECMDB context results
 

context = pd.read_csv(
    "ipapy2_ecmdb_context_candidates.csv"
)

print(
    "Step 3 candidate rows:",
    len(context)
)


 
# 2. Keep only ECMDB-supported candidates
 

supported = context[
    context["ecmdb_match"] == True
].copy()

print(
    "ECMDB-supported candidate rows:",
    len(supported)
)


 
# 3. Check KEGG IDs
 

supported["kegg_id"] = (
    supported["kegg_id"]
    .astype("string")
    .str.strip()
)

supported_kegg = set(
    supported[
        "kegg_id"
    ]
    .dropna()
)

print(
    "Unique ECMDB-supported KEGG metabolites:",
    len(supported_kegg)
)


 
# 4. Load biochemical connections
 

connections = pd.read_csv(
    "ipaPy2/DB/allBIO_reactions.csv"
)

connections = connections[
    ["X0", "X1"]
].dropna()

connections["X0"] = (
    connections["X0"]
    .astype(str)
    .str.strip()
)

connections["X1"] = (
    connections["X1"]
    .astype(str)
    .str.strip()
)

print(
    "Total biochemical connections:",
    len(connections)
)


 
# 5. Find connections between ECMDB-supported metabolites
 

connected_pairs = []

for _, row in connections.iterrows():

    a = row["X0"]
    b = row["X1"]

    if (
        a in supported_kegg
        and b in supported_kegg
    ):

        connected_pairs.append(
            (a, b)
        )


connected_pairs = sorted(
    set(connected_pairs)
)


print(
    "Connections between ECMDB-supported metabolites:",
    len(connected_pairs)
)


 
# 6. Identify KEGG compounds with biochemical support
 

connected_kegg = set()

for a, b in connected_pairs:

    connected_kegg.add(a)
    connected_kegg.add(b)


 
# 7. Mark candidate annotations
 

supported["biochemical_context"] = (
    supported["kegg_id"]
    .isin(connected_kegg)
)


 
# 8. Save candidate-level results
 

supported.to_csv(
    "ipapy2_ecmdb_biochemical_candidates.csv",
    index=False
)


 
# 9. Feature-level summary
 

feature_summary = (
    supported
    .groupby("feature_id")
    .agg(
        n_ecmdb_candidates=(
            "ipapy2_id",
            "count"
        ),
        n_biochemically_supported=(
            "biochemical_context",
            "sum"
        ),
        biochemical_context=(
            "biochemical_context",
            "any"
        )
    )
    .reset_index()
)


feature_summary.to_csv(
    "ipapy2_biochemical_feature_summary.csv",
    index=False
)


 
# 10. Convert KEGG pairs into readable metabolite pairs
 

kegg_to_name = (
    supported[
        ["kegg_id", "name_ecmdb"]
    ]
    .dropna()
    .drop_duplicates(
        subset=["kegg_id"]
    )
    .set_index("kegg_id")[
        "name_ecmdb"
    ]
    .to_dict()
)


pair_table = pd.DataFrame(
    connected_pairs,
    columns=[
        "kegg_id_1",
        "kegg_id_2"
    ]
)

pair_table["name_1"] = (
    pair_table["kegg_id_1"]
    .map(kegg_to_name)
)

pair_table["name_2"] = (
    pair_table["kegg_id_2"]
    .map(kegg_to_name)
)


pair_table.to_csv(
    "ipapy2_ecmdb_biochemical_connections.csv",
    index=False
)


 
# 11. Summary
 

print()
print("STEP 4 COMPLETE")
print("----------------")

print(
    "ECMDB-supported KEGG metabolites:",
    len(supported_kegg)
)

print(
    "Biochemical connections between "
    "supported metabolites:",
    len(connected_pairs)
)

print(
    "Features with biochemical context:",
    feature_summary[
        "biochemical_context"
    ].sum()
)

print()
print("Saved:")
print(
    "  ipapy2_ecmdb_biochemical_candidates.csv"
)
print(
    "  ipapy2_biochemical_feature_summary.csv"
)
print(
    "  ipapy2_ecmdb_biochemical_connections.csv"
)