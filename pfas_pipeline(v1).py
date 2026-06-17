"""
PFAS Bioaccumulation Data Pipeline
===================================
Sources: ECOTOX (CSV exports) + EPA CompTox (bulk download or API)

Outputs
-------
  pfas_bioaccumulation_dataset.csv   – cleaned, ML-ready dataset
  pfas_gap_heatmap.png               – data gap heatmap
  feature_importance.png             – Random Forest feature importances
  model_predictions.png              – predicted vs actual log concentrations
  cross_species_validation.png       – leave-one-species-out performance

Usage
-----
  python pfas_pipeline.py

Adjust the INPUT PATHS section below before running.
"""

# ─────────────────────────────────────────────
# IMPORTS
# ─────────────────────────────────────────────
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score, LeaveOneGroupOut
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import mean_squared_error, r2_score

warnings.filterwarnings("ignore")


# ═══════════════════════════════════════════════════════════════════════════════
# INPUT PATHS  ← edit these before running
# ═══════════════════════════════════════════════════════════════════════════════
ECOTOX_EXPORT_DIR  = "/Users/navink.admin/Desktop/ecotox_exports/"
COMPTOX_SNAPSHOT   = "/Users/navink.admin/Desktop/comptox_snapshot.csv"
OUTPUT_DIR         = "/Users/navink.admin/Desktop/"


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 — CompTox chemical feature table
# ═══════════════════════════════════════════════════════════════════════════════

# PFAS of interest (Tier 1 → Tier 3)
PFAS_CASRN = {
    # Tier 1
    "PFOS":    "1763-23-1",
    "PFOA":    "335-67-1",
    "PFHxS":   "355-46-4",
    "PFNA":    "375-95-1",
    "PFBS":    "375-73-5",
    # Tier 2
    "PFDA":    "335-76-2",
    "PFUnDA":  "2058-94-8",
    "PFDoDA":  "307-55-1",
    "PFHpA":   "375-85-9",
    "PFHxA":   "307-24-4",
    # Tier 3 – emerging
    "GenX":    "13252-13-6",
    "ADONA":   "958445-44-8",
    "F53B":    "73606-19-6",
}

# Manually curated chemical features (extend with CompTox data below)
PFAS_FEATURES = pd.DataFrame([
    # name,      casrn,           class,        chain, MW,     logKow
    ("PFOS",    "1763-23-1",   "Sulfonate",    8,  500.1,  5.26),
    ("PFOA",    "335-67-1",    "Carboxylate",  8,  414.1,  5.30),
    ("PFHxS",   "355-46-4",    "Sulfonate",    6,  400.1,  4.14),
    ("PFNA",    "375-95-1",    "Carboxylate",  9,  464.1,  6.05),
    ("PFBS",    "375-73-5",    "Sulfonate",    4,  300.1,  1.82),
    ("PFDA",    "335-76-2",    "Carboxylate", 10,  514.1,  6.83),
    ("PFUnDA",  "2058-94-8",   "Carboxylate", 11,  564.1,  7.59),
    ("PFDoDA",  "307-55-1",    "Carboxylate", 12,  614.1,  8.35),
    ("PFHpA",   "375-85-9",    "Carboxylate",  7,  364.1,  4.55),
    ("PFHxA",   "307-24-4",    "Carboxylate",  6,  314.1,  3.77),
    ("GenX",    "13252-13-6",  "Carboxylate",  6,  330.1,  2.50),
    ("ADONA",   "958445-44-8", "Carboxylate",  8,  380.1,  2.80),
    ("F53B",    "73606-19-6",  "Sulfonate",    6,  570.1,  4.00),
], columns=["PFAS_Name", "CASRN", "PFAS_Class", "Chain_Length", "MW", "LogKow"])


def load_comptox(path):
    """
    Merge CompTox bulk table (if available) into the feature table.
    Falls back silently to the manually curated table above.
    """
    try:
        ct = pd.read_csv(path, low_memory=False)
        ct = ct[ct["CASRN"].isin(PFAS_CASRN.values())]
        # Keep only columns we care about; rename to standard names
        col_map = {
            "PREFERRED_NAME":       "CompTox_Name",
            "MOLECULAR_WEIGHT":     "MW_ct",
            "LOGKOW":               "LogKow_ct",
        }
        ct = ct.rename(columns={k: v for k, v in col_map.items() if k in ct.columns})
        merged = PFAS_FEATURES.merge(ct[["CASRN"] + [v for v in col_map.values() if v in ct.columns]],
                                     on="CASRN", how="left")
        # Prefer CompTox values when available
        if "MW_ct" in merged:
            merged["MW"] = merged["MW_ct"].fillna(merged["MW"])
        if "LogKow_ct" in merged:
            merged["LogKow"] = merged["LogKow_ct"].fillna(merged["LogKow"])
        print(f"[CompTox] Merged {len(ct)} records from {path}")
        return merged[PFAS_FEATURES.columns]
    except FileNotFoundError:
        print(f"[CompTox] '{path}' not found — using curated feature table.")
        return PFAS_FEATURES


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 — ECOTOX ingestion
# ═══════════════════════════════════════════════════════════════════════════════

# Column aliases across different ECOTOX export formats
ECOTOX_COL_ALIASES = {
    "CASRN":             ["CAS Number", "CASRN", "CAS_Number", "cas_number"],
    "Chemical Name":     ["Chemical Name", "chemical_name", "Compound"],
    "Species":           ["Species Common Name", "Species", "species", "Common Name"],
    "Organism Group":    ["Species Group", "Organism Group", "organism_group", "Tax Group"],
    "Tissue":            ["Response Site", "Tissue", "tissue", "Body Part", "Matrix"],
    "Exposure Route":    ["Exposure Type", "Exposure Route", "exposure_route", "Route"],
    "Exposure Duration": ["Observed Duration Mean (Days)", "Exposure Duration Mean", "Exposure Duration", "Duration Mean", "duration"],
    "Duration Unit":     ["Observed Duration Units (Days)", "Exposure Duration Unit", "Duration Unit", "duration_unit"],
    "Result Value":      ["Conc 1 Mean (Standardized)", "Result Value Mean", "Result Value", "Concentration", "result_value"],
    "Result Unit":       ["Conc 1 Units (Standardized)", "Result Unit", "result_unit", "Conc Unit"],
    "BCF":               ["BCF 1 Value", "BCF", "bcf"],
    "Study Year":        ["Publication Year", "Study Year", "year"],
    "DOI":               ["DOI", "doi", "Reference DOI"],
}


def resolve_columns(df, alias_map):
    """Rename df columns to standard names using the alias map."""
    rename = {}
    for standard, aliases in alias_map.items():
        for alias in aliases:
            if alias in df.columns and standard not in df.columns:
                rename[alias] = standard
                break
    return df.rename(columns=rename)


def load_ecotox(export_dir):
    """Load and concatenate all ECOTOX CSVs/XLSXs in export_dir."""
    files = glob.glob(f"{export_dir}*.csv") + glob.glob(f"{export_dir}*.xlsx")
    if not files:
        print(f"[ECOTOX] No CSV/XLSX files found in '{export_dir}'. "
              "Place your ECOTOX web-export files there and re-run.")
        return pd.DataFrame(columns=list(ECOTOX_COL_ALIASES.keys()))

    frames = []
    for f in files:
        try:
            df = pd.read_excel(f) if f.endswith(".xlsx") else pd.read_csv(f, low_memory=False)
            df = resolve_columns(df, ECOTOX_COL_ALIASES)
            frames.append(df)
            print(f"[ECOTOX] Loaded {len(df):,} rows from {f}")
        except Exception as e:
            print(f"[ECOTOX] Could not load {f}: {e}")

    combined = pd.concat(frames, ignore_index=True)
    # Keep only columns we need (drop anything extra)
    keep = [c for c in ECOTOX_COL_ALIASES.keys() if c in combined.columns]
    return combined[keep].copy()


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Harmonization layer
# ═══════════════════════════════════════════════════════════════════════════════

UNIT_TO_NG_G = {
    # mass/mass (dry or wet weight)
    "ng/g":  1,
    "ug/g":  1_000,
    "µg/g":  1_000,
    "mg/kg": 1_000,
    "mg/g":  1_000_000,
    "pg/g":  0.001,
    # mass/L (aqueous; approximate 1 g ≈ 1 mL)
    "ng/L":  0.001,
    "ug/L":  1,
    "µg/L":  1,
    "mg/L":  1_000,
    "ai mg/l":  1_000,
    "ug/l":  1,
    "ai mg/kg bdwt": 1_000,
    "ai mg/kg food": 1_000,
}

DURATION_TO_DAYS = {
    "d": 1, "day": 1, "days": 1,
    "h": 1/24, "hr": 1/24, "hour": 1/24, "hours": 1/24,
    "wk": 7, "week": 7, "weeks": 7,
    "mo": 30.4, "month": 30.4, "months": 30.4,
}

SPECIES_GROUP_MAP = {
    # Fish
    "rainbow trout":    "Fish",
    "zebrafish":        "Fish",
    "fathead minnow":   "Fish",
    "common carp":      "Fish",
    "carp":             "Fish",
    "atlantic salmon":  "Fish",
    "salmon":           "Fish",
    "medaka":           "Fish",
    # Mammals
    "mouse":            "Mammal",
    "rat":              "Mammal",
    "rhesus monkey":    "Mammal",
    "monkey":           "Mammal",
    "cynomolgus monkey":"Mammal",
    # Plants
    "lettuce":          "Plant",
    "wheat":            "Plant",
    "corn":             "Plant",
    "maize":            "Plant",
    "rice":             "Plant",
    "soybean":          "Plant",
    # Human biomonitoring
    "human":            "Human",
}

TROPHIC_LEVEL = {
    "Plant":  1,
    "Fish":   3,
    "Mammal": 4,
    "Human":  5,
    "Other":  2,
}

AQUATIC = {"Fish": 1, "Plant": 0, "Mammal": 0, "Human": 0, "Other": 0}


def convert_concentration(row):
    """Convert Result Value to ng/g; return NaN if unknown unit."""
    unit = str(row.get("Result Unit", "")).strip().lower()
    factor = UNIT_TO_NG_G.get(unit, None)
    if factor is None:
        return np.nan
    return row["Result Value"] * factor


def convert_duration(row):
    """Convert exposure duration to days."""
    try:
        val = float(row.get("Exposure Duration", np.nan))
    except (ValueError, TypeError):
        return np.nan
    unit = str(row.get("Duration Unit", "d")).strip().lower()
    factor = DURATION_TO_DAYS.get(unit, None)
    if factor is None:
        return np.nan
    return val * factor


def map_species_group(species_name):
    key = str(species_name).lower().strip()
    for fragment, group in SPECIES_GROUP_MAP.items():
        if fragment in key:
            return group
    return "Other"


def harmonize(df):
    """Apply all cleaning and unit conversion steps."""
    df = df.copy()

    if "Result Value" not in df.columns and "Conc 1 Mean (Standardized)" in df.columns:
        df["Result Value"] = df["Conc 1 Mean (Standardized)"]
    if "Result Unit" not in df.columns and "Conc 1 Units (Standardized)" in df.columns:
        df["Result Unit"] = df["Conc 1 Units (Standardized)"]

    # 1. Numeric coercion
    df["Result Value"] = pd.to_numeric(df["Result Value"], errors="coerce")
    df["Exposure Duration"] = pd.to_numeric(df["Exposure Duration"], errors="coerce")

    # 2. Unit conversion
    df["Concentration_ng_g"] = df.apply(convert_concentration, axis=1)
    df["Duration_days"]      = df.apply(convert_duration, axis=1)

    # 3. Drop rows with no usable concentration
    df = df.dropna(subset=["Concentration_ng_g"])
    df = df[df["Concentration_ng_g"] > 0]

    # 4. Log-transform target
    df["log_concentration"] = np.log10(df["Concentration_ng_g"])

    # 5. Species group
    df["Species_Group"] = df["Species"].apply(map_species_group)

    # 6. Trophic level + aquatic flag
    df["Trophic_Level"] = df["Species_Group"].map(TROPHIC_LEVEL).fillna(2)
    df["Is_Aquatic"]    = df["Species_Group"].map(AQUATIC).fillna(0)

    # 7. Exposure route encoding
    df["Exposure Route"] = df["Exposure Route"].fillna("Unknown")

    return df


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Merge + feature engineering
# ═══════════════════════════════════════════════════════════════════════════════

def build_dataset(ecotox_df, chem_df):
    """Merge ECOTOX observations with CompTox chemical features."""
    def format_casrn(x):
        s = str(int(float(x))) if str(x).replace('.','').isdigit() else str(x)
        if '-' not in s and len(s) >= 5:
            return s[:-3] + '-' + s[-3:-1] + '-' + s[-1]
        return s
    ecotox_df["CASRN"] = ecotox_df["CASRN"].apply(format_casrn)
    chem_df["CASRN"] = chem_df["CASRN"].astype(str)
    merged = ecotox_df.merge(chem_df, on="CASRN", how="left")

    # One-hot encode species group
    for grp in ["Fish", "Mammal", "Plant", "Human"]:
        merged[f"is_{grp.lower()}"] = (merged["Species_Group"] == grp).astype(int)

    # Encode exposure route
    le_route = LabelEncoder()
    merged["Route_encoded"] = le_route.fit_transform(merged["Exposure Route"].fillna("Unknown"))

    return merged, le_route


def build_ml_ready(df):
    """Return X, y, and group labels for leave-one-species-out CV."""
    feature_cols = [
        "Chain_Length",
        "MW",
        "LogKow",
        "Duration_days",
        "Route_encoded",
        "Trophic_Level",
        "Is_Aquatic",
        "is_fish",
        "is_mammal",
        "is_plant",
        "is_human",
    ]
    available = [c for c in feature_cols if c in df.columns]
    core = ["Chain_Length", "MW", "LogKow", "log_concentration", "Species_Group"]
    ml = df[available + ["log_concentration", "Species_Group"]].dropna(subset=core)

    X = ml[available]
    y = ml["log_concentration"]
    groups = ml["Species_Group"]
    return X, y, groups


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Output CSV
# ═══════════════════════════════════════════════════════════════════════════════

OUTPUT_COLS = [
    "PFAS_Name", "CASRN", "PFAS_Class", "Chain_Length", "MW", "LogKow",
    "Species", "Species_Group", "Trophic_Level", "Is_Aquatic",
    "Tissue", "Exposure Route", "Duration_days",
    "Concentration_ng_g", "log_concentration",
    "Study_Year", "DOI", "Source",
]


def save_dataset(df, path):
    available = [c for c in OUTPUT_COLS if c in df.columns]
    df[available].to_csv(path, index=False)
    print(f"[Output] Dataset saved → {path}  ({len(df):,} rows)")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Data gap analysis
# ═══════════════════════════════════════════════════════════════════════════════

def plot_gap_heatmap(df, path):
    """Heatmap: PFAS (rows) × Species Group (cols), count of observations."""
    pfas_col = "PFAS_Name" if "PFAS_Name" in df.columns else "Chemical Name"
    if pfas_col not in df.columns or "Species_Group" not in df.columns:
        print("[Gap] Required columns missing — skipping heatmap.")
        return

    gap = (df.groupby([pfas_col, "Species_Group"])
             .size()
             .unstack(fill_value=0))

    fig, ax = plt.subplots(figsize=(max(8, gap.shape[1] * 1.5),
                                    max(5, gap.shape[0] * 0.6)))
    sns.heatmap(
        gap, annot=True, fmt="d", cmap="YlOrRd_r",
        linewidths=0.5, ax=ax,
        cbar_kws={"label": "# observations"},
    )
    ax.set_title("PFAS Bioaccumulation — Data Gaps by Species Group", fontsize=13)
    ax.set_xlabel("Species Group")
    ax.set_ylabel("PFAS")
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] Gap heatmap → {path}")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 7 — Random Forest model + visualisations
# ═══════════════════════════════════════════════════════════════════════════════

def train_rf(X, y):
    rf = RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1)
    rf.fit(X, y)
    return rf


def plot_feature_importance(rf, feature_names, path):
    importances = pd.Series(rf.feature_importances_, index=feature_names).sort_values()
    fig, ax = plt.subplots(figsize=(7, max(4, len(feature_names) * 0.4)))
    importances.plot(kind="barh", color="steelblue", ax=ax)
    ax.set_title("Random Forest — Feature Importances")
    ax.set_xlabel("Importance")
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] Feature importance → {path}")


def plot_predictions(rf, X, y, path):
    y_pred = rf.predict(X)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(y, y_pred, alpha=0.4, edgecolors="steelblue", facecolors="none", s=40)
    lims = [min(y.min(), y_pred.min()) - 0.5, max(y.max(), y_pred.max()) + 0.5]
    ax.plot(lims, lims, "r--", lw=1)
    r2  = r2_score(y, y_pred)
    rmse = mean_squared_error(y, y_pred) ** 0.5
    ax.set_title(f"Predicted vs Actual  (R²={r2:.3f}, RMSE={rmse:.3f})")
    ax.set_xlabel("Actual log₁₀ concentration (ng/g)")
    ax.set_ylabel("Predicted log₁₀ concentration (ng/g)")
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] Predictions plot → {path}")


def plot_cross_species_validation(X, y, groups, path):
    """Leave-one-species-out cross-validation."""
    logo = LeaveOneGroupOut()
    rf   = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1)

    results = {}
    for train_idx, test_idx in logo.split(X, y, groups):
        left_out = groups.iloc[test_idx[0]]
        rf.fit(X.iloc[train_idx], y.iloc[train_idx])
        preds = rf.predict(X.iloc[test_idx])
        r2    = r2_score(y.iloc[test_idx], preds) if len(test_idx) > 1 else np.nan
        rmse  = mean_squared_error(y.iloc[test_idx], preds) ** 0.5
        results[left_out] = {"R²": r2, "RMSE": rmse, "n": len(test_idx)}

    res_df = pd.DataFrame(results).T.sort_values("R²")

    fig, axes = plt.subplots(1, 2, figsize=(10, max(4, len(res_df) * 0.6)))
    res_df["R²"].plot(kind="barh", ax=axes[0], color="steelblue")
    axes[0].axvline(0, color="red", lw=0.8)
    axes[0].set_title("Leave-one-species-out R²")
    axes[0].set_xlabel("R²")

    res_df["RMSE"].plot(kind="barh", ax=axes[1], color="coral")
    axes[1].set_title("Leave-one-species-out RMSE")
    axes[1].set_xlabel("RMSE (log₁₀ ng/g)")

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] Cross-species validation → {path}")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("\n" + "═" * 60)
    print("  PFAS Bioaccumulation Pipeline")
    print("═" * 60 + "\n")

    # 1. Chemical features
    chem_df = load_comptox(COMPTOX_SNAPSHOT)
    print(f"[CompTox] {len(chem_df)} PFAS in chemical feature table\n")

    # 2. ECOTOX data
    raw_ecotox = load_ecotox(ECOTOX_EXPORT_DIR)
    if raw_ecotox.empty:
        print("\n[INFO] No ECOTOX data loaded. "
              "The pipeline still runs with the chemical feature table.\n"
              "       To build the full dataset:\n"
              "         1. Go to https://cfpub.epa.gov/ecotox/\n"
              "         2. Search each PFAS name → export CSV\n"
              "         3. Place CSVs in the folder: {ECOTOX_EXPORT_DIR}\n"
              "         4. Re-run this script.\n")
        return

    print(f"[ECOTOX] Raw rows: {len(raw_ecotox):,}\n")

    # 3. Harmonize
    clean = harmonize(raw_ecotox)
    print(f"[Harmonize] Rows after cleaning: {len(clean):,}\n")

    # 4. Merge + feature engineering
    dataset, _ = build_dataset(clean, chem_df)

    # Propagate PFAS_Name from chem_df if not already in ECOTOX export
    if "PFAS_Name" not in dataset.columns and "PFAS_Name" in chem_df.columns:
        dataset = dataset.merge(chem_df[["CASRN", "PFAS_Name"]], on="CASRN", how="left")

    # 5. Save ML-ready CSV
    save_dataset(dataset, OUTPUT_DIR + "pfas_bioaccumulation_dataset.csv")

    # 6. Gap heatmap
    plot_gap_heatmap(dataset, OUTPUT_DIR + "pfas_gap_heatmap.png")

    # 7. ML
    X, y, groups = build_ml_ready(dataset)
    print(f"[ML] Training set: {len(X):,} rows × {X.shape[1]} features\n")

    if len(X) < 30:
        print("[ML] Too few rows for meaningful modelling — "
              "add more ECOTOX exports and re-run.")
        return

    rf = train_rf(X, y)
    plot_feature_importance(rf, X.columns.tolist(),
                            OUTPUT_DIR + "feature_importance.png")
    plot_predictions(rf, X, y,
                     OUTPUT_DIR + "model_predictions.png")

    unique_groups = groups.unique()
    if len(unique_groups) >= 2:
        plot_cross_species_validation(X, y, groups,
                                      OUTPUT_DIR + "cross_species_validation.png")
    else:
        print("[CV] Need ≥2 species groups for leave-one-out — skipping.")

    print("\n" + "═" * 60)
    print("  Pipeline complete.")
    print("═" * 60)


if __name__ == "__main__":
    main()

#Add main pipeline script
