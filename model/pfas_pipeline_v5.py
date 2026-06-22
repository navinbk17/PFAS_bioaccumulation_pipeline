"""
PFAS Bioaccumulation Data Pipeline v4.0
========================================
Sources: ECOTOX + EPA CompTox + CDC NHANES

Outputs
-------
  pfas_bioaccumulation_dataset.csv
  pfas_gap_heatmap.png
  feature_importance.png
  model_predictions.png
  cross_species_validation.png
  bcf_feature_importance.png
  bcf_predictions.png
  linear_coefficients.png
"""

import glob
import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import mean_squared_error, r2_score

warnings.filterwarnings("ignore")

ECOTOX_EXPORT_DIR = "/Users/navink.admin/Desktop/ecotox_exports/"
COMPTOX_SNAPSHOT  = "/Users/navink.admin/Desktop/comptox_snapshot.csv"
NHANES_PATH       = "/Users/navink.admin/Desktop/nhanes_pfas_processed.csv"
OUTPUT_DIR        = "/Users/navink.admin/Desktop/"

PFAS_FEATURES = pd.DataFrame([
    ("PFOS",   "1763-23-1",  "Sulfonate",    8,  500.1, 5.26),
    ("PFOA",   "335-67-1",   "Carboxylate",  8,  414.1, 5.30),
    ("PFHxS",  "355-46-4",   "Sulfonate",    6,  400.1, 4.14),
    ("PFNA",   "375-95-1",   "Carboxylate",  9,  464.1, 6.05),
    ("PFBS",   "375-73-5",   "Sulfonate",    4,  300.1, 1.82),
    ("PFDA",   "335-76-2",   "Carboxylate", 10,  514.1, 6.83),
    ("PFUnDA", "2058-94-8",  "Carboxylate", 11,  564.1, 7.59),
    ("PFDoDA", "307-55-1",   "Carboxylate", 12,  614.1, 8.35),
    ("PFHpA",  "375-85-9",   "Carboxylate",  7,  364.1, 4.55),
    ("PFHxA",  "307-24-4",   "Carboxylate",  6,  314.1, 3.77),
    ("GenX",   "13252-13-6", "Carboxylate",  6,  330.1, 2.50),
    ("ADONA",  "958445-44-8","Carboxylate",  8,  380.1, 2.80),
    ("F53B",   "73606-19-6", "Sulfonate",    6,  570.1, 4.00),
], columns=["PFAS_Name", "CASRN", "PFAS_Class", "Chain_Length", "MW", "LogKow"])

ECOTOX_COL_ALIASES = {
    "CASRN":             ["CAS Number", "CASRN"],
    "Chemical Name":     ["Chemical Name", "Compound"],
    "Species":           ["Species Common Name", "Species", "Common Name"],
    "Organism Group":    ["Species Group", "Organism Group"],
    "Tissue":            ["Response Site", "Tissue", "Body Part"],
    "Exposure Route":    ["Exposure Type", "Exposure Route", "Route"],
    "Exposure Duration": ["Observed Duration Mean (Days)", "Exposure Duration"],
    "Duration Unit":     ["Observed Duration Units (Days)", "Duration Unit"],
    "Result Value":      ["Conc 1 Mean (Author)", "Conc 1 Mean (Standardized)", "Result Value"],
    "Result Unit":       ["Conc 1 Units (Author)", "Conc 1 Units (Standardized)", "Result Unit"],
    "BCF":               ["BCF 1 Value", "BCF"],
    "Study Year":        ["Publication Year", "Study Year"],
    "DOI":               ["DOI"],
}

UNIT_TO_NG_G = {
    "ng/g": 1, "ug/g": 1_000, "µg/g": 1_000, "mg/kg": 1_000,
    "mg/g": 1_000_000, "pg/g": 0.001, "ng/l": 0.001, "ug/l": 1,
    "µg/l": 1, "mg/l": 1_000, "ai mg/l": 1_000, "ai mg/kg bdwt": 1_000,
    "ai mg/kg food": 1_000, "mg/kg soil": 1_000, "mg/g soil": 1_000_000,
    "ng/g soil": 1, "ug/kg soil": 1, "ug/g egg": 1_000,
    "mg/kg dry soil": 1_000, "mg/kg egg": 1_000, "mg/kg diet": 1_000,
    "ppb dry wt": 1, "ppm": 1_000, "ug/kg dry soil": 1,
    "ng/g dw soil": 1, "ug/ml": 1,
}

DURATION_TO_DAYS = {
    "d": 1, "day": 1, "days": 1, "h": 1/24, "hr": 1/24,
    "wk": 7, "week": 7, "weeks": 7, "mo": 30.4, "month": 30.4,
}

SPECIES_GROUP_MAP = {
    "rainbow trout": "Fish", "zebrafish": "Fish", "fathead minnow": "Fish",
    "common carp": "Fish", "carp": "Fish", "salmon": "Fish", "medaka": "Fish",
    "mouse": "Mammal", "rat": "Mammal", "monkey": "Mammal",
    "lettuce": "Plant", "wheat": "Plant", "corn": "Plant",
    "maize": "Plant", "rice": "Plant", "soybean": "Plant",
    "human": "Human",
}

TROPHIC_LEVEL = {"Plant": 1, "Fish": 3, "Mammal": 4, "Human": 5, "Other": 2}
AQUATIC       = {"Fish": 1, "Plant": 0, "Mammal": 0, "Human": 0, "Other": 0}

CONC_FEATURES = [
    "Chain_Length", "MW", "LogKow", "Trophic_Level", "Is_Aquatic",
    "is_fish", "is_mammal", "is_plant", "is_human", "PFAS_Class_encoded", "Duration_days",
]
BCF_FEATURES = [
    "Chain_Length", "MW", "LogKow", "Trophic_Level", "Is_Aquatic",
    "is_fish", "is_mammal", "is_plant", "is_human", "PFAS_Class_encoded",
]


def load_comptox(path):
    try:
        ct = pd.read_csv(path, low_memory=False)
        print(f"[CompTox] Loaded {len(ct)} records from {path}")
    except FileNotFoundError:
        print(f"[CompTox] '{path}' not found — using curated table.")
    return PFAS_FEATURES


def resolve_columns(df, alias_map):
    rename = {}
    for standard, aliases in alias_map.items():
        for alias in aliases:
            if alias in df.columns and standard not in df.columns:
                rename[alias] = standard
                break
    return df.rename(columns=rename)


def load_ecotox(export_dir):
    files = glob.glob(f"{export_dir}*.csv") + glob.glob(f"{export_dir}*.xlsx")
    if not files:
        print(f"[ECOTOX] No files found in '{export_dir}'.")
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
    keep = [c for c in ECOTOX_COL_ALIASES.keys() if c in combined.columns]
    return combined[keep].copy()


def harmonize(df):
    df = df.copy()
    df["Result Value"] = pd.to_numeric(df["Result Value"], errors="coerce")
    df["Exposure Duration"] = pd.to_numeric(df["Exposure Duration"], errors="coerce")

    def conc(row):
        unit = str(row.get("Result Unit", "")).strip().lower()
        f = UNIT_TO_NG_G.get(unit)
        return row["Result Value"] * f if f else np.nan

    def dur(row):
        try:
            val = float(row.get("Exposure Duration", np.nan))
        except:
            return np.nan
        unit = str(row.get("Duration Unit", "d")).strip().lower()
        f = DURATION_TO_DAYS.get(unit)
        return val * f if f else np.nan

    df["Concentration_ng_g"] = df.apply(conc, axis=1)
    df["Duration_days"] = df.apply(dur, axis=1)
    df = df[df["Concentration_ng_g"].notna() & (df["Concentration_ng_g"] > 0)]
    df["log_concentration"] = np.log10(df["Concentration_ng_g"])

    def sg(name):
        key = str(name).lower().strip()
        for frag, grp in SPECIES_GROUP_MAP.items():
            if frag in key:
                return grp
        return "Other"

    df["Species_Group"] = df["Species"].apply(sg)
    df["Trophic_Level"] = df["Species_Group"].map(TROPHIC_LEVEL).fillna(2)
    df["Is_Aquatic"] = df["Species_Group"].map(AQUATIC).fillna(0)
    if "BCF" in df.columns:
        df["BCF"] = pd.to_numeric(df["BCF"], errors="coerce")
        df["log_BCF"] = df["BCF"].apply(lambda x: np.log10(x) if pd.notna(x) and x > 0 else np.nan)
    df["Exposure Route"] = df["Exposure Route"].fillna("Unknown")
    return df


def format_casrn(x):
    s = str(int(float(x))) if str(x).replace('.', '').isdigit() else str(x)
    if '-' not in s and len(s) >= 5:
        return s[:-3] + '-' + s[-3:-1] + '-' + s[-1]
    return s


def build_dataset(ecotox_df, chem_df):
    ecotox_df = ecotox_df.copy()
    ecotox_df["CASRN"] = ecotox_df["CASRN"].apply(format_casrn)
    chem_df["CASRN"] = chem_df["CASRN"].astype(str)
    merged = ecotox_df.merge(chem_df, on="CASRN", how="left")
    for grp in ["Fish", "Mammal", "Plant", "Human"]:
        merged[f"is_{grp.lower()}"] = (merged["Species_Group"] == grp).astype(int)
    merged["PFAS_Class_encoded"] = merged["PFAS_Class"].map({"Sulfonate": 0, "Carboxylate": 1}).fillna(-1)
    le = LabelEncoder()
    merged["Route_encoded"] = le.fit_transform(merged["Exposure Route"].fillna("Unknown"))
    return merged


def plot_gap_heatmap(df, path):
    pfas_col = "PFAS_Name" if "PFAS_Name" in df.columns else "Chemical Name"
    if pfas_col not in df.columns:
        return
    gap = df.groupby([pfas_col, "Species_Group"]).size().unstack(fill_value=0)
    fig, ax = plt.subplots(figsize=(max(8, gap.shape[1] * 1.5), max(5, gap.shape[0] * 0.6)))
    sns.heatmap(gap, annot=True, fmt="d", cmap="YlOrRd_r", linewidths=0.5, ax=ax,
                cbar_kws={"label": "# observations"})
    ax.set_title("PFAS Bioaccumulation — Data Gaps by Species Group", fontsize=13)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] Gap heatmap → {path}")


def train_rf(X, y):
    rf = RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1)
    rf.fit(X, y)
    return rf


def plot_importance(importances, path, title, color="steelblue"):
    importances = importances.sort_values()
    fig, ax = plt.subplots(figsize=(7, max(4, len(importances) * 0.4)))
    importances.plot(kind="barh", color=color, ax=ax)
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}")


def plot_predictions(model, X, y, path, xlabel, ylabel):
    y_pred = model.predict(X)
    r2 = r2_score(y, y_pred)
    rmse = mean_squared_error(y, y_pred) ** 0.5
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(y, y_pred, alpha=0.4, edgecolors="steelblue", facecolors="none", s=40)
    lims = [min(y.min(), y_pred.min()) - 0.5, max(y.max(), y_pred.max()) + 0.5]
    ax.plot(lims, lims, "r--", lw=1)
    ax.set_title(f"Predicted vs Actual  (R²={r2:.3f}, RMSE={rmse:.3f})")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}  R²={r2:.3f}")
    return r2


def plot_cross_species(X, y, groups, path):
    logo = LeaveOneGroupOut()
    rf = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1)
    results = {}
    for train_idx, test_idx in logo.split(X, y, groups):
        left_out = groups.iloc[test_idx[0]]
        rf.fit(X.iloc[train_idx], y.iloc[train_idx])
        preds = rf.predict(X.iloc[test_idx])
        r2 = r2_score(y.iloc[test_idx], preds) if len(test_idx) > 1 else np.nan
        rmse = mean_squared_error(y.iloc[test_idx], preds) ** 0.5
        results[left_out] = {"R²": r2, "RMSE": rmse}
    res_df = pd.DataFrame(results).T.sort_values("R²")
    fig, axes = plt.subplots(1, 2, figsize=(10, max(4, len(res_df) * 0.6)))
    res_df["R²"].plot(kind="barh", ax=axes[0], color="steelblue")
    axes[0].axvline(0, color="red", lw=0.8)
    axes[0].set_title("Leave-one-species-out R²")
    res_df["RMSE"].plot(kind="barh", ax=axes[1], color="coral")
    axes[1].set_title("Leave-one-species-out RMSE")
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}")


def main():
    print("\n" + "═" * 60)
    print("  PFAS Bioaccumulation Pipeline v5.0")
    print("═" * 60 + "\n")

    chem_df = load_comptox(COMPTOX_SNAPSHOT)
    print(f"[CompTox] {len(chem_df)} PFAS in feature table\n")

    raw_ecotox = load_ecotox(ECOTOX_EXPORT_DIR)
    if raw_ecotox.empty:
        print("[INFO] No ECOTOX data. Add xlsx files and re-run.")
        return
    print(f"[ECOTOX] Raw rows: {len(raw_ecotox):,}\n")

    clean = harmonize(raw_ecotox)
    print(f"[Harmonize] Rows after cleaning: {len(clean):,}\n")

    if os.path.exists(NHANES_PATH):
        nhanes = pd.read_csv(NHANES_PATH)
        clean = pd.concat([clean, nhanes], ignore_index=True)
        print("[NHANES] Added " + str(len(nhanes)) + " human blood serum rows")
        print("[Combined] Total rows: " + str(len(clean)) + "\n")

    dataset = build_dataset(clean, chem_df)
    if "PFAS_Name" not in dataset.columns:
        dataset = dataset.merge(chem_df[["CASRN", "PFAS_Name"]], on="CASRN", how="left")

    OUTPUT_COLS = [
        "PFAS_Name", "CASRN", "PFAS_Class", "Chain_Length", "MW", "LogKow",
        "Species", "Species_Group", "Trophic_Level", "Is_Aquatic",
        "is_fish", "is_mammal", "is_plant", "is_human", "PFAS_Class_encoded",
        "Tissue", "Exposure Route", "Duration_days",
        "Concentration_ng_g", "log_concentration", "BCF", "log_BCF",
    ]
    available = [c for c in OUTPUT_COLS if c in dataset.columns]
    dataset[available].to_csv(OUTPUT_DIR + "pfas_bioaccumulation_dataset.csv", index=False)
    print(f"[Output] Dataset saved → {OUTPUT_DIR}pfas_bioaccumulation_dataset.csv  ({len(dataset):,} rows)")

    plot_gap_heatmap(dataset, OUTPUT_DIR + "pfas_gap_heatmap.png")

    # ── CONCENTRATION MODEL ─────────────────────────────────────────
    core = ["Chain_Length", "MW", "LogKow", "log_concentration", "Species_Group"]
    avail_conc = [c for c in CONC_FEATURES if c in dataset.columns]
    ml = dataset[avail_conc + ["log_concentration", "Species_Group"]].dropna(subset=core)
    X, y, groups = ml[avail_conc], ml["log_concentration"], ml["Species_Group"]
    print(f"\n[RF Conc] Training set: {len(X):,} rows × {X.shape[1]} features")

    rf = train_rf(X, y)
    plot_importance(pd.Series(rf.feature_importances_, index=avail_conc),
                    OUTPUT_DIR + "feature_importance.png",
                    "Random Forest — Feature Importances (Concentration)")
    r2_rf = plot_predictions(rf, X, y, OUTPUT_DIR + "model_predictions.png",
                             "Actual log₁₀ concentration (ng/g)",
                             "Predicted log₁₀ concentration (ng/g)")
    if groups.nunique() >= 2:
        plot_cross_species(X, y, groups, OUTPUT_DIR + "cross_species_validation.png")

    # ── BCF MODEL ───────────────────────────────────────────────────
    avail_bcf = [c for c in BCF_FEATURES if c in dataset.columns]
    ml_bcf = dataset[avail_bcf + ["log_BCF"]].dropna()
    print(f"\n[RF BCF] Training set: {len(ml_bcf):,} rows × {len(avail_bcf)} features")
    if len(ml_bcf) >= 30:
        X_bcf, y_bcf = ml_bcf[avail_bcf], ml_bcf["log_BCF"]
        rf_bcf = train_rf(X_bcf, y_bcf)
        plot_importance(pd.Series(rf_bcf.feature_importances_, index=avail_bcf),
                        OUTPUT_DIR + "bcf_feature_importance.png",
                        "Random Forest — Feature Importances (BCF)")
        r2_bcf = plot_predictions(rf_bcf, X_bcf, y_bcf, OUTPUT_DIR + "bcf_predictions.png",
                                  "Actual log₁₀ BCF", "Predicted log₁₀ BCF")

    # ── LINEAR REGRESSION BASELINE ──────────────────────────────────
    print(f"\n[Linear] Training baseline model...")
    X_lr = X.copy()
    for col in X_lr.columns:
        X_lr[col] = X_lr[col].fillna(X_lr[col].median())
    X_lr = X_lr.fillna(0)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_lr)
    lr = LinearRegression()
    lr.fit(X_scaled, y)
    y_pred_lr = lr.predict(X_scaled)
    r2_lr = r2_score(y, y_pred_lr)
    rmse_lr = mean_squared_error(y, y_pred_lr) ** 0.5
    print(f"[Linear] R²={r2_lr:.3f}  RMSE={rmse_lr:.3f}")

    coef = pd.Series(lr.coef_, index=avail_conc).sort_values()
    colors = ["coral" if v < 0 else "steelblue" for v in coef]
    fig, ax = plt.subplots(figsize=(7, max(4, len(avail_conc) * 0.4)))
    coef.plot(kind="barh", color=colors, ax=ax)
    ax.set_title("Linear Regression — Feature Coefficients (standardized)")
    ax.set_xlabel("Coefficient")
    ax.axvline(0, color="black", lw=0.8)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR + "linear_coefficients.png", dpi=150)
    plt.close()
    print(f"[Output] → {OUTPUT_DIR}linear_coefficients.png")

    print("\n── Model Comparison ──────────────────────────────")
    print(f"  Random Forest (concentration): R²={r2_rf:.3f}")
    print(f"  Linear Regression (baseline):  R²={r2_lr:.3f}")
    if len(ml_bcf) >= 30:
        print(f"  Random Forest (BCF):           R²={r2_bcf:.3f}")

    print("\n" + "═" * 60)
    print("  Pipeline complete.")
    print("═" * 60)


if __name__ == "__main__":
    main()