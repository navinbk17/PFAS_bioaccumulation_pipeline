"""
PFAS Bioaccumulation Data Pipeline v6.0
========================================
Sources: ECOTOX + EPA CompTox + CDC NHANES (2015-2016, 2017-2018)

Changes in v6.0
---------------
- Proper 80/20 train/test split (headline metrics are held-out only)
- Per-species-group metrics + group mean baseline
- Fish-only model (chemistry features only)
- Plant-only model (chemistry features only)
- Chain length vs BCF scatter plot
- NHANES time trend plot
- In-sample metrics removed from summary

Outputs
-------
  pfas_bioaccumulation_dataset.csv
  pfas_gap_heatmap.png
  feature_importance.png
  model_predictions.png
  cross_species_validation.png
  per_group_metrics.png
  bcf_feature_importance.png
  bcf_predictions.png
  linear_coefficients.png
  fish_feature_importance.png
  fish_predictions.png
  plant_feature_importance.png
  chain_length_bcf_scatter.png
  nhanes_time_trend.png
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
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import LeaveOneGroupOut, train_test_split
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
CHEM_ONLY_FEATURES = ["Chain_Length", "MW", "LogKow", "PFAS_Class_encoded", "Duration_days"]


def load_comptox(path):
    try:
        pd.read_csv(path, low_memory=False)
        print(f"[CompTox] Loaded from {path}")
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


def train_rf(X, y):
    rf = RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1)
    rf.fit(X, y)
    return rf


def plot_importance(importances, path, title, color="steelblue"):
    importances = importances.sort_values()
    fig, ax = plt.subplots(figsize=(7, max(4, len(importances) * 0.45)))
    importances.plot(kind="barh", color=color, ax=ax)
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}")


def plot_predictions_test(y_test, y_pred, path, r2, rmse, title_prefix="RF"):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(y_test, y_pred, alpha=0.4, edgecolors="steelblue", facecolors="none", s=40)
    lims = [min(y_test.min(), y_pred.min()) - 0.5, max(y_test.max(), y_pred.max()) + 0.5]
    ax.plot(lims, lims, "r--", lw=1)
    ax.set_title(f"{title_prefix} — Held-out Test Set\n(R²={r2:.3f}, RMSE={rmse:.3f})")
    ax.set_xlabel("Actual log10 concentration (ng/g)")
    ax.set_ylabel("Predicted log10 concentration (ng/g)")
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}  R²={r2:.3f}")


def run_group_model(dataset, group_name, color, out_prefix):
    """Train a chemistry-only model on a single species group."""
    print(f"\n[{group_name} Model] Training...")
    sub = dataset[dataset["Species_Group"] == group_name].copy()
    avail = [c for c in CHEM_ONLY_FEATURES if c in sub.columns]
    ml = sub[avail + ["log_concentration", "PFAS_Name"]].dropna(
        subset=["Chain_Length", "MW", "LogKow", "log_concentration"])
    print(f"[{group_name} Model] {len(ml):,} rows across {ml['PFAS_Name'].nunique()} PFAS")

    if len(ml) < 20:
        print(f"[{group_name} Model] Too few rows — skipping.")
        return

    X_g, y_g = ml[avail], ml["log_concentration"]
    X_tr, X_te, y_tr, y_te = train_test_split(X_g, y_g, test_size=0.2, random_state=42)
    rf_g = train_rf(X_tr, y_tr)
    y_pred = rf_g.predict(X_te)
    r2_g   = r2_score(y_te, y_pred)
    rmse_g = mean_squared_error(y_te, y_pred) ** 0.5

    # PFAS mean baseline
    pfas_means = y_tr.groupby(ml.loc[y_tr.index, "PFAS_Name"]).mean()
    y_base = ml.loc[y_te.index, "PFAS_Name"].map(pfas_means).fillna(y_tr.mean())
    r2_base = r2_score(y_te, y_base)

    print(f"[{group_name} Model] Held-out R²={r2_g:.3f}  RMSE={rmse_g:.3f}")
    print(f"[{group_name} Model] PFAS mean baseline R²={r2_base:.3f}")
    print(f"[{group_name} Model] Chemistry gain: {r2_g - r2_base:+.3f}")

    plot_importance(pd.Series(rf_g.feature_importances_, index=avail),
                    OUTPUT_DIR + f"{out_prefix}_feature_importance.png",
                    f"{group_name}-Only Model — Feature Importances\n(chemistry features only)",
                    color=color)
    plot_predictions_test(y_te, y_pred,
                          OUTPUT_DIR + f"{out_prefix}_predictions.png",
                          r2_g, rmse_g, f"{group_name}-Only Model")


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


def plot_per_group_metrics(results_df, path):
    fig, axes = plt.subplots(1, 2, figsize=(12, max(4, len(results_df) * 0.7)))
    x = range(len(results_df))
    width = 0.35
    rf_colors   = ["steelblue" if v >= 0 else "lightcoral" for v in results_df["RF R²"]]
    base_colors = ["darkorange" if v >= 0 else "lightyellow" for v in results_df["Baseline R²"]]
    axes[0].barh([i + width/2 for i in x], results_df["RF R²"],       width, color=rf_colors,   label="Random Forest")
    axes[0].barh([i - width/2 for i in x], results_df["Baseline R²"], width, color=base_colors, label="Group Mean Baseline")
    axes[0].set_yticks(list(x))
    axes[0].set_yticklabels(results_df["Group"])
    axes[0].axvline(0, color="black", lw=0.8)
    axes[0].set_title("Held-out R² by Species Group")
    axes[0].set_xlabel("R²")
    axes[0].legend()
    gain_colors = ["steelblue" if v > 0 else "coral" for v in results_df["Gain"]]
    axes[1].barh(results_df["Group"], results_df["Gain"], color=gain_colors)
    axes[1].axvline(0, color="black", lw=0.8)
    axes[1].set_title("RF Gain over Group Mean Baseline")
    axes[1].set_xlabel("R² gain")
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}")


def plot_cross_species(X, y, groups, path):
    logo = LeaveOneGroupOut()
    rf = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1)
    results = {}
    for train_idx, test_idx in logo.split(X, y, groups):
        left_out = groups.iloc[test_idx[0]]
        rf.fit(X.iloc[train_idx], y.iloc[train_idx])
        preds = rf.predict(X.iloc[test_idx])
        r2   = r2_score(y.iloc[test_idx], preds) if len(test_idx) > 1 else np.nan
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


def plot_chain_length_bcf(df, path):
    bcf = df.dropna(subset=["Chain_Length", "log_BCF", "PFAS_Class"])
    if len(bcf) < 5:
        print("[Chain-BCF] Not enough BCF data — skipping.")
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = {"Sulfonate": "steelblue", "Carboxylate": "coral"}
    for cls, grp in bcf.groupby("PFAS_Class"):
        ax.scatter(grp["Chain_Length"], grp["log_BCF"],
                   color=colors.get(cls, "gray"), alpha=0.6, s=50,
                   label=cls, edgecolors="white", linewidth=0.5)
    z = np.polyfit(bcf["Chain_Length"], bcf["log_BCF"], 1)
    p = np.poly1d(z)
    x_line = np.linspace(bcf["Chain_Length"].min(), bcf["Chain_Length"].max(), 100)
    ax.plot(x_line, p(x_line), "k--", lw=1.5, label=f"Trend (slope={z[0]:.2f})")
    ax.set_xlabel("Fluorinated Carbon Chain Length")
    ax.set_ylabel("log10 BCF")
    ax.set_title("Chain Length vs Bioconcentration Factor by PFAS Class")
    ax.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}")


def plot_nhanes_trend(nhanes_path, path):
    if not os.path.exists(nhanes_path):
        return
    nhanes = pd.read_csv(nhanes_path)
    if "Study_Year" not in nhanes.columns:
        return
    trend = nhanes.groupby(["Study_Year", "PFAS_Name"])["Concentration_ng_g"].median().unstack()
    if trend.shape[0] < 2:
        print("[NHANES Trend] Only one year — skipping trend plot.")
        return
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = ["steelblue","coral","forestgreen","purple","darkorange","teal"]
    for i, pfas in enumerate(trend.columns):
        ax.plot(trend.index, trend[pfas], marker="o", linewidth=2,
                color=colors[i % len(colors)], label=pfas)
    ax.set_xlabel("NHANES Survey Year")
    ax.set_ylabel("Median Blood Serum Concentration (ng/g)")
    ax.set_title("PFAS Trends in Human Blood — CDC NHANES 2015-2018")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.set_xticks(sorted(trend.index.tolist()))
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}")


def main():
    print("\n" + "═" * 60)
    print("  PFAS Bioaccumulation Pipeline v6.0")
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
    avail_out = [c for c in OUTPUT_COLS if c in dataset.columns]
    dataset[avail_out].to_csv(OUTPUT_DIR + "pfas_bioaccumulation_dataset.csv", index=False)
    print(f"[Output] Dataset saved ({len(dataset):,} rows)\n")

    plot_gap_heatmap(dataset, OUTPUT_DIR + "pfas_gap_heatmap.png")
    plot_chain_length_bcf(dataset, OUTPUT_DIR + "chain_length_bcf_scatter.png")
    plot_nhanes_trend(NHANES_PATH, OUTPUT_DIR + "nhanes_time_trend.png")

    # ── MAIN CONCENTRATION MODEL ────────────────────────────────────
    core = ["Chain_Length", "MW", "LogKow", "log_concentration", "Species_Group"]
    avail_conc = [c for c in CONC_FEATURES if c in dataset.columns]
    ml = dataset[avail_conc + ["log_concentration", "Species_Group"]].dropna(subset=core)
    X = ml[avail_conc]
    y = ml["log_concentration"]
    groups = ml["Species_Group"]

    print(f"[RF Conc] Total ML rows: {len(X):,} × {X.shape[1]} features")
    X_train, X_test, y_train, y_test, grp_train, grp_test = train_test_split(
        X, y, groups, test_size=0.2, random_state=42)
    print(f"[Split] Train: {len(X_train):,} | Test: {len(X_test):,}\n")

    rf = train_rf(X_train, y_train)
    y_pred_test = rf.predict(X_test)
    r2_test  = r2_score(y_test, y_pred_test)
    rmse_test = mean_squared_error(y_test, y_pred_test) ** 0.5

    group_means = y_train.groupby(grp_train).mean()
    y_pred_baseline = grp_test.map(group_means).fillna(y_train.mean())
    r2_baseline = r2_score(y_test, y_pred_baseline)

    print(f"[RF]       Held-out R²={r2_test:.3f}  RMSE={rmse_test:.3f}")
    print(f"[Baseline] Group mean R²={r2_baseline:.3f}")
    print(f"[Gain]     RF adds {r2_test - r2_baseline:+.3f} beyond group mean\n")

    print("[Per-group held-out metrics]")
    test_df = pd.DataFrame({
        "y_true": y_test.values, "y_pred_rf": y_pred_test,
        "y_pred_base": y_pred_baseline.values, "group": grp_test.values
    })
    group_results = []
    for grp in sorted(test_df["group"].unique()):
        sub = test_df[test_df["group"] == grp]
        if len(sub) < 5:
            continue
        r2_g   = r2_score(sub["y_true"], sub["y_pred_rf"])
        r2_b_g = r2_score(sub["y_true"], sub["y_pred_base"])
        gain   = r2_g - r2_b_g
        group_results.append({"Group": grp, "n": len(sub), "RF R²": r2_g,
                               "Baseline R²": r2_b_g, "Gain": gain})
        print(f"  {grp:10s} n={len(sub):5d}  RF R²={r2_g:+.3f}  Base R²={r2_b_g:+.3f}  Gain={gain:+.3f}")

    grp_df = pd.DataFrame(group_results).sort_values("RF R²")
    plot_per_group_metrics(grp_df, OUTPUT_DIR + "per_group_metrics.png")
    plot_importance(pd.Series(rf.feature_importances_, index=avail_conc),
                    OUTPUT_DIR + "feature_importance.png",
                    "Random Forest — Feature Importances (Concentration)")
    plot_predictions_test(y_test, y_pred_test, OUTPUT_DIR + "model_predictions.png",
                          r2_test, rmse_test, "Random Forest Concentration")
    if grp_train.nunique() >= 2:
        plot_cross_species(X_train, y_train, grp_train, OUTPUT_DIR + "cross_species_validation.png")

    # ── BCF MODEL ───────────────────────────────────────────────────
    avail_bcf = [c for c in BCF_FEATURES if c in dataset.columns]
    ml_bcf = dataset[avail_bcf + ["log_BCF"]].dropna()
    print(f"\n[RF BCF] {len(ml_bcf):,} rows × {len(avail_bcf)} features")
    r2_bcf = rmse_bcf = None
    if len(ml_bcf) >= 30:
        X_bcf_tr, X_bcf_te, y_bcf_tr, y_bcf_te = train_test_split(
            ml_bcf[avail_bcf], ml_bcf["log_BCF"], test_size=0.2, random_state=42)
        rf_bcf = train_rf(X_bcf_tr, y_bcf_tr)
        y_bcf_pred = rf_bcf.predict(X_bcf_te)
        r2_bcf  = r2_score(y_bcf_te, y_bcf_pred)
        rmse_bcf = mean_squared_error(y_bcf_te, y_bcf_pred) ** 0.5
        print(f"[RF BCF] Held-out R²={r2_bcf:.3f}  RMSE={rmse_bcf:.3f}")
        plot_importance(pd.Series(rf_bcf.feature_importances_, index=avail_bcf),
                        OUTPUT_DIR + "bcf_feature_importance.png",
                        "Random Forest — Feature Importances (BCF)")
        plot_predictions_test(y_bcf_te, y_bcf_pred, OUTPUT_DIR + "bcf_predictions.png",
                              r2_bcf, rmse_bcf, "Random Forest BCF")

    # ── LINEAR REGRESSION BASELINE ──────────────────────────────────
    print(f"\n[Linear] Training baseline...")
    X_lr_tr = X_train.copy().apply(lambda c: c.fillna(c.median())).fillna(0)
    X_lr_te = X_test.copy().apply(lambda c: c.fillna(X_train[c.name].median())).fillna(0)
    scaler = StandardScaler()
    lr = LinearRegression()
    lr.fit(scaler.fit_transform(X_lr_tr), y_train)
    y_lr_pred = lr.predict(scaler.transform(X_lr_te))
    r2_lr   = r2_score(y_test, y_lr_pred)
    rmse_lr = mean_squared_error(y_test, y_lr_pred) ** 0.5
    print(f"[Linear] Held-out R²={r2_lr:.3f}  RMSE={rmse_lr:.3f}")
    coef = pd.Series(lr.coef_, index=avail_conc).sort_values()
    colors = ["coral" if v < 0 else "steelblue" for v in coef]
    fig, ax = plt.subplots(figsize=(7, max(4, len(avail_conc) * 0.4)))
    coef.plot(kind="barh", color=colors, ax=ax)
    ax.set_title("Linear Regression — Standardized Coefficients")
    ax.axvline(0, color="black", lw=0.8)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR + "linear_coefficients.png", dpi=150)
    plt.close()
    print(f"[Output] → {OUTPUT_DIR}linear_coefficients.png")

    # ── PER-GROUP CHEMISTRY MODELS ──────────────────────────────────
    run_group_model(dataset, "Fish",  "teal",        "fish")
    run_group_model(dataset, "Plant", "forestgreen", "plant")

    # ── SUMMARY ─────────────────────────────────────────────────────
    print("\n" + "═" * 60)
    print("  Held-out Test Results (80/20 split)")
    print("═" * 60)
    print(f"  RF (concentration):   R²={r2_test:.3f}  RMSE={rmse_test:.3f}")
    print(f"  Linear (baseline):    R²={r2_lr:.3f}  RMSE={rmse_lr:.3f}")
    if r2_bcf is not None:
        print(f"  RF (BCF):             R²={r2_bcf:.3f}  RMSE={rmse_bcf:.3f}")
    print(f"  Group mean baseline:  R²={r2_baseline:.3f}")
    print(f"  RF gain over baseline:{r2_test - r2_baseline:+.3f}")
    print("═" * 60)
    print("  Pipeline complete — v6.0")
    print("═" * 60 + "\n")


if __name__ == "__main__":
    main()