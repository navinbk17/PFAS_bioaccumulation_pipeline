"""
PFAS Bioaccumulation Data Pipeline v10.0
=========================================
Sources: ECOTOX + EPA CompTox + CDC NHANES (2015-2016, 2017-2018)

Changes in v10.0 — Apparent half-life estimation
---------------------------------------------------
- compute_apparent_half_life(): one-compartment first-order elimination
  estimate per PFAS from the two available NHANES survey-wave medians.
  Closed-form two-point solve (k = ln(C1/C2)/Δt, t½ = ln(2)/k) — no
  curve-fitting library needed since only 2 time points exist.
- Reports an "apparent population half-life," explicitly NOT a clinical/
  longitudinal elimination rate: cross-sectional NHANES medians conflate
  true biological elimination with declining population exposure and
  cohort effects between waves. This caveat is surfaced in the docstring,
  the chart title, and the console output rather than left implicit.
- Published literature half-lives (ATSDR Toxicological Profile for
  Perfluoroalkyls, 2021) included as a reference overlay for sanity-
  checking, not as ground truth the model is scored against.
- Compounds with only one NHANES wave, or with rising (non-decreasing)
  concentration between waves, are explicitly flagged rather than
  silently dropped or forced into a misleading negative half-life.
- New output: nhanes_half_life.png — bar chart of NHANES apparent
  half-life vs. published literature half-life, per PFAS.

Changes in v9.0 — Per-PFAS models
-----------------------------------
- run_per_pfas_models(): trains a SEPARATE RF model for each PFAS with
  enough data (n >= MIN_ROWS_FOR_OWN_MODEL), and falls back to reporting
  the existing global RF model's held-out performance broken down by PFAS
  for compounds too sparse to support their own model.
- Per-PFAS summary table: which compounds are reliably predictable
  (own model, decent R²) vs. data-starved (falls back to global model,
  or insufficient data even for that).
- New output: per_pfas_r2.png — bar chart of R² per PFAS, colour-coded by
  whether it got its own model or used the global model fallback.
- Per-PFAS results feed directly into the simulator: a compound with its
  own well-performing model gets a tighter, compound-specific prediction;
  a data-starved compound flags lower confidence to the user.

Changes in v8.0 — Prediction confidence intervals
--------------------------------------------------
- Residual-calibrated (NOT tree-variance) prediction intervals — tree
  variance was found to under-cover (~2% vs 80% target); fixed with a
  proper 3-way Fit/Calibration/Test split.
- Per-species-group interval calibration via apply_group_intervals().
- New outputs: prediction_intervals.png, interval_coverage.png

Changes in v7.0 — Three high-impact leakage/bias fixes
-------------------------------------------------------
FIX 1 - is_human/is_fish/is_mammal/is_plant removed from CONC_FEATURES
FIX 2 - Duration_days removed from CONC_FEATURES (NHANES null-leakage)
FIX 3 - Stratified train/test split by Species_Group

Changes in v6.x (folded into v7.0)
------------------------------------
- XGBoost concentration + BCF models
- BCF chemistry-only features, per-PFAS baseline, BCF 1/2/3 coalescing
- NHANES 2015-2016 data, fish/plant chemistry models

Outputs
-------
  pfas_bioaccumulation_dataset.csv
  pfas_gap_heatmap.png
  feature_importance.png
  model_predictions.png
  prediction_intervals.png
  interval_coverage.png
  per_pfas_r2.png
  nhanes_half_life.png            <- NEW
  cross_species_validation.png
  per_group_metrics.png
  bcf_feature_importance.png
  bcf_predictions.png
  bcf_xgb_feature_importance.png
  bcf_xgb_predictions.png
  linear_coefficients.png
  xgboost_feature_importance.png
  xgboost_predictions.png
  model_comparison.png
  fish_feature_importance.png
  fish_predictions.png
  plant_feature_importance.png
  plant_predictions.png
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
from xgboost import XGBRegressor

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
    # FIX 1: is_fish/is_mammal/is_plant/is_human removed — encode data source
    #         identity rather than biology; breaks simulator generalisability.
    #         Trophic_Level + Is_Aquatic carry the same biological signal cleanly.
    # FIX 2: Duration_days removed — null for all NHANES rows, acting as an
    #         implicit human-row flag (second leakage path).
    "Chain_Length", "MW", "LogKow", "Trophic_Level", "Is_Aquatic", "PFAS_Class_encoded",
]
BCF_FEATURES = [
    # Chemistry-only: BCF already normalizes for exposure concentration,
    # so species/trophic flags add noise rather than signal.
    # Duration_days intentionally excluded — BCF rows rarely have it populated.
    "Chain_Length", "MW", "LogKow", "PFAS_Class_encoded",
]
CHEM_ONLY_FEATURES = ["Chain_Length", "MW", "LogKow", "PFAS_Class_encoded", "Duration_days"]

# v9.0: minimum rows a PFAS needs to get its own dedicated model.
# Below this, we report the global RF model's performance on that PFAS instead.
MIN_ROWS_FOR_OWN_MODEL = 60

# v10.0: published literature apparent elimination half-lives (years), human
# serum. Sources: ATSDR Toxicological Profile for Perfluoroalkyls (2021);
# occupational cohort studies. Reference overlay only -- see
# compute_apparent_half_life() docstring for the cross-sectional-vs-true-
# elimination caveat.
LITERATURE_HALF_LIFE_YEARS = {
    "PFOA": 3.5, "PFOS": 5.4, "PFHxS": 8.5, "PFNA": 2.5,
    "PFDA": 4.5, "PFUnDA": 4.0, "PFBS": 0.08,
}
NHANES_CYCLE_MIDPOINT = {"2015-2016": 2015.5, "2017-2018": 2017.5}


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
    # Coalesce BCF 1/2/3 Value columns — ECOTOX exports all three but
    # resolve_columns only picks BCF 1 Value as "BCF". Fill nulls from 2 and 3.
    for bcf_col in ["BCF 2 Value", "BCF 3 Value"]:
        if bcf_col in df.columns and "BCF" in df.columns:
            df["BCF"] = df["BCF"].fillna(pd.to_numeric(df[bcf_col], errors="coerce"))
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


def train_xgb(X, y):
    xgb = XGBRegressor(
        n_estimators=300,
        learning_rate=0.05,
        max_depth=6,
        subsample=0.8,
        colsample_bytree=0.8,
        random_state=42,
        n_jobs=-1,
        verbosity=0,
    )
    xgb.fit(X, y)
    return xgb


def predict_with_intervals(rf, X, y_true=None, alpha_80=0.10, alpha_95=0.025):
    """
    Compute point predictions + prediction intervals.

    IMPORTANT (v8.0 fix): Raw RF tree-to-tree variance dramatically
    UNDERSTATES true prediction uncertainty — trees trained on overlapping
    bootstrap samples of the same data tend to agree with each other even
    when they're all wrong, so naive tree-spread intervals can show <5%
    coverage against an 80% target. This function instead calibrates
    interval width from the empirical distribution of held-out residuals
    (point prediction vs actual), which is the standard, properly-calibrated
    approach for tree ensembles without specialized quantile forests.

    Parameters
    ----------
    rf       : fitted RandomForestRegressor
    X        : feature matrix to predict on
    y_true   : true values for X, used ONLY to calibrate residual spread.
               Pass the same held-out test y_true used to generate X's
               predictions — this is what makes intervals honest about
               real-world error, not just model self-agreement.
    alpha_80, alpha_95 : tail probabilities for 80%/95% intervals

    Returns
    -------
    dict with mean, lower/upper at 80%/95%, and interval widths.
    If y_true is provided, intervals are residual-calibrated (correct).
    If y_true is omitted, falls back to tree-variance (NOT recommended —
    will under-cover; kept only for cases with no held-out labels available).
    """
    point_pred = rf.predict(X)

    if y_true is not None:
        # Residual-calibrated intervals (correct approach)
        residuals = np.asarray(y_true) - point_pred
        # Symmetric quantiles of the residual distribution
        r_lo_80, r_hi_80 = np.percentile(residuals, [alpha_80 * 100, (1 - alpha_80) * 100])
        r_lo_95, r_hi_95 = np.percentile(residuals, [alpha_95 * 100, (1 - alpha_95) * 100])
        return {
            "mean":     point_pred,
            "std":      np.full_like(point_pred, residuals.std()),
            "lower_80": point_pred + r_lo_80,
            "upper_80": point_pred + r_hi_80,
            "lower_95": point_pred + r_lo_95,
            "upper_95": point_pred + r_hi_95,
            "width_80": np.full_like(point_pred, r_hi_80 - r_lo_80),
            "width_95": np.full_like(point_pred, r_hi_95 - r_lo_95),
        }

    # Fallback: tree-variance (documented as miscalibrated — under-covers)
    tree_preds = np.array([tree.predict(X) for tree in rf.estimators_])
    return {
        "mean":     tree_preds.mean(axis=0),
        "std":      tree_preds.std(axis=0),
        "lower_80": np.percentile(tree_preds, alpha_80 * 100,       axis=0),
        "upper_80": np.percentile(tree_preds, (1 - alpha_80) * 100, axis=0),
        "lower_95": np.percentile(tree_preds, alpha_95 * 100,       axis=0),
        "upper_95": np.percentile(tree_preds, (1 - alpha_95) * 100, axis=0),
        "width_80": np.percentile(tree_preds, (1 - alpha_80) * 100, axis=0)
                  - np.percentile(tree_preds, alpha_80 * 100,       axis=0),
        "width_95": np.percentile(tree_preds, (1 - alpha_95) * 100, axis=0)
                  - np.percentile(tree_preds, alpha_95 * 100,       axis=0),
    }


def predict_with_intervals_per_group(rf, X_cal, y_cal, groups_cal, alpha_80=0.10, alpha_95=0.025):
    """
    Calibrate residual spread per species group on a CALIBRATION set that
    the model has not seen during training. Returns only the calibration
    dict (lower/upper residual offsets per group) — apply it to new
    predictions with apply_group_intervals().

    Using a disjoint calibration fold (not the training set, not the final
    test set) is what makes the resulting coverage numbers honest rather
    than circular.
    """
    point_pred = rf.predict(X_cal)
    residuals_all = np.asarray(y_cal) - point_pred
    grp_arr = np.asarray(groups_cal)

    group_calibration = {}
    for grp in np.unique(grp_arr):
        mask = grp_arr == grp
        n = mask.sum()
        # Too few points to calibrate reliably — fall back to global residuals
        res = residuals_all[mask] if n >= 10 else residuals_all
        r_lo_80, r_hi_80 = np.percentile(res, [alpha_80 * 100, (1 - alpha_80) * 100])
        r_lo_95, r_hi_95 = np.percentile(res, [alpha_95 * 100, (1 - alpha_95) * 100])
        group_calibration[grp] = {
            "lower_80": r_lo_80, "upper_80": r_hi_80,
            "lower_95": r_lo_95, "upper_95": r_hi_95, "n_calibration": n,
        }
    # Fallback entry for any group seen at test/predict time but absent in calibration
    res_all = residuals_all
    r_lo_80, r_hi_80 = np.percentile(res_all, [alpha_80 * 100, (1 - alpha_80) * 100])
    r_lo_95, r_hi_95 = np.percentile(res_all, [alpha_95 * 100, (1 - alpha_95) * 100])
    group_calibration["_default"] = {
        "lower_80": r_lo_80, "upper_80": r_hi_80,
        "lower_95": r_lo_95, "upper_95": r_hi_95, "n_calibration": len(res_all),
    }
    return group_calibration


def apply_group_intervals(rf, X, groups, group_calibration):
    """
    Apply pre-calibrated per-group residual offsets to new predictions.
    This is the function the simulator calls at inference time: given a
    new PFAS's features and its species group, return point estimate +
    calibrated 80%/95% interval using that group's calibration (falling
    back to '_default' if the group wasn't seen during calibration).
    """
    point_pred = rf.predict(X)
    grp_arr = np.asarray(groups)

    lower_80 = np.zeros_like(point_pred)
    upper_80 = np.zeros_like(point_pred)
    lower_95 = np.zeros_like(point_pred)
    upper_95 = np.zeros_like(point_pred)
    width_80 = np.zeros_like(point_pred)
    width_95 = np.zeros_like(point_pred)

    for i, grp in enumerate(grp_arr):
        cal = group_calibration.get(grp, group_calibration["_default"])
        lower_80[i] = point_pred[i] + cal["lower_80"]
        upper_80[i] = point_pred[i] + cal["upper_80"]
        lower_95[i] = point_pred[i] + cal["lower_95"]
        upper_95[i] = point_pred[i] + cal["upper_95"]
        width_80[i] = cal["upper_80"] - cal["lower_80"]
        width_95[i] = cal["upper_95"] - cal["lower_95"]

    return {
        "mean": point_pred, "lower_80": lower_80, "upper_80": upper_80,
        "lower_95": lower_95, "upper_95": upper_95,
        "width_80": width_80, "width_95": width_95,
    }


def plot_prediction_intervals(y_test, intervals, groups, path):
    """
    Plot predicted vs actual with 80% and 95% shaded intervals.
    Points are coloured by species group. Sorted by predicted value
    so intervals are readable as a ribbon chart.
    """
    mean  = intervals["mean"]
    lo80  = intervals["lower_80"]
    hi80  = intervals["upper_80"]
    lo95  = intervals["lower_95"]
    hi95  = intervals["upper_95"]

    # Sort by predicted mean for ribbon readability
    order  = np.argsort(mean)
    x_axis = np.arange(len(order))

    group_colors = {
        "Human": "steelblue", "Fish": "teal",
        "Plant": "forestgreen", "Mammal": "purple", "Other": "gray"
    }
    grp_arr = np.array(groups)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # ── Left panel: ribbon chart ────────────────────────────────────
    ax = axes[0]
    ax.fill_between(x_axis, lo95[order], hi95[order],
                    alpha=0.20, color="steelblue", label="95% interval")
    ax.fill_between(x_axis, lo80[order], hi80[order],
                    alpha=0.35, color="steelblue", label="80% interval")
    ax.plot(x_axis, mean[order], color="steelblue", lw=1.2, label="Predicted mean")

    for grp, col in group_colors.items():
        mask = grp_arr[order] == grp
        if mask.any():
            ax.scatter(x_axis[mask],
                       np.array(y_test.values)[order][mask],
                       color=col, s=12, alpha=0.6, label=grp, zorder=3)

    ax.set_xlabel("Samples (sorted by predicted value)")
    ax.set_ylabel("log10 concentration (ng/g)")
    ax.set_title("Prediction Intervals — RF Concentration Model\n(sorted by predicted mean)")
    ax.legend(fontsize=8, ncol=2)

    # ── Right panel: interval width by species group ─────────────────
    ax2 = axes[1]
    width_df = pd.DataFrame({
        "group": grp_arr,
        "width_80": intervals["width_80"],
        "width_95": intervals["width_95"],
    })
    grp_width = width_df.groupby("group")[["width_80", "width_95"]].mean().sort_values("width_95")
    x2 = np.arange(len(grp_width))
    ax2.barh([i + 0.2 for i in x2], grp_width["width_95"], 0.35,
             color="steelblue", alpha=0.7, label="95% interval width")
    ax2.barh([i - 0.2 for i in x2], grp_width["width_80"], 0.35,
             color="coral", alpha=0.7, label="80% interval width")
    ax2.set_yticks(list(x2))
    ax2.set_yticklabels(grp_width.index)
    ax2.set_xlabel("Mean interval width (log10 ng/g)")
    ax2.set_title("Mean Prediction Interval Width\nby Species Group")
    ax2.legend()

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}")


def plot_interval_coverage(y_test, intervals, groups, path):
    """
    Coverage diagnostic: what fraction of true values fall within each
    interval? A well-calibrated 80% interval should cover ~80% of points.
    Plots overall and per-species-group coverage.
    """
    y_arr   = y_test.values
    in_80   = (y_arr >= intervals["lower_80"]) & (y_arr <= intervals["upper_80"])
    in_95   = (y_arr >= intervals["lower_95"]) & (y_arr <= intervals["upper_95"])
    grp_arr = np.array(groups)

    rows = []
    for grp in sorted(set(grp_arr)):
        mask = grp_arr == grp
        if mask.sum() < 5:
            continue
        rows.append({
            "Group":        grp,
            "n":            mask.sum(),
            "Coverage 80%": in_80[mask].mean() * 100,
            "Coverage 95%": in_95[mask].mean() * 100,
        })
    rows.append({
        "Group":        "Overall",
        "n":            len(y_arr),
        "Coverage 80%": in_80.mean() * 100,
        "Coverage 95%": in_95.mean() * 100,
    })
    cov_df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(9, max(4, len(cov_df) * 0.7)))
    x = np.arange(len(cov_df))
    w = 0.35
    ax.barh([i + w/2 for i in x], cov_df["Coverage 95%"], w,
            color="steelblue", alpha=0.8, label="95% interval coverage")
    ax.barh([i - w/2 for i in x], cov_df["Coverage 80%"], w,
            color="coral",     alpha=0.8, label="80% interval coverage")
    ax.axvline(95, color="steelblue", lw=1.2, ls="--", alpha=0.6, label="Target 95%")
    ax.axvline(80, color="coral",     lw=1.2, ls="--", alpha=0.6, label="Target 80%")
    ax.set_yticks(list(x))
    ax.set_yticklabels(
        [f"{r['Group']} (n={r['n']})" for _, r in cov_df.iterrows()]
    )
    ax.set_xlabel("Coverage (%)")
    ax.set_title("Prediction Interval Coverage by Species Group\n"
                 "(dashed = target; bars > target = conservative intervals)")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 110)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}")

    # Print coverage table
    print("\n[Coverage] Prediction interval calibration:")
    print(f"  {'Group':<12} {'n':>6}  {'80% cov':>8}  {'95% cov':>8}")
    print(f"  {'-'*12}  {'-'*6}  {'-'*8}  {'-'*8}")
    for _, row in cov_df.iterrows():
        flag_80 = "✓" if abs(row["Coverage 80%"] - 80) < 10 else "!"
        flag_95 = "✓" if abs(row["Coverage 95%"] - 95) < 5  else "!"
        print(f"  {row['Group']:<12} {int(row['n']):>6}  "
              f"{row['Coverage 80%']:>7.1f}% {flag_80}  "
              f"{row['Coverage 95%']:>7.1f}% {flag_95}")
    return cov_df


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


def compute_apparent_half_life(nhanes_df, pfas_col="PFAS_Name",
                                year_col="Study_Year",
                                conc_col="Concentration_ng_g"):
    """
    Apparent population half-life per PFAS from NHANES wave medians.

    CAVEAT: cross-sectional, not longitudinal. Conflates true elimination
    with declining population exposure and cohort effects between survey
    waves. Treat as a literature sanity-check, not a clinical/regulatory
    elimination rate. With only 2 time points, k is solved in closed form
    and no uncertainty interval is derivable from the data alone.
    """
    waves = sorted(nhanes_df[year_col].dropna().unique())
    rows = []
    if len(waves) < 2:
        print(f"[Half-Life] Only {len(waves)} NHANES wave(s) present — need >=2. Skipping.")
        return pd.DataFrame(rows)

    wave1, wave2 = waves[0], waves[-1]
    t1 = NHANES_CYCLE_MIDPOINT.get(str(wave1))
    t2 = NHANES_CYCLE_MIDPOINT.get(str(wave2))
    if t1 is None or t2 is None:
        print(f"[Half-Life] Unrecognized wave labels ({wave1}, {wave2}); assuming 2-year spacing.")
        t1, t2 = 0.0, 2.0
    dt = t2 - t1

    med = nhanes_df.groupby([pfas_col, year_col])[conc_col].median().unstack()

    for pfas in med.index:
        c1 = med.loc[pfas].get(wave1, np.nan)
        c2 = med.loc[pfas].get(wave2, np.nan)
        lit = LITERATURE_HALF_LIFE_YEARS.get(pfas, np.nan)

        if pd.isna(c1) or pd.isna(c2) or c1 <= 0 or c2 <= 0:
            rows.append({"PFAS_Name": pfas, "C_wave1": c1, "C_wave2": c2, "dt_years": dt,
                        "k_per_year": np.nan, "half_life_years": np.nan,
                        "literature_half_life_years": lit, "pct_diff_vs_literature": np.nan,
                        "status": "insufficient_temporal_data"})
            continue

        k = np.log(c1 / c2) / dt
        if k <= 0:
            rows.append({"PFAS_Name": pfas, "C_wave1": c1, "C_wave2": c2, "dt_years": dt,
                        "k_per_year": k, "half_life_years": np.nan,
                        "literature_half_life_years": lit, "pct_diff_vs_literature": np.nan,
                        "status": "non_decreasing_between_waves"})
            continue

        t_half = np.log(2) / k
        pct_diff = (t_half - lit) / lit * 100 if pd.notna(lit) else np.nan
        rows.append({"PFAS_Name": pfas, "C_wave1": c1, "C_wave2": c2, "dt_years": dt,
                    "k_per_year": k, "half_life_years": t_half,
                    "literature_half_life_years": lit, "pct_diff_vs_literature": pct_diff,
                    "status": "estimated"})

    return pd.DataFrame(rows)


def plot_half_life_comparison(hl_df, path):
    """Bar chart: NHANES apparent half-life vs. published literature half-life per PFAS."""
    plot_df = hl_df[hl_df["status"] == "estimated"].copy()
    flagged = hl_df[hl_df["status"] != "estimated"]

    if plot_df.empty:
        print("[Half-Life] No PFAS had a valid 2-wave decay estimate — skipping plot.")
        if not flagged.empty:
            print(f"[Half-Life] Flagged (no estimate): {', '.join(flagged['PFAS_Name'].tolist())}")
        return

    plot_df = plot_df.sort_values("half_life_years")
    y = np.arange(len(plot_df))
    fig, ax = plt.subplots(figsize=(8, max(4, len(plot_df) * 0.6)))
    ax.barh(y - 0.18, plot_df["half_life_years"], height=0.35, color="steelblue",
            label="NHANES apparent half-life\n(cross-sectional, 2 waves)")
    ax.barh(y + 0.18, plot_df["literature_half_life_years"], height=0.35, color="coral",
            label="Published literature half-life\n(longitudinal cohort studies)")
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["PFAS_Name"])
    ax.set_xlabel("Half-life (years)")
    ax.set_title("Apparent Population Half-Life vs. Literature\n"
                "(NHANES estimate is cross-sectional, not a clinical elimination rate)")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}")

    print("\n[Half-Life] Apparent population half-life vs. literature:")
    print(f"  {'PFAS':<8} {'NHANES t½':>10}  {'Lit. t½':>9}  {'%diff':>8}  status")
    print(f"  {'-'*8}  {'-'*10}  {'-'*9}  {'-'*8}")
    for _, r in hl_df.iterrows():
        if r["status"] == "estimated":
            lit_str = f"{r['literature_half_life_years']:.2f}y" if pd.notna(r['literature_half_life_years']) else "n/a"
            diff_str = f"{r['pct_diff_vs_literature']:+.1f}%" if pd.notna(r['pct_diff_vs_literature']) else "n/a"
            print(f"  {r['PFAS_Name']:<8} {r['half_life_years']:>9.2f}y  {lit_str:>9}  {diff_str:>8}  estimated")
        else:
            print(f"  {r['PFAS_Name']:<8} {'—':>10}  {'—':>9}  {'—':>8}  {r['status']}")
    if not flagged.empty:
        print(f"\n[Half-Life] Note: {len(flagged)} compound(s) could not be estimated "
              f"({', '.join(flagged['PFAS_Name'].tolist())}).")


def run_per_pfas_models(dataset, avail_conc, global_rf, X_test_global, y_test_global,
                        grp_test_global, pfas_test_global, min_rows=MIN_ROWS_FOR_OWN_MODEL):
    """
    For each PFAS compound, either:
      (a) train a dedicated RF model if it has enough rows (own model), or
      (b) fall back to reporting the GLOBAL RF model's held-out performance
          restricted to that PFAS's test rows (global fallback).

    This directly answers: "which compounds can we predict reliably, and
    which are too data-starved to trust?" — the question the simulator
    needs answered before showing a user a prediction for a specific PFAS.

    Returns a results DataFrame with columns:
      PFAS_Name, n_total, n_test, Model_Type ('own'/'global_fallback'/'insufficient'),
      R2, RMSE
    """
    print("\n[Per-PFAS Models] Evaluating predictability by compound...")
    core_cols = ["Chain_Length", "MW", "LogKow", "log_concentration"]
    results = []

    # v9.0 fix: iterate over the FULL curated PFAS table, not just compounds
    # present in the merged dataset, so zero-data compounds are explicitly
    # reported as a data gap rather than silently omitted.
    all_pfas_names = sorted(PFAS_FEATURES["PFAS_Name"].unique())

    for pfas_name in all_pfas_names:
        sub = dataset[dataset["PFAS_Name"] == pfas_name]
        avail = [c for c in avail_conc if c in sub.columns]
        ml = sub[avail + ["log_concentration"]].dropna(subset=core_cols + avail) if len(sub) else sub
        n_total = len(ml)

        if n_total == 0:
            results.append({"PFAS_Name": pfas_name, "n_total": 0, "n_test": 0,
                            "Model_Type": "no_data", "R2": np.nan, "RMSE": np.nan})
            print(f"  {pfas_name:10s} n=    0  [no data in dataset — zero ECOTOX/NHANES rows]")
            continue

        if n_total >= min_rows:
            # Own dedicated model — has enough data to train and test independently
            X_p, y_p = ml[avail], ml["log_concentration"]
            X_p_tr, X_p_te, y_p_tr, y_p_te = train_test_split(
                X_p, y_p, test_size=0.2, random_state=42)
            model_p = train_rf(X_p_tr, y_p_tr)
            y_p_pred = model_p.predict(X_p_te)
            r2_p   = r2_score(y_p_te, y_p_pred)
            rmse_p = mean_squared_error(y_p_te, y_p_pred) ** 0.5
            results.append({"PFAS_Name": pfas_name, "n_total": n_total, "n_test": len(X_p_te),
                            "Model_Type": "own", "R2": r2_p, "RMSE": rmse_p})
            print(f"  {pfas_name:10s} n={n_total:5d}  [own model]      "
                  f"R²={r2_p:+.3f}  RMSE={rmse_p:.3f}")
        else:
            # Fall back to global model's performance on this PFAS's test rows
            mask = (pfas_test_global == pfas_name)
            n_test = mask.sum()
            if n_test >= 5:
                y_sub_true = y_test_global[mask]
                y_sub_pred = global_rf.predict(X_test_global[mask])
                r2_p   = r2_score(y_sub_true, y_sub_pred)
                rmse_p = mean_squared_error(y_sub_true, y_sub_pred) ** 0.5
                results.append({"PFAS_Name": pfas_name, "n_total": n_total, "n_test": int(n_test),
                                "Model_Type": "global_fallback", "R2": r2_p, "RMSE": rmse_p})
                print(f"  {pfas_name:10s} n={n_total:5d}  [global fallback] "
                      f"R²={r2_p:+.3f}  RMSE={rmse_p:.3f}  (n_test={n_test})")
            else:
                results.append({"PFAS_Name": pfas_name, "n_total": n_total, "n_test": int(n_test),
                                "Model_Type": "insufficient", "R2": np.nan, "RMSE": np.nan})
                print(f"  {pfas_name:10s} n={n_total:5d}  [insufficient data — no reliable estimate]")

    return pd.DataFrame(results)


def plot_per_pfas_r2(results_df, path):
    """
    Bar chart of held-out R² per PFAS, colour-coded by whether the
    compound got its own dedicated model or fell back to the global
    model's performance on its subset.
    """
    plot_df = results_df.dropna(subset=["R2"]).sort_values("R2")
    if plot_df.empty:
        print("[Per-PFAS] No PFAS had enough data to compute R² — skipping plot.")
        return

    colors = plot_df["Model_Type"].map({
        "own": "steelblue", "global_fallback": "coral"
    }).fillna("gray")

    fig, ax = plt.subplots(figsize=(8, max(4, len(plot_df) * 0.5)))
    bars = ax.barh(plot_df["PFAS_Name"], plot_df["R2"], color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Held-out R²")
    ax.set_title("Per-PFAS Predictability\n(blue = dedicated model, coral = global model fallback)")

    for i, (_, row) in enumerate(plot_df.iterrows()):
        ax.text(row["R2"] + (0.01 if row["R2"] >= 0 else -0.01), i,
                f"n={int(row['n_total'])}", va="center",
                ha="left" if row["R2"] >= 0 else "right", fontsize=8)

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}")


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


def plot_model_comparison(results, path):
    """Bar chart comparing all models side by side on R² and RMSE."""
    names = [r["Model"] for r in results]
    r2s   = [r["R²"]   for r in results]
    rmses = [r["RMSE"] for r in results]

    fig, axes = plt.subplots(1, 2, figsize=(12, max(4, len(names) * 0.6)))
    colors_r2   = ["steelblue" if v >= 0 else "lightcoral" for v in r2s]
    colors_rmse = ["coral"] * len(rmses)

    axes[0].barh(names, r2s, color=colors_r2)
    axes[0].axvline(0, color="black", lw=0.8)
    axes[0].set_title("Model Comparison — Held-out R²")
    axes[0].set_xlabel("R²")
    for i, v in enumerate(r2s):
        axes[0].text(max(v, 0) + 0.005, i, f"{v:.3f}", va="center", fontsize=9)

    axes[1].barh(names, rmses, color=colors_rmse)
    axes[1].set_title("Model Comparison — Held-out RMSE")
    axes[1].set_xlabel("RMSE (log10 ng/g)")
    for i, v in enumerate(rmses):
        axes[1].text(v + 0.005, i, f"{v:.3f}", va="center", fontsize=9)

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}")


def main():
    print("\n" + "═" * 60)
    print("  PFAS Bioaccumulation Pipeline v10.0")
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

    # ── APPARENT HALF-LIFE (v10.0) ──────────────────────────────────
    hl_df = pd.DataFrame()
    if os.path.exists(NHANES_PATH):
        nhanes_for_hl = pd.read_csv(NHANES_PATH)
        if "Study_Year" in nhanes_for_hl.columns:
            hl_df = compute_apparent_half_life(nhanes_for_hl)
            if not hl_df.empty:
                plot_half_life_comparison(hl_df, OUTPUT_DIR + "nhanes_half_life.png")
        else:
            print("[Half-Life] NHANES file has no 'Study_Year' column — skipping.")

    # ── SHARED TRAIN/TEST SPLIT ─────────────────────────────────────
    core = ["Chain_Length", "MW", "LogKow", "log_concentration", "Species_Group"]
    avail_conc = [c for c in CONC_FEATURES if c in dataset.columns]
    ml = dataset[avail_conc + ["log_concentration", "Species_Group", "PFAS_Name"]].dropna(subset=core)
    X = ml[avail_conc]
    y = ml["log_concentration"]
    groups = ml["Species_Group"]
    pfas_id = ml["PFAS_Name"]  # v9.0: carried through split for per-PFAS breakdown

    print(f"[Models] Total ML rows: {len(X):,} × {X.shape[1]} features")
    # FIX 3: Stratified split by Species_Group — ensures proportional
    # representation of all groups in train and test regardless of random seed.
    X_train, X_test, y_train, y_test, grp_train, grp_test, pfas_train, pfas_test = train_test_split(
        X, y, groups, pfas_id, test_size=0.2, random_state=42, stratify=groups)
    print(f"[Split]  Train: {len(X_train):,} | Test: {len(X_test):,}")
    print(f"[Split]  Test group counts: {grp_test.value_counts().to_dict()}\n")

    # v8.0: carve a calibration set out of training data for honest interval
    # calibration. Calibrating on the same set you measure coverage on is
    # circular (trivially hits the target) — calibration set must be disjoint
    # from both the model's training fold and the final test fold.
    X_fit, X_cal, y_fit, y_cal, grp_fit, grp_cal = train_test_split(
        X_train, y_train, grp_train, test_size=0.2, random_state=7, stratify=grp_train)
    print(f"[Split]  Fit: {len(X_fit):,} | Calibration: {len(X_cal):,} | Test: {len(X_test):,}\n")

    # Group mean baseline
    group_means = y_train.groupby(grp_train).mean()
    y_pred_baseline = grp_test.map(group_means).fillna(y_train.mean())
    r2_baseline  = r2_score(y_test, y_pred_baseline)
    rmse_baseline = mean_squared_error(y_test, y_pred_baseline) ** 0.5

    # ── RANDOM FOREST ───────────────────────────────────────────────
    print("[RF] Training (on fit fold, holding out calibration fold)...")
    rf = train_rf(X_fit, y_fit)
    y_pred_rf = rf.predict(X_test)
    r2_rf   = r2_score(y_test, y_pred_rf)
    rmse_rf = mean_squared_error(y_test, y_pred_rf) ** 0.5
    print(f"[RF]       Held-out R²={r2_rf:.3f}  RMSE={rmse_rf:.3f}")

    # ── CONFIDENCE INTERVALS ────────────────────────────────────────
    # Calibrate residual spread on X_cal (unseen by rf during training),
    # then APPLY that calibrated width to predictions on X_test. This is
    # the honest 3-way split: fit / calibrate / evaluate are all disjoint.
    print("[CI]  Calibrating prediction intervals on held-out calibration fold...")
    group_cal = predict_with_intervals_per_group(rf, X_cal, y_cal, grp_cal)
    intervals = apply_group_intervals(rf, X_test, grp_test, group_cal)
    mean_width_80 = intervals["width_80"].mean()
    mean_width_95 = intervals["width_95"].mean()
    print(f"[CI]  Mean 80% interval width: {mean_width_80:.3f} log10 ng/g")
    print(f"[CI]  Mean 95% interval width: {mean_width_95:.3f} log10 ng/g")
    print(f"[CI]  Per-group calibration (n used to calibrate each group):")
    for grp, cal in group_cal.items():
        print(f"    {grp:10s} n_cal={cal['n_calibration']:5d}  "
              f"80% width={cal['upper_80']-cal['lower_80']:.3f}  "
              f"95% width={cal['upper_95']-cal['lower_95']:.3f}")
    plot_prediction_intervals(y_test, intervals, grp_test,
                              OUTPUT_DIR + "prediction_intervals.png")
    cov_df = plot_interval_coverage(y_test, intervals, grp_test,
                                    OUTPUT_DIR + "interval_coverage.png")

    plot_importance(pd.Series(rf.feature_importances_, index=avail_conc),
                    OUTPUT_DIR + "feature_importance.png",
                    "Random Forest — Feature Importances (Concentration)")
    plot_predictions_test(y_test, y_pred_rf, OUTPUT_DIR + "model_predictions.png",
                          r2_rf, rmse_rf, "Random Forest Concentration")

    # ── XGBOOST ─────────────────────────────────────────────────────
    print("\n[XGB] Training...")
    X_train_xgb = X_train.fillna(X_train.median())
    X_test_xgb  = X_test.fillna(X_train.median())
    xgb = train_xgb(X_train_xgb, y_train)
    y_pred_xgb = xgb.predict(X_test_xgb)
    r2_xgb   = r2_score(y_test, y_pred_xgb)
    rmse_xgb = mean_squared_error(y_test, y_pred_xgb) ** 0.5
    print(f"[XGB]      Held-out R²={r2_xgb:.3f}  RMSE={rmse_xgb:.3f}")

    xgb_importances = pd.Series(xgb.feature_importances_, index=avail_conc)
    plot_importance(xgb_importances,
                    OUTPUT_DIR + "xgboost_feature_importance.png",
                    "XGBoost — Feature Importances (Concentration)",
                    color="darkorange")
    plot_predictions_test(y_test, y_pred_xgb, OUTPUT_DIR + "xgboost_predictions.png",
                          r2_xgb, rmse_xgb, "XGBoost Concentration")

    # ── LINEAR REGRESSION ───────────────────────────────────────────
    print("\n[Linear] Training...")
    X_lr_tr = X_train.copy().apply(lambda c: c.fillna(c.median())).fillna(0)
    X_lr_te = X_test.copy().apply(lambda c: c.fillna(X_train[c.name].median())).fillna(0)
    scaler = StandardScaler()
    lr = LinearRegression()
    lr.fit(scaler.fit_transform(X_lr_tr), y_train)
    y_lr_pred = lr.predict(scaler.transform(X_lr_te))
    r2_lr   = r2_score(y_test, y_lr_pred)
    rmse_lr = mean_squared_error(y_test, y_lr_pred) ** 0.5
    print(f"[Linear]   Held-out R²={r2_lr:.3f}  RMSE={rmse_lr:.3f}")

    coef = pd.Series(lr.coef_, index=avail_conc).sort_values()
    colors_lr = ["coral" if v < 0 else "steelblue" for v in coef]
    fig, ax = plt.subplots(figsize=(7, max(4, len(avail_conc) * 0.4)))
    coef.plot(kind="barh", color=colors_lr, ax=ax)
    ax.set_title("Linear Regression — Standardized Coefficients")
    ax.axvline(0, color="black", lw=0.8)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR + "linear_coefficients.png", dpi=150)
    plt.close()
    print(f"[Output] → {OUTPUT_DIR}linear_coefficients.png")

    # ── PER-GROUP METRICS ───────────────────────────────────────────
    print("\n[Per-group held-out metrics]")
    test_df = pd.DataFrame({
        "y_true": y_test.values, "y_pred_rf": y_pred_rf,
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
    if grp_train.nunique() >= 2:
        plot_cross_species(X_train, y_train, grp_train, OUTPUT_DIR + "cross_species_validation.png")

    # ── BCF MODELS ──────────────────────────────────────────────────
    avail_bcf = [c for c in BCF_FEATURES if c in dataset.columns]
    ml_bcf = dataset[avail_bcf + ["log_BCF", "PFAS_Name"]].dropna(subset=avail_bcf + ["log_BCF"])
    print(f"\n[BCF Models] {len(ml_bcf):,} rows × {len(avail_bcf)} features")
    print(f"[BCF Models] Features used: {avail_bcf}")
    r2_bcf = rmse_bcf = r2_bcf_xgb = rmse_bcf_xgb = r2_bcf_baseline = rmse_bcf_baseline = None
    if len(ml_bcf) >= 30:
        X_bcf = ml_bcf[avail_bcf]
        y_bcf = ml_bcf["log_BCF"]
        pfas_bcf = ml_bcf["PFAS_Name"]

        X_bcf_tr, X_bcf_te, y_bcf_tr, y_bcf_te, pfas_bcf_tr, pfas_bcf_te = train_test_split(
            X_bcf, y_bcf, pfas_bcf, test_size=0.2, random_state=42)

        # Per-PFAS mean baseline (appropriate for BCF since it's already species-normalized)
        pfas_means_bcf = y_bcf_tr.groupby(pfas_bcf_tr).mean()
        y_bcf_baseline = pfas_bcf_te.map(pfas_means_bcf).fillna(y_bcf_tr.mean())
        r2_bcf_baseline   = r2_score(y_bcf_te, y_bcf_baseline)
        rmse_bcf_baseline = mean_squared_error(y_bcf_te, y_bcf_baseline) ** 0.5
        print(f"[BCF Base] Per-PFAS mean baseline R²={r2_bcf_baseline:.3f}  RMSE={rmse_bcf_baseline:.3f}")

        # RF BCF
        rf_bcf = train_rf(X_bcf_tr, y_bcf_tr)
        y_bcf_pred = rf_bcf.predict(X_bcf_te)
        r2_bcf   = r2_score(y_bcf_te, y_bcf_pred)
        rmse_bcf = mean_squared_error(y_bcf_te, y_bcf_pred) ** 0.5
        print(f"[RF BCF]   Held-out R²={r2_bcf:.3f}  RMSE={rmse_bcf:.3f}  Gain={r2_bcf - r2_bcf_baseline:+.3f}")
        plot_importance(pd.Series(rf_bcf.feature_importances_, index=avail_bcf),
                        OUTPUT_DIR + "bcf_feature_importance.png",
                        "Random Forest — Feature Importances (BCF, chemistry-only)")
        plot_predictions_test(y_bcf_te, y_bcf_pred, OUTPUT_DIR + "bcf_predictions.png",
                              r2_bcf, rmse_bcf, "Random Forest BCF")

        # XGBoost BCF
        xgb_bcf = train_xgb(X_bcf_tr.fillna(X_bcf_tr.median()), y_bcf_tr)
        y_bcf_xgb_pred = xgb_bcf.predict(X_bcf_te.fillna(X_bcf_tr.median()))
        r2_bcf_xgb   = r2_score(y_bcf_te, y_bcf_xgb_pred)
        rmse_bcf_xgb = mean_squared_error(y_bcf_te, y_bcf_xgb_pred) ** 0.5
        print(f"[XGB BCF]  Held-out R²={r2_bcf_xgb:.3f}  RMSE={rmse_bcf_xgb:.3f}  Gain={r2_bcf_xgb - r2_bcf_baseline:+.3f}")
        plot_importance(pd.Series(xgb_bcf.feature_importances_, index=avail_bcf),
                        OUTPUT_DIR + "bcf_xgb_feature_importance.png",
                        "XGBoost — Feature Importances (BCF, chemistry-only)",
                        color="darkorange")
        plot_predictions_test(y_bcf_te, y_bcf_xgb_pred, OUTPUT_DIR + "bcf_xgb_predictions.png",
                              r2_bcf_xgb, rmse_bcf_xgb, "XGBoost BCF")

    # ── PER-GROUP CHEMISTRY MODELS ──────────────────────────────────
    run_group_model(dataset, "Fish",  "teal",        "fish")
    run_group_model(dataset, "Plant", "forestgreen", "plant")

    # ── PER-PFAS MODELS ──────────────────────────────────────────────
    per_pfas_results = run_per_pfas_models(
        dataset, avail_conc, rf, X_test, y_test, grp_test, pfas_test)
    plot_per_pfas_r2(per_pfas_results, OUTPUT_DIR + "per_pfas_r2.png")

    # ── MODEL COMPARISON CHART ──────────────────────────────────────
    comparison_results = [
        {"Model": "Group Mean Baseline", "R²": r2_baseline,  "RMSE": rmse_baseline},
        {"Model": "Linear Regression",   "R²": r2_lr,        "RMSE": rmse_lr},
        {"Model": "Random Forest",        "R²": r2_rf,        "RMSE": rmse_rf},
        {"Model": "XGBoost",              "R²": r2_xgb,       "RMSE": rmse_xgb},
    ]
    if r2_bcf is not None:
        comparison_results += [
            {"Model": "BCF Per-PFAS Baseline", "R²": r2_bcf_baseline, "RMSE": rmse_bcf_baseline},
            {"Model": "RF (BCF)",               "R²": r2_bcf,          "RMSE": rmse_bcf},
            {"Model": "XGBoost (BCF)",          "R²": r2_bcf_xgb,      "RMSE": rmse_bcf_xgb},
        ]
    plot_model_comparison(comparison_results, OUTPUT_DIR + "model_comparison.png")

    # ── SUMMARY ─────────────────────────────────────────────────────
    print("\n" + "═" * 60)
    print("  Held-out Test Results (80/20 stratified split) — v10.0")
    print("═" * 60)
    print(f"  {'Model':<30} {'R²':>7}  {'RMSE':>7}")
    print(f"  {'-'*30}  {'-'*7}  {'-'*7}")
    for r in comparison_results:
        print(f"  {r['Model']:<30} {r['R²']:>7.3f}  {r['RMSE']:>7.3f}")
    print(f"\n  [Concentration]")
    print(f"    XGBoost gain over RF:        {r2_xgb - r2_rf:+.3f}")
    print(f"    RF gain over baseline:       {r2_rf - r2_baseline:+.3f}")
    print(f"    XGBoost gain over baseline:  {r2_xgb - r2_baseline:+.3f}")
    if r2_bcf is not None:
        print(f"\n  [BCF — chemistry-only features]")
        print(f"    RF gain over per-PFAS baseline:      {r2_bcf - r2_bcf_baseline:+.3f}")
        print(f"    XGBoost gain over per-PFAS baseline: {r2_bcf_xgb - r2_bcf_baseline:+.3f}")
    print(f"\n  [Prediction Intervals — RF Concentration]")
    print(f"    Mean 80% interval width: {mean_width_80:.3f} log10 ng/g")
    print(f"    Mean 95% interval width: {mean_width_95:.3f} log10 ng/g")
    overall_cov = cov_df[cov_df["Group"] == "Overall"].iloc[0]
    print(f"    Overall 80% coverage:    {overall_cov['Coverage 80%']:.1f}%  (target: 80%)")
    print(f"    Overall 95% coverage:    {overall_cov['Coverage 95%']:.1f}%  (target: 95%)")
    if not hl_df.empty:
        print(f"\n  [Apparent Half-Life — cross-sectional NHANES 2-wave estimate]")
        estimated = hl_df[hl_df["status"] == "estimated"]
        print(f"    Compounds with a valid 2-wave estimate: {len(estimated)}")
        flagged_hl = hl_df[hl_df["status"] != "estimated"]
        if len(flagged_hl) > 0:
            print(f"    Compounds flagged (insufficient/non-decreasing data): "
                  f"{', '.join(flagged_hl['PFAS_Name'].tolist())}")
        if not estimated.empty:
            closest = estimated.loc[estimated["pct_diff_vs_literature"].abs().idxmin()] \
                if estimated["pct_diff_vs_literature"].notna().any() else None
            if closest is not None:
                print(f"    Closest match to literature: {closest['PFAS_Name']} "
                      f"({closest['pct_diff_vs_literature']:+.1f}% vs. published value)")
    print(f"\n  [Per-PFAS Predictability]")
    own_models = per_pfas_results[per_pfas_results["Model_Type"] == "own"]
    fallback   = per_pfas_results[per_pfas_results["Model_Type"] == "global_fallback"]
    insuff     = per_pfas_results[per_pfas_results["Model_Type"] == "insufficient"]
    no_data    = per_pfas_results[per_pfas_results["Model_Type"] == "no_data"]
    print(f"    Compounds with own dedicated model: {len(own_models)} "
          f"(n>={MIN_ROWS_FOR_OWN_MODEL} rows)")
    print(f"    Compounds using global model fallback: {len(fallback)}")
    print(f"    Compounds with insufficient data: {len(insuff)}")
    print(f"    Compounds with zero data (curated but unmeasured): {len(no_data)}")
    if len(no_data) > 0:
        print(f"      → {', '.join(no_data['PFAS_Name'].tolist())}")
    if not own_models.empty:
        best = own_models.loc[own_models["R2"].idxmax()]
        worst = own_models.loc[own_models["R2"].idxmin()]
        print(f"    Most predictable:  {best['PFAS_Name']} (R²={best['R2']:.3f}, n={int(best['n_total'])})")
        print(f"    Least predictable: {worst['PFAS_Name']} (R²={worst['R2']:.3f}, n={int(worst['n_total'])})")
    print("═" * 60)
    print("  Pipeline complete — v10.0")
    print("═" * 60 + "\n")


if __name__ == "__main__":
    main()
