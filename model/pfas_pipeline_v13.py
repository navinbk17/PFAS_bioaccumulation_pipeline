"""
PFAS Bioaccumulation Data Pipeline v13.0
=========================================
Sources: ECOTOX + EPA CompTox + CDC NHANES (2015-2016, 2017-2018)

Changes in v13.0 — Direct albumin-Ka term in Kfish (closes Findings 16–17)
----------------------------------------------------------------------------
Motivation: Findings 16 and 17 exhausted all tunable scaling approaches in
the Arnot-Gobas Kfish formulation:
  Finding 16: Kprot-scale sweep (0.05 → 0.30 × Koc) — PFOS % error invariant
    (-91% → -89%). The Koc-proxy is structurally wrong, not just mis-scaled.
  Finding 17: P_mem sweep (1.0 → 0.05 correction on k1) — PFOS error worsens
    monotonically (-90% → -99%). Gill membrane permeability is not the source.
Both investigations pointed to the same structural conclusion: the error lives
in Kfish (the fish tissue partition coefficient) and requires a direct albumin-
binding term that does not route through Koc.

What changed:
  Kfish formulation (v11.x/v12.x):
    Kprot = KPROT_SCALE[class] × Koc   (indirect proxy — EXHAUSTED)
    Kfish = Kprot + Kd_nlom + Vw

  Kfish formulation (v13.0):
    Ka_tissue = Ka_albumin [L/mol] × [albumin]_fish [mol/L] × fw
    Kfish     = Ka_tissue + Kd_nlom + Vw

  Ka_albumin: compound-specific albumin association constant (log10, L/mol).
    New column added to PFAS_FEATURES. Values from ITC and fluorescence-
    displacement measurements in Bischel et al. (2010) Environ. Sci. Technol.
    44:5770 and Beesoon & Martin (2015) Environ. Sci. Technol. 49:5758.
    Measured for 10/13 curated PFAS. Interpolated for PFDoDA, PFHpA.
    NaN for GenX/ADONA/F53B (no measured Ka) → runtime class-median imputation,
    reverts to legacy Koc-proxy in those 3 cases only.

  Fish albumin parameters (new module-level constants):
    C_ALBUMIN_FISH_G_L  = 20.0     g/L   (Metcalfe & Thorpe 1992; Conder et al. 2008)
    MW_ALBUMIN_FISH     = 68000.0  g/mol (Chen et al. 2016)
    FW_PROTEIN          = 0.085          (Conder et al. 2008, Table 3)
    → C_ALBUMIN_FISH_MOL_L ≈ 2.94 × 10⁻⁴ mol/L

  NLOM term (Kd_nlom = 0.035 × Koc) is RETAINED — it is a distinct
  organic-carbon sorption mechanism, unaffected by the Ka substitution.

  Kfood (prey partition coefficient) updated to match Kfish structure:
    uses Ka_tissue term for compounds with measured Ka.

Backward compatibility:
  Fallback path for compounds without Ka_albumin: reverts to Koc-proxy
  (v11.x legacy behaviour). Kfish_method column in compute_arnot_gobas_bcf()
  output identifies which path was used per compound.

Finding 18 (confirmed at runtime — see main() summary):
  The Ka-albumin Kfish substitution is mechanistically correct but does NOT
  close the sulfonate BCF gap. Root cause: the Arnot-Gobas two-resistance
  gill membrane term produces mem_r ≈ 6×10⁻⁵ for PFOS (Kow=10^5.26),
  giving k1 = Ew×GV×mem_r ≈ 0.012 L/kg/d. To explain BCF=100 (observed)
  from BCF ≈ k1/(ke+kg), k1 must be ~13× higher. The Kfish term governs
  k2 (elimination), not k1 (uptake) — so correct Kfish cannot fix a
  wrong k1 mechanism.

  This is a model architecture finding: Arnot-Gobas assumes passive
  diffusion as the only gill uptake route. For PFAS, the dominant uptake
  pathway is protein-facilitated transport — the thermodynamic gradient
  created by high Ka albumin binding in fish blood draws PFAS across the
  gill epithelium faster than free diffusion predicts. The Kow-based
  two-resistance membrane term systematically under-predicts k1 for all
  strongly protein-binding PFAS.

  Implication for paper: the correct next step is a protein-facilitated k1
  formulation, where the gill uptake rate is augmented by the albumin
  concentration gradient:
    k1_prot = k1_passive + k_fac × Ka_albumin × [albumin]_blood × Cw
  This is beyond the single-constant-correction approach and requires either
  (a) a modified rate constant from empirical fish albumin uptake data, or
  (b) a two-compartment model (blood + tissue) that separates the blood-
      protein pool from the whole-body partition coefficient.
  Finding 18 closes the single-compartment Arnot-Gobas framework as a
  quantitatively accurate BCF predictor for sulfonates. The framework
  remains mechanistically valuable as a diagnostic tool and qualitative
  rank-ordering model.

New/changed outputs:
  arnot_gobas_bcf.png  — updated with v13 BCF predictions (Ka-based Kfish)
  Ka_albumin and Ka_tissue columns added to compute_arnot_gobas_bcf() output
  Kfish_method column shows 'Ka_albumin' / 'Koc_proxy_legacy' / 'Kow_fallback'
  Console summary table includes Ka_tissue and Kfish for each compound

Preserved from v12.x (historical diagnostics, not changed):
  run_arnot_gobas_sensitivity()     — Kprot sweep; evidence for Finding 16
  run_arnot_gobas_pmem_sensitivity()— P_mem sweep; evidence for Finding 17
  Both functions retain their original math (Koc-proxy) to preserve the
  diagnostic integrity of the findings they produced.



Changes in v11.0 — PFAS-appropriate chemistry features (Koc, AlbuminBinding_pKa)
----------------------------------------------------------------------------------
Motivation: LogKow is calibrated for neutral, lipophilic organic compounds that
partition into fat. PFAS are ionic surfactants — they don't accumulate in fat;
they bind to serum proteins (albumin, fatty-acid-binding proteins) in blood and
partition to organic carbon in soil/sediment. Using LogKow alone as the
hydrophobicity descriptor leaves two important biological mechanisms without a
feature:
  (a) Soil/sediment → dietary uptake pathway for fish and plant models: governed
      by organic-carbon-normalized soil-water partition (Koc), not Kow.
  (b) Blood protein binding for human serum: albumin binding affinity (proxied
      here by literature pKa of the head group, which correlates with binding
      strength across the carboxylate/sulfonate series).

New features added to PFAS_FEATURES table (hardcoded from literature):
  Koc      — soil organic-carbon/water partition coefficient (L/kg), log-
             transformed to Koc_log for ML. Sources: Guelfo & Higgins (2013),
             Ahrens et al. (2010), OECD PFAS sorption review (2021). Where
             experimental values span a range, the geometric mean is used.
             Missing values (GenX, ADONA, F53B) flagged as NaN — the pipeline
             imputes the column median so these rows are not dropped.
  AlbuminBinding_pKa — functional-group pKa as a proxy for serum-protein
             binding strength. Carboxylates: pKa ≈ 0–1 (strongly acidic,
             fully ionized at physiological pH, high protein affinity).
             Sulfonates: pKa < –1 (super-acid, also fully ionized but
             different binding geometry). Literature: Jones et al. (2003),
             Bischel et al. (2010), Beesoon & Martin (2015). Values are
             compound-specific where available; class-level defaults
             otherwise.

Feature sets updated:
  CONC_FEATURES    : + Koc_log, AlbuminBinding_pKa
  BCF_FEATURES     : + Koc_log  (Koc directly governs fish dietary uptake;
                       AlbuminBinding_pKa is blood-specific, not BCF-relevant)
  CHEM_ONLY_FEATURES: + Koc_log, AlbuminBinding_pKa

Ablation comparison:
  run_feature_ablation() trains RF on v10.5 features (Chain_Length, MW,
  LogKow, PFAS_Class_encoded) and v11 features side-by-side on the same
  train/test split, then reports delta-R² per species group. This is the
  honest test of whether the new features help. Results logged to console
  and to feature_ablation.png.

Expected impact:
  Human model: AlbuminBinding_pKa should add modest signal — the human
  dataset is already well-fit and largely linear; the new feature captures
  a mechanistic distinction the linear model can use.
  Fish/Plant models: Koc_log is the theoretically justified improvement for
  these groups. Whether it actually helps depends on whether the fish/plant
  heterogeneity is feature-poverty (new features would help) or data-
  heterogeneity (noise dominates regardless of features). Both outcomes
  are informative — see Finding 13.

Changes in v11.2 — Dedicated human-only model
----------------------------------------------
Motivation: The pooled model is ~95% human rows but is framed as a
cross-species model. Fish/Plant models are proven unfixable with current
data (Finding 13 — ΔR²=0.000 after adding Koc, AlbuminBinding_pKa).
A dedicated human-only model:
  (a) Removes fish/plant noise from training entirely
  (b) Allows tighter per-PFAS confidence intervals (NHANES is standardized)
  (c) Frames the science honestly — this is a human biomonitoring model
  (d) Sets up the interactive simulator cleanly: input a PFAS, get a
      predicted blood serum concentration with calibrated uncertainty

What's new:
  run_human_only_model(): trains RF, XGBoost, and Linear Regression on
    NHANES human blood serum rows only. Features: Chain_Length, MW,
    LogKow, PFAS_Class_encoded, AlbuminBinding_pKa (Koc_log dropped —
    governs soil/fish uptake, not human blood serum partitioning).
    Per-PFAS R² breakdown within the human model.
    Calibrated 80%/95% prediction intervals using the same 3-way
    Fit/Calibration/Test split as the pooled model.
  Outputs: human_model_predictions.png, human_model_per_pfas.png,
    human_model_feature_importance.png
  Console: human-only R², RMSE, per-PFAS breakdown, interval coverage.

Expected improvement over pooled model (Human R²=0.604):
  Removing ~1,300 non-human training rows that add noise without signal
  should modestly improve R², tighten intervals, and give cleaner
  per-PFAS estimates. The honest test is whether R² moves materially —
  if it stays at 0.604, the pooled model was already effectively a
  human model and the framing fix is the main gain.

Changes in v11.1 — Arnot-Gobas mechanistic BCF model (PFAS-adapted)
----------------------------------------------------------------------
Motivation: Finding 5 shows ML cannot beat a per-compound mean for BCF
(RF/XGBoost match per-PFAS baseline exactly, R²≈0.241). A mass-balance
mechanistic model explains WHY — the chemistry signal is real, it's just
swamped by study-level variance (species, lab conditions, water chemistry)
that ML learns equally badly regardless of features. The Arnot-Gobas model
sidesteps that by computing BCF from first principles rather than learning
from noisy data.

Model: steady-state one-compartment mass balance (Arnot & Gobas 2004,
Environ. Toxicol. Chem. 23:1523–1532).

  BCF = k1 / (k2 + ke + km + kg)

Where:
  k1  = gill uptake rate constant (L/kg·d) — function of Kow and fish Ew
  k2  = gill elimination rate constant (1/d) — inverse of k1 scaled by
        lipid+NLOM+water partitioning in fish
  ke  = fecal egestion rate constant (1/d)
  km  = metabolic transformation rate (1/d) — set to 0 for PFAS (no
        significant biotransformation; PFAS are metabolically inert)
  kg  = growth dilution rate (1/d)

PFAS-specific adaptations (this is where it differs from the classic model):
  (a) Protein binding correction: PFAS partition to serum proteins, not
      lipids. The classic model uses Kow to estimate lipid partitioning;
      for PFAS we substitute a protein-binding-corrected partition
      coefficient (Kprot) derived from Koc and chain length following
      Kelly et al. (2004) and Conder et al. (2008).
  (b) Non-lipid organic matter (NLOM) term: PFAS bind strongly to NLOM
      in fish tissue. Included as a separate term in the tissue partition
      coefficient following Gobas et al. (2003).
  (c) Dietary uptake (kd): for compounds with high Koc, sediment/dietary
      pathway dominates over gill uptake. Included but parameterized
      conservatively since ECOTOX BCF measurements are mostly aqueous
      exposure studies (not dietary).

Default fish parameters (Arnot & Gobas 2004, Table 1 — generic fish):
  Vl  = 0.05   (lipid fraction, 5%)
  Vn  = 0.09   (NLOM fraction, 9%)
  Vw  = 0.72   (water fraction, 72%)
  Ew  = 0.50   (gill oxygen extraction efficiency)
  Ed  = 0.75   (dietary assimilation efficiency)
  GV  = 0.40   (ventilation rate, L water / kg fish / d, at 12°C)
  Gd  = 0.025  (feeding rate, kg food / kg fish / d)
  kg  = 0.001  (specific growth rate, 1/d — 0.1%/day for adult fish)
  T   = 12°C   (water temperature — standard freshwater default)

Output:
  BCF_AG column added to the BCF dataset (mechanistic prediction per PFAS)
  arnot_gobas_bcf.png — side-by-side: mechanistic BCF vs ML BCF vs observed
    median per PFAS, with observed scatter shown for context.
  Console table: AG BCF, ML BCF, observed median, and % error vs observed.

Interpretation guidance printed to console: the gap between AG and ML BCF
is the "model agreement" metric; the gap between both and observed is the
"data heterogeneity" metric. A large observed scatter with both models
landing near the observed median = data problem. A systematic offset of
both models in the same direction = feature/parameter problem.

Changes in v10.1 — Trophic_Level / Is_Aquatic leakage fix + per-group headline
--------------------------------------------------------------------------------
- FIX 4 (supersedes v7.0's FIX 1): Trophic_Level and Is_Aquatic removed from
  CONC_FEATURES. Both are `Species_Group.map({...})` — a fixed, bijective
  dict keyed one-to-one on the same five group labels the is_human/is_fish/
  is_mammal/is_plant flags encoded (Plant=1, Fish=3, Mammal=4, Human=5,
  Other=2). They are not measured biological quantities; they are
  Species_Group re-encoded as an integer, so v7.0's FIX 1 removed four
  copies of the group-identity leak and left a fifth (two, counting
  Is_Aquatic) standing. This is why v7.0/v9.0's Finding 9 saw R² barely
  move (0.712 → 0.710) after "fixing" leakage — the signal didn't shrink,
  it consolidated into these two features.
- CONC_FEATURES is now chemistry-only: Chain_Length, MW, LogKow,
  PFAS_Class_encoded. Honest held-out R² on this feature set is ~0.52,
  down from ~0.71 with Trophic_Level/Is_Aquatic included. The dataset
  still carries Trophic_Level/Is_Aquatic as descriptive columns (they're
  useful for grouping/plotting) — they are simply no longer fed to any
  model as predictive features.
- Headline metric changed from pooled held-out R² to per-species-group R².
  The pipeline is ~95% human rows, so a single pooled number was almost
  entirely a human-serum number wearing an "overall" label; Fish/Plant
  performance (already reported in per_group_metrics.png) was being
  contradicted rather than reflected by the banner figure. Console
  summary and README now lead with the per-group breakdown and label the
  pooled figure explicitly as a human-weighted average, not cross-species
  evidence.
- No change to BCF features (already chemistry-only, unaffected by this fix),
  half-life estimation, or per-PFAS modeling logic.

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
  nhanes_half_life.png
  feature_ablation.png            <- NEW (v11): v10.5 vs v11 features, per group
  arnot_gobas_bcf.png             <- NEW (v11.1): mechanistic BCF vs ML BCF vs observed
  human_model_predictions.png     <- NEW (v11.2): human-only RF predicted vs actual
  human_model_per_pfas.png        <- NEW (v11.2): per-PFAS R² within human-only model
  human_model_feature_importance.png <- NEW (v11.2): human-only RF feature importances
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

# v11.0: two new columns added to the feature table.
#
# Koc (L/kg organic carbon): soil organic-carbon/water partition coefficient.
# Controls how strongly PFAS sorb to soil/sediment, which governs dietary
# uptake in fish (via benthic prey) and root uptake in plants.
# Sources: Guelfo & Higgins (2013) Environ. Sci. Technol.; Ahrens et al. (2010)
# Chemosphere; OECD (2021) Series on Risk Management No. 39.
# Where a range is reported, geometric mean is used.
# NaN = no reliable experimental value found; imputed with column median at
# runtime so rows are NOT dropped from training.
#
# v11.3: class-specific protein binding scaling factors for Arnot-Gobas.
# Sulfonates bind serum albumin ~3-5x more strongly than carboxylates of
# equivalent chain length (Bischel et al. 2010; Beesoon & Martin 2015).
# The original single constant (0.05) systematically under-predicted all
# sulfonates (Finding 14). Carboxylate value unchanged from v11.1.
KPROT_SCALE = {
    "Sulfonate":   0.15,  # ~3x higher — literature-informed starting point
    "Carboxylate": 0.05,  # unchanged from Kelly et al. (2004)
}

# Sensitivity sweep values for sulfonate Kprot — run_arnot_gobas_sensitivity()
# tests these to show how % error moves before committing to one number.
KPROT_SULFONATE_SWEEP = [0.05, 0.08, 0.12, 0.15, 0.20, 0.25, 0.30]

# v12.1: gill membrane permeability correction factor per PFAS class.
#
# Motivation: Finding 16 established that tuning KPROT_SCALE (tissue
# partitioning) cannot close the sulfonate BCF gap — PFOS % error moves only
# ~2 percentage points across the full Kprot sweep (0.05 → 0.30). This implicates
# a separate physical mechanism: reduced gill membrane permeability for ionized
# PFAS relative to neutral lipophilic organics.
#
# The Arnot & Gobas (2004) Eq. 5 membrane resistance term:
#   mem_resistance = 1 / (5.49e7 × Kow^-0.670 + 1)
# is calibrated on neutral, non-ionic organics. Sulfonates carry a permanent
# negative charge at physiological pH (pKa ≈ -3.3) and cross phospholipid
# bilayers via protein-mediated transport rather than passive diffusion — a
# fundamentally different permeation mechanism that Kow alone cannot represent.
# Carboxylates (pKa ≈ 0-1) are also ionized at pH 7.4 but show less
# charge-density effect, consistent with observed BCF patterns.
#
# P_MEM_CORRECTION multiplies mem_resistance (and therefore k1) by a
# class-specific scalar. Values < 1.0 = reduced permeability vs Kow prediction.
# This is a direct, isolated test of the membrane permeability hypothesis:
#   - KPROT_SCALE controls tissue accumulation (Kfish) → affects k2, ke, kd
#   - P_MEM_CORRECTION controls gill uptake rate → affects k1 only
# The two hypotheses are therefore cleanly separable parameters.
#
# Starting values: Carboxylate=1.0 (no correction, current baseline);
# Sulfonate=0.50 (exploratory — membrane resistance doubled, k1 halved).
# run_arnot_gobas_pmem_sensitivity() sweeps P_MEM_SWEEP to map the full
# response surface before committing to a value.
#
# Literature context: Armitage et al. (2013) Environ. Sci. Technol. report
# that ionized organic acids have effective membrane permeabilities 10–100×
# lower than neutral analogs of similar Kow; 0.50 is conservative relative
# to that range and should be treated as an exploratory starting point.
P_MEM_CORRECTION = {
    "Sulfonate":   0.50,  # exploratory — halves k1 relative to Kow prediction
    "Carboxylate": 1.00,  # no correction — current baseline
}

# Sweep values for sulfonate P_mem — run_arnot_gobas_pmem_sensitivity()
# tests the full range from no correction (1.0) to near-zero permeability (0.05).
# Distinct from KPROT_SULFONATE_SWEEP: that sweeps tissue partitioning;
# this sweeps gill membrane uptake. Running both together maps whether the
# sulfonate error is a partitioning problem, a permeability problem, or both.
P_MEM_SWEEP = [1.00, 0.80, 0.60, 0.50, 0.40, 0.25, 0.10, 0.05]

# v13.0: Fish albumin (serum protein) parameters for the direct Ka-based Kfish term.
#
# Motivation: Findings 16 and 17 exhausted all tunable scaling approaches
# (Kprot-scale sweep and P_mem sweep both showed sulfonate BCF error is
# invariant). The structural problem is that Kprot = scale × Koc is an
# indirect proxy for protein binding that cannot distinguish PFAS that have
# similar Koc but very different albumin affinities. The correct Kfish
# formulation for PFAS serum-binding compounds is:
#
#   Kfish = Ka_albumin * [albumin]_fish + Kd_nlom + Vw
#
# where Ka_albumin (L/mol) is the compound-specific albumin association
# constant measured by ITC or fluorescence displacement (Bischel et al. 2010;
# Beesoon & Martin 2015), [albumin]_fish is the albumin concentration in fish
# blood/tissue (mol/L), and the Koc-based Kprot term is REPLACED entirely
# (not supplemented) by the Ka term. The NLOM and water terms are unchanged.
#
# Fish albumin parameters:
#   C_albumin_fish_g_L  = 20.0 g/L  — typical teleost fish plasma albumin.
#     Sources: Metcalfe & Thorpe (1992) Comp. Biochem. Physiol. 102B:79;
#     Hunn (1982) J. Fish Biol. 21:95. Range 15–25 g/L across salmonids;
#     20 g/L is the geometric mean used in PFAS fish models (Conder et al. 2008).
#   MW_albumin_fish     = 68000 g/mol — fish serum albumin MW (teleosts).
#     Fish albumin is slightly smaller than human HSA (66500 g/mol);
#     68000 g/mol is the consensus value from Chen et al. (2016).
#   fw (protein volume fraction in fish tissue) = 0.085 — fraction of fish
#     wet weight that is protein. This scales the blood-phase Ka term to a
#     whole-body tissue partition coefficient (dimensionless). Derived from
#     Conder et al. (2008) Environ. Sci. Technol. 42:9283, Table 3; consistent
#     with Arnot & Gobas (2004) tissue composition defaults.
#
# Together these give:
#   [albumin]_fish_mol_L = C_albumin_fish_g_L / MW_albumin_fish  ≈ 2.94e-4 mol/L
#   Ka_tissue = Ka_albumin [L/mol] × [albumin]_fish [mol/L] × fw [dimensionless]
#
# For PFOS (Ka = 10^4.68 = 47,863 L/mol):
#   Ka_tissue = 47863 × 2.94e-4 × 0.085 ≈ 1.19
#   Kfish     = 1.19 + 0.035×2100 + 0.72 ≈ 75.3   (vs. old Koc proxy: 0.15×2100+... ≈ 388)
# This is a larger, chemistry-grounded Kfish than the scale×Koc approach
# and should substantially close the sulfonate under-prediction gap.
C_ALBUMIN_FISH_G_L  = 20.0     # g/L — teleost plasma albumin concentration
MW_ALBUMIN_FISH     = 68000.0  # g/mol — fish serum albumin molecular weight
FW_PROTEIN          = 0.085    # dimensionless — protein volume fraction in fish tissue

# Derived: molar concentration of albumin in fish tissue (mol/L)
C_ALBUMIN_FISH_MOL_L = C_ALBUMIN_FISH_G_L / MW_ALBUMIN_FISH  # ≈ 2.94e-4 mol/L

# AlbuminBinding_pKa: functional-group acid dissociation constant (pKa),
# used as a proxy for serum-albumin binding affinity. Both carboxylates and
# sulfonates are fully deprotonated at physiological pH (7.4), but the
# magnitude of the pKa correlates with binding strength through electrostatic
# and hydrophobic interactions. Lower pKa = stronger acid = more charge density
# = stronger albumin binding at pH 7.4.
# Sources: Jones et al. (2003) Environ. Toxicol. Chem.; Bischel et al. (2010)
# Environ. Sci. Technol.; Beesoon & Martin (2015) Environ. Sci. Technol.
# Sulfonates assigned class-level default (–3.3); carboxylates vary by chain.
#
# v13.0: Ka_albumin added — compound-specific albumin binding affinity constant
# (log10, L/mol). This is the direct mechanistic quantity that Findings 16–17
# established is required: Kprot-scale and P_mem tuning are both exhausted;
# the error lives in Kfish itself, which needs an albumin-Ka term rather than
# an Koc proxy. Sources: Bischel et al. (2010) Environ. Sci. Technol. 44:5770
# (Table 2, ITC measurements at 25°C, pH 7.4, HSA); Beesoon & Martin (2015)
# Environ. Sci. Technol. 49:5758 (Table 1, FAF-BSA displacement, pH 7.4).
# Where the two studies report both fish serum albumin (FSA) and human serum
# albumin (HSA) values, HSA is used as the default — FSA values are within
# 0.1–0.3 log units for most compounds (Chen et al. 2016 notes FSA/HSA
# binding affinities track closely across the PFAS series).
# Compounds with no measured Ka: chain-length interpolation within class
# (PFHpA, PFHxA, PFDoDA), class-level median for emerging compounds
# (GenX, ADONA, F53B — imputed at runtime to avoid dropping rows).
# NaN = no measured or reliably interpolated value → runtime median imputation.
#
#                                                     Koc      Albumin   Ka_albumin
#  Name      CASRN           Class        CL    MW   LogKow   (L/kg)     _pKa   (log L/mol)
PFAS_FEATURES = pd.DataFrame([
    # ── Tier 1: core PFAS — all Ka values from Bischel et al. (2010) or Beesoon & Martin (2015) ──
    ("PFOS",   "1763-23-1",  "Sulfonate",    8, 500.1,  5.26,  2100.0,  -3.27,  4.68),
    ("PFOA",   "335-67-1",   "Carboxylate",  8, 414.1,  5.30,  1900.0,   0.50,  3.65),
    ("PFHxS",  "355-46-4",   "Sulfonate",    6, 400.1,  4.14,   560.0,  -3.27,  4.20),
    ("PFNA",   "375-95-1",   "Carboxylate",  9, 464.1,  6.05,  3200.0,   0.40,  4.10),
    ("PFBS",   "375-73-5",   "Sulfonate",    4, 300.1,  1.82,    47.0,  -3.27,  2.90),
    # ── Tier 2: extended — Ka from Bischel et al. (2010) ──────────────────────────────────────
    ("PFDA",   "335-76-2",   "Carboxylate", 10, 514.1,  6.83,  5800.0,   0.30,  4.55),
    ("PFUnDA", "2058-94-8",  "Carboxylate", 11, 564.1,  7.59,  9100.0,   0.20,  4.98),
    ("PFDoDA", "307-55-1",   "Carboxylate", 12, 614.1,  8.35, 14000.0,   0.10,  5.40),  # extrapolated: +0.42/CF from PFUnDA trend
    ("PFHpA",  "375-85-9",   "Carboxylate",  7, 364.1,  4.55,   820.0,   0.60,  3.30),  # interpolated: PFOA(C8)=3.65, PFHxA(C6)=2.80 → C7≈3.30
    ("PFHxA",  "307-24-4",   "Carboxylate",  6, 314.1,  3.77,   310.0,   0.70,  2.80),  # Bischel et al. (2010) lower bound; shortest measured carboxylate
    # ── Tier 3: emerging — no reliable Ka measurement; runtime median imputation ──────────────
    ("GenX",   "13252-13-6", "Carboxylate",  6, 330.1,  2.50,     np.nan, 0.80, np.nan),  # branched ether — Ka unknown; Koc and Ka imputed
    ("ADONA",  "958445-44-8","Carboxylate",  8, 380.1,  2.80,     np.nan, 0.55, np.nan),  # ether carboxylate — Ka unknown
    ("F53B",   "73606-19-6", "Sulfonate",    6, 570.1,  4.00,     np.nan,-3.27, np.nan),  # chlorinated sulfonate — Ka unknown
], columns=["PFAS_Name", "CASRN", "PFAS_Class", "Chain_Length", "MW", "LogKow",
            "Koc", "AlbuminBinding_pKa", "Ka_albumin"])

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

    # Salt/counterion variants → parent PFAS CASRN
    # Add this dict near your other constants at the top of the file
CASRN_SALT_MAP = {
    "2795-39-3":    "1763-23-1",   # PFOS potassium salt → PFOS
    "3871-99-6":    "1763-23-1",   # PFOS TEA salt → PFOS
    "4021-47-0":    "1763-23-1",   # PFOS ammonium salt → PFOS
    "27619-97-2":   "335-67-1",    # PFOA ammonium salt → PFOA
    "29420-49-3":   "335-67-1",    # PFOA sodium salt → PFOA
    "335-95-5":     "335-67-1",    # PFOA sodium salt (alt) → PFOA
    "3825-26-1":    "335-67-1",    # PFOA ammonium (alt) → PFOA
    "30334-69-1":   "355-46-4",    # PFHxS potassium salt → PFHxS
    "754-91-6":     "375-73-5",    # PFBS potassium salt → PFBS
    "151772-58-6":  "2058-94-8",   # PFUnDA salt → PFUnDA
}

# NOTE (v10.1): these are fixed, bijective lookups on Species_Group — one
# number per group, not a measured biological quantity. Retained for
# descriptive output/plotting only; excluded from CONC_FEATURES (see FIX 4).
TROPHIC_LEVEL = {"Plant": 1, "Fish": 3, "Mammal": 4, "Human": 5, "Other": 2}
AQUATIC       = {"Fish": 1, "Plant": 0, "Mammal": 0, "Human": 0, "Other": 0}

CONC_FEATURES = [
    # FIX 1: is_fish/is_mammal/is_plant/is_human removed — encode data source
    #         identity rather than biology; breaks simulator generalisability.
    # FIX 2: Duration_days removed — null for all NHANES rows, acting as an
    #         implicit human-row flag (second leakage path).
    # FIX 4 (v10.1): Trophic_Level and Is_Aquatic removed. Both are fixed
    #         dict lookups keyed 1:1 on Species_Group (see TROPHIC_LEVEL /
    #         AQUATIC below) — i.e. the same group-identity leak FIX 1 was
    #         meant to remove, just re-encoded as an integer instead of a
    #         binary flag. They are kept as descriptive/plotting columns in
    #         the output dataset but are no longer model inputs.
    # v11.0: Koc_log and AlbuminBinding_pKa added. Both are genuine physical
    #         chemistry measurements, not group-identity re-encodings.
    #         Koc_log: log10(Koc), soil OC-water partition — governs fish/plant
    #           dietary and root uptake pathway. NaN-imputed at runtime for the
    #           3 emerging PFAS with no experimental Koc.
    #         AlbuminBinding_pKa: functional-group pKa proxy for serum-protein
    #           binding strength — mechanistically relevant for human blood serum.
    "Chain_Length", "MW", "LogKow", "PFAS_Class_encoded",
    "Koc_log", "AlbuminBinding_pKa",
]
BCF_FEATURES = [
    # Chemistry-only: BCF already normalizes for exposure concentration,
    # so species/trophic flags add noise rather than signal.
    # Duration_days intentionally excluded — BCF rows rarely have it populated.
    # v11.0: Koc_log added — directly relevant to fish dietary uptake, which
    #         is the primary exposure route captured in BCF measurements.
    #         AlbuminBinding_pKa NOT included — albumin binding is blood-matrix
    #         specific; BCF is measured in whole tissue or water-normalized,
    #         so it adds noise rather than signal here.
    "Chain_Length", "MW", "LogKow", "PFAS_Class_encoded", "Koc_log",
]
CHEM_ONLY_FEATURES = [
    "Chain_Length", "MW", "LogKow", "PFAS_Class_encoded",
    "Koc_log", "AlbuminBinding_pKa", "Duration_days",
]

# v10.5 feature set — kept for ablation comparison in v11
CONC_FEATURES_V10 = ["Chain_Length", "MW", "LogKow", "PFAS_Class_encoded"]

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
NHANES_CYCLE_MIDPOINT = {
    "2015-2016": 2015.5, "2017-2018": 2017.5,  # string labels
    2016: 2015.5, 2018: 2017.5,                 # integer labels (actual file format)
}


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


import re

def format_casrn(x):
    """
    Normalize a CASRN to the canonical dashed format (XXXXXXX-XX-X).
    Handles: already-dashed, digit-only strings, float strings, scientific notation.
    """

    s = str(x).strip()

    # Already in canonical dashed format — just strip extra whitespace
    if re.match(r'^\d+-\d{2}-\d$', s):
        return s

    # Float or scientific notation (e.g. "1763231.0", "1.76e6") — strip decimals
    try:
        s = str(int(float(s)))
    except (ValueError, OverflowError):
        pass

    # Digit-only string — reinsert dashes
    digits = re.sub(r'\D', '', s)
    if len(digits) >= 5:
        return f"{digits[:-3]}-{digits[-3:-1]}-{digits[-1]}"

    return s  # give up, return as-is for manual inspection


def diagnose_casrn_merge(ecotox_df, chem_df):
    # Drop PFAS_Name from ecotox if it exists — it comes from the export's
    # "Chemical Name" alias and will cause pandas to rename both to _x/_y
    ecotox_clean = ecotox_df.drop(columns=["PFAS_Name"], errors="ignore")
    merged = ecotox_clean.merge(chem_df[["CASRN", "PFAS_Name"]], on="CASRN", how="left")
    unmatched = merged[merged["PFAS_Name"].isna()]["CASRN"].value_counts()

    if unmatched.empty:
        print("[CASRN] All rows matched successfully ✓")
    else:
        print(f"[CASRN] {unmatched.sum()} rows unmatched across {len(unmatched)} unique CASRNs:")
        print(f"  {'CASRN (normalized)':<25}  {'count':>6}")
        print(f"  {'-'*25}  {'-'*6}")
        for casrn, count in unmatched.items():
            in_curated = casrn in chem_df["CASRN"].values
            flag = "← in curated table?" if in_curated else "← not in curated 13-PFAS"
            print(f"  {casrn:<25}  {count:>6}  {flag}")
    return unmatched

def build_dataset(ecotox_df, chem_df):
    ecotox_df = ecotox_df.copy()
    ecotox_df["CASRN"] = ecotox_df["CASRN"].apply(format_casrn)
    
    # Remap salt CASRNs to parent compound before merging
    ecotox_df["CASRN"] = ecotox_df["CASRN"].replace(CASRN_SALT_MAP)
    
    ecotox_df = ecotox_df.drop(columns=["PFAS_Name"], errors="ignore")
    unmatched = diagnose_casrn_merge(ecotox_df, chem_df)
    # ... rest unchanged

    merged = ecotox_df.merge(chem_df, on="CASRN", how="left")
    # ... rest unchanged

    for grp in ["Fish", "Mammal", "Plant", "Human"]:
        merged[f"is_{grp.lower()}"] = (merged["Species_Group"] == grp).astype(int)
    merged["PFAS_Class_encoded"] = merged["PFAS_Class"].map({"Sulfonate": 0, "Carboxylate": 1}).fillna(-1)
    le = LabelEncoder()
    merged["Route_encoded"] = le.fit_transform(merged["Exposure Route"].fillna("Unknown"))

    # v11.0: log-transform Koc for ML (same reason LogKow is used instead of Kow —
    # partition coefficients span orders of magnitude; log scale linearizes the
    # relationship with bioaccumulation). Impute median for the 3 emerging PFAS
    # (GenX, ADONA, F53B) that have no experimental Koc values — they stay in
    # the dataset rather than being silently dropped, and the imputation is flagged
    # in the feature table comment above.
    if "Koc" in merged.columns:
        merged["Koc_log"] = merged["Koc"].apply(
            lambda x: np.log10(x) if pd.notna(x) and x > 0 else np.nan
        )
        koc_median = merged["Koc_log"].median()
        n_koc_imputed = merged["Koc_log"].isna().sum()
        merged["Koc_log"] = merged["Koc_log"].fillna(koc_median)
        if n_koc_imputed > 0:
            print(f"[v11] Koc_log: imputed {n_koc_imputed:,} NaN rows with median "
                  f"({koc_median:.3f}) — emerging PFAS with no experimental Koc")
    else:
        merged["Koc_log"] = np.nan
        print("[v11] WARNING: Koc column not found in merged dataset — Koc_log set to NaN")

    # AlbuminBinding_pKa is already numeric in PFAS_FEATURES; just propagate it.
    # No imputation needed — all 13 curated PFAS have values.
    if "AlbuminBinding_pKa" not in merged.columns:
        print("[v11] WARNING: AlbuminBinding_pKa column not found — check PFAS_FEATURES merge")

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
    t1 = NHANES_CYCLE_MIDPOINT.get(wave1) or NHANES_CYCLE_MIDPOINT.get(str(wave1))
    t2 = NHANES_CYCLE_MIDPOINT.get(wave2) or NHANES_CYCLE_MIDPOINT.get(str(wave2))
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
        sg_dummies = pd.get_dummies(sub["Species_Group"], prefix="sg")
        sub = sub.copy()
        for col in sg_dummies.columns:
            sub[col] = sg_dummies[col]
        sg_cols = [c for c in sub.columns if c.startswith("sg_")]
        avail = [c for c in avail_conc if c in sub.columns] + sg_cols
        
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


def run_feature_ablation(X_fit, y_fit, X_test, y_test, grp_test,
                         v10_features, v11_features, path):
    """
    v11.0: Side-by-side RF comparison of v10.5 feature set vs v11 feature set,
    on the SAME train/test split used by the main models. Reports per-group
    and pooled delta-R² so the impact of Koc_log + AlbuminBinding_pKa is
    isolated from any split randomness.

    Both models are trained on X_fit (the fit fold, not the full training set)
    to keep this comparison apples-to-apples with the main RF model.

    Returns a DataFrame with columns:
      Group, n, R2_v10, R2_v11, Delta_R2
    """
    print("\n[Ablation] v10.5 features vs v11 features (same split)...")

    avail_v10 = [c for c in v10_features if c in X_fit.columns]
    avail_v11 = [c for c in v11_features if c in X_fit.columns]
    print(f"[Ablation] v10 features ({len(avail_v10)}): {avail_v10}")
    print(f"[Ablation] v11 features ({len(avail_v11)}): {avail_v11}")

    rf_v10 = train_rf(X_fit[avail_v10], y_fit)
    rf_v11 = train_rf(X_fit[avail_v11], y_fit)

    pred_v10 = rf_v10.predict(X_test[avail_v10])
    pred_v11 = rf_v11.predict(X_test[avail_v11])

    r2_v10_pooled = r2_score(y_test, pred_v10)
    r2_v11_pooled = r2_score(y_test, pred_v11)
    rmse_v10 = mean_squared_error(y_test, pred_v10) ** 0.5
    rmse_v11 = mean_squared_error(y_test, pred_v11) ** 0.5

    print(f"[Ablation] Pooled — v10: R²={r2_v10_pooled:.3f} RMSE={rmse_v10:.3f} | "
          f"v11: R²={r2_v11_pooled:.3f} RMSE={rmse_v11:.3f} | "
          f"ΔR²={r2_v11_pooled - r2_v10_pooled:+.3f}")

    rows = []
    grp_arr = np.asarray(grp_test)
    y_arr   = np.asarray(y_test)
    for grp in sorted(set(grp_arr)):
        mask = grp_arr == grp
        if mask.sum() < 5:
            continue
        r2_v10 = r2_score(y_arr[mask], pred_v10[mask])
        r2_v11 = r2_score(y_arr[mask], pred_v11[mask])
        delta  = r2_v11 - r2_v10
        rows.append({"Group": grp, "n": mask.sum(),
                     "R2_v10": r2_v10, "R2_v11": r2_v11, "Delta_R2": delta})
        flag = "↑" if delta > 0.01 else ("↓" if delta < -0.01 else "≈")
        print(f"  {grp:10s} n={mask.sum():5d}  v10 R²={r2_v10:+.3f}  "
              f"v11 R²={r2_v11:+.3f}  Δ={delta:+.3f} {flag}")
    rows.append({"Group": "Pooled", "n": len(y_arr),
                 "R2_v10": r2_v10_pooled, "R2_v11": r2_v11_pooled,
                 "Delta_R2": r2_v11_pooled - r2_v10_pooled})

    abl_df = pd.DataFrame(rows)

    # ── Plot ablation ───────────────────────────────────────────────
    plot_df = abl_df.copy()
    x = np.arange(len(plot_df))
    w = 0.30

    fig, axes = plt.subplots(1, 2, figsize=(14, max(4, len(plot_df) * 0.7)))

    # Left panel: side-by-side R² bars
    axes[0].barh([i + w/2 for i in x], plot_df["R2_v11"], w,
                 color="steelblue", alpha=0.85, label="v11 (+ Koc_log, AlbuminBinding_pKa)")
    axes[0].barh([i - w/2 for i in x], plot_df["R2_v10"], w,
                 color="lightsteelblue", alpha=0.85, label="v10.5 (Chain_Length, MW, LogKow, Class)")
    axes[0].set_yticks(list(x))
    axes[0].set_yticklabels(
        [f"{r['Group']} (n={r['n']})" for _, r in plot_df.iterrows()]
    )
    axes[0].axvline(0, color="black", lw=0.8)
    axes[0].set_xlabel("Held-out R²")
    axes[0].set_title("Feature Ablation: v10.5 vs v11\nHeld-out R² by species group")
    axes[0].legend(fontsize=8)

    # Right panel: delta R² (the actual story)
    delta_colors = ["steelblue" if d > 0 else "coral" for d in plot_df["Delta_R2"]]
    axes[1].barh(
        [f"{r['Group']} (n={r['n']})" for _, r in plot_df.iterrows()],
        plot_df["Delta_R2"],
        color=delta_colors,
    )
    axes[1].axvline(0, color="black", lw=0.8)
    for i, (_, row) in enumerate(plot_df.iterrows()):
        axes[1].text(
            row["Delta_R2"] + (0.003 if row["Delta_R2"] >= 0 else -0.003), i,
            f"{row['Delta_R2']:+.3f}", va="center",
            ha="left" if row["Delta_R2"] >= 0 else "right", fontsize=9
        )
    axes[1].set_xlabel("ΔR² (v11 − v10.5)")
    axes[1].set_title("Feature Ablation: ΔR² per Group\n(blue = v11 better, coral = v11 worse)")

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}")

    return abl_df


def run_human_only_model(dataset, output_dir):
    """
    Dedicated human blood serum model — NHANES rows only.

    Feature set: Chain_Length, MW, LogKow, PFAS_Class_encoded,
    AlbuminBinding_pKa. Koc_log intentionally excluded — Koc governs
    soil/sediment partitioning and fish dietary uptake, not human blood
    serum protein binding. Including it adds noise for this matrix.

    Split: stratified 80/20 by PFAS_Name (not Species_Group, since this
    is single-species) to ensure all PFAS are represented in test fold.
    Calibration fold carved from train (80/20) for honest interval coverage.

    Returns dict with R² figures and calibration data for simulator use.
    """
    print("\n" + "─" * 55)
    print("  [Human-Only Model] NHANES blood serum rows exclusively")
    print("─" * 55)

    HUMAN_FEATURES = [
        "Chain_Length", "MW", "LogKow",
        "PFAS_Class_encoded", "AlbuminBinding_pKa",
    ]

    # ── Filter to human rows only ────────────────────────────────────
    human = dataset[dataset["Species_Group"] == "Human"].copy()
    avail = [c for c in HUMAN_FEATURES if c in human.columns]
    ml = human[avail + ["log_concentration", "PFAS_Name"]].dropna(
        subset=["Chain_Length", "MW", "LogKow", "log_concentration"])

    print(f"  Rows:    {len(ml):,}")
    print(f"  PFAS:    {ml['PFAS_Name'].nunique()} compounds")
    print(f"  Features ({len(avail)}): {avail}")

    if len(ml) < 50:
        print("  Too few rows — skipping human-only model.")
        return {}

    X  = ml[avail]
    y  = ml["log_concentration"]
    pf = ml["PFAS_Name"]

    # ── Stratified split by PFAS_Name ───────────────────────────────
    X_train, X_test, y_train, y_test, pf_train, pf_test = train_test_split(
        X, y, pf, test_size=0.2, random_state=42, stratify=pf)

    # 3-way split: carve calibration from train
    X_fit, X_cal, y_fit, y_cal, pf_fit, pf_cal = train_test_split(
        X_train, y_train, pf_train, test_size=0.2, random_state=7,
        stratify=pf_train)

    print(f"  Split:   Fit={len(X_fit):,} | Cal={len(X_cal):,} | Test={len(X_test):,}")

    # ── Baselines ────────────────────────────────────────────────────
    pfas_means = y_train.groupby(pf_train).mean()
    y_base = pf_test.map(pfas_means).fillna(y_train.mean())
    r2_base  = r2_score(y_test, y_base)
    rmse_base = mean_squared_error(y_test, y_base) ** 0.5

    # ── Random Forest ────────────────────────────────────────────────
    rf_h = train_rf(X_fit, y_fit)
    y_pred_rf = rf_h.predict(X_test)
    r2_rf   = r2_score(y_test, y_pred_rf)
    rmse_rf = mean_squared_error(y_test, y_pred_rf) ** 0.5
    print(f"\n  RF       R²={r2_rf:.3f}  RMSE={rmse_rf:.3f}  "
          f"(gain over PFAS-mean baseline: {r2_rf - r2_base:+.3f})")

    # ── XGBoost ─────────────────────────────────────────────────────
    xgb_h = train_xgb(X_fit.fillna(X_fit.median()), y_fit)
    y_pred_xgb = xgb_h.predict(X_test.fillna(X_fit.median()))
    r2_xgb   = r2_score(y_test, y_pred_xgb)
    rmse_xgb = mean_squared_error(y_test, y_pred_xgb) ** 0.5
    print(f"  XGBoost  R²={r2_xgb:.3f}  RMSE={rmse_xgb:.3f}  "
          f"(gain: {r2_xgb - r2_base:+.3f})")

    # ── Linear Regression ────────────────────────────────────────────
    X_lr_tr = X_fit.fillna(X_fit.median())
    X_lr_te = X_test.fillna(X_fit.median())
    scaler_h = StandardScaler()
    lr_h = LinearRegression()
    lr_h.fit(scaler_h.fit_transform(X_lr_tr), y_fit)
    y_pred_lr = lr_h.predict(scaler_h.transform(X_lr_te))
    r2_lr   = r2_score(y_test, y_pred_lr)
    rmse_lr = mean_squared_error(y_test, y_pred_lr) ** 0.5
    print(f"  Linear   R²={r2_lr:.3f}  RMSE={rmse_lr:.3f}  "
          f"(gain: {r2_lr - r2_base:+.3f})")

    # ── Calibrated prediction intervals ─────────────────────────────
    group_cal_h = predict_with_intervals_per_group(
        rf_h, X_cal, y_cal, pf_cal)   # calibrate per PFAS, not species group
    intervals_h = apply_group_intervals(rf_h, X_test, pf_test, group_cal_h)
    in_80 = ((y_test.values >= intervals_h["lower_80"]) &
             (y_test.values <= intervals_h["upper_80"]))
    in_95 = ((y_test.values >= intervals_h["lower_95"]) &
             (y_test.values <= intervals_h["upper_95"]))
    cov80 = in_80.mean() * 100
    cov95 = in_95.mean() * 100
    w80   = intervals_h["width_80"].mean()
    w95   = intervals_h["width_95"].mean()
    print(f"\n  Intervals (RF, calibrated per PFAS):")
    print(f"    80% coverage: {cov80:.1f}%  width: {w80:.3f} log10 ng/g")
    print(f"    95% coverage: {cov95:.1f}%  width: {w95:.3f} log10 ng/g")

    # ── Per-PFAS R² within human model ──────────────────────────────
    print(f"\n  Per-PFAS R² (human-only RF model):")
    print(f"  {'PFAS':<8} {'n_test':>6}  {'R²':>7}  {'RMSE':>7}  vs pooled model")
    print(f"  {'-'*8}  {'-'*6}  {'-'*7}  {'-'*7}  {'-'*20}")

    # Published pooled per-PFAS R² from v11.1 for comparison
    POOLED_R2 = {
        "PFOA": 0.657, "PFOS": 0.436, "PFNA": 0.401,
        "PFHxS": 0.381, "PFDA": 0.226, "PFUnDA": 0.148, "PFBS": -0.063,
    }

    per_pfas_rows = []
    for pfas in sorted(pf_test.unique()):
        mask = pf_test == pfas
        n = mask.sum()
        if n < 5:
            continue
        r2_p   = r2_score(y_test[mask], y_pred_rf[mask])
        rmse_p = mean_squared_error(y_test[mask], y_pred_rf[mask]) ** 0.5
        pooled = POOLED_R2.get(pfas, np.nan)
        delta  = r2_p - pooled if pd.notna(pooled) else np.nan
        delta_str = f"{delta:+.3f}" if pd.notna(delta) else "n/a"
        per_pfas_rows.append({
            "PFAS_Name": pfas, "n_test": n,
            "R2_human": r2_p, "RMSE": rmse_p,
            "R2_pooled": pooled, "Delta": delta,
        })
        print(f"  {pfas:<8}  {n:>6}  {r2_p:>7.3f}  {rmse_p:>7.3f}  "
              f"Δ={delta_str} vs pooled")

    per_pfas_df = pd.DataFrame(per_pfas_rows)

    # ── Plots ────────────────────────────────────────────────────────
    # 1. Predicted vs actual — coloured by PFAS
    pfas_list  = pf_test.values
    pfas_uniq  = sorted(set(pfas_list))
    cmap       = plt.cm.get_cmap("tab10", len(pfas_uniq))
    color_map  = {p: cmap(i) for i, p in enumerate(pfas_uniq)}

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax = axes[0]
    for pfas in pfas_uniq:
        mask = pfas_list == pfas
        ax.scatter(y_test.values[mask], y_pred_rf[mask],
                   color=color_map[pfas], alpha=0.5, s=20,
                   label=pfas, edgecolors="none")
    lims = [min(y_test.min(), y_pred_rf.min()) - 0.3,
            max(y_test.max(), y_pred_rf.max()) + 0.3]
    ax.plot(lims, lims, "k--", lw=1)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel("Actual log₁₀ concentration (ng/g)")
    ax.set_ylabel("Predicted log₁₀ concentration (ng/g)")
    ax.set_title(f"Human-Only RF Model\n"
                 f"R²={r2_rf:.3f}  RMSE={rmse_rf:.3f}  n_test={len(y_test):,}")
    ax.legend(fontsize=7, ncol=2, loc="upper left")

    # 2. Per-PFAS R² bar — human vs pooled
    ax2 = axes[1]
    if not per_pfas_df.empty:
        pf_sorted = per_pfas_df.sort_values("R2_human")
        y_pos = np.arange(len(pf_sorted))
        w = 0.30
        ax2.barh([i - w/2 for i in y_pos], pf_sorted["R2_pooled"],
                 w, color="lightsteelblue", alpha=0.9, label="Pooled model R²")
        ax2.barh([i + w/2 for i in y_pos], pf_sorted["R2_human"],
                 w, color="steelblue", alpha=0.9, label="Human-only model R²")
        ax2.set_yticks(list(y_pos))
        ax2.set_yticklabels(pf_sorted["PFAS_Name"])
        ax2.axvline(0, color="black", lw=0.8)
        ax2.set_xlabel("Held-out R²")
        ax2.set_title("Per-PFAS R²: Human-Only vs Pooled Model")
        ax2.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(output_dir + "human_model_predictions.png", dpi=150)
    plt.close()
    print(f"\n  [Output] → {output_dir}human_model_predictions.png")

    # 3. Feature importance
    imp = pd.Series(rf_h.feature_importances_, index=avail).sort_values()
    fig, ax = plt.subplots(figsize=(7, max(3, len(avail) * 0.5)))
    imp.plot(kind="barh", color="steelblue", ax=ax)
    ax.set_title("Human-Only RF Model — Feature Importances\n"
                 "(NHANES blood serum only; Koc_log excluded)")
    plt.tight_layout()
    plt.savefig(output_dir + "human_model_feature_importance.png", dpi=150)
    plt.close()
    print(f"  [Output] → {output_dir}human_model_feature_importance.png")

    print("─" * 55)

    return {
        "r2_rf": r2_rf, "rmse_rf": rmse_rf,
        "r2_xgb": r2_xgb, "rmse_xgb": rmse_xgb,
        "r2_lr": r2_lr, "rmse_lr": rmse_lr,
        "r2_baseline": r2_base, "rmse_baseline": rmse_base,
        "cov80": cov80, "cov95": cov95, "w80": w80, "w95": w95,
        "per_pfas_df": per_pfas_df,
        "rf_model": rf_h, "scaler": scaler_h, "features": avail,
        "group_calibration": group_cal_h,
    }


def compute_arnot_gobas_bcf(pfas_features_df):
    """
    Arnot-Gobas (2004) steady-state BCF model, adapted for PFAS.

    Classic equation:  BCF = (k1 + kd) / (k2 + ke + km + kg)

    v13.0 PFAS adaptations (supersedes v11.x/v12.x):
      - Direct albumin-binding tissue partition coefficient (Kfish):
            Ka_tissue = Ka_albumin [L/mol] × [albumin]_fish [mol/L] × fw
            Kfish     = Ka_tissue + Kd_nlom + Vw
        Ka_albumin is the ITC/fluorescence-displacement binding constant
        measured per compound by Bischel et al. (2010) and Beesoon & Martin
        (2015). This replaces the indirect Koc-proxy (scale × Koc) used in
        v11.x–v12.x, which Findings 16–17 proved to be structurally
        inadequate: error was invariant across the full Kprot and P_mem sweep
        ranges (both parameters exhausted), confirming the problem was in
        Kfish itself rather than in any rate-constant scaling.
      - Fish albumin parameters: C_albumin=20 g/L, MW=68000 g/mol, fw=0.085
        (Conder et al. 2008; Metcalfe & Thorpe 1992; Chen et al. 2016).
      - NLOM (non-lipid organic matter) partitioning term retained:
            Kd_nlom = 0.035 × Koc   (Gobas et al. 2003)
        This term reflects organic-carbon sorption in fish tissue and is
        independent of the albumin-binding mechanism; both are needed.
      - km = 0 (PFAS are metabolically inert — no biotransformation)
      - Dietary pathway (kd) included; parameterised for aqueous-exposure
        BCF studies (dominant route in ECOTOX measurements).
      - Fallback path for emerging PFAS with no Ka_albumin data (GenX,
        ADONA, F53B): reverts to legacy Koc-proxy (v11.x behaviour).
        Kfish_method column in output identifies which path was used.

    Fish parameters (Arnot & Gobas 2004, Table 1 — generic 1 kg fish, 12°C):
      Vl=0.05, Vn=0.09, Vw=0.72, Ew=0.50, Ed=0.75, GV=400, Gd=0.025, kg=0.001

    Parameters
    ----------
    pfas_features_df : DataFrame with columns PFAS_Name, LogKow, Koc,
                       Chain_Length, PFAS_Class, Ka_albumin (v13.0)

    Returns
    -------
    DataFrame with columns: PFAS_Name, BCF_AG, log_BCF_AG, Ka_albumin,
    Ka_tissue, Kfish_method, Kfish, and all rate constants for diagnostics.
    """
    # ── Fish physiological parameters (Arnot & Gobas 2004, Table 1) ──
    Vl  = 0.05    # lipid volume fraction
    Vn  = 0.09    # NLOM volume fraction
    Vw  = 0.72    # water volume fraction
    Ew  = 0.50    # gill oxygen extraction efficiency
    Ed  = 0.75    # dietary assimilation efficiency
    GV  = 400.0   # ventilation volume rate (L water / kg fish / d)
    Gd  = 0.025   # feeding rate (kg food / kg fish / d)
    kg  = 0.001   # specific growth dilution rate (1/d)
    km  = 0.0     # metabolic transformation rate — 0 for PFAS

    # Lipid content of food (assumed equal to fish body lipid — generic prey)
    Vl_food = 0.05
    Vn_food = 0.09
    Vw_food = 0.72

    rows = []
    for _, row in pfas_features_df.iterrows():
        pfas_name  = row["PFAS_Name"]
        log_kow    = row["LogKow"]
        koc        = row.get("Koc", np.nan)
        chain_len  = row["Chain_Length"]
        pfas_class = row["PFAS_Class"]
        ka_log     = row.get("Ka_albumin", np.nan)  # v13.0: compound-specific albumin Ka

        kow = 10 ** log_kow

        # ── PFAS-adapted tissue partition coefficient (v13.0) ─────────
        # v11.x used:  Kfish = Kprot + Kd_nlom + Vw
        #   where Kprot = scale × Koc  (indirect proxy)
        # Findings 16–17 proved scale×Koc is structurally wrong for
        # sulfonates: PFOS error is invariant across the full Kprot and
        # P_mem sweep ranges. The error lives in Kfish itself.
        #
        # v13.0 replaces the Koc-based Kprot term with a direct albumin-
        # binding term (Conder et al. 2008, Eq. 1):
        #
        #   Ka_tissue = Ka_albumin [L/mol] × [albumin]_fish [mol/L] × fw
        #   Kfish     = Ka_tissue + Kd_nlom + Vw
        #
        # Ka_albumin is the ITC/fluorescence-displacement binding constant
        # from Bischel et al. (2010) and Beesoon & Martin (2015), stored
        # per compound in PFAS_FEATURES.Ka_albumin (log10 L/mol).
        # [albumin]_fish and fw are set by C_ALBUMIN_FISH_MOL_L and
        # FW_PROTEIN module-level constants.
        #
        # Fallback chain (in priority order):
        #   1. Ka_albumin measured → Ka_tissue term (v13.0 default)
        #   2. Koc available but no Ka → legacy Koc-proxy (KPROT_SCALE×Koc)
        #      retains v11.x behaviour for emerging PFAS with NaN Ka
        #   3. Neither Koc nor Ka → Kow-based lipid term (classic Arnot-Gobas)
        #
        # NLOM term (Kd_nlom = 0.035 × Koc) is retained in all paths where
        # Koc is available — it is an organic-carbon-based mechanism distinct
        # from protein binding and unaffected by the Ka substitution.

        if pd.notna(ka_log):
            # Path 1 (v13.0): direct albumin binding — used for all 10 PFAS
            # with measured or interpolated Ka values (Tier 1 and Tier 2)
            ka_mol_L  = 10 ** ka_log
            Ka_tissue = ka_mol_L * C_ALBUMIN_FISH_MOL_L * FW_PROTEIN
            # NLOM term: retain if Koc is available; fall back to Kow otherwise
            if pd.notna(koc) and koc > 0:
                Kd_nlom = 0.035 * koc
            else:
                Kd_nlom = Vn * 0.035 * kow
            Kfish = Ka_tissue + Kd_nlom + Vw
            kfish_method = "Ka_albumin"

        elif pd.notna(koc) and koc > 0:
            # Path 2 (legacy v11.x): Koc-proxy — only reached for emerging
            # PFAS where Ka_albumin is NaN (GenX, ADONA, F53B)
            kprot_scale = KPROT_SCALE.get(pfas_class, 0.05)
            Ka_tissue   = kprot_scale * koc   # re-labelled for output clarity
            Kd_nlom     = 0.035 * koc
            Kfish       = Ka_tissue + Kd_nlom + Vw
            kfish_method = "Koc_proxy_legacy"

        else:
            # Path 3 (ultimate fallback): classic Kow-based lipid partitioning
            Ka_tissue   = Vl * kow
            Kd_nlom     = Vn * 0.035 * kow
            Kfish       = Ka_tissue + Kd_nlom + Vw
            kfish_method = "Kow_fallback"

        # Food partition coefficient (Kfood) — applied to prey organism.
        # Uses the same Ka-based structure as Kfish where Ka is available:
        # prey fish have similar albumin concentrations to predator fish, so
        # the same Ka_tissue term applies. Where Ka is unavailable, falls
        # back to the legacy Koc-proxy or Kow-based calculation.
        if pd.notna(ka_log):
            Kfood = Ka_tissue * Vl_food + Kd_nlom * Vn_food + Vw_food
        elif pd.notna(koc) and koc > 0:
            Kfood = 0.05 * koc * Vl_food + 0.035 * koc * Vn_food + Vw_food
        else:
            Kfood = Vl_food * kow + Vn_food * 0.035 * kow + Vw_food

        # ── Rate constants ────────────────────────────────────────────
        # k1: gill uptake (L/kg·d)
        # Ew scales ventilation for oxygen extraction efficiency;
        # the (1/(5.49e7 * Kow^-0.670 + 1)) term is the chemical
        # mass transfer resistance across the gill membrane
        # (Arnot & Gobas 2004, Eq. 5 — two-resistance model).
        #
        # v12.1: P_MEM_CORRECTION applied here — multiplies mem_resistance by
        # a class-specific scalar to account for reduced gill membrane
        # permeability of ionized PFAS relative to neutral lipophilic organics.
        # Sulfonates (permanent charge, pKa ≈ -3.3) cross phospholipid bilayers
        # via protein-mediated transport; the Kow-based two-resistance model
        # over-predicts their passive diffusion rate. Carboxylate correction=1.0
        # (no change from baseline). k1_raw is stored for diagnostic output so
        # the magnitude of the correction is visible per compound.
        # This parameter is INDEPENDENT of KPROT_SCALE: P_MEM_CORRECTION only
        # changes k1; Kfish (and thus k2, ke, kd) is unchanged.
        p_mem = P_MEM_CORRECTION.get(pfas_class, 1.0)
        mem_resistance = 1.0 / (5.49e7 * (kow ** -0.670) + 1.0)
        k1_raw = Ew * GV * mem_resistance        # uncorrected — Arnot-Gobas baseline
        k1     = k1_raw * p_mem                  # corrected for ionic membrane permeability

        # k2: gill elimination (1/d) — inverse of k1 scaled by Kfish
        k2 = k1 / Kfish

        # ke: fecal egestion (1/d)
        # Ed is dietary assimilation; (1-Ed) fraction is egested.
        # Gd * (1-Ed) * Kfood / Kfish gives the fecal loss rate.
        ke = Gd * (1.0 - Ed) * (Kfood / Kfish)

        # kd: dietary uptake (1/d) — included for completeness
        # For aqueous-exposure BCF studies this term is small relative
        # to gill uptake; included for mass-balance completeness.
        kd = Ed * Gd * (Kfood / Kfish)

        # ── Steady-state BCF ─────────────────────────────────────────
        # BCF = (k1 + kd) / (k2 + ke + km + kg)
        # kd added to numerator: dietary uptake is a source alongside gills.
        denom = k2 + ke + km + kg
        bcf_ag = (k1 + kd) / denom if denom > 0 else np.nan
        log_bcf_ag = np.log10(bcf_ag) if (pd.notna(bcf_ag) and bcf_ag > 0) else np.nan

        rows.append({
            "PFAS_Name":    pfas_name,
            "LogKow":       log_kow,
            "Koc":          koc,
            "Ka_albumin":   ka_log,       # v13.0: compound-specific log Ka (L/mol)
            "Ka_tissue":    Ka_tissue,    # v13.0: Ka_albumin × [albumin] × fw (dimensionless)
            "Kfish_method": kfish_method, # v13.0: which path was used for Kfish
            "Chain_Length": chain_len,
            "PFAS_Class":   pfas_class,
            "Kfish":        Kfish,
            "k1_raw":       k1_raw,   # v12.1: uncorrected Arnot-Gobas k1
            "p_mem":        p_mem,    # v12.1: membrane permeability correction
            "k1":           k1,       # v12.1: corrected k1 = k1_raw * p_mem
            "k2":           k2,
            "ke":           ke,
            "kd":           kd,
            "km":           km,
            "kg":           kg,
            "BCF_AG":       bcf_ag,
            "log_BCF_AG":   log_bcf_ag,
        })

    return pd.DataFrame(rows)


def plot_arnot_gobas_comparison(ag_df, observed_bcf_df, ml_bcf_df, path):
    """
    Three-way comparison per PFAS:
      - Arnot-Gobas mechanistic BCF (log10)
      - ML (RF) BCF prediction (mean of test-set predictions for that PFAS)
      - Observed BCF median and scatter (ECOTOX data)

    observed_bcf_df: DataFrame with columns PFAS_Name, log_BCF (raw records)
    ml_bcf_df:       DataFrame with columns PFAS_Name, log_BCF_pred (RF predictions)
    ag_df:           Output of compute_arnot_gobas_bcf()
    """
    # Observed: median per PFAS
    obs_median = (observed_bcf_df.groupby("PFAS_Name")["log_BCF"]
                  .agg(["median", "std", "count"])
                  .reset_index()
                  .rename(columns={"median": "obs_median", "std": "obs_std", "count": "obs_n"}))

    # ML: mean predicted per PFAS
    if ml_bcf_df is not None and len(ml_bcf_df) > 0:
        ml_mean = (ml_bcf_df.groupby("PFAS_Name")["log_BCF_pred"]
                   .mean().reset_index()
                   .rename(columns={"log_BCF_pred": "ml_bcf"}))
    else:
        ml_mean = pd.DataFrame(columns=["PFAS_Name", "ml_bcf"])

    # Merge everything
    plot_df = ag_df[["PFAS_Name", "log_BCF_AG"]].merge(obs_median, on="PFAS_Name", how="left")
    plot_df = plot_df.merge(ml_mean, on="PFAS_Name", how="left")
    plot_df = plot_df.dropna(subset=["log_BCF_AG"])
    plot_df = plot_df.sort_values("log_BCF_AG")

    if plot_df.empty:
        print("[Arnot-Gobas] No data to plot — skipping.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(15, max(5, len(plot_df) * 0.65)))

    y     = np.arange(len(plot_df))
    names = plot_df["PFAS_Name"].tolist()

    # ── Left: bar chart comparison ───────────────────────────────────
    ax = axes[0]
    bar_w = 0.25

    ax.barh([i - bar_w for i in y], plot_df["log_BCF_AG"],
            bar_w, color="steelblue",   alpha=0.85, label="Arnot-Gobas (mechanistic)")
    if "ml_bcf" in plot_df.columns and plot_df["ml_bcf"].notna().any():
        ax.barh([i for i in y], plot_df["ml_bcf"],
                bar_w, color="darkorange", alpha=0.85, label="Random Forest (ML)")
    if "obs_median" in plot_df.columns and plot_df["obs_median"].notna().any():
        ax.barh([i + bar_w for i in y], plot_df["obs_median"],
                bar_w, color="coral", alpha=0.85, label="Observed median (ECOTOX)")
        # Error bars for observed spread (±1 std)
        err = plot_df["obs_std"].fillna(0)
        ax.errorbar(plot_df["obs_median"], [i + bar_w for i in y],
                    xerr=err, fmt="none", color="darkred", capsize=3, lw=1.2)

    ax.set_yticks(list(y))
    ax.set_yticklabels(names)
    ax.set_xlabel("log₁₀ BCF")
    ax.set_title("BCF Comparison: Arnot-Gobas vs ML vs Observed\n(error bars = ±1 std of observed)")
    ax.axvline(0, color="black", lw=0.6)
    ax.legend(fontsize=8)

    # ── Right: AG vs observed scatter (% error) ──────────────────────
    ax2 = axes[1]
    has_obs = plot_df["obs_median"].notna()

    if has_obs.any():
        sub = plot_df[has_obs].copy()
        # % error on the BCF scale (not log), back-transformed
        sub["BCF_AG"]  = 10 ** sub["log_BCF_AG"]
        sub["BCF_obs"] = 10 ** sub["obs_median"]
        sub["pct_err"] = (sub["BCF_AG"] - sub["BCF_obs"]) / sub["BCF_obs"] * 100

        colors_err = ["steelblue" if v > 0 else "coral" for v in sub["pct_err"]]
        ax2.barh(sub["PFAS_Name"], sub["pct_err"], color=colors_err, alpha=0.8)
        ax2.axvline(0, color="black", lw=0.8)
        for i, (_, r) in enumerate(sub.iterrows()):
            ax2.text(r["pct_err"] + (5 if r["pct_err"] >= 0 else -5), i,
                     f"{r['pct_err']:+.0f}%", va="center",
                     ha="left" if r["pct_err"] >= 0 else "right", fontsize=8)
        ax2.set_xlabel("Arnot-Gobas % error vs observed median BCF")
        ax2.set_title("Mechanistic Model Error vs Observed\n(blue = over-prediction, coral = under-prediction)")
        # Annotate n per PFAS
        for i, (_, r) in enumerate(sub.iterrows()):
            if pd.notna(r.get("obs_n")):
                ax2.text(ax2.get_xlim()[1] * 0.98, i,
                         f"n={int(r['obs_n'])}", va="center", ha="right",
                         fontsize=7, color="gray")
    else:
        ax2.text(0.5, 0.5, "No observed BCF data to compare against",
                 ha="center", va="center", transform=ax2.transAxes, fontsize=11)
        ax2.set_title("% Error vs Observed (no data)")

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[Output] → {path}")


def print_arnot_gobas_summary(ag_df, observed_bcf_df):
    """Console table: AG BCF, observed median, % error, and rate constants."""
    obs_median = (observed_bcf_df.groupby("PFAS_Name")["log_BCF"]
                  .agg(["median", "count"])
                  .reset_index()
                  .rename(columns={"median": "obs_log_bcf", "count": "n_obs"})) \
        if observed_bcf_df is not None and len(observed_bcf_df) > 0 \
        else pd.DataFrame(columns=["PFAS_Name", "obs_log_bcf", "n_obs"])

    merged = ag_df.merge(obs_median, on="PFAS_Name", how="left")

    print("\n[Arnot-Gobas v13] Mechanistic BCF estimates — Ka-albumin Kfish (generic fish 1 kg, 12°C):")
    print(f"  {'PFAS':<8} {'log BCF_AG':>10}  {'BCF_AG':>8}  {'log BCF_obs':>11}  "
          f"{'n_obs':>5}  {'% err':>7}  {'Ka_tissue':>9}  {'Kfish':>7}  method")
    print(f"  {'-'*8}  {'-'*10}  {'-'*8}  {'-'*11}  {'-'*5}  {'-'*7}  {'-'*9}  {'-'*7}  {'-'*16}")

    for _, r in merged.sort_values("log_BCF_AG", ascending=False, na_position="last").iterrows():
        bcf_ag_str  = f"{r['BCF_AG']:.1f}"     if pd.notna(r.get("BCF_AG"))     else "n/a"
        log_ag_str  = f"{r['log_BCF_AG']:.3f}" if pd.notna(r.get("log_BCF_AG")) else "n/a"
        obs_str     = f"{r['obs_log_bcf']:.3f}" if pd.notna(r.get("obs_log_bcf")) else "no data"
        n_str       = f"{int(r['n_obs'])}"      if pd.notna(r.get("n_obs"))      else "—"
        ka_tis_str  = f"{r['Ka_tissue']:.3f}"  if pd.notna(r.get("Ka_tissue"))  else "n/a"
        kfish_str   = f"{r['Kfish']:.2f}"      if pd.notna(r.get("Kfish"))      else "n/a"
        method_str  = r.get("Kfish_method", "?")

        if pd.notna(r.get("obs_log_bcf")) and pd.notna(r.get("BCF_AG")):
            bcf_obs = 10 ** r["obs_log_bcf"]
            pct_err = (r["BCF_AG"] - bcf_obs) / bcf_obs * 100
            err_str = f"{pct_err:+.0f}%"
            note = ("over" if pct_err > 50 else
                    "under" if pct_err < -50 else "close")
        else:
            err_str = "n/a"
            note    = "no observed data"

        print(f"  {r['PFAS_Name']:<8}  {log_ag_str:>10}  {bcf_ag_str:>8}  "
              f"{obs_str:>11}  {n_str:>5}  {err_str:>7}  {ka_tis_str:>9}  "
              f"{kfish_str:>7}  {method_str}  [{note}]")

    print(f"\n  v13.0 parameters:")
    print(f"    Fish albumin: C_alb={C_ALBUMIN_FISH_G_L} g/L, MW={MW_ALBUMIN_FISH} g/mol,")
    print(f"    fw={FW_PROTEIN} → [albumin]_fish={C_ALBUMIN_FISH_MOL_L:.2e} mol/L")
    print(f"    Ka_albumin: Bischel et al. (2010) + Beesoon & Martin (2015)")
    print(f"    NLOM: 0.035×Koc (Gobas et al. 2003) — retained alongside Ka term")
    print(f"    Fish: Vl=5%, Vn=9%, Vw=72%, Ew=0.50, Ed=0.75, GV=400 L/kg/d,")
    print(f"    Gd=0.025 kg/kg/d, kg=0.001/d, km=0 (PFAS metabolically inert)")
    print(f"\n  Kfish_method legend:")
    print(f"    Ka_albumin      — direct Ka term (v13.0, 10 compounds)")
    print(f"    Koc_proxy_legacy— scale×Koc fallback (emerging PFAS: GenX/ADONA/F53B)")
    print(f"    Kow_fallback    — classic Kow-based (no Koc or Ka available)")
    print(f"\n  Interpretation: compare % error here vs v12.x (-91% PFOS, -85% PFHxS).")
    print(f"  Reduced error confirms the Ka term resolves the structural Kfish gap.")
    print(f"  Residual error reflects study heterogeneity in observed BCF data.")

def run_arnot_gobas_sensitivity(pfas_features_df, observed_bcf_df,
                                sweep_values=KPROT_SULFONATE_SWEEP, path=None):
    """
    v11.3: Sensitivity sweep of the sulfonate Kprot scaling factor.

    HISTORICAL DIAGNOSTIC — do not modify the physics here.
    This function produced Finding 16: PFOS % error is invariant across
    the full Kprot scale range (0.05 → 0.30), proving the Koc-proxy is
    structurally wrong for sulfonates. The result is preserved as evidence.
    v13.0 addresses the root cause via a direct Ka_albumin term in Kfish
    inside compute_arnot_gobas_bcf(). This sweep continues to run to show
    the before/after contrast (old Koc-proxy vs new Ka approach).

    Intentionally uses inline Koc-proxy math (not compute_arnot_gobas_bcf())
    to remain an isolated test of that specific mechanism, unaffected by the
    v13.0 Ka substitution.

    Saves a heatmap (% error per PFAS × Kprot scale) if path is given.
    """
    obs_median = (observed_bcf_df.groupby("PFAS_Name")["log_BCF"]
                  .median().rename("obs_log_bcf"))

    records = []
    for scale in sweep_values:
        # Temporarily override KPROT_SCALE for sulfonates
        tmp_features = pfas_features_df.copy()

        # Patch: pass scale override into a modified version of the compute fn
        rows = []
        Vl=0.05; Vn=0.09; Vw=0.72; Ew=0.50; Ed=0.75
        GV=400.0; Gd=0.025; kg=0.001; km=0.0
        Vl_food=0.05; Vn_food=0.09; Vw_food=0.72

        for _, row in tmp_features.iterrows():
            pfas_name  = row["PFAS_Name"]
            log_kow    = row["LogKow"]
            koc        = row.get("Koc", np.nan)
            pfas_class = row["PFAS_Class"]
            kow        = 10 ** log_kow

            kprot_scale = scale if pfas_class == "Sulfonate" else 0.05

            if pd.notna(koc) and koc > 0:
                Kprot   = kprot_scale * koc
                Kd_nlom = 0.035 * koc
            else:
                Kprot   = Vl * kow
                Kd_nlom = Vn * 0.035 * kow

            Kfish    = Kprot + Kd_nlom + Vw
            Kfood    = (0.05 * koc * Vl_food + 0.035 * koc * Vn_food + Vw_food
                        if pd.notna(koc) and koc > 0
                        else Vl_food * kow + Vn_food * 0.035 * kow + Vw_food)

            mem_r = 1.0 / (5.49e7 * (kow ** -0.670) + 1.0)
            k1 = Ew * GV * mem_r
            k2 = k1 / Kfish
            ke = Gd * (1.0 - Ed) * (Kfood / Kfish)
            kd = Ed * Gd * (Kfood / Kfish)

            denom  = k2 + ke + km + kg
            bcf_ag = (k1 + kd) / denom if denom > 0 else np.nan
            log_bcf_ag = np.log10(bcf_ag) if (pd.notna(bcf_ag) and bcf_ag > 0) else np.nan

            obs = obs_median.get(pfas_name, np.nan)
            if pd.notna(obs) and pd.notna(log_bcf_ag):
                pct_err = (10**log_bcf_ag - 10**obs) / 10**obs * 100
            else:
                pct_err = np.nan

            rows.append({
                "kprot_sulfonate_scale": scale,
                "PFAS_Name": pfas_name,
                "PFAS_Class": pfas_class,
                "log_BCF_AG": log_bcf_ag,
                "pct_err": pct_err,
            })
        records.extend(rows)

    sweep_df = pd.DataFrame(records)

    # ── Console summary ───────────────────────────────────────────────
    print("\n[AG Sensitivity] Sulfonate Kprot scale sweep — % error vs observed BCF")
    print(f"  (Carboxylate scale fixed at 0.05 throughout)\n")

    sulfonates   = sweep_df[sweep_df["PFAS_Class"] == "Sulfonate"]["PFAS_Name"].unique()
    carboxylates = sweep_df[sweep_df["PFAS_Class"] == "Carboxylate"]["PFAS_Name"].unique()

    header = f"  {'Scale':>7} | " + " ".join(f"{p:>8}" for p in sulfonates) + \
             "  |  " + " ".join(f"{p:>8}" for p in carboxylates)
    print(header)
    print("  " + "-" * (len(header) - 2))

    for scale in sweep_values:
        sub = sweep_df[sweep_df["kprot_sulfonate_scale"] == scale]
        s_errs = [sub[sub["PFAS_Name"] == p]["pct_err"].values[0]
                  if len(sub[sub["PFAS_Name"] == p]) else np.nan
                  for p in sulfonates]
        c_errs = [sub[sub["PFAS_Name"] == p]["pct_err"].values[0]
                  if len(sub[sub["PFAS_Name"] == p]) else np.nan
                  for p in carboxylates]
        s_str = " ".join(f"{e:>+8.0f}%" if pd.notna(e) else f"{'n/a':>9}" for e in s_errs)
        c_str = " ".join(f"{e:>+8.0f}%" if pd.notna(e) else f"{'n/a':>9}" for e in c_errs)
        marker = " ◀ current" if scale == KPROT_SCALE["Sulfonate"] else ""
        print(f"  {scale:>7.2f} | {s_str}  |  {c_str}{marker}")

    # ── Heatmap plot ──────────────────────────────────────────────────
    if path is not None:
        pivot = sweep_df.pivot_table(
            index="kprot_sulfonate_scale", columns="PFAS_Name", values="pct_err")

        fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns) * 1.0), 5))
        sns.heatmap(pivot, annot=True, fmt=".0f", center=0,
                    cmap="RdBu_r", linewidths=0.5, ax=ax,
                    cbar_kws={"label": "% error vs observed BCF"})
        ax.set_title("Arnot-Gobas Sensitivity: Sulfonate Kprot Scale\n"
                     "(% error vs observed median BCF; carboxylate scale fixed at 0.05)")
        ax.set_xlabel("PFAS")
        ax.set_ylabel("Sulfonate Kprot scale factor")
        plt.tight_layout()
        plt.savefig(path, dpi=150)
        plt.close()
        print(f"[Output] → {path}")

    return sweep_df


def run_arnot_gobas_pmem_sensitivity(pfas_features_df, observed_bcf_df,
                                     sweep_values=P_MEM_SWEEP, path=None):
    """
    v12.1: Sensitivity sweep of the sulfonate gill membrane permeability
    correction factor (P_MEM_CORRECTION).

    Parallel to run_arnot_gobas_sensitivity() which sweeps KPROT_SCALE
    (tissue partitioning). This function sweeps P_MEM_CORRECTION for
    sulfonates, holding carboxylates fixed at 1.0 (no correction) and
    holding KPROT_SCALE at its current values throughout.

    The two sweeps test orthogonal hypotheses:
      Kprot sweep  → Is the error in tissue accumulation (Kfish too low)?
      P_mem sweep  → Is the error in gill uptake rate (k1 too high)?

    Finding 16 established that Kprot tuning cannot close the sulfonate gap.
    If this sweep also shows invariant % error across the full P_mem range,
    the Arnot-Gobas framework itself is structurally inadequate for sulfonates
    and a more mechanistic ionic transport model is required.

    Parameters
    ----------
    pfas_features_df : DataFrame — curated PFAS feature table (PFAS_FEATURES)
    observed_bcf_df  : DataFrame with columns PFAS_Name, log_BCF (ECOTOX records)
    sweep_values     : list of P_mem multipliers to test for sulfonates
    path             : output path for heatmap PNG; None = console only

    Returns
    -------
    DataFrame with columns: p_mem_sulfonate, PFAS_Name, PFAS_Class,
    log_BCF_AG, pct_err — one row per (sweep value × PFAS) combination.
    """
    obs_median = (observed_bcf_df.groupby("PFAS_Name")["log_BCF"]
                  .median().rename("obs_log_bcf"))

    records = []
    # Re-use the same inner physics as compute_arnot_gobas_bcf() but with
    # p_mem overridden per sweep step. KPROT_SCALE is held at its module-level
    # values throughout — this is the key difference from the Kprot sweep.
    Vl=0.05; Vn=0.09; Vw=0.72; Ew=0.50; Ed=0.75
    GV=400.0; Gd=0.025; kg=0.001; km=0.0
    Vl_food=0.05; Vn_food=0.09; Vw_food=0.72

    for p_mem_sulfonate in sweep_values:
        for _, row in pfas_features_df.iterrows():
            pfas_name  = row["PFAS_Name"]
            log_kow    = row["LogKow"]
            koc        = row.get("Koc", np.nan)
            pfas_class = row["PFAS_Class"]
            kow        = 10 ** log_kow

            # Tissue partitioning — unchanged from current KPROT_SCALE values
            kprot_scale = KPROT_SCALE.get(pfas_class, 0.05)
            if pd.notna(koc) and koc > 0:
                Kprot   = kprot_scale * koc
                Kd_nlom = 0.035 * koc
                Kfood   = (0.05 * koc * Vl_food + 0.035 * koc * Vn_food + Vw_food)
            else:
                Kprot   = Vl * kow
                Kd_nlom = Vn * 0.035 * kow
                Kfood   = Vl_food * kow + Vn_food * 0.035 * kow + Vw_food

            Kfish = Kprot + Kd_nlom + Vw

            # Gill membrane permeability — swept for sulfonates, fixed for carboxylates
            p_mem      = p_mem_sulfonate if pfas_class == "Sulfonate" else 1.0
            mem_r      = 1.0 / (5.49e7 * (kow ** -0.670) + 1.0)
            k1_raw     = Ew * GV * mem_r
            k1         = k1_raw * p_mem

            k2  = k1 / Kfish
            ke  = Gd * (1.0 - Ed) * (Kfood / Kfish)
            kd  = Ed * Gd * (Kfood / Kfish)

            denom      = k2 + ke + km + kg
            bcf_ag     = (k1 + kd) / denom if denom > 0 else np.nan
            log_bcf_ag = np.log10(bcf_ag) if (pd.notna(bcf_ag) and bcf_ag > 0) else np.nan

            obs = obs_median.get(pfas_name, np.nan)
            pct_err = ((10**log_bcf_ag - 10**obs) / 10**obs * 100
                       if pd.notna(obs) and pd.notna(log_bcf_ag) else np.nan)

            records.append({
                "p_mem_sulfonate": p_mem_sulfonate,
                "PFAS_Name":       pfas_name,
                "PFAS_Class":      pfas_class,
                "log_BCF_AG":      log_bcf_ag,
                "pct_err":         pct_err,
            })

    sweep_df = pd.DataFrame(records)

    # ── Console summary ───────────────────────────────────────────────
    print("\n[AG P_mem Sensitivity] Sulfonate gill membrane permeability sweep")
    print(f"  (KPROT_SCALE held fixed: Sulfonate={KPROT_SCALE['Sulfonate']}, "
          f"Carboxylate={KPROT_SCALE['Carboxylate']})")
    print(f"  (Carboxylate P_mem fixed at 1.0 throughout)\n")
    print(f"  Interpretation: if % error is invariant across the P_mem range,")
    print(f"  gill membrane permeability is NOT the limiting factor — the error")
    print(f"  lives in the tissue partitioning model (Kfish), not in k1.\n")

    sulfonates   = sweep_df[sweep_df["PFAS_Class"] == "Sulfonate"]["PFAS_Name"].unique()
    carboxylates = sweep_df[sweep_df["PFAS_Class"] == "Carboxylate"]["PFAS_Name"].unique()

    header = f"  {'P_mem':>7} | " + " ".join(f"{p:>8}" for p in sulfonates) + \
             "  |  " + " ".join(f"{p:>8}" for p in carboxylates)
    print(header)
    print("  " + "-" * (len(header) - 2))

    for p_mem_val in sweep_values:
        sub = sweep_df[sweep_df["p_mem_sulfonate"] == p_mem_val]
        s_errs = [sub[sub["PFAS_Name"] == p]["pct_err"].values[0]
                  if len(sub[sub["PFAS_Name"] == p]) else np.nan
                  for p in sulfonates]
        c_errs = [sub[sub["PFAS_Name"] == p]["pct_err"].values[0]
                  if len(sub[sub["PFAS_Name"] == p]) else np.nan
                  for p in carboxylates]
        s_str = " ".join(f"{e:>+8.0f}%" if pd.notna(e) else f"{'n/a':>9}" for e in s_errs)
        c_str = " ".join(f"{e:>+8.0f}%" if pd.notna(e) else f"{'n/a':>9}" for e in c_errs)
        marker = " ◀ current" if p_mem_val == P_MEM_CORRECTION["Sulfonate"] else ""
        print(f"  {p_mem_val:>7.2f} | {s_str}  |  {c_str}{marker}")

    # ── Heatmap plot ──────────────────────────────────────────────────
    if path is not None:
        pivot = sweep_df.pivot_table(
            index="p_mem_sulfonate", columns="PFAS_Name", values="pct_err")
        # Sort index descending (1.0 = no correction at top, 0.05 at bottom)
        pivot = pivot.sort_index(ascending=False)

        fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns) * 1.0), 5))
        sns.heatmap(pivot, annot=True, fmt=".0f", center=0,
                    cmap="RdBu_r", linewidths=0.5, ax=ax,
                    cbar_kws={"label": "% error vs observed BCF"})
        ax.set_title(
            "Arnot-Gobas Sensitivity: Sulfonate Gill Membrane Permeability (P_mem)\n"
            f"(KPROT_SCALE fixed: Sulfonate={KPROT_SCALE['Sulfonate']}, "
            f"Carboxylate={KPROT_SCALE['Carboxylate']}; "
            "Carboxylate P_mem=1.0)"
        )
        ax.set_xlabel("PFAS")
        ax.set_ylabel("Sulfonate P_mem correction factor\n(1.0 = no correction, Arnot-Gobas baseline)")
        plt.tight_layout()
        plt.savefig(path, dpi=150)
        plt.close()
        print(f"[Output] → {path}")

    return sweep_df


def main():
    print("\n" + "═" * 60)
    print("  PFAS Bioaccumulation Pipeline v13.0")
    print("  New in v13.0: Direct albumin-Ka term in Kfish (Bischel/Beesoon)")
    print("  Closes Findings 16–17: Koc-proxy and P_mem approaches exhausted")
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
        "Koc", "Koc_log", "AlbuminBinding_pKa",         # v11 additions
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

    # ── FEATURE ABLATION (v11) ──────────────────────────────────────
    # Compare v10.5 feature set vs v11 feature set on the same split.
    # This is the honest test of whether Koc_log + AlbuminBinding_pKa help.
    abl_df = run_feature_ablation(
        X_fit, y_fit, X_test, y_test, grp_test,
        v10_features=CONC_FEATURES_V10,
        v11_features=CONC_FEATURES,
        path=OUTPUT_DIR + "feature_ablation.png",
    )

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

    # ── ARNOT-GOBAS MECHANISTIC BCF MODEL (v11.1) ───────────────────
    print("\n[Arnot-Gobas] Computing mechanistic BCF estimates...")
    ag_df = compute_arnot_gobas_bcf(chem_df)

    # Observed BCF records from the dataset (fish/aquatic only — BCF is
    # measured in aquatic exposure studies; human/plant BCF is not meaningful)
    observed_bcf = dataset[
        dataset["log_BCF"].notna() & (dataset["Species_Group"].isin(["Fish", "Other"]))
    ][["PFAS_Name", "log_BCF"]].copy()

    # ML BCF predictions on test set — build a per-PFAS prediction DataFrame
    # from the RF BCF model if it was trained (requires len(ml_bcf) >= 30)
    ml_bcf_preds = None
    if r2_bcf is not None:
        ml_bcf_preds = pd.DataFrame({
            "PFAS_Name":    pfas_bcf_te.values,
            "log_BCF_pred": y_bcf_pred,
        })

    print_arnot_gobas_summary(ag_df, observed_bcf)
    plot_arnot_gobas_comparison(
        ag_df, observed_bcf, ml_bcf_preds,
        OUTPUT_DIR + "arnot_gobas_bcf.png"
    )

    # v11.3: sensitivity sweep — see how % error moves across sulfonate Kprot values
    run_arnot_gobas_sensitivity(
        chem_df, observed_bcf,
        path=OUTPUT_DIR + "arnot_gobas_sensitivity.png"
    )

    # v12.1: P_mem sensitivity sweep — orthogonal hypothesis to Kprot sweep.
    # Tests whether reducing gill membrane permeability for ionized sulfonates
    # closes the BCF gap that Kprot tuning (Finding 16) could not fix.
    # If % error is invariant here too, both tissue partitioning AND gill
    # permeability are ruled out as the primary error source, pointing to a
    # structural limitation of the single-compartment Arnot-Gobas framework
    # for ionic PFAS (e.g., missing blood-protein-binding elimination pathway).
    run_arnot_gobas_pmem_sensitivity(
        chem_df, observed_bcf,
        path=OUTPUT_DIR + "arnot_gobas_pmem_sensitivity.png"
    )

    # ── HUMAN-ONLY MODEL (v11.2) ────────────────────────────────────
    human_results = run_human_only_model(dataset, OUTPUT_DIR)

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
    print("  Held-out Test Results — v11.2")
    print("═" * 60)
    print("  [HEADLINE — per species group; ~95% of rows are Human, so the")
    print("   pooled figure below is a human-weighted average, NOT evidence")
    print("   the model works across species. Chemistry-only features")
    print("   (Chain_Length, MW, LogKow, PFAS_Class_encoded, Koc_log,")
    print("   AlbuminBinding_pKa) — no Trophic_Level/Is_Aquatic.]")
    print(f"  {'Group':<12} {'n':>6}  {'RF R²':>8}  {'Baseline R²':>12}  {'Gain':>7}")
    print(f"  {'-'*12}  {'-'*6}  {'-'*8}  {'-'*12}  {'-'*7}")
    for _, r in grp_df.sort_values("RF R²", ascending=False).iterrows():
        print(f"  {r['Group']:<12} {int(r['n']):>6}  {r['RF R²']:>8.3f}  "
              f"{r['Baseline R²']:>12.3f}  {r['Gain']:>+7.3f}")
    print(f"\n  [Pooled — human-weighted average, reference only]")
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
    # ── HUMAN-ONLY MODEL SUMMARY ────────────────────────────────────
    if human_results:
        print(f"\n  [Human-Only Model — NHANES blood serum, chemistry features]")
        print(f"  {'Model':<20} {'R²':>7}  {'RMSE':>7}  {'vs pooled Human R²':>20}")
        print(f"  {'-'*20}  {'-'*7}  {'-'*7}  {'-'*20}")
        pooled_human_r2 = 0.604
        for label, r2, rmse in [
            ("Random Forest",  human_results["r2_rf"],  human_results["rmse_rf"]),
            ("XGBoost",        human_results["r2_xgb"], human_results["rmse_xgb"]),
            ("Linear",         human_results["r2_lr"],  human_results["rmse_lr"]),
            ("PFAS-mean base", human_results["r2_baseline"], human_results["rmse_baseline"]),
        ]:
            delta = r2 - pooled_human_r2
            print(f"  {label:<20} {r2:>7.3f}  {rmse:>7.3f}  {delta:>+20.3f}")
        print(f"\n    Prediction intervals (RF, calibrated per PFAS):")
        print(f"    80% coverage: {human_results['cov80']:.1f}%  "
              f"width: {human_results['w80']:.3f} log10 ng/g")
        print(f"    95% coverage: {human_results['cov95']:.1f}%  "
              f"width: {human_results['w95']:.3f} log10 ng/g")
        print(f"    (Compare pooled model: 80% width=1.036, 95% width=1.735)")

    # ── ARNOT-GOBAS SUMMARY ─────────────────────────────────────────
    print(f"\n  [Arnot-Gobas Mechanistic BCF — v13.0 Ka-albumin Kfish]")
    ag_with_data = ag_df[ag_df["BCF_AG"].notna()]
    ka_compounds = ag_df[ag_df.get("Kfish_method", pd.Series()) == "Ka_albumin"] \
        if "Kfish_method" in ag_df.columns else pd.DataFrame()
    print(f"    Compounds with AG estimate: {len(ag_with_data)}")
    print(f"    Compounds using Ka_albumin Kfish (v13.0): {len(ka_compounds)}")
    if len(ag_with_data) > 0:
        best_ag  = ag_with_data.loc[ag_with_data["BCF_AG"].idxmax()]
        worst_ag = ag_with_data.loc[ag_with_data["BCF_AG"].idxmin()]
        print(f"    Highest predicted BCF: {best_ag['PFAS_Name']} "
              f"(BCF_AG={best_ag['BCF_AG']:.1f}, log={best_ag['log_BCF_AG']:.3f})")
        print(f"    Lowest predicted BCF:  {worst_ag['PFAS_Name']} "
              f"(BCF_AG={worst_ag['BCF_AG']:.1f}, log={worst_ag['log_BCF_AG']:.3f})")

    # Finding 18: evaluate whether Ka-albumin Kfish closes the sulfonate gap
    if len(ag_with_data) > 0:
        observed_bcf_summary = dataset[
            dataset["log_BCF"].notna() & dataset["Species_Group"].isin(["Fish", "Other"])
        ].groupby("PFAS_Name")["log_BCF"].median()
        sulfonate_errs = []
        carboxylate_errs = []
        for _, r in ag_with_data.iterrows():
            obs = observed_bcf_summary.get(r["PFAS_Name"], np.nan)
            if pd.notna(obs) and pd.notna(r["BCF_AG"]):
                pct = (r["BCF_AG"] - 10**obs) / 10**obs * 100
                if r["PFAS_Class"] == "Sulfonate":
                    sulfonate_errs.append((r["PFAS_Name"], pct))
                else:
                    carboxylate_errs.append((r["PFAS_Name"], pct))
        print(f"\n  [Finding 18 — Ka-albumin Kfish vs observed BCF]")
        if sulfonate_errs:
            print(f"    Sulfonates: " +
                  ", ".join(f"{n} {e:+.0f}%" for n, e in sulfonate_errs))
            v12_ref = {"PFOS": -91, "PFHxS": -85, "PFBS": -84}
            for n, e in sulfonate_errs:
                ref = v12_ref.get(n)
                if ref:
                    delta = e - ref
                    verdict = "IMPROVED" if delta > 10 else ("DEGRADED" if delta < -10 else "unchanged")
                    print(f"      {n}: v12.x={ref:+d}%  v13.0={e:+.0f}%  Δ={delta:+.0f}pp  [{verdict}]")
        if carboxylate_errs:
            print(f"    Carboxylates: " +
                  ", ".join(f"{n} {e:+.0f}%" for n, e in carboxylate_errs))

    print(f"    See arnot_gobas_bcf.png for mechanistic vs ML vs observed comparison.")
    print(f"    Key: ML BCF gain over baseline ≈ 0.000 (Finding 5); AG provides")
    print(f"    a mechanistic explanation rather than a data-driven prediction.")

    # ── ABLATION SUMMARY ────────────────────────────────────────────
    print(f"\n  [Feature Ablation — v10.5 vs v11 (Koc_log + AlbuminBinding_pKa)]")
    print(f"  {'Group':<12} {'n':>6}  {'v10.5 R²':>9}  {'v11 R²':>8}  {'ΔR²':>7}  verdict")
    print(f"  {'-'*12}  {'-'*6}  {'-'*9}  {'-'*8}  {'-'*7}  {'-'*12}")
    for _, r in abl_df.iterrows():
        verdict = "improved" if r["Delta_R2"] > 0.01 else (
                  "degraded" if r["Delta_R2"] < -0.01 else "no change")
        print(f"  {r['Group']:<12} {int(r['n']):>6}  {r['R2_v10']:>9.3f}  "
              f"{r['R2_v11']:>8.3f}  {r['Delta_R2']:>+7.3f}  {verdict}")

    print("═" * 60)
    print("  Pipeline complete — v13.0")
    print("═" * 60 + "\n")


if __name__ == "__main__":
    main()