# PFAS Bioaccumulation Research Pipeline v12.0

A reproducible, multi-source data pipeline for studying PFAS bioaccumulation in human populations and aquatic/terrestrial species. Integrates EPA ECOTOX biological exposure data, EPA CompTox chemical properties, and CDC NHANES human biomonitoring data to build a machine-learning-ready dataset, identify critical data gaps, predict bioaccumulation from chemical structure alone, and model BCF from first-principles mass balance — with calibrated uncertainty and per-compound confidence.

**Current dataset: 25,056 observations | 13 curated PFAS (7 individually modelable) | 5 species groups | 5 ML models + Arnot-Gobas mechanistic BCF | Best Human R²=0.658 (human-only model) | Calibrated 80%/95% prediction intervals | Apparent half-life estimated for 4/6 modelable PFAS | Sulfonate Kprot investigation closed**

---

## Table of Contents
- [Why This Research Matters](#why-this-research-matters)
- [Key Findings](#key-findings)
- [Version History](#version-history)
- [Outputs](#outputs)
- [Dataset Schema](#dataset-schema)
- [PFAS Chemicals](#pfas-chemicals)
- [Setup](#setup)
- [Usage](#usage)
- [Data Sources](#data-sources)
- [Pipeline Architecture](#pipeline-architecture)
- [Model Results](#model-results)
- [Data Gaps](#data-gaps)
- [Roadmap](#roadmap)
- [How to Add Data](#how-to-add-data)

---

## Why This Research Matters

### PFAS Are Everywhere — And They Don't Leave
Per- and polyfluoroalkyl substances (PFAS) are a class of over 12,000 synthetic chemicals used in non-stick cookware, food packaging, firefighting foam, waterproof clothing, and hundreds of industrial applications. They are called "forever chemicals" for a reason: the carbon-fluorine bond is one of the strongest in chemistry. PFAS do not break down in the environment. They do not break down in the human body.

They accumulate.

### The Food Chain Problem
When PFAS enter an ecosystem — through industrial discharge, agricultural runoff, or contaminated groundwater — they are absorbed by plants and small organisms at the base of the food chain. As larger animals eat smaller ones, PFAS concentrations multiply at each trophic level. This process, called biomagnification, means a fish at the top of an aquatic food chain can carry concentrations thousands of times higher than the water it swims in.

Humans sit at the top of the food chain.

### What the Numbers Say
- PFAS have been detected in the blood of **97% of Americans**
- The EPA has set drinking water limits for PFAS at **4 parts per trillion** — so low it required new analytical methods to measure
- PFAS exposure has been linked to thyroid disease, immune suppression, certain cancers, reproductive harm, and developmental delays in children
- Our pipeline finds median PFOS levels of **2.83 ng/g** in human blood serum from CDC data — in people with no known occupational exposure

### The Scientific Gap We're Addressing
Despite this, our understanding of how PFAS move through ecosystems remains deeply fragmented. Data is scattered across hundreds of studies, measured in inconsistent units, tested on different species, and reported under different conditions. No single database cleanly maps PFAS bioaccumulation from soil → plant → fish → mammal → human.

**That gap is what this project addresses.**

---

## Key Findings

### Finding 1 — Trophic level is the strongest predictor of PFAS accumulation
Across 25,056 observations, trophic level explains more variance in PFAS tissue concentration than any chemical property. This directly confirms biomagnification: the higher you are in the food chain, the more PFAS you accumulate. Humans at trophic level 5 show the highest and most consistent concentrations. Note: trophic level was removed as a model feature in v10.5 because it encodes species group identity rather than a measured biological quantity — but it remains a valid descriptive finding.

### Finding 2 — Within human blood serum, PFAS accumulation is largely linear with chemistry
After removing all group-identity leakage, Linear Regression (R²=0.649) and RF/XGBoost (R²=0.658) perform nearly identically on the human-only model — confirming the relationship between chain length, LogKow, and blood concentration is genuinely near-linear. Longer-chain and more hydrophobic PFAS accumulate more in human blood, and a linear model captures most of that signal.

### Finding 3 — Human blood levels are predictable; environmental data is not
Human blood measurements collected under standardized NHANES laboratory protocols are far more consistent than heterogeneous environmental studies. The human-only model achieves R²=0.658 — but cannot transfer that knowledge to fish or plant predictions, which remain worse than a simple mean baseline.

### Finding 4 — Environmental data cannot predict human exposure
Leave-one-species-out validation shows that a model trained on fish, plant, and mammal data performs dramatically worse than baseline when predicting human blood levels. Two distinct problems are entangled: unit and matrix incompatibility (blood serum vs. tissue vs. water concentrations), and genuine biological non-transferability. Both must be solved before environmental data can inform human exposure estimates.

### Finding 5 — BCF cannot be predicted from chemical structure alone
Bioconcentration Factor (BCF) normalizes tissue concentration by exposure concentration. With chemistry-only features, both RF and XGBoost match a per-compound mean baseline exactly (gain ≈ 0.000). Study-level variance from species physiology, lab conditions, and water chemistry dominates any chemistry signal. This is a data gap finding, not a model failure.

### Finding 6 — Data gaps are severe, systematic, and chemically biased
PFHxS is detected in the blood of nearly every American — yet has zero fish and zero mammal tissue records in ECOTOX. Five of the 13 curated PFAS (GenX, ADONA, F53B, PFDoDA, PFHxA) have zero measured records anywhere despite being fully characterized in the chemical feature table. This pattern repeats across short-chain and emerging PFAS.

### Finding 7 — Mammalian bioaccumulation data is nearly unusable
ECOTOX mammal records for PFAS are almost entirely dose-response studies rather than tissue residue measurements. Only 4 mammalian tissue records exist across all 8 PFAS. Cross-species prediction from environmental mammals to humans is currently impossible.

### Finding 8 — Environmental-only prediction remains extremely difficult
Species-specific chemistry-only models produce:

| Model | Held-Out R² |
|---|---|
| Fish-only | -0.008 |
| Plant-only | -0.018 |

Both perform worse than a PFAS-mean baseline, confirming the environmental bioaccumulation dataset is too sparse and heterogeneous for reliable chemistry-based prediction.

### Finding 9 — Leakage was more significant than initially reported
The v7.0 leakage fix removed `is_human/fish/mammal/plant` flags but left `Trophic_Level` and `Is_Aquatic` — both fixed dict lookups on `Species_Group` — in the feature set. Removing all group-identity features in v10.5 drops pooled R² from 0.710 to 0.490. The ~0.22 R² that disappeared was the model learning which measurement protocol was used (NHANES vs. ECOTOX), not anything about PFAS chemistry.

### Finding 10 — Naive uncertainty estimates from tree ensembles are dangerously overconfident
Raw Random Forest tree-to-tree variance produced 80% intervals that covered only 2.0% of true values. After implementing a proper three-way Fit/Calibration/Test split with residual-based calibration, overall coverage corrected to 79–80% (80% target) and 94–95% (95% target) across all species groups.

### Finding 11 — Predictive reliability varies enormously across PFAS compounds
Of 13 curated PFAS, only 7 have enough data for a dedicated model. PFOA is the most predictable (R²=0.657, n=4,594). PFBS performs worse than predicting the mean (R²=-0.063) despite having 83 records. Five compounds have zero measured records entirely.

### Finding 12 — Apparent NHANES half-lives are inflated by ongoing population exposure
One-compartment elimination estimates from two NHANES survey waves dramatically exceed published clinical values: PFOS 42.3 years vs. 5.4 years published (+683%), PFOA 18.7 years vs. 3.5 years (+435%), PFHxS 15.9 years vs. 8.5 years (+87%), PFNA 3.4 years vs. 2.5 years (+37%). This is not a measurement error — it is the expected signature of ongoing population exposure contaminating a cross-sectional estimate. The magnitude of inflation scales with true half-life, making it a proxy signal for how much ongoing exposure remains for a given compound.

### Finding 13 — Fish and plant model failures are a data problem, not a feature problem
Adding PFAS-appropriate chemistry descriptors (Koc, AlbuminBinding_pKa) produced ΔR²=0.000 across every species group. The fish and plant model failures cannot be fixed with better features — the data is too sparse and heterogeneous for any chemistry descriptor to extract a reliable signal. More and better standardized environmental field measurements are required.

### Finding 14 — Arnot-Gobas mechanistic model reveals sulfonate-specific bioaccumulation
The PFAS-adapted Arnot-Gobas mass balance model predicts BCF accurately for mid-chain carboxylates (PFOA +5%, PFDA +2%) but systematically under-predicts all sulfonates (PFOS -91%, PFHxS -85%, PFBS -77%). This systematic bias is not random — it points to a specific mechanism: sulfonates bind serum albumin more strongly than carboxylates of equivalent chain length, and the standard protein-binding scaling factor (Kprot = 0.05 × Koc) is too conservative for sulfonates.

### Finding 15 — Within human blood serum, PFAS identity dominates over chemical structure
The dedicated human-only model (NHANES rows exclusively) achieves R²=0.658 overall — but per-PFAS R² is near zero for every individual compound (PFOA R²=-0.000, PFOS R²=-0.001, etc.). The PFAS-mean baseline also hits R²=0.658, matching RF and XGBoost exactly. The model is learning which PFAS it is, not anything about the chemistry within a compound class. Between-compound variance is large; within-compound chemistry signal is near zero at current NHANES sample sizes.

### Finding 16 — Sulfonate Kprot scaling cannot be fixed by tuning a single constant
A systematic sweep of the sulfonate Kprot scale factor across 7 values (0.05 to 0.30) shows PFOS error is essentially invariant: -91% at scale=0.05, -89% at scale=0.30. PFHxS stays at -85% throughout. The Koc-based protein binding term is structurally inadequate for sulfonates — not just under-scaled. This closes the simple-tuning approach and establishes that a direct albumin-binding-affinity term (based on measured Ka values from Bischel et al. 2010, Beesoon & Martin 2015) is the correct next mechanistic step, not further constant adjustment. The sensitivity heatmap (`arnot_gobas_sensitivity.png`) provides the quantitative evidence for this conclusion.

---

## Version History

### v12.0 (current) — July 2026
- **Sulfonate Kprot sensitivity sweep:** `run_arnot_gobas_sensitivity()` tests 7 Kprot scale values (0.05–0.30) for sulfonates, holding carboxylate scale fixed at 0.05. Sweep runs at pipeline startup automatically.
- **Class-specific `KPROT_SCALE` dict:** Sulfonate=0.15, Carboxylate=0.05 — implemented as a named constant replacing the previous hardcoded value, making future per-class tuning explicit.
- **Finding 16:** Sweep proves Kprot scaling is structurally inadequate for sulfonates — PFOS % error moves only 2 percentage points across the full sweep range. Simple constant tuning is a dead end; a direct albumin-binding-affinity term is the required next step.
- **Kprot investigation closed** as a standalone tuning task. Finding 16 documented. Roadmap updated accordingly.
- New output: `arnot_gobas_sensitivity.png` — heatmap of % error per PFAS × Kprot scale factor.

### v11.2 — July 2026
- **Human-only model:** `run_human_only_model()` trains RF, XGBoost, and Linear Regression on NHANES blood serum rows exclusively. Features: Chain_Length, MW, LogKow, PFAS_Class_encoded, AlbuminBinding_pKa (Koc_log excluded — governs soil/fish uptake, not blood serum). Overall R² improves 0.604 → 0.658; interval width narrows from 1.036 → 0.858 log10 ng/g at 80%.
- **Finding 15:** Human-only model reveals the pooled model gain was driven by between-compound variance, not within-compound chemistry — PFAS identity dominates, per-PFAS R² ≈ 0.000 for all compounds.
- New outputs: `human_model_predictions.png`, `human_model_feature_importance.png`

### v11.1 — July 2026
- **Arnot-Gobas mechanistic BCF model:** PFAS-adapted steady-state mass balance (Arnot & Gobas 2004). Replaces Kow-based lipid partitioning with Kprot = 0.05 × Koc (Kelly et al. 2004) and NLOM = 0.035 × Koc (Gobas et al. 2003). Sets km = 0 (PFAS metabolically inert). Includes both gill (k1) and dietary (kd) uptake pathways.
- **Finding 14:** Mechanistic model accurately predicts PFOA (+5%) and PFDA (+2%) but systematically under-predicts all sulfonates (PFOS -91%, PFHxS -85%), revealing sulfonate-specific protein binding not captured by Koc alone.
- New output: `arnot_gobas_bcf.png`

### v11.0 — July 2026
- **PFAS-appropriate chemistry features:** Koc and AlbuminBinding_pKa added to PFAS_FEATURES table and ML feature sets.
- **Feature ablation:** `run_feature_ablation()` compares v10.5 vs v11 features — ΔR²=0.000 across all groups.
- **Finding 13:** New features added zero predictive value — confirmed fish/plant failure is a data problem, not a feature problem.
- New output: `feature_ablation.png`

### v10.5 — July 2026
- CASRN salt mapping: recovered 396/401 previously dropped ECOTOX rows via `CASRN_SALT_MAP`
- Leakage fix (FIX 4): removed `Trophic_Level` and `Is_Aquatic` from model features
- Headline metric changed from pooled R² to per-species-group R²
- PFHxS per-PFAS R² improved 0.172 → 0.381 from recovered salt rows

### v10.0 — July 2026
- Apparent half-life estimation: `compute_apparent_half_life()` from two NHANES wave medians
- Finding 12: NHANES apparent half-lives dramatically exceed literature values

### v9.0 — June–July 2026
- Per-PFAS models: dedicated RF for each PFAS with n≥60 rows (7 compounds)

### v8.0 — June 2026
- Residual-calibrated prediction intervals replacing naive tree-variance
- Three-way Fit/Calibration/Test split; coverage corrected to 80.6%/94.8%

### v7.0 — June–July 2026
- XGBoost models added; leakage fixes for is_human/fish/mammal/plant flags and Duration_days; stratified split by Species_Group; NHANES 2015–2016 added

### v6.0 — July 2026
- Proper 80/20 held-out validation; 13-PFAS feature table; NHANES 2017–2018

### v5.0–v1.0 — April–June 2026
- v5.0: Linear Regression baseline; v4.0: CDC NHANES 2017–2018; v3.0: BCF model; v2.x: terrestrial data, Tier 2 PFAS; v1.0: PFOS only (411 rows)

---

## Outputs

| File | Description |
|---|---|
| `pfas_bioaccumulation_dataset.csv` | 25,056 rows, fully cleaned and ML-ready |
| `pfas_gap_heatmap.png` | Observations per PFAS × species group |
| `feature_importance.png` | RF concentration model feature importances |
| `model_predictions.png` | RF concentration predicted vs actual (pooled R²=0.490) |
| `human_model_predictions.png` | Human-only RF predicted vs actual + per-PFAS R² comparison |
| `human_model_feature_importance.png` | Human-only RF feature importances |
| `arnot_gobas_bcf.png` | Mechanistic BCF vs ML BCF vs observed median per PFAS |
| `arnot_gobas_sensitivity.png` | **NEW (v12)** Heatmap: % error per PFAS × sulfonate Kprot scale factor |
| `feature_ablation.png` | v10.5 vs v11 features ΔR² per species group |
| `prediction_intervals.png` | RF prediction ribbon (80%/95% intervals) by species group |
| `interval_coverage.png` | Calibration diagnostic — actual vs target coverage by species group |
| `per_pfas_r2.png` | Held-out R² per individual PFAS compound |
| `per_group_metrics.png` | Held-out R² and baseline comparison by species group |
| `cross_species_validation.png` | Leave-one-species-out R² and RMSE |
| `bcf_feature_importance.png` | RF BCF model feature importances |
| `bcf_predictions.png` | RF BCF predicted vs actual |
| `bcf_xgb_predictions.png` | XGBoost BCF predicted vs actual |
| `linear_coefficients.png` | Linear regression standardized coefficients |
| `xgboost_predictions.png` | XGBoost concentration predicted vs actual |
| `model_comparison.png` | Side-by-side R² and RMSE across all models |
| `fish_predictions.png` | Fish-only chemistry model (R²=-0.008) |
| `plant_predictions.png` | Plant-only chemistry model (R²=-0.018) |
| `chain_length_bcf_scatter.png` | Chain length vs BCF by PFAS class |
| `nhanes_time_trend.png` | NHANES PFAS blood concentration trends 2015–2018 |
| `nhanes_half_life.png` | NHANES apparent population half-life vs published literature |

---

## Dataset Schema

| Column | Type | Description |
|---|---|---|
| `PFAS_Name` | string | Common name (PFOS, PFOA, etc.) |
| `CASRN` | string | CAS Registry Number |
| `PFAS_Class` | string | Sulfonate or Carboxylate |
| `PFAS_Class_encoded` | int | 0=Sulfonate, 1=Carboxylate |
| `Chain_Length` | int | Fluorinated carbon chain length |
| `MW` | float | Molecular weight (g/mol) |
| `LogKow` | float | Octanol-water partition coefficient |
| `Koc` | float | Soil organic-carbon/water partition coefficient (L/kg) |
| `Koc_log` | float | log10(Koc) — ML feature |
| `AlbuminBinding_pKa` | float | Functional-group pKa proxy for serum protein binding |
| `Species` | string | Common species name |
| `Species_Group` | string | Fish / Mammal / Plant / Human / Other |
| `Trophic_Level` | int | 1 (plant) → 5 (human) — descriptive only, not a model feature |
| `Is_Aquatic` | int | 1=aquatic, 0=terrestrial — descriptive only, not a model feature |
| `Tissue` | string | Tissue or response site |
| `Exposure Route` | string | Exposure type |
| `Duration_days` | float | Exposure duration in days (retained in dataset; excluded from ML) |
| `Concentration_ng_g` | float | Measured tissue concentration (ng/g) |
| `log_concentration` | float | log₁₀ concentration — ML target 1 |
| `BCF` | float | Bioconcentration factor |
| `log_BCF` | float | log₁₀ BCF — ML target 2 |
| `Source` | string | CDC NHANES / ECOTOX |

---

## PFAS Chemicals

### Tier 1 — Core
| PFAS | CASRN | Class | Chain | MW | LogKow | Koc (L/kg) |
|---|---|---|---|---|---|---|
| PFOS | 1763-23-1 | Sulfonate | 8 | 500.1 | 5.26 | 2100 |
| PFOA | 335-67-1 | Carboxylate | 8 | 414.1 | 5.30 | 1900 |
| PFHxS | 355-46-4 | Sulfonate | 6 | 400.1 | 4.14 | 560 |
| PFNA | 375-95-1 | Carboxylate | 9 | 464.1 | 6.05 | 3200 |
| PFBS | 375-73-5 | Sulfonate | 4 | 300.1 | 1.82 | 47 |

### Tier 2 — Extended
| PFAS | CASRN | Class | Chain | MW | LogKow | Koc (L/kg) |
|---|---|---|---|---|---|---|
| PFDA | 335-76-2 | Carboxylate | 10 | 514.1 | 6.83 | 5800 |
| PFUnDA | 2058-94-8 | Carboxylate | 11 | 564.1 | 7.59 | 9100 |
| PFDoDA | 307-55-1 | Carboxylate | 12 | 614.1 | 8.35 | 14000 |
| PFHpA | 375-85-9 | Carboxylate | 7 | 364.1 | 4.55 | 820 |
| PFHxA | 307-24-4 | Carboxylate | 6 | 314.1 | 3.77 | 310 |

### Tier 3 — Emerging (zero measured records)
| PFAS | CASRN | Class | Chain | MW | LogKow | Koc (L/kg) |
|---|---|---|---|---|---|---|
| GenX (HFPO-DA) | 13252-13-6 | Carboxylate | 6 | 330.1 | 2.50 | — |
| ADONA | 958445-44-8 | Carboxylate | 8 | 380.1 | 2.80 | — |
| F53B | 73606-19-6 | Sulfonate | 6 | 570.1 | 4.00 | — |

---

## Setup

```bash
pip3 install pandas numpy matplotlib seaborn scikit-learn openpyxl requests xgboost
```

On macOS, XGBoost requires OpenMP:
```bash
brew install libomp
```

---

## Usage

### Set paths in `pfas_pipeline_v12.py`
```python
ECOTOX_EXPORT_DIR = "/path/to/ecotox_exports/"
COMPTOX_SNAPSHOT  = "/path/to/comptox_snapshot.csv"
NHANES_PATH       = "/path/to/nhanes_pfas_processed.csv"
OUTPUT_DIR        = "/path/to/outputs/"
```

### Run
```bash
python3 pfas_pipeline_v12.py
```

### Required files
| File | Where to get it |
|---|---|
| ECOTOX xlsx exports | https://cfpub.epa.gov/ecotox/ — search each PFAS, export XLSX |
| `nhanes_pfas_processed.csv` | Included in this repo |
| `comptox_snapshot.csv` | Optional — https://comptox.epa.gov/dashboard/batch-search |

---

## Data Sources

| Source | URL | What it provides |
|---|---|---|
| EPA ECOTOX | https://cfpub.epa.gov/ecotox/ | Species, tissue, concentration, BCF |
| EPA CompTox | https://comptox.epa.gov/dashboard/batch-search | MW, LogKow, chemical properties |
| CDC NHANES 2015–2016 & 2017–2018 | https://wwwn.cdc.gov/nchs/nhanes/ | Human blood serum PFAS levels |

---

## Pipeline Architecture

```
EPA CompTox          EPA ECOTOX              CDC NHANES
(chemical traits)  (species + tissue)     (human blood serum)
      │                   │                      │
      ▼                   ▼                      │
Chemical Feature    ECOTOX Ingestion             │
Table (13 PFAS)     (18 xlsx files)              │
      │                   │                      │
      │            Harmonization                 │
      │        (units → ng/g, taxonomy,          │
      │         BCF coalescing, CASRN            │
      │         salt remapping)                  │
      │                   │                      │
      │                   └──────────┬───────────┘
      │                              │
      │                       Combined Dataset
      └──────────────────────────────┘
                             │
                      Merged Dataset
                      (25,056 rows)
                             │
              ┌──────────────┼──────────────────┐
              ▼              ▼                  ▼
       Pooled Model   Human-Only Model   Arnot-Gobas BCF
       (all species)  (NHANES only)      (mechanistic)
              │              │                  │
     Stratified Split   Stratified Split   Kprot Sensitivity
     by Species_Group   by PFAS_Name       Sweep (v12) ←NEW
              │
     Fit / Cal / Test (3-way)
              │
     RF | XGBoost | Linear
              │
     Residual-calibrated
     prediction intervals
     (per species group)
              │
     ┌────────┼──────────────┬──────────────┐
     ▼        ▼              ▼              ▼
  Gap      Per-Group     BCF Models    Per-PFAS
  Heatmap  Metrics       RF + XGB      Models (7)
```

---

## Model Results

### Headline — Per-Group Performance (v12.0, chemistry-only features)

| Group | n | RF R² | Baseline R² | Gain |
|---|---|---|---|---|
| Human | 4,707 | **0.604** | 0.000 | +0.604 |
| Other | 198 | -1.059 | -0.001 | -1.057 |
| Fish | 49 | -1.680 | -0.012 | -1.668 |
| Plant | 55 | -8.429 | -0.005 | -8.424 |

Human blood serum is reliably predictable from chemistry alone. Fish and Plant are not — environmental data heterogeneity dominates any chemistry signal. Mammal excluded (only 2 test rows).

### Human-Only Model (v11.2)

| Model | R² | RMSE | vs Pooled Human R² |
|---|---|---|---|
| Random Forest | **0.658** | 0.360 | +0.054 |
| XGBoost | 0.658 | 0.360 | +0.054 |
| Linear Regression | 0.649 | 0.365 | +0.045 |
| PFAS-mean Baseline | 0.658 | 0.360 | +0.054 |

Note: RF matches the PFAS-mean baseline exactly — the model is learning per-compound averages, not chemistry within compounds (Finding 15). Interval width narrows vs pooled model: 80% width 0.858 vs 1.036 log10 ng/g.

### Feature Set (v12.0)

**Pooled model (6 features):**

| Feature | Description |
|---|---|
| `Chain_Length` | Fluorinated carbon chain length |
| `MW` | Molecular weight (g/mol) |
| `LogKow` | Octanol-water partition coefficient |
| `PFAS_Class_encoded` | 0=Sulfonate, 1=Carboxylate |
| `Koc_log` | log10(Koc) — soil OC/water partition |
| `AlbuminBinding_pKa` | Functional-group pKa proxy for protein binding |

**Human-only model (5 features):** same minus `Koc_log` (irrelevant for blood serum partitioning).

### Pooled Held-Out Results (v12.0, human-weighted)

| Model | Pooled R² | RMSE |
|---|---|---|
| Random Forest | 0.490 | 0.636 |
| XGBoost | 0.490 | 0.636 |
| Linear Regression | 0.479 | 0.642 |
| Group Mean Baseline | 0.413 | 0.682 |
| RF (BCF) | 0.240 | 0.877 |
| XGBoost (BCF) | 0.240 | 0.877 |
| BCF Per-PFAS Baseline | 0.241 | 0.877 |

### Prediction Intervals (v12.0)

| Group | 80% width | 95% width | 80% coverage | 95% coverage |
|---|---|---|---|---|
| Human | 0.900 | 1.530 | 79.3% ✓ | 94.6% ✓ |
| Fish | 2.891 | 5.088 | 79.6% ✓ | 93.9% ✓ |
| Plant | 1.293 | 1.859 | 78.2% ✓ | 96.4% ✓ |
| Other | 3.753 | 5.723 | 80.3% ✓ | 96.0% ✓ |
| **Overall** | **1.036** | **1.735** | **79.3% ✓** | **94.6% ✓** |

### Per-PFAS Predictability (v12.0)

| PFAS | n | Model Type | Held-Out R² |
|---|---|---|---|
| PFOA | 4,594 | own model | **0.657** |
| PFNA | 3,978 | own model | 0.401 |
| PFOS | 4,496 | own model | 0.436 |
| PFHxS | 4,004 | own model | 0.381 |
| PFDA | 3,949 | own model | 0.226 |
| PFUnDA | 3,946 | own model | 0.148 |
| PFBS | 83 | own model | -0.063 (worse than baseline) |
| PFHpA | 1 | insufficient data | — |
| GenX | 0 | no data | — |
| ADONA | 0 | no data | — |
| F53B | 0 | no data | — |
| PFDoDA | 0 | no data | — |
| PFHxA | 0 | no data | — |

### Arnot-Gobas Mechanistic BCF (v11.1) + Kprot Sensitivity (v12.0)

PFAS-adapted steady-state mass balance, generic fish (1 kg, 12°C). Kprot = class-specific (Carboxylate: 0.05 × Koc; Sulfonate: 0.15 × Koc), NLOM = 0.035 × Koc (Gobas et al. 2003), km = 0.

| PFAS | log BCF_AG | BCF_AG | log BCF_obs | % error |
|---|---|---|---|---|
| PFOA | 0.973 | 9.4 | 0.953 | +5% |
| PFDA | 1.910 | 81.4 | 1.903 | +2% |
| PFNA | 1.429 | 26.9 | 1.863 | -63% |
| PFOS | 1.016 | 10.4 | 2.000 | -90% |
| PFHxS | 0.359 | 2.3 | 1.187 | -85% |
| PFBS | 0.086 | 1.2 | 0.874 | -84% |
| PFUnDA | 2.348 | 222.7 | 2.980 | -77% |

Carboxylates (PFOA, PFDA) predicted accurately. Sulfonates remain systematically under-predicted across all tested Kprot scale values (Finding 16) — simple scaling is a dead end. A direct albumin-binding-affinity term is the required next mechanistic step.

---

## Data Gaps

| PFAS | Fish | Human | Mammal | Plant | Other |
|---|---|---|---|---|---|
| PFOS | 22 | 1,929 | 2 | 130 | 121 |
| PFOA | 150 | 1,929 | 2 | 123 | 352 |
| PFNA | 16 | 1,929 | 0 | 0 | 40 |
| PFHxS | 0 | 1,929 | 0 | 0 | 46 |
| PFDA | 0 | 1,929 | 0 | 0 | 27 |
| PFUnDA | 2 | 1,929 | 0 | 0 | 17 |
| PFBS | 6 | 0 | 0 | 0 | 66 |
| PFHpA | 0 | 0 | 0 | 0 | 1 |

**Critical gap: PFHxS has 1,929 human blood measurements and zero fish records.**

---

## Roadmap

### Phase 1 — Data Expansion ✅ Complete
- [x] PFAS class as ML feature
- [x] BCF as second target variable
- [x] Terrestrial species data
- [x] Tier 2 PFAS
- [x] CDC NHANES human blood serum data (2015–2016 and 2017–2018)

### Phase 2 — Better Models ✅ Complete
- [x] Linear Regression, Random Forest, XGBoost concentration and BCF models
- [x] Held-out validation framework (stratified 80/20 split)
- [x] Leakage audit and fix (all group-identity features removed)
- [x] Residual-calibrated prediction intervals (3-way split)
- [x] Per-PFAS models (7 own models)
- [x] Apparent NHANES half-life estimation (Finding 12)
- [x] CASRN salt mapping — recovered 396/401 previously dropped rows
- [x] Per-group R² as headline metric
- [x] PFAS-appropriate chemistry features: Koc, AlbuminBinding_pKa (Finding 13)
- [x] Feature ablation: confirmed ΔR²=0.000 — data problem not feature problem
- [x] Arnot-Gobas PFAS-adapted mechanistic BCF model (Finding 14)
- [x] Human-only model — NHANES rows exclusively (Finding 15)

### Phase 3 — Outputs & Extensions (In Progress)
- [x] Sulfonate Kprot sensitivity sweep — `run_arnot_gobas_sensitivity()` with 7 scale values; Finding 16 closes simple tuning as a viable approach
- [ ] Treatment removal efficiency module — GAC/AER/RO removal as function of chain length and PFAS class; explains why short-chain GenX replacements are harder to remove
- [ ] Albumin-binding-affinity term for sulfonates — replace `Kprot = scale × Koc` with a direct Ka-based partition coefficient for sulfonates (Bischel et al. 2010, Beesoon & Martin 2015); the correct mechanistic follow-up to Finding 16
- [ ] Longitudinal half-life validation — individual-level cohort data to resolve exposure-vs-elimination conflation (Finding 12 follow-up)
- [ ] Interactive bioaccumulation simulator (per-PFAS models + calibrated intervals)
- [ ] Streamlit interactive dashboard
- [ ] Research poster
- [ ] Full methods + results report
- [ ] Published GitHub repository

---

## How to Add Data

### New ECOTOX exports
```bash
# Download from https://cfpub.epa.gov/ecotox/
# Filter: Effect → Accumulation if >10,000 results
# Export as XLSX, then:
mv ~/Downloads/ECOTOX-*.xlsx /path/to/ecotox_exports/
python3 pfas_pipeline_v12.py
```

### New PFAS chemicals
Add to `PFAS_FEATURES` in `pfas_pipeline_v12.py`:
```python
("PFDA", "335-76-2", "Carboxylate", 10, 514.1, 6.83, 5800.0, 0.30),
# (Name, CASRN, Class, Chain_Length, MW, LogKow, Koc, AlbuminBinding_pKa)
```

---

## Citation
- EPA ECOTOX Knowledgebase: https://cfpub.epa.gov/ecotox/
- EPA CompTox Dashboard: https://comptox.epa.gov/dashboard/
- CDC NHANES 2015–2016 and 2017–2018: https://wwwn.cdc.gov/nchs/nhanes/
- Arnot & Gobas (2004) Environ. Toxicol. Chem. 23:1523–1532
- Kelly et al. (2004) Environ. Sci. Technol.
- Gobas et al. (2003) Environ. Sci. Technol.
- Guelfo & Higgins (2013) Environ. Sci. Technol.
- Bischel et al. (2010) Environ. Sci. Technol.
- Beesoon & Martin (2015) Environ. Sci. Technol.
- ATSDR Toxicological Profile for Perfluoroalkyls (2021)

---

## Author
PFAS Environmental Informatics Research Project (2026)
