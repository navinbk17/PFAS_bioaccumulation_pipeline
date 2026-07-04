# PFAS Bioaccumulation Research Pipeline v10.5
A reproducible, multi-source data pipeline for studying PFAS bioaccumulation across aquatic and terrestrial species and human populations. Integrates EPA ECOTOX biological exposure data, EPA CompTox chemical properties, and CDC NHANES human biomonitoring data to build a machine-learning-ready dataset, identify critical data gaps, and predict bioaccumulation from chemical structure alone — with calibrated uncertainty and per-compound confidence.
**Current dataset: 25,056 observations | 13 curated PFAS (7 individually modelable) | 5 species groups | 5 ML models | Best Human R²=0.604 (chemistry-only features) | Calibrated 80%/95% prediction intervals | Apparent half-life estimated for 4/6 modelable PFAS**

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
- [How to Push to GitHub](#how-to-push-to-github)
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
**That gap is what this project addresses.

---
## Key Findings
### Finding 1 — Trophic level is the strongest predictor of PFAS accumulation
Across 25,056 observations, trophic level (where an organism sits in the food chain) explains more variance in PFAS tissue concentration than any chemical property. This directly confirms biomagnification: the higher you are in the food chain, the more PFAS you accumulate. Humans at trophic level 5 show the highest and most consistent concentrations.
### Finding 2 — Within human blood serum, PFAS accumulation is largely linear with chemistry
After removing group-identity leakage (v10.5), Linear Regression (R²=0.471) and XGBoost (R²=0.490) remain close on the chemistry-only feature set — but this narrow gap should not be over-interpreted. Both models are predicting primarily within the human cluster (~95% of rows), where the relationship between chain length, LogKow, and blood concentration is genuinely near-linear. The earlier claim that this validated classical bioaccumulation theory across species was an artifact of the trophic level leakage — predicting group means is an inherently linear problem regardless of the underlying biology. The honest interpretation is that within standardized human serum data, longer-chain and more hydrophobic PFAS accumulate more, and a linear model captures most of that signal.
### Finding 3 — Human blood levels are predictable; environmental data is not
When CDC NHANES human biomonitoring data was introduced into the pipeline, model performance increased dramatically relative to environmental-only datasets. Human blood measurements collected under standardized laboratory protocols are far more consistent than heterogeneous environmental studies. The model learned human PFAS patterns well — but cannot transfer that knowledge to fish or plant predictions.
### Finding 4 — Environmental data cannot predict human exposure (cross-species performance far worse than baseline)
Leave-one-species-out validation shows that a model trained on fish, plant, and mammal data performs dramatically worse than baseline when predicting human blood levels. This is a critical finding — but two distinct problems are entangled in that number. First, human measurements are blood serum (ng/mL) while fish are tissue (ng/g) and some inputs are water concentrations; these are different physical quantities, so part of the failure measures unit and matrix incompatibility rather than purely biological non-transferability. Second, even after accounting for that, the environmental-to-human knowledge gap is real: measurement protocols, exposure routes, and biological matrices differ so fundamentally that a model trained on environmental residue studies cannot reliably predict standardized clinical biomonitoring. Both problems need to be solved — not just the biology — before environmental data can inform human exposure estimates.
### Finding 5 — BCF cannot be predicted from chemical structure alone
Bioconcentration Factor (BCF) normalizes tissue concentration by exposure concentration. With chemistry-only features and a proper per-PFAS baseline, both RF and XGBoost achieve R²≈0.331 — essentially identical to simply knowing the average BCF for that compound. Study-level variance (species physiology, lab conditions, water chemistry) dominates any chemistry signal. This is a data gap finding, not a model failure.
### Finding 6 — Data gaps are severe, systematic, and chemically biased
PFHxS is detected in the blood of nearly every American — yet has zero fish and zero mammal tissue records in ECOTOX. We know it's in human blood but have almost no data on how it gets there through the food chain. This pattern repeats across short-chain and emerging PFAS.
### Finding 7 — Mammalian bioaccumulation data is nearly unusable
ECOTOX mammal records for PFAS are almost entirely dose-response studies rather than tissue residue measurements. Only 4 mammalian tissue records exist across all 8 PFAS. Cross-species prediction from environmental mammals to humans is currently impossible.
### Finding 8 — Environmental-only prediction remains extremely difficult
To test whether PFAS chemistry alone can explain accumulation patterns, species-specific models were trained using only chemical descriptors (chain length, molecular weight, LogKow, PFAS class, and exposure duration).
Results were poor:
| Model | Held-Out R² |
|--------|--------|
| Fish-only | 0.113 |
| Plant-only | 0.003 |
For fish, a simple PFAS-average baseline actually outperformed the chemistry-only model. This suggests current environmental bioaccumulation datasets remain too sparse and heterogeneous for reliable within-species prediction. More standardized residue measurements are needed before chemistry-based ecological forecasting becomes feasible.
### Finding 9 — Leakage was more significant than initially reported; chemistry-only R² is the honest number
The v7.0 leakage fix (removing `is_human`, `is_fish`, `is_mammal`, `is_plant`, `Duration_days`) appeared to have minimal impact (R² 0.712 → 0.710) because `Trophic_Level` and `Is_Aquatic` — both fixed dict lookups on `Species_Group` — remained in the feature set, consolidating the same group-identity signal into two integer columns. Removing all group-identity features in v10.5 drops pooled R² to 0.490, with Human R²=0.604. The ~0.16 R² that disappeared was not biological signal — it was the model learning which measurement protocol was used (NHANES vs. ECOTOX), not anything about PFAS chemistry. The current chemistry-only feature set (Chain_Length, MW, LogKow, PFAS_Class_encoded) is fully generalisable to new PFAS predictions.
### Finding 10 — Naive uncertainty estimates from tree ensembles are dangerously overconfident
Random Forest tree-to-tree prediction variance — a commonly used shortcut for uncertainty — produced 80% intervals that covered only 2.0% of true values and 95% intervals that covered only 3.3%. Trees trained on overlapping bootstrap samples of the same data agree with each other even when collectively wrong, making this approach unsafe for any decision-facing tool. After implementing a proper three-way Fit/Calibration/Test split with residual-based calibration, overall coverage corrected to 80.6% (80% target) and 94.8% (95% target). Species groups differ sharply in uncertainty: Human predictions carry tight intervals (0.82 log10 ng/g at 80%) reflecting thousands of standardized NHANES samples, while Fish, Plant, and especially uncategorized "Other" records carry intervals 2–4x wider, honestly reflecting how little reliable data exists for those groups.
### Finding 11 — Predictive reliability varies enormously across PFAS compounds
Of the 13 chemically-characterized PFAS in the curated feature table, only 6 (PFOA, PFOS, PFNA, PFDA, PFHxS, PFUnDA) have sufficient monitoring data to support an individually-trained, meaningfully predictive model. PFOA is the most predictable (R²=0.649, n=4,549) and also the most data-rich. PFBS, despite having 72 records — enough to clear the modeling threshold — performs *worse than predicting the mean* (R²=-0.140), indicating its short-chain environmental behavior is not well captured by chain-length/LogKow chemistry alone at current sample sizes. PFHpA has essentially no usable data (n=1). Five chemically-characterized compounds — GenX, ADONA, F53B, PFDoDA, and PFHxA — have **zero** measured records in either ECOTOX or NHANES despite being fully characterized in the chemical feature table. In short: less than half of the PFAS this pipeline can theoretically predict for actually have the data to back up that prediction.

### Finding 12 — Apparent NHANES half-lives are inflated by ongoing population exposure, not a measurement error
A one-compartment first-order elimination estimate (`compute_apparent_half_life()`, v10.0) was computed from the two available NHANES survey-wave medians (2015–2016, 2017–2018) and compared against published longitudinal cohort half-lives. The apparent half-lives came out dramatically longer than the literature values: PFOS 42.3 years vs. a published 5.4 years (+683%), PFOA 18.7 years vs. 3.5 years (+435%), PFHxS 15.9 years vs. 8.5 years (+87%), and PFNA 3.4 years vs. 2.5 years (+37%). PFDA and PFUnDA showed non-decreasing concentration between waves and could not be estimated at all.

This is not a bug — it's the expected signature of a cross-sectional population estimate contaminated by ongoing exposure. A true elimination half-life assumes no new intake during the measurement window; NHANES respondents keep being exposed to PFAS (drinking water, food packaging, etc.) throughout the survey period, so population blood levels decline far more slowly than any individual's true clearance rate would predict. Critically, the size of the distortion scales with the true half-life itself — PFOS and PFOA (longest true half-lives, most persistent) show the largest inflation, while PFNA (shortest true half-life) shows the smallest — because a slower-clearing chemical gives ongoing exposure more time to mask the decline within a fixed 2-year window. PFDA and PFUnDA, both long-chain and both essentially flat between waves, represent the limiting case where ongoing exposure and elimination roughly cancel out.

This reframes the apparent-vs-literature gap as a signal rather than noise: the magnitude of inflation is a rough proxy for how much ongoing population-level exposure remains for a given compound, independent of any explicit exposure measurement. It also reinforces Finding 4 from a different angle — cross-sectional human biomonitoring data cannot cleanly separate "how the body handles this chemical" from "how exposed the population currently is," the same environmental-to-human disconnect that shows up in the cross-species prediction failure, arrived at here through kinetics rather than cross-species ML.

---
## Version History

### v10.5 (current) — July 2026
- **CASRN salt mapping:** added `CASRN_SALT_MAP` to remap salt/counterion variants (e.g. PFOS potassium salt `2795-39-3`) to their parent compound CASRN before the ECOTOX merge. Previously, 401 rows were silently dropped because ECOTOX stores these under different CASRNs than the parent compounds in the curated feature table.
- **CASRN diagnostic:** added `diagnose_casrn_merge()` — runs before every merge and prints a table of unmatched CASRNs with counts, so future data gaps are visible rather than silent. Unmatched rows drop from 401 → 5 (2 genuinely unknown CASRNs, 3 rows total — not worth mapping).
- **Leakage fix (FIX 4):** removed `Trophic_Level` and `Is_Aquatic` from `CONC_FEATURES`. Both are fixed dict lookups keyed 1:1 on `Species_Group` — the same group-identity leak FIX 1 (v7.0) removed for the binary flags, just re-encoded as integers. Feature set is now chemistry-only: `Chain_Length`, `MW`, `LogKow`, `PFAS_Class_encoded`.
- **Headline metric change:** pooled held-out R² replaced by per-species-group R² as the headline figure. The dataset is ~95% human rows, so a single pooled number was almost entirely a human-serum number mislabeled as "overall." Per-group breakdown (Human R²=0.604, Fish R²=-1.68, Plant R²=-8.4) is now the honest lead.
- **NHANES wave label fix:** `NHANES_CYCLE_MIDPOINT` now includes both string keys (`"2015-2016"`) and integer keys (`2016`) to match the actual `Study_Year` column format in the NHANES file. Previously the lookup silently fell back to assuming 2-year spacing — the half-life numbers were correct by coincidence, but the warning `Unrecognized wave labels (2016, 2018)` has been eliminated.
- **Key results after fixes:** Human R²=0.604 (chemistry-only, honest); Fish/Plant remain unpredictable from chemistry alone — confirms environmental heterogeneity, not a model bug. All 5 species groups now hit calibration coverage targets (79–80% at 80% target). PFHxS per-PFAS R² improved 0.172 → 0.381 from recovered salt rows.

### v10.0 — July 2026
- Apparent half-life estimation: `compute_apparent_half_life()` — one-compartment first-order elimination per PFAS from two NHANES wave medians (2015–2016, 2017–2018). Closed-form two-point solve; no curve-fitting needed.
- Apparent half-lives dramatically exceed literature values (PFOS +683%, PFOA +435%, PFHxS +87%, PFNA +37%) — documented as Finding 12, a signature of ongoing population exposure contaminating a cross-sectional estimate, not a code defect.
- PFDA and PFUnDA flagged `non_decreasing_between_waves` rather than forced into a negative half-life.
- **Caveat:** this is a cross-sectional apparent population half-life, not a clinical elimination rate. Treat as a literature sanity-check and exposure-trend signal only.
- New output: `nhanes_half_life.png`

### v9.0 — June–July 2026
- Per-PFAS models: dedicated RF model for each PFAS with n≥60 rows (6 compounds). PFBS gets its own model but underperforms the mean baseline (R²=-0.140) — flagged, not hidden. Zero-data compounds (GenX, ADONA, F53B, PFDoDA, PFHxA) explicitly reported rather than silently omitted.
- New output: `per_pfas_r2.png`

### v8.0 — June 2026
- Residual-calibrated prediction intervals replacing naive tree-variance (which covered only 2% of values at an 80% target). Three-way Fit/Calibration/Test split; calibration computed separately per species group. Final coverage: 80.6%/94.8%.
- New outputs: `prediction_intervals.png`, `interval_coverage.png`

### v7.0 — June–July 2026
- XGBoost concentration + BCF models added.
- Leakage fixes: removed `is_human/fish/mammal/plant` flags and `Duration_days` from concentration features (both encoded data-source identity, not biology).
- Stratified train/test split by `Species_Group`.
- BCF coalescing: reads `BCF 1/2/3 Value`, recovering 321 previously missed rows.
- NHANES 2015–2016 added; dataset grew from 13,098 → 25,056 rows.

### v6.0 — July 2026
- Proper 80/20 held-out validation; group-mean baseline; NHANES 2017–2018; 13-PFAS feature table.

### v5.0 — June 2026
- Linear Regression baseline added. Finding: RF barely outperforms Linear — relationship is largely linear.

### v4.0 — June 2026
- CDC NHANES 2017–2018 added (11,574 rows). Dataset grew 1,524 → 13,098 rows. Discovered cross-species prediction failure (human R²= -8 in leave-one-out).

### v3.0 — June 2026
- BCF as second ML target; PFAS class (Sulfonate/Carboxylate) as feature.

### v2.1 — May 2026
- Terrestrial ECOTOX exports; Tier 2 PFAS; Route_encoded leakage removed.

### v2.0 — May 2026
- Aquatic exports for PFNA, PFHxS, PFBS, PFOA; first gap heatmap.

### v1.0 — April 2026
- Initial pipeline, PFOS only (411 rows), first Random Forest model.

---
## Outputs
| File | Description |
|---|---|
| `pfas_bioaccumulation_dataset.csv` | 25,056 rows, fully cleaned and ML-ready |
| `pfas_gap_heatmap.png` | Observations per PFAS × species group |
| `feature_importance.png` | RF concentration model feature importances |
| `model_predictions.png` | RF concentration predicted vs actual (Held-Out R²=0.705) |
| `prediction_intervals.png` | RF prediction ribbon (80%/95% intervals) + interval width by species group |
| `interval_coverage.png` | Calibration diagnostic — actual vs target coverage by species group |
| `per_pfas_r2.png` | Held-out R² per individual PFAS compound, colour-coded by model type |
| `cross_species_validation.png` | Leave-one-species-out R² and RMSE |
| `bcf_feature_importance.png` | RF BCF model feature importances (chemistry-only) |
| `bcf_predictions.png` | RF BCF predicted vs actual (Held-Out R²=0.330) |
| `bcf_xgb_feature_importance.png` | XGBoost BCF feature importances |
| `bcf_xgb_predictions.png` | XGBoost BCF predicted vs actual (Held-Out R²=0.332) |
| `linear_coefficients.png` | Linear regression standardized coefficients (Held-Out R²=0.690) |
| `xgboost_feature_importance.png` | XGBoost concentration model feature importances |
| `xgboost_predictions.png` | XGBoost concentration predicted vs actual (Held-Out R²=0.710) |
| `model_comparison.png` | Side-by-side R² and RMSE across all models |
| `per_group_metrics.png` | Held-out R² and baseline comparison by species group |
| `fish_feature_importance.png` | Fish-only chemistry model feature importance |
| `fish_predictions.png` | Fish-only model predictions (R²=0.113) |
| `plant_feature_importance.png` | Plant-only chemistry model feature importance |
| `plant_predictions.png` | Plant-only model predictions (R²=0.003) |
| `chain_length_bcf_scatter.png` | Chain length vs BCF relationship by PFAS class |
| `nhanes_time_trend.png` | NHANES PFAS blood concentration trends 2015–2018 |
| `nhanes_half_life.png` | NHANES apparent population half-life vs. published literature half-life, per PFAS (v10.0) |

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
| `Species` | string | Common species name |
| `Species_Group` | string | Fish / Mammal / Plant / Human / Other |
| `Trophic_Level` | int | 1 (plant) → 5 (human) |
| `Is_Aquatic` | int | 1=aquatic, 0=terrestrial |
| `Tissue` | string | Tissue or response site |
| `Exposure Route` | string | Exposure type |
| `Duration_days` | float | Exposure duration in days (retained in dataset; removed from ML features in v7.0) |
| `Concentration_ng_g` | float | Measured tissue concentration (ng/g) |
| `log_concentration` | float | Log₁₀ concentration — ML target 1 |
| `BCF` | float | Bioconcentration factor (coalesced from BCF 1/2/3 Value) |
| `log_BCF` | float | Log₁₀ BCF — ML target 2 |
| `Source` | string | CDC NHANES / ECOTOX |

---
## PFAS Chemicals
### Tier 1 — Core
| PFAS | CASRN | Class | Chain | MW | LogKow |
|---|---|---|---|---|---|
| PFOS | 1763-23-1 | Sulfonate | 8 | 500.1 | 5.26 |
| PFOA | 335-67-1 | Carboxylate | 8 | 414.1 | 5.30 |
| PFHxS | 355-46-4 | Sulfonate | 6 | 400.1 | 4.14 |
| PFNA | 375-95-1 | Carboxylate | 9 | 464.1 | 6.05 |
| PFBS | 375-73-5 | Sulfonate | 4 | 300.1 | 1.82 |
### Tier 2 — Extended
| PFAS | CASRN | Class | Chain | MW | LogKow |
|---|---|---|---|---|---|
| PFDA | 335-76-2 | Carboxylate | 10 | 514.1 | 6.83 |
| PFUnDA | 2058-94-8 | Carboxylate | 11 | 564.1 | 7.59 |
| PFDoDA | 307-55-1 | Carboxylate | 12 | 614.1 | 8.35 |
| PFHpA | 375-85-9 | Carboxylate | 7 | 364.1 | 4.55 |
| PFHxA | 307-24-4 | Carboxylate | 6 | 314.1 | 3.77 |
### Tier 3 — Emerging
| PFAS | CASRN | Class | Chain | MW | LogKow |
|---|---|---|---|---|---|
| GenX (HFPO-DA) | 13252-13-6 | Carboxylate | 6 | 330.1 | 2.50 |
| ADONA | 958445-44-8 | Carboxylate | 8 | 380.1 | 2.80 |
| F53B | 73606-19-6 | Sulfonate | 6 | 570.1 | 4.00 |

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
### Set paths in `pfas_pipeline_v10_1.py`
```python
ECOTOX_EXPORT_DIR = "/path/to/ecotox_exports/"
COMPTOX_SNAPSHOT  = "/path/to/comptox_snapshot.csv"
NHANES_PATH       = "/path/to/nhanes_pfas_processed.csv"
OUTPUT_DIR        = "/path/to/outputs/"
```
### Run
```bash
python3 pfas_pipeline_v10_1.py
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
| CDC NHANES 2015-2016 & 2017-2018 | https://wwwn.cdc.gov/nchs/nhanes/ | Human blood serum PFAS levels |

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
      │         BCF 1/2/3 coalescing)            │
      │                   │                      │
      │                   └──────────┬───────────┘
      │                              │
      │                       Combined Dataset
      └──────────────────────────────┘
                             │
                      Merged Dataset
                      (25,056 rows)
                             │
              Stratified Split (by Species_Group)
                    │                    │
              Train (80%)            Test (20%, held out)
                    │
          Fit/Calibration Split (80/20)
              │              │
          Fit fold     Calibration fold
              │              │
              ▼              │
        RF / XGBoost         │
              │              ▼
              │      Residual Calibration
              │      (per species group)
              └──────┬───────┘
                     ▼
          Apply calibrated intervals
            to held-out Test fold
                     │
     ┌───────────────┼────────────────────┬────────────────────┐
     ▼               ▼                    ▼                    ▼
Gap Heatmap   Per-Group Metrics    RF BCF / XGB BCF    Per-PFAS Models
                     │                                  (6 own models,
                     ▼                                   1 fails baseline,
              Linear Baseline                            1 insufficient,
                     │                                    5 zero-data)
            Fish / Plant Chemistry
                  Models
```

---
## Model Results
| Version | Rows | Best R² | Key change |
|---|---|---|---|
| v1.0 | 411 | — | PFOS only |
| v2.0 | 697 | 0.279* | 5 PFAS aquatic |
| v2.1 | 1,273 | 0.110 | Terrestrial, leakage fixed |
| v3.0 | 1,524 | 0.176 | BCF model added |
| v4.0 | 13,098 | 0.742 | NHANES added |
| v5.0 | 13,098 | 0.742 | Linear baseline added |
| v6.0 | 25,056 | 0.712* | Held-out validation |
| v7.0 | 25,056 | 0.710* | Leakage fixed, XGBoost, stratified split |
| v8.0 | 25,056 | 0.705 | Calibrated prediction intervals (honest 80.6%/94.8% coverage) |
| v9.0 | 25,056 | 0.710 | Per-PFAS models — 6/13 compounds individually reliable |
| v10.0 | 25,056 | 0.710* | Apparent half-life estimation — Finding 12 |
| v10.5 | 25,056 | **0.604 (Human)** | CASRN salt fix, leakage fix (Trophic_Level/Is_Aquatic), per-group headline |

*v2.0 R²=0.279 included Route_encoded data leakage.
*v6.0 R²=0.712 and v7.0/v9.0/v10.0 R²=0.710 included Trophic_Level and Is_Aquatic — group-identity features, not measured biology (see FIX 4, v10.5). Not comparable to v10.5.
*v10.5 headline is Human R²=0.604 on chemistry-only features (Chain_Length, MW, LogKow, PFAS_Class_encoded). Pooled R²=0.490 but is ~95% human-weighted and no longer reported as the primary metric.

### Feature Set (4 features, v10.5)
Chemistry-only — all are genuine chemical properties, not group-identity re-encodings:

| Feature | Description |
|---|---|
| `Chain_Length` | Fluorinated carbon chain length |
| `MW` | Molecular weight (g/mol) |
| `LogKow` | Octanol-water partition coefficient |
| `PFAS_Class_encoded` | 0=Sulfonate, 1=Carboxylate |

Note: `Trophic_Level` and `Is_Aquatic` are retained in the output dataset as descriptive columns but removed from model inputs in v10.5 (FIX 4) — both are fixed dict lookups on `Species_Group`, not measured biological quantities.

### Held-Out Test Results (80/20 stratified split, v10.5)
| Model | Pooled R² | RMSE | Note |
|---------|---------|---------|---------|
| Random Forest | 0.490 | 0.636 | ~95% human-weighted |
| XGBoost | 0.490 | 0.636 | ~95% human-weighted |
| Linear Regression | 0.471 | 0.647 | |
| Group Mean Baseline | 0.413 | 0.682 | |
| RF (BCF) | 0.240 | 0.877 | matches baseline — data gap |
| XGBoost (BCF) | 0.240 | 0.877 | matches baseline — data gap |
| BCF Per-PFAS Baseline | 0.241 | 0.877 | |

Pooled R² is human-weighted and no longer the primary metric. See per-group results below.

### Per-Group Held-Out Performance (v10.5 — primary headline)
| Group | n | RF R² | Baseline R² | Gain |
|---------|---------|---------|---------|---------|
| Human | 4,707 | **0.604** | 0.000 | +0.604 |
| Other | 198 | -1.059 | -0.001 | -1.057 |
| Fish | 49 | -1.680 | -0.012 | -1.668 |
| Plant | 55 | -8.429 | -0.005 | -8.424 |

Human blood serum is reliably predictable from chemistry alone. Fish and Plant are not — environmental study heterogeneity dominates any chemistry signal. This is the central finding, not a model failure. Mammal excluded (only 2 test rows).

### Confidence Intervals (v10.5)
Residual-calibrated per species group on a disjoint calibration fold. All groups now hit coverage targets after the CASRN salt fix recovered additional fish records.

| Group | n (calibration) | 80% width | 95% width | 80% coverage | 95% coverage |
|---|---|---|---|---|---|
| Human | 3,765 | 0.900 | 1.530 | 79.3% ✓ | 94.6% ✓ |
| Other | 158 | 3.753 | 5.723 | 80.3% ✓ | 96.5% ✓ |
| Plant | 44 | 1.293 | 1.859 | 78.2% ✓ | 96.4% ✓ |
| Fish | 39 | 2.891 | 5.088 | 79.6% ✓ | 93.9% ✓ |
| **Overall** | 4,008 | **1.036** | **1.735** | **79.3% ✓** | **94.7% ✓** |

### Per-PFAS Predictability (v10.5)
| PFAS | n (total) | Model Type | Held-Out R² |
|---|---|---|---|
| PFOA | 4,594 | own model | **0.657** |
| PFNA | 3,978 | own model | 0.401 |
| PFOS | 4,496 | own model | 0.436 |
| PFHxS | 4,004 | own model | 0.381 |
| PFDA | 3,949 | own model | 0.226 |
| PFUnDA | 3,946 | own model | 0.147 |
| PFBS | 83 | own model | -0.063 (worse than baseline) |
| PFHpA | 1 | insufficient data | — |
| GenX | 0 | no data | — |
| ADONA | 0 | no data | — |
| F53B | 0 | no data | — |
| PFDoDA | 0 | no data | — |
| PFHxA | 0 | no data | — |

7 of 13 PFAS now have enough data for a dedicated model (up from 6 in v9.0 — PFBS crossed the threshold via recovered salt rows). PFHxS improved most dramatically: R² 0.172 → 0.381 from the 36 recovered PFHxS potassium salt records. PFBS remains a cautionary case — more rows didn't improve predictability.

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

**Resolved (v10.5):** 401 previously unmatched ECOTOX rows recovered via `CASRN_SALT_MAP`, which remaps salt/counterion CASRNs to their parent compound. 5 rows (2 CASRNs) remain unmatched — confirmed outside the curated 13-PFAS scope.

---
## Roadmap
### Phase 1 — Data Expansion ✅ Complete
- [x] PFAS class as ML feature
- [x] BCF as second target variable
- [x] Terrestrial species data
- [x] Tier 2 PFAS
- [x] CDC NHANES human blood serum data (2015–2016 and 2017–2018)

### Phase 2 — Better Models ✅ Complete
- [x] Linear Regression baseline
- [x] Random Forest concentration model
- [x] Random Forest BCF model
- [x] XGBoost concentration model
- [x] XGBoost BCF model
- [x] Model comparison summary and chart
- [x] Per-species-group evaluation
- [x] Fish-only chemistry model
- [x] Plant-only chemistry model
- [x] Held-out validation framework (stratified 80/20 split)
- [x] NHANES trend analysis
- [x] Leakage audit and fix (is_human, Duration_days, is_* flags)
- [x] BCF feature set corrected (chemistry-only, proper per-PFAS baseline)
- [x] BCF data recovery (BCF 1/2/3 Value coalescing)
- [x] Prediction confidence intervals (residual-calibrated, per species group)
- [x] Per-PFAS models (6 own models, 1 fails baseline, 1 insufficient, 5 zero-data)
- [x] Apparent NHANES half-life estimation vs. published literature (Finding 12)
- [x] CASRN salt mapping — recovered 396/401 previously dropped ECOTOX rows (v10.5)
- [x] Leakage fix: Trophic_Level and Is_Aquatic removed from model features (v10.5)
- [x] Per-group R² as headline metric replacing misleading pooled figure (v10.5)

### Phase 3 — Mechanistic Modeling & Outputs (July–August 2026)
**Mechanistic / chemical-engineering additions** — the current pipeline is purely statistical
(ML on tabular features); Phase 3 adds physically-grounded models that explain *why* the
statistical findings look the way they do, rather than only reporting that they do:
- [ ] **Arnot-Gobas / fugacity bioaccumulation model for BCF** — Finding 5 shows BCF can't be
      predicted from chemistry-only ML (RF/XGBoost match a per-compound-mean baseline exactly,
      R²≈0.331). A mass-balance mechanistic model (gill/dietary uptake vs. respiratory/fecal
      egestion and growth dilution, parameterized by Kow) is the standard chemE alternative —
      but needs a PFAS-adapted formulation, since PFAS behave as protein-binding surfactants
      rather than the neutral-organic lipid-partitioning chemicals the classic model assumes.
      Goal: report mechanistic-model BCF alongside ML BCF as a direct comparison, and use the
      gap between them to explain *why* chemistry-only ML underperforms.
- [ ] **PFAS-appropriate physical chemistry features** — replace/augment LogKow (a weak
      descriptor for PFAS, which don't partition like classic neutral organics) with features
      that reflect actual PFAS environmental behavior: soil-water partition coefficient (Koc)
      for the terrestrial/plant pathway, air-water partition coefficient (Henry's constant) for
      volatilization, a serum-albumin protein-binding proxy, and critical micelle concentration
      (PFAS are surfactants). Test whether this measurably improves the currently very weak
      chemistry-only fish (R²=0.113) and plant (R²=0.003) models — a real, testable hypothesis
      that should be attempted before investing further in the Arnot-Gobas model above, since a
      better mechanistic BCF model is only as good as the chemistry features feeding it.
- [ ] **Engineered treatment removal-efficiency module** — literature-correlation estimates of
      PFAS removal by GAC adsorption, anion exchange resin, and reverse osmosis as a function of
      chain length and PFAS class. Connects the pipeline's existing chain-length/class features
      to a practical "what do we do about it" outcome, and explains *why* short-chain
      replacements like GenX are an emerging problem (harder to remove by the treatment
      technologies designed around longer-chain legacy PFAS). Lower research novelty than the
      two items above, but highest value for a general/policy audience.
- [ ] Longitudinal (individual-level, not cross-sectional) half-life validation, if a suitable
      public cohort dataset can be identified — would resolve the exposure-vs-elimination
      conflation flagged in Finding 12 rather than just quantifying it

**Deliverables:**
- [ ] Interactive bioaccumulation simulator (feeds from per-PFAS models + mechanistic BCF model)
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
python3 pfas_pipeline_v10_1.py
```
### New PFAS chemicals
Add to `PFAS_FEATURES` in `pfas_pipeline_v10_1.py`:
```python
("PFDA", "335-76-2", "Carboxylate", 10, 514.1, 6.83),
# (Name, CASRN, Class, Chain_Length, MW, LogKow)
```

---
## Citation
- EPA ECOTOX Knowledgebase: https://cfpub.epa.gov/ecotox/
- EPA CompTox Dashboard: https://comptox.epa.gov/dashboard/
- CDC NHANES 2015–2016 and 2017–2018: https://wwwn.cdc.gov/nchs/nhanes/

---
## Author
PFAS Environmental Informatics Research Project (2026)
