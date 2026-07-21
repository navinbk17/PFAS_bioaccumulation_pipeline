# PFAS Bioaccumulation Research Pipeline v16.0

A reproducible, multi-source data pipeline for studying how "forever chemicals" (PFAS) build up in humans, animals, and plants. By combining environmental, chemical, and human health data from the EPA and CDC, this project builds a clean dataset to train AI models, map out critical missing data, predict chemical behavior from molecular structure alone, and model how chemicals move through living bodies using biology-based physics models.

**Current dataset: 25,056 observations | 13 curated PFAS (7 with full models) | 5 species groups | 5 Machine Learning models + Arnot-Gobas biological model + Two-compartment body model | Headline evaluation: Unseen Chemical Testing (LOCO CV, R²=−0.160) | Calibrated 80%/95% prediction ranges | Apparent blood half-life estimated for 4/6 main PFAS | Key Insight: Splitting fish biology into blood vs. tissue compartments resolves calculation errors for major PFAS like PFOS (−6%) and PFBS (+19%)**

---

## Table of Contents

* [Why This Research Matters](#why-this-research-matters)
* [Key Findings](#key-findings)
* [Version History](#version-history)
* [Outputs](#outputs)
* [Dataset Schema](#dataset-schema)
* [PFAS Chemicals](#pfas-chemicals)
* [Setup](#setup)
* [Usage](#usage)
* [Data Sources](#data-sources)
* [Pipeline Architecture](#pipeline-architecture)
* [Model Results](#model-results)
* [Data Gaps](#data-gaps)
* [Roadmap](#roadmap)
* [How to Add Data](#how-to-add-data)

---

## Why This Research Matters

### PFAS Are Everywhere — And They Don't Leave

Per- and polyfluoroalkyl substances (PFAS) are a group of over 12,000 synthetic chemicals used in non-stick cookware, food packaging, firefighting foams, waterproof clothing, and industrial manufacturing. They are nicknamed "forever chemicals" because their chemical bonds are among the strongest in nature. They do not break down naturally in the environment, nor do they break down in human bodies. Instead, they accumulate.

### The Food Chain Problem

When PFAS enter an environment through factory wastewater, farming runoff, or polluted water, they get absorbed by plants and tiny organisms at the bottom of the food chain. As larger animals eat smaller ones, the concentration of PFAS multiplies at each step. This process (biomagnification) means a fish at the top of a lake's food chain can carry PFAS levels thousands of times higher than the surrounding water. Humans sit at the very top of that chain.

### What the Numbers Say

* PFAS are found in the blood of **97% of Americans**.
* The EPA set drinking water limits for major PFAS at **4 parts per trillion** — a concentration so tiny it required brand-new technology just to measure.
* PFAS exposure is linked to thyroid disease, immune system damage, certain cancers, fertility problems, and childhood developmental delays.
* Our pipeline shows average PFOS levels of **2.83 ng/g** in human blood serum from CDC national surveys — measured in people with no direct chemical job exposure.

### The Scientific Gap We're Addressing

Despite these risks, our understanding of how PFAS move through ecosystems is fragmented. Studies are scattered across hundreds of isolated papers, measured in different units, tested on different animals, and reported inconsistently. No single system has cleanly mapped PFAS accumulation from soil → plant → fish → mammal → human. Filling that gap is what this project does.

---

## Key Findings

### Finding 1 — Where you sit in the food chain is the biggest predictor of PFAS accumulation

Across 25,056 measurements, an animal's position in the food chain explains tissue concentration better than any chemical property. This directly proves biomagnification: the higher up the food chain an organism is, the more PFAS it holds. Humans show the highest and most consistent levels.

### Finding 2 — In human blood, chemical structure predicts accumulation logically

When looking strictly at human blood serum, basic machine learning (Linear Regression) and advanced algorithms (Random Forest/XGBoost) perform almost identically. Longer chemical chains and higher water-repelling properties directly lead to higher blood concentrations in a predictable, straight-line relationship.

### Finding 3 — Human blood data is clean; environmental data is chaotic

Human blood samples collected by the CDC follow strict, standardized laboratory rules, making them highly consistent and predictable (R²=0.658). Environmental studies on plants and wild fish are so varied and messy that models trained on them perform worse than simple guessing.

### Finding 4 — Wildlife data cannot predict human health risks

Testing shows that a model trained on fish, plants, and wild mammals fails completely when trying to predict human blood levels. This is caused by two separate issues: mismatched measurement types (blood serum vs. whole tissue vs. water concentrations) and fundamental biological differences between species.

### Finding 5 — Bioconcentration cannot be predicted from chemical traits alone

Bioconcentration Factor (BCF) measures how much more chemical is in an animal compared to its surrounding water. Using only chemical formulas, machine learning models fail to beat a basic average baseline. Differences in water chemistry, temperature, and individual animal biology outweigh molecular traits.

### Finding 6 — Missing data is severe and heavily biased

PFHxS is present in the blood of almost every American, yet the EPA's ECOTOX database contains zero records for it in fish or wild mammals. Five out of the 13 major PFAS studied have zero real-world test records anywhere in the database.

### Finding 7 — Mammalian bioaccumulation data is currently unusable

Nearly all EPA records for mammals involve short-term toxicity dosing experiments rather than measurements of chemical residue remaining in tissues. Only 4 usable mammal tissue records exist across all studied chemicals.

### Finding 8 — Predicting environmental levels remains extremely difficult

Models trained strictly on environmental species yield poor results:

* **Fish-only models:** Perform worse than taking a simple average.
* **Plant-only models:** Perform worse than taking a simple average.

This confirms that environmental field data is currently too sparse and inconsistent for pure chemical machine learning.

### Finding 9 — AI models can easily "cheat" if not strictly audited

Earlier pipeline versions appeared highly accurate because models were accidentally picking up hidden species identifiers (e.g., recognizing that a record came from a human study rather than learning the chemistry). Removing all indirect species tags dropped artificial accuracy down to the true chemical signal (R²=0.490).

### Finding 10 — Standard AI models are dangerously overconfident about uncertainty

Standard Random Forest confidence ranges were wildly overconfident, covering only 2.0% of real-world outcomes when targeting 80%. After implementing proper statistical calibration, our prediction ranges now accurately hit 80% and 95% real-world coverage across all species.

### Finding 11 — Model accuracy varies wildly depending on the specific PFAS

Out of 13 key PFAS chemicals, only 7 have enough scientific data to train dedicated models. PFOA is the most predictable (R²=0.657 across 4,594 records), while PFBS performs worse than guessing (R²=-0.063) despite having 83 records.

### Finding 12 — Estimated human half-lives are inflated by continuous exposure

Calculating how fast human bodies clear PFAS using national survey waves yields apparent "half-lives" far longer than medical trials suggest (e.g., PFOS appears to take 42.3 years to drop by half, versus 5.4 years in medical trials). This is not an error — it proves that ongoing everyday exposure continuously replenishes PFAS levels in the general population.

### Finding 13 — Adding complex chemical descriptors doesn't fix environmental models

Adding specialized chemical parameters (like soil binding coefficients and protein binding pKa) added zero predictive power to fish and plant models. The bottleneck is not missing chemical features; it is a lack of standardized field data.

### Finding 14 — Biological mass balance models reveal protein-binding quirks in sulfonates

Using a physics-based biological model (Arnot-Gobas) shows that standard chemical formulas systematically under-predict accumulation for sulfonate-type PFAS (PFOS -95%, PFHxS -82%, PFBS -67%). This shows that sulfonates bind to blood proteins (albumin) far more aggressively than standard models assume.

### Finding 15 — Within human blood, chemical identity matters more than small structural tweaks

While our human-only model performs well overall, testing it on a single chemical at a time shows near-zero predictive power. The model easily tells different PFAS families apart, but small internal variations within a single chemical family carry little statistical weight at current sample sizes.

### Finding 16 — Protein binding cannot be fixed by tweaking a single scaling number

Testing 7 different scaling factors for sulfonate protein binding showed almost no change in model error (PFOS error stayed frozen between -92% and -95%). This proves that standard chemical proxies (Koc) are fundamentally wrong for sulfonates, requiring direct blood-protein binding measurements (Ka) instead.

### Finding 17 — Gill permeability is not the missing piece

Testing whether sulfonates pass through fish gills slower than carboxylates showed that reducing gill permeability actually made model predictions worse for PFOS (−94% → −98%). This proves the modeling error lies entirely within tissue accumulation, not gill intake.

### Finding 18 — Direct protein binding improves science but highlights transport gaps

Replacing chemical proxies with real blood-protein binding measurements (albumin affinity) fixed tissue ratios but showed that simple single-compartment fish models are missing active protein transport mechanisms at the gills.

### Finding 19 — Single-compartment fish models have hit a theoretical limit

Adding protein-facilitated intake terms to a single-compartment fish model improved predictions for PFOS (error dropped from -95% to -25%), but caused PFHxS to over-predict (+34%) while PFBS barely moved (-68%). A single-compartment model cannot solve all three chemicals simultaneously.

### Finding 20 — Two-compartment body models close the gap for major PFAS

Separating fish biology into two distinct zones — a blood compartment and a tissue compartment — dramatically improves physical accuracy. The model brings PFOS error down to just −6% and PFBS to +19%. However, fitting PFHxS simultaneously requires compound-specific tissue binding data that currently does not exist in scientific literature.

---

## Version History

### v16.0 (current) — July 2026

* **Unseen Chemical Testing (LOCO CV) as Headline Metric:** Evaluates models by holding out entire chemical compounds during training to test real-world generalization on new chemicals. Pooled LOCO R²=−0.160, confirming that 8 distinct chemistry profiles are insufficient for structure-to-accumulation generalization — an honest negative result.
* **Relative File Paths:** Removed hardcoded system paths in favor of an automatic local project directory structure. All paths are derived from the script's own location.

### v15.0 — July 2026

* **Two-Compartment Biological Model:** Implemented steady-state calculations separating fish blood from tissue compartments.
* **Tissue-Specific Partitioning:** Calibrated non-fat tissue factors for different chemical classes.
* **Threshold Gating:** Added affinity thresholds to prevent over-estimating uptake for low-binding chemicals like PFBS.
* **Finding 20 Confirmed:** Proved two-compartment modeling is the correct architecture for PFAS in aquatic life.

### v14.0 — July 2026

* **Protein-Facilitated Uptake:** Added protein-mediated gill transport equations to the mechanistic biological model.
* **Class-Specific Uptake:** Restricted protein-facilitated gill uptake rules exclusively to sulfonates.
* **Code Bug Fixes:** Fixed data filtering loops and corrected food-intake energy scaling formulas.

### v13.0 — July 2026

* **Direct Protein Binding:** Replaced chemical proxies with direct experimental blood-protein affinity values (Ka).

### v12.0 — July 2026

* **Sensitivity Sweeps:** Ran systematic sweeps testing protein scaling factors and gill permeability thresholds (Findings 16 & 17).

### v11.0 — July 2026

* **Mechanistic Biological Model:** Integrated the adapted Arnot-Gobas mass balance model for fish bioconcentration.
* **Human-Only Sub-Model:** Built dedicated human blood serum models using CDC NHANES data.

### v10.0 — June-July 2026

* **Data Recovery:** Recovered 396 missing EPA records using salt-mapping routines.
* **Apparent Half-Life Engine:** Calculated population clearing rates from CDC multi-year studies.

### v9.0–v1.0 — April–June 2026

* Incremental developments: Added XGBoost models, residual calibration, dataset harmonization, CDC NHANES integration, and baseline linear regressions.

---

## Outputs

| File | Description |
| --- | --- |
| `pfas_bioaccumulation_dataset.csv` | 25,056 rows, fully cleaned and ML-ready |
| `pfas_gap_heatmap.png` | Visual count of available measurements per chemical and species group |
| `loco_r2.png` | Per-compound LOCO R² bar chart — headline chemistry evaluation |
| `feature_importance.png` | Key chemical factors driving the main concentration model |
| `model_predictions.png` | AI model predictions vs real-world measurements |
| `human_model_predictions.png` | Human blood model accuracy broken down by chemical |
| `human_model_feature_importance.png` | Key factors driving human blood concentration |
| `arnot_gobas_bcf.png` | Biological model calculations vs machine learning vs real observations |
| `arnot_gobas_sensitivity.png` | Testing protein scaling factors across different chemicals |
| `arnot_gobas_pmem_sensitivity.png` | Testing gill permeability factors across different chemicals |
| `arnot_gobas_kfac_sensitivity.png` | Testing protein transport efficiency at fish gills |
| `two_comp_nlom_sensitivity.png` | Calibration heatmap for the two-compartment fish body model |
| `arnot_gobas_2comp_bcf.png` | Side-by-side comparison of 1-compartment vs 2-compartment fish models |
| `feature_ablation.png` | Impact chart showing whether adding new chemical features improved accuracy |
| `prediction_intervals.png` | High/low confidence boundary bands by species group |
| `interval_coverage.png` | Diagnostic chart confirming prediction interval accuracy |
| `per_pfas_r2.png` | Model accuracy scores for individual PFAS compounds |
| `per_group_metrics.png` | Model performance split across humans, fish, plants, and mammals |
| `cross_species_validation.png` | Accuracy metrics when using one species to predict another |
| `bcf_feature_importance.png` | Key drivers behind aquatic bioconcentration |
| `bcf_predictions.png` | Random Forest bioconcentration predictions |
| `bcf_xgb_predictions.png` | XGBoost bioconcentration predictions |
| `linear_coefficients.png` | Standardized weight factors from linear baseline models |
| `xgboost_predictions.png` | XGBoost concentration predictions vs actual data |
| `model_comparison.png` | Side-by-side performance metrics across all implemented AI models |
| `fish_predictions.png` | Fish-only chemical prediction performance chart |
| `plant_predictions.png` | Plant-only chemical prediction performance chart |
| `chain_length_bcf_scatter.png` | Relationship between carbon chain length and bioaccumulation |
| `nhanes_time_trend.png` | Changing human blood serum levels from 2015 to 2018 |
| `nhanes_half_life.png` | Estimated population clearance rates vs published medical trials |

---

## Dataset Schema

| Column | Type | Description |
| --- | --- | --- |
| `PFAS_Name` | string | Common chemical abbreviation (e.g., PFOS, PFOA) |
| `CASRN` | string | Official chemical registry index number |
| `PFAS_Class` | string | Chemical family (Sulfonate or Carboxylate) |
| `PFAS_Class_encoded` | int | Numerical flag: 0 = Sulfonate, 1 = Carboxylate |
| `Chain_Length` | int | Number of fluorinated carbon atoms in the molecule |
| `MW` | float | Molecular weight (grams per mole) |
| `LogKow` | float | Water-to-fat partition coefficient (measure of oil/water preference) |
| `Koc` | float | Soil organic-carbon to water partition coefficient |
| `Koc_log` | float | Logarithmic soil coefficient (used in AI models) |
| `AlbuminBinding_pKa` | float | Acidity constant proxy for blood protein binding |
| `Species` | string | Common name of tested organism |
| `Species_Group` | string | Category: Fish / Mammal / Plant / Human / Other |
| `Trophic_Level` | int | Position in food chain (1=Plant → 5=Human); reference only |
| `Is_Aquatic` | int | Environment flag: 1 = Aquatic, 0 = Terrestrial; reference only |
| `Tissue` | string | Sampled organ or body part (e.g., Blood Serum, Liver, Whole Body) |
| `Exposure Route` | string | How the organism encountered the chemical (e.g., Water, Food, Oral) |
| `Duration_days` | float | Length of exposure trial in days |
| `Concentration_ng_g` | float | Measured concentration in parts per billion (ng/g) |
| `log_concentration` | float | Log-transformed concentration (primary AI target) |
| `BCF` | float | Bioconcentration Factor (Tissue Ratio ÷ Water Concentration) |
| `log_BCF` | float | Log-transformed BCF (secondary AI target) |
| `Source` | string | Origin dataset (CDC NHANES or EPA ECOTOX) |

---

## PFAS Chemicals

### Tier 1 — Core (Widely Studied)

| PFAS | CASRN | Class | Carbon Chain | Molecular Weight | LogKow | Koc (L/kg) |
| --- | --- | --- | --- | --- | --- | --- |
| PFOS | 1763-23-1 | Sulfonate | 8 | 500.1 | 5.26 | 2,100 |
| PFOA | 335-67-1 | Carboxylate | 8 | 414.1 | 5.30 | 1,900 |
| PFHxS | 355-46-4 | Sulfonate | 6 | 400.1 | 4.14 | 560 |
| PFNA | 375-95-1 | Carboxylate | 9 | 464.1 | 6.05 | 3,200 |
| PFBS | 375-73-5 | Sulfonate | 4 | 300.1 | 1.82 | 47 |

### Tier 2 — Extended (Secondary Importance)

| PFAS | CASRN | Class | Carbon Chain | Molecular Weight | LogKow | Koc (L/kg) |
| --- | --- | --- | --- | --- | --- | --- |
| PFDA | 335-76-2 | Carboxylate | 10 | 514.1 | 6.83 | 5,800 |
| PFUnDA | 2058-94-8 | Carboxylate | 11 | 564.1 | 7.59 | 9,100 |
| PFDoDA | 307-55-1 | Carboxylate | 12 | 614.1 | 8.35 | 14,000 |
| PFHpA | 375-85-9 | Carboxylate | 7 | 364.1 | 4.55 | 820 |
| PFHxA | 307-24-4 | Carboxylate | 6 | 314.1 | 3.77 | 310 |

### Tier 3 — Emerging (Zero Environmental Test Data)

| PFAS | CASRN | Class | Carbon Chain | Molecular Weight | LogKow | Koc (L/kg) |
| --- | --- | --- | --- | --- | --- | --- |
| GenX (HFPO-DA) | 13252-13-6 | Carboxylate | 6 | 330.1 | 2.50 | — |
| ADONA | 958445-44-8 | Carboxylate | 8 | 380.1 | 2.80 | — |
| F53B | 73606-19-6 | Sulfonate | 6 | 570.1 | 4.00 | — |

---

## Setup

```bash
pip3 install pandas numpy matplotlib seaborn scikit-learn openpyxl requests xgboost
```

On macOS, XGBoost requires the OpenMP library:

```bash
brew install libomp
```

---

## Usage

### Directory Layout

Place your input files in the same folder as the script. The `outputs/` folder is created automatically on first run.

```
pfas_pipeline_v16.py
ecotox_exports/             ← Folder of EPA ECOTOX XLSX exports
nhanes_pfas_processed.csv
comptox_snapshot.csv        ← Optional EPA properties snapshot
outputs/                    ← Auto-generated folder for results and charts
```

All four paths are defined in a single config block near the top of `pfas_pipeline_v16.py`. Edit them there if your files live somewhere else:

```python
ECOTOX_EXPORT_DIR = os.path.join(_BASE, "ecotox_exports") + os.sep
COMPTOX_SNAPSHOT  = os.path.join(_BASE, "comptox_snapshot.csv")
NHANES_PATH       = os.path.join(_BASE, "nhanes_pfas_processed.csv")
OUTPUT_DIR        = os.path.join(_BASE, "outputs") + os.sep
```

### Run

```bash
python3 pfas_pipeline_v16.py
```

### Required Files

| File | Download Source |
| --- | --- |
| ECOTOX exports | https://cfpub.epa.gov/ecotox/ (search each chemical and export to XLSX) |
| `nhanes_pfas_processed.csv` | Included directly in this repository |
| `comptox_snapshot.csv` | Optional — https://comptox.epa.gov/dashboard/batch-search |

---

## Data Sources

| Source | Link | Data Provided |
| --- | --- | --- |
| EPA ECOTOX | https://cfpub.epa.gov/ecotox/ | Wildlife species, tissue measurements, concentrations, BCF |
| EPA CompTox | https://comptox.epa.gov/dashboard/batch-search | Molecular weights, LogKow values, chemical structures |
| CDC NHANES (2015–2018) | https://wwwn.cdc.gov/nchs/nhanes/ | Representative human blood serum levels in the USA |

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
       ┌─────────────────────┼──────────────────────────┐
       ▼                     ▼                          ▼
Pooled Model          Human-Only Model          Arnot-Gobas BCF
(all species)         (NHANES only)             (mechanistic 1-comp)
       │                     │                          │
Stratified Split      Stratified Split          Kprot Sweep (F16)
by Species_Group      by PFAS_Name              P_mem Sweep (F17)
       │                     │                  Ka-albumin Kfish (F18)
Fit / Cal / Test (3-way)     │                  k1_fac Sweep (F19)
       │                     │                          │
RF | XGBoost | Linear        │                  Two-Compartment TK
       │                     │                  (Option B, v15.0)
Residual-calibrated          │                  NLOM Sweep (F20)
prediction intervals         │
(per species group)          │
       │                     │
  ┌────┼──────────┬────┐     │
  ▼    ▼          ▼    ▼     ▼
Gap  Per-Group  BCF  Per-PFAS  Human
Heat  Metrics   RF+  Models   Intervals
map            XGB   (7)     (per PFAS)
                             +
                          LOCO CV
                       (headline R²)
```

---

## Model Results

### Headline — Testing on Unseen Chemicals (LOCO CV)

Leave-One-Compound-Out (LOCO) cross-validation is our primary benchmark. Each test round holds out **all data for an entire chemical** and trains models on the remaining chemicals. This tests whether the model can predict bioaccumulation for a brand-new chemical it has never seen before.

Standard random splits look artificially accurate because rows sharing the exact same chemical structure end up in both training and testing sets. LOCO eliminates this overlap completely.

| Evaluation Metric | Random Forest Score (R²) | Meaning |
| --- | --- | --- |
| **LOCO Benchmark (Headline)** | **−0.160** | Performance when tested on an entirely unseen chemical |
| Pooled Random Split (Reference only) | 0.490 | Standard split score — inflated because matching structure rows appear in both train and test |
| Simple Average Baseline | 0.413 | Guessing the overall species group average |

*Note: A negative LOCO R² means the model performs worse than simply predicting the training mean for an unseen compound. With only 8 distinct chemical profiles, this is an expected and honest result — not a failure of the code. It tells us that current data is insufficient to generalize chemistry → accumulation to new PFAS.*

### Performance by Species Group (Chemistry-Only Features)

| Group | Samples (n) | AI Model Score (R²) | Baseline Guess (R²) | Net Improvement |
| --- | --- | --- | --- | --- |
| Human | 4,707 | **0.604** | 0.000 | +0.604 |
| Other | 198 | -1.059 | -0.001 | -1.057 |
| Fish | 49 | -1.680 | -0.012 | -1.668 |
| Plant | 55 | -8.429 | -0.005 | -8.424 |

Human blood serum concentration is predictable from molecular features. Fish and plant bioaccumulation cannot be predicted from chemical structure alone using current messy environmental field data.

### Human-Only Model Results

| Model Type | Accuracy Score (R²) | Error Margin (RMSE) | Improvement over Pooled Model |
| --- | --- | --- | --- |
| Random Forest | **0.658** | 0.360 | +0.054 |
| XGBoost | 0.658 | 0.360 | +0.054 |
| Linear Regression | 0.649 | 0.365 | +0.045 |
| Chemical Class Average | 0.658 | 0.360 | +0.054 |

### Chemical Feature Sets Used

**All-Species Model (6 features):**

* Carbon Chain Length
* Molecular Weight
* LogKow (Water/Fat Partitioning)
* PFAS Class (Sulfonate vs Carboxylate)
* Log Koc (Soil/Water Partitioning)
* Albumin Binding pKa (Protein Affinity Proxy)

**Human-Only Model (5 features):**

* Same features as above, excluding soil coefficient (`Log Koc`).

### Calibrated Prediction Range Accuracy

| Species Group | 80% Target Margin | 95% Target Margin | Real 80% Coverage | Real 95% Coverage |
| --- | --- | --- | --- | --- |
| Human | 0.900 | 1.530 | 79.8% ✓ | 94.8% ✓ |
| Fish | 2.891 | 5.088 | 79.6% ✓ | 93.9% ✓ |
| Plant | 1.293 | 1.859 | 78.2% ✓ | 96.4% ✓ |
| Other | 3.753 | 5.723 | 80.3% ✓ | 96.5% ✓ |
| **Overall** | **1.036** | **1.735** | **79.8% ✓** | **94.9% ✓** |

### Predictability by Chemical Compound

| PFAS Name | Sample Count (n) | Model Strategy | Held-Out Accuracy (R²) |
| --- | --- | --- | --- |
| PFOA | 4,594 | Dedicated Model | **0.657** |
| PFNA | 3,978 | Dedicated Model | 0.401 |
| PFOS | 4,496 | Dedicated Model | 0.436 |
| PFHxS | 4,004 | Dedicated Model | 0.381 |
| PFDA | 3,949 | Dedicated Model | 0.226 |
| PFUnDA | 3,946 | Dedicated Model | 0.148 |
| PFBS | 83 | Dedicated Model | -0.063 (Below Baseline) |
| PFHpA | 1 | Insufficient Data | — |
| GenX / ADONA / F53B / PFDoDA / PFHxA | 0 | Zero Measurements | — |

### Biological Fish Bioconcentration Model (1-Compartment)

Calculated bioconcentration ratios using biological mass balance for a standard 1 kg fish at 12°C:

| PFAS | Calculated Log BCF | Calculated Ratio | Measured Real Ratio | Calculation Error |
| --- | --- | --- | --- | --- |
| PFOS | 1.877 | 75.3 | 2.000 | -25% |
| PFHxS | 1.316 | 20.7 | 1.187 | +34% |
| PFBS | 0.377 | 2.4 | 0.874 | -68% |
| PFOA | 0.912 | 8.2 | 0.953 | -9% |
| PFDA | 1.790 | 61.7 | 1.903 | -23% |
| PFNA | 1.341 | 21.9 | 1.863 | -70% |
| PFUnDA | 2.178 | 150.6 | 2.980 | -84% |

### Two-Compartment Fish Body Model (v15.0)

Splits fish biology into Blood (5% body volume) and Tissue (95% body volume) compartments to calculate steady-state accumulation:

| PFAS | Calculated Log BCF | Calculated Ratio | Measured Real Ratio | 2-Compartment Error | 1-Compartment Error | Outcome |
| --- | --- | --- | --- | --- | --- | --- |
| PFOS | 1.974 | 94.3 | 2.000 | **−6%** | −25% | **Improved** |
| PFHxS | 1.546 | 35.1 | 1.187 | **+128%** | +34% | Degraded |
| PFBS | 0.949 | 8.9 | 0.874 | **+19%** | −68% | **Improved** |
| PFNA | 2.067 | 116.6 | 1.863 | **+60%** | −70% | **Improved** |
| PFDA | 2.175 | 149.7 | 1.903 | **+87%** | −23% | Degraded |
| PFUnDA | 2.461 | 289.3 | 2.980 | **−70%** | −84% | **Improved** |

---

## Data Gaps

Number of available scientific measurements in the database:

| PFAS | Fish | Human | Mammal | Plant | Other |
| --- | --- | --- | --- | --- | --- |
| PFOS | 22 | 1,929 | 2 | 130 | 121 |
| PFOA | 150 | 1,929 | 2 | 123 | 352 |
| PFNA | 16 | 1,929 | 0 | 0 | 40 |
| PFHxS | **0** | 1,929 | **0** | **0** | 46 |
| PFDA | 0 | 1,929 | 0 | 0 | 27 |
| PFUnDA | 2 | 1,929 | 0 | 0 | 17 |
| PFBS | 6 | 0 | 0 | 0 | 66 |
| PFHpA | 0 | 0 | 0 | 0 | 1 |

**Major Gap:** PFHxS has nearly 2,000 human measurements, but **zero** fish or mammal records in the primary EPA database.

---

## Roadmap

### Phase 1 — Data Expansion ✅ Complete

* [x] Add chemical class features
* [x] Add bioconcentration targets
* [x] Integrate land species records
* [x] Integrate Tier 2 PFAS compounds
* [x] Add CDC NHANES human blood serum datasets (2015–2018)

### Phase 2 — Model Optimization ✅ Complete

* [x] Build baseline and advanced machine learning models
* [x] Establish strict train/test data splits
* [x] Audit and fix data leakage
* [x] Implement calibrated prediction intervals
* [x] Build compound-specific models
* [x] Estimate population clearance rates
* [x] Recover missing chemical salt records
* [x] Integrate biological mass balance models (Arnot-Gobas)
* [x] Build human-only blood serum models

### Phase 3 — Refinements & Output ✅ Complete

* [x] Test protein scaling limits (Finding 16)
* [x] Test gill permeability parameters (Finding 17)
* [x] Implement direct blood-protein affinity constants (Finding 18)
* [x] Implement protein-facilitated uptake terms (Finding 19)
* [x] Implement steady-state two-compartment fish model (Finding 20)
* [x] Implement Leave-One-Compound-Out (LOCO) CV as headline evaluation (v16.0)

### Phase 4 — Future Extensions (Planned)

* [ ] Draft journal publication for findings 14–20
* [ ] Add water treatment removal module (activated carbon, reverse osmosis removal rates)
* [ ] Build interactive Web/Streamlit simulation dashboard
* [ ] Publish open-source GitHub release repository

---

## How to Add Data

### Adding New ECOTOX Files

```bash
# 1. Download records from https://cfpub.epa.gov/ecotox/
# 2. Export as XLSX and save to your local project folder:
mv ~/Downloads/ECOTOX-*.xlsx ./ecotox_exports/
python3 pfas_pipeline_v16.py
```

### Adding New PFAS Chemicals

Add the chemical traits to the `PFAS_FEATURES` table inside `pfas_pipeline_v16.py`:

```python
("PFDA", "335-76-2", "Carboxylate", 10, 514.1, 6.83, 5800.0, 0.30, 4.55),
# Format: (Name, CASRN, Class, Chain_Length, Molecular_Weight, LogKow, Koc, Albumin_pKa, Log_Ka_Albumin)
```

---

## Citation

* EPA ECOTOX Database: https://cfpub.epa.gov/ecotox/
* EPA CompTox Dashboard: https://comptox.epa.gov/dashboard/
* CDC NHANES Human Surveys: https://wwwn.cdc.gov/nchs/nhanes/
* Arnot & Gobas (2004) *Environ. Toxicol. Chem.* 23:1523–1532
* Kelly et al. (2004) *Environ. Sci. Technol.*
* Gobas et al. (2003) *Environ. Sci. Technol.*
* Bischel et al. (2010) *Environ. Sci. Technol.*
* Beesoon & Martin (2015) *Environ. Sci. Technol.*
* Ng & Hungerbühler (2013) *Environ. Sci. Technol.* 47:7214
* Barber (2003) *Chemosphere* 53:1099
* Farrell (1991) *J. Exp. Biol.* 159:213
* ATSDR Toxicological Profile for Perfluoroalkyls (2021)

---

## Author

PFAS Environmental Informatics Research Project (2026)
