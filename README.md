# PFAS Bioaccumulation Research Pipeline v4.0

> *"PFAS chemicals don't break down. They move up food chains, accumulate in tissues, and end up in human blood. This pipeline maps where they go — and where our scientific understanding runs out."*

A reproducible, multi-source data pipeline for studying PFAS bioaccumulation across aquatic and terrestrial species and human populations. Integrates EPA ECOTOX biological exposure data, EPA CompTox chemical properties, and CDC NHANES human biomonitoring data to build a machine-learning-ready dataset, identify critical data gaps, and predict bioaccumulation from chemical structure alone.

**Current dataset: 13,098 observations | 8 PFAS | 5 species groups | 2 ML models | R²=0.742**

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

Despite this, our understanding of *how* PFAS move through ecosystems remains deeply fragmented. Data is scattered across hundreds of studies, measured in inconsistent units, tested on different species, and reported under different conditions. No single database cleanly maps PFAS bioaccumulation from soil → plant → fish → mammal → human.

**That gap is what this project addresses.**

---

## Key Findings

### Finding 1 — Trophic level is the strongest predictor of PFAS accumulation

Across 13,098 observations, trophic level (where an organism sits in the food chain) explains more variance in PFAS tissue concentration than any chemical property. This directly confirms biomagnification: the higher you are in the food chain, the more PFAS you accumulate. Humans at trophic level 5 show the highest and most consistent concentrations.

### Finding 2 — Human blood levels are predictable; environmental data is not

When CDC NHANES human biomonitoring data is included, model R² improves from 0.174 to 0.742. This is because human blood measurements (2,133 people, same lab, same protocol) are far more consistent than heterogeneous environmental studies. The model learned human blood PFAS patterns extremely well — but cannot transfer that knowledge to fish or plant predictions.

### Finding 3 — Environmental data cannot predict human exposure (cross-species R²= -8)

Leave-one-species-out validation shows that a model trained on fish, plant, and mammal data achieves R²= -8 when predicting human blood levels. This is a critical finding: **we cannot use environmental bioaccumulation data to predict human exposure with current available data.** The disconnect between environmental measurements and human biomonitoring is the central unsolved problem in PFAS risk assessment.

### Finding 4 — BCF is a better normalized target than raw concentration (R²=0.408 vs 0.174)

Bioconcentration Factor (BCF) normalizes tissue concentration by water concentration, removing a major source of cross-study noise. Our Random Forest predicts BCF more than twice as accurately as raw tissue concentration, confirming that study-level variation in exposure dose — not chemistry — drives most of the variance in raw concentration data.

### Finding 5 — Data gaps are severe, systematic, and chemically biased

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

**PFHxS is detected in the blood of nearly every American — yet has zero fish and zero mammal tissue records in ECOTOX.** We know it's in human blood but have almost no data on how it gets there through the food chain.

### Finding 6 — Mammalian bioaccumulation data is nearly unusable

ECOTOX mammal records for PFAS are almost entirely dose-response studies (mg/kg/day administered) rather than tissue residue measurements (ng/g accumulated). Only 4 mammalian tissue residue records exist across all 8 PFAS. This makes cross-species prediction from environmental mammals to humans currently impossible.

---

## Version History

### v4.0 (current) — June 2026
- Added CDC NHANES 2017-2018 human blood serum data (11,574 rows, 6 PFAS, 2,133 people)
- Dataset grew from 1,524 → 13,098 rows (+760%)
- Concentration model R² improved from 0.174 → 0.742
- Human column now populated in gap heatmap
- NHANES data loaded post-harmonization (pre-calculated concentrations)
- Discovered critical cross-species prediction failure (human R²= -8 in leave-one-out)
- Full 5-species gap analysis now possible

### v3.0 — June 2026
- Added BCF as second ML target (R²=0.408)
- Added PFAS class (Sulfonate/Carboxylate) as feature
- Clean pipeline rewrite — both models in single integrated script
- Dataset: 1,524 rows, R²=0.174
- All 7 outputs generated in one command

### v2.1 — May 2026
- Added terrestrial ECOTOX exports
- Added Tier 2 PFAS (PFDA, PFUnDA, PFHpA)
- Dataset grew from 697 → 1,273 rows
- Plant species group appeared for first time
- Removed Route_encoded (identified as data leakage)

### v2.0 — May 2026
- Added aquatic exports for PFNA, PFHxS, PFBS, PFOA
- First gap heatmap and cross-species validation
- R²=0.279

### v1.0 — April 2026
- Initial pipeline, PFOS only (411 rows)
- First Random Forest model

---

## Outputs

| File | Description |
|---|---|
| `pfas_bioaccumulation_dataset.csv` | 13,098 rows, fully cleaned and ML-ready |
| `pfas_gap_heatmap.png` | Observations per PFAS × species group |
| `feature_importance.png` | Concentration model feature importances |
| `model_predictions.png` | Concentration predicted vs actual (R²=0.742) |
| `cross_species_validation.png` | Leave-one-species-out R² and RMSE |
| `bcf_feature_importance.png` | BCF model feature importances |
| `bcf_predictions.png` | BCF predicted vs actual (R²=0.408) |

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
| `is_fish` | int | Binary species flag |
| `is_mammal` | int | Binary species flag |
| `is_plant` | int | Binary species flag |
| `is_human` | int | Binary species flag |
| `Tissue` | string | Tissue or response site |
| `Exposure Route` | string | Exposure type |
| `Duration_days` | float | Exposure duration in days |
| `Concentration_ng_g` | float | Measured tissue concentration (ng/g) |
| `log_concentration` | float | Log₁₀ concentration — ML target 1 |
| `BCF` | float | Bioconcentration factor |
| `log_BCF` | float | Log₁₀ BCF — ML target 2 |
| `Source` | string | CDC NHANES / ECOTOX |

---

## PFAS Chemicals

### Tier 1 — Core (fully included)
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

### Tier 3 — Emerging (not yet included)
| PFAS | CASRN | Class |
|---|---|---|
| GenX (HFPO-DA) | 13252-13-6 | Carboxylate |
| ADONA | 958445-44-8 | Carboxylate |
| F53B | 73606-19-6 | Sulfonate |

---

## Setup

```bash
pip3 install pandas numpy matplotlib seaborn scikit-learn openpyxl requests
```

---

## Usage

### 1. Set paths in `pfas_pipeline_v4.py`
```python
ECOTOX_EXPORT_DIR = "/path/to/ecotox_exports/"
COMPTOX_SNAPSHOT  = "/path/to/comptox_snapshot.csv"
OUTPUT_DIR        = "/path/to/outputs/"
```

### 2. Run
```bash
python3 pfas_pipeline_v4.py
```

### Required files
| File | Where to get it |
|---|---|
| ECOTOX xlsx exports | https://cfpub.epa.gov/ecotox/ — search each PFAS, export as XLSX |
| `nhanes_pfas_processed.csv` | Included in this repo — generated from CDC NHANES 2017-2018 |
| `comptox_snapshot.csv` | Optional — https://comptox.epa.gov/dashboard/batch-search |

---

## Data Sources

| Source | URL | What it provides |
|---|---|---|
| EPA ECOTOX | https://cfpub.epa.gov/ecotox/ | Species, tissue, concentration, BCF |
| EPA CompTox | https://comptox.epa.gov/dashboard/batch-search | MW, LogKow, chemical properties |
| CDC NHANES 2017-2018 | https://wwwn.cdc.gov/nchs/nhanes/ | Human blood serum PFAS levels |

---

## Pipeline Architecture

```
EPA CompTox              EPA ECOTOX               CDC NHANES
(chemical traits)    (species + tissue data)   (human blood serum)
      │                      │                        │
      ▼                      ▼                        │
Chemical Feature       ECOTOX Ingestion               │
Table (13 PFAS)        (18 xlsx files)                │
      │                      │                        │
      │               Harmonization                   │
      │           (units → ng/g, taxonomy)            │
      │                      │                        │
      │                      └──────────┬─────────────┘
      │                                 │
      │                          Combined Dataset
      └─────────────────────────────────┘
                                │
                         Merged Dataset
                        (13,098 rows)
                                │
              ┌─────────────────┼─────────────────┐
              ▼                 ▼                  ▼
         Gap Heatmap    Conc RF Model        BCF RF Model
         (5 species     (R²=0.742)           (R²=0.408)
          × 8 PFAS)
```

---

## Model Results

| Version | Rows | R² Conc | R² BCF | Key change |
|---|---|---|---|---|
| v1.0 | 411 | — | — | PFOS only |
| v2.0 | 697 | 0.279* | — | 5 PFAS aquatic |
| v2.1 | 1,273 | 0.110 | — | Terrestrial added, leakage fixed |
| v3.0 | 1,524 | 0.176 | 0.408 | BCF model added |
| v4.0 | 13,098 | 0.742 | 0.408 | NHANES human data added |

*R²=0.279 included Route_encoded data leakage

### Feature importance (concentration model, v4.0)
1. Trophic Level (0.55) — food chain position
2. is_human (0.15) — human vs environmental
3. LogKow (0.13) — hydrophobicity
4. Chain Length (0.12) — carbon chain length
5. MW (0.05) — molecular weight

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

### Phase 1 — Data Expansion (June 2026)
- [x] PFAS class as ML feature
- [x] BCF as second target variable
- [x] Terrestrial species data
- [x] Tier 2 PFAS
- [x] CDC NHANES human blood serum data
- [ ] Real CompTox physicochemical properties
- [ ] Real trophic levels per species

### Phase 2 — Better Models (July 2026)
- [ ] XGBoost / Gradient Boosting
- [ ] Linear Regression baseline
- [ ] Per-PFAS models
- [ ] Prediction confidence intervals
- [ ] BCF cross-species validation

### Phase 3 — Outputs (August 2026)
- [ ] Streamlit interactive dashboard
- [ ] Research poster
- [ ] Full methods + results report
- [ ] Published GitHub repository

---

## How to Add Data

### New ECOTOX exports
1. Go to https://cfpub.epa.gov/ecotox/
2. Search chemical → filter **Effect → Accumulation** if >10,000 results
3. Export as XLSX → move to `ecotox_exports/`
4. Run pipeline — all files picked up automatically

### New PFAS chemicals
Add to `PFAS_FEATURES` in `pfas_pipeline_v4.py`:
```python
("PFDA", "335-76-2", "Carboxylate", 10, 514.1, 6.83),
# (Name, CASRN, Class, Chain_Length, MW, LogKow)
```

### New unit conversions
Add to `UNIT_TO_NG_G` in `pfas_pipeline_v4.py`:
```python
"your_unit": conversion_factor_to_ng_per_g,
```

---

## How to Push to GitHub

### First time setup
```bash
brew install git
git config --global user.name "Your Name"
git config --global user.email "your@email.com"

# Go to github.com → New Repository → name "pfas-bioaccumulation"

cd /Users/navink.admin/Desktop
git init
git remote add origin https://github.com/YOUR_USERNAME/pfas-bioaccumulation.git

# Create .gitignore
echo "ecotox_exports/
comptox_snapshot.csv
pfas_bioaccumulation_dataset.csv
PFAS_J.XPT
__pycache__/
.DS_Store" > .gitignore

# Add files
git add pfas_pipeline_v4.py
git add README.md
git add nhanes_pfas_processed.csv
git add pfas_gap_heatmap.png feature_importance.png model_predictions.png
git add cross_species_validation.png bcf_feature_importance.png bcf_predictions.png

git commit -m "v4.0: NHANES integration, 13098 rows, R2=0.742"
git branch -M main
git push -u origin main
```

### Every update
```bash
git add pfas_pipeline_v4.py README.md
git commit -m "describe what changed"
git push
```

### What goes on GitHub
| Include | Exclude |
|---|---|
| `pfas_pipeline_v4.py` | `ecotox_exports/` folder |
| `README.md` | `PFAS_J.XPT` |
| `nhanes_pfas_processed.csv` | `pfas_bioaccumulation_dataset.csv` |
| All 6 PNG plots | `comptox_snapshot.csv` |

---

## Citation

- EPA ECOTOX Knowledgebase: https://cfpub.epa.gov/ecotox/
- EPA CompTox Dashboard: https://comptox.epa.gov/dashboard/
- CDC NHANES 2017-2018 PFAS Data: https://wwwn.cdc.gov/nchs/nhanes/

---

## Author

PFAS Environmental Informatics Research Project, 2026
