# PFAS Bioaccumulation Research Pipeline

A reproducible, multi-source data pipeline for studying PFAS bioaccumulation across aquatic and terrestrial species. Integrates EPA ECOTOX biological exposure data with EPA CompTox chemical properties to build a machine-learning-ready dataset and identify critical data gaps.

**Current dataset: 1,454 observations | 8 PFAS | 4 species groups | 1,056 ML-ready rows**

---

## Table of Contents

- [Project Overview](#project-overview)
- [Key Findings](#key-findings)
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
- [How to Contribute / Add Data](#how-to-contribute--add-data)

---

## Project Overview

Per- and polyfluoroalkyl substances (PFAS) are persistent synthetic chemicals that accumulate in biological tissues across food webs. Despite regulatory concern, bioaccumulation data remains fragmented across species, tissues, and study designs.

This pipeline:
1. Ingests raw ECOTOX export files (aquatic + terrestrial)
2. Harmonizes units, species taxonomy, and exposure metadata
3. Merges with CompTox chemical property data
4. Outputs a clean ML-ready dataset
5. Trains a Random Forest to predict log tissue concentration
6. Generates data gap heatmaps and cross-species validation plots

---

## Key Findings

**Finding 1 — Chemistry drives bioaccumulation**
Molecular weight (MW), hydrophobicity (LogKow), and chain length are the top 3 predictors of tissue concentration — consistent with established bioaccumulation theory.

**Finding 2 — Data gaps are severe and systematic**
- PFHxS has zero fish and zero mammal tissue concentration records
- Mammalian bioaccumulation data is almost entirely dose-response studies, not tissue residue measurements
- Emerging PFAS (GenX, ADONA, F53B) have no usable bioaccumulation data
- Short-chain PFAS (PFBS, PFHxS) are severely underrepresented

**Finding 3 — Study design dominates over chemistry**
Exposure route explains more variance than chemical properties alone, suggesting inconsistent study designs limit cross-study prediction.

---

## Outputs

| File | Description |
|---|---|
| `pfas_bioaccumulation_dataset.csv` | 1,454 rows, cleaned and ML-ready |
| `pfas_gap_heatmap.png` | Observations per PFAS × species group |
| `feature_importance.png` | Random Forest feature importances |
| `model_predictions.png` | Predicted vs actual log₁₀ concentration |
| `cross_species_validation.png` | Leave-one-species-out R² and RMSE |

---

## Dataset Schema

Each row is one bioaccumulation observation.

| Column | Type | Description |
|---|---|---|
| `PFAS_Name` | string | Common name (e.g. PFOS, PFOA) |
| `CASRN` | string | CAS Registry Number |
| `PFAS_Class` | string | Sulfonate or Carboxylate |
| `Chain_Length` | int | Fluorinated carbon chain length |
| `MW` | float | Molecular weight (g/mol) |
| `LogKow` | float | Octanol-water partition coefficient |
| `Species` | string | Common species name |
| `Species_Group` | string | Fish / Mammal / Plant / Human / Other |
| `Trophic_Level` | int | 1 (plant) → 5 (human) |
| `Is_Aquatic` | int | 1 = aquatic, 0 = terrestrial |
| `Tissue` | string | Tissue or response site |
| `Exposure Route` | string | Exposure type (Static, Flow-through, etc.) |
| `Duration_days` | float | Exposure duration in days |
| `Concentration_ng_g` | float | Measured concentration (ng/g) |
| `log_concentration` | float | Log₁₀ concentration — ML target variable |

---

## PFAS Chemicals

### Tier 1 — Core (fully included)
| PFAS | CASRN | Class | Chain Length | MW | LogKow |
|---|---|---|---|---|---|
| PFOS | 1763-23-1 | Sulfonate | 8 | 500.1 | 5.26 |
| PFOA | 335-67-1 | Carboxylate | 8 | 414.1 | 5.30 |
| PFHxS | 355-46-4 | Sulfonate | 6 | 400.1 | 4.14 |
| PFNA | 375-95-1 | Carboxylate | 9 | 464.1 | 6.05 |
| PFBS | 375-73-5 | Sulfonate | 4 | 300.1 | 1.82 |

### Tier 2 — Extended (partially included)
| PFAS | CASRN | Class | Chain Length | MW | LogKow |
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

### Requirements

```
Python 3.8+
pandas
numpy
matplotlib
seaborn
scikit-learn
openpyxl
requests
```

Install:
```bash
pip3 install pandas numpy matplotlib seaborn scikit-learn openpyxl requests
```

### Folder Structure

```
pfas-pipeline/
├── pfas_pipeline.py              ← main script
├── README.md                     ← this file
├── comptox_snapshot.csv          ← CompTox bulk download (add manually)
├── ecotox_exports/               ← ECOTOX xlsx exports (add manually)
│   ├── ECOTOX_PFOS.xlsx
│   ├── ECOTOX-Aquatic-Export_*.xlsx
│   └── ECOTOX-Terrestrial-Export_*.xlsx
└── outputs/                      ← generated automatically
    ├── pfas_bioaccumulation_dataset.csv
    ├── pfas_gap_heatmap.png
    ├── feature_importance.png
    ├── model_predictions.png
    └── cross_species_validation.png
```

---

## Usage

### 1. Set paths in `pfas_pipeline.py`

```python
ECOTOX_EXPORT_DIR = "/path/to/ecotox_exports/"
COMPTOX_SNAPSHOT  = "/path/to/comptox_snapshot.csv"
OUTPUT_DIR        = "/path/to/outputs/"
```

### 2. Run

```bash
python3 pfas_pipeline.py
```

---

## Data Sources

| Source | URL | What it provides |
|---|---|---|
| EPA ECOTOX | https://cfpub.epa.gov/ecotox/ | Species, tissue, concentration data |
| EPA CompTox | https://comptox.epa.gov/dashboard/batch-search | MW, LogKow, chemical properties |
| CDC NHANES | https://www.cdc.gov/nchs/nhanes/ | Human blood serum PFAS levels (planned) |

---

## Pipeline Architecture

```
CompTox (chemical traits)
        │
        ▼
Chemical Master Table
(PFAS name, CASRN, MW, LogKow, Chain Length, Class)
        │
        ▼
ECOTOX Data Pull
(species, tissue, concentration, exposure)
        │
        ▼
Harmonization Layer
(unit conversion → ng/g, duration → days, species mapping)
        │
        ▼
Joined ML-Ready Dataset (1,454 rows)
        │
   ┌────┴────┐
   ▼         ▼
Gap       Random Forest
Analysis  (log concentration prediction)
Heatmap   Feature Importance
          Cross-species Validation
```

---

## Model Results

| Version | Rows | Features | R² | RMSE |
|---|---|---|---|---|
| v1 — Aquatic only | 697 | 11 | 0.110 | 1.294 |
| v2 — + Terrestrial | 1,273 | 10 | 0.176 | 1.270 |
| v3 — + Tier 2 PFAS | 1,454 | 10 | 0.176 | 1.270 |

### Feature importance ranking
1. Molecular Weight (MW)
2. LogKow
3. Chain Length
4. Trophic Level
5. Is Aquatic

---

## Data Gaps

| PFAS | Fish | Mammal | Plant | Other |
|---|---|---|---|---|
| PFOA | 150 | 0 | 104 | 326 |
| PFOS | 22 | 1 | 114 | 118 |
| PFBS | 6 | 0 | 0 | 66 |
| PFHxS | 0 | 0 | 0 | 46 |
| PFNA | 16 | 0 | 0 | 40 |
| PFDA | 0 | 0 | 0 | 27 |
| PFHpA | 0 | 0 | 0 | 1 |

**Key gap: PFHxS has zero fish and zero mammal tissue records.**

---

## Roadmap

### Phase 1 — Data Expansion (June)
- [ ] Add PFAS class (Sulfonate/Carboxylate) as ML feature
- [ ] Add BCF as second target variable
- [ ] Integrate CDC NHANES human blood serum data
- [ ] Add real trophic levels per species
- [ ] Add CompTox real physicochemical properties

### Phase 2 — Better Models (July)
- [ ] Add XGBoost / Gradient Boosting
- [ ] Add Linear Regression baseline
- [ ] Train per-PFAS models
- [ ] Add prediction confidence intervals

### Phase 3 — Outputs (August)
- [ ] Build Streamlit interactive dashboard
- [ ] Design research poster
- [ ] Write full methods + results report
- [ ] Publish GitHub repository

---

## How to Contribute / Add Data

### Adding new ECOTOX exports

1. Go to https://cfpub.epa.gov/ecotox/
2. Search a chemical name (e.g. PFDA)
3. Filter by **Effect → Accumulation** if results > 10,000
4. Click **Export as** → XLSX
5. Move the file into `ecotox_exports/`:
```bash
mv ~/Downloads/ECOTOX-*.xlsx /path/to/ecotox_exports/
```
6. Run the pipeline — it picks up all files automatically:
```bash
python3 pfas_pipeline.py
```

### Adding new PFAS chemicals

Open `pfas_pipeline.py` and add to the `PFAS_FEATURES` table:

```python
("PFDA", "335-76-2", "Carboxylate", 10, 514.1, 6.83),
```

Fields: `(Name, CASRN, Class, Chain_Length, MW, LogKow)`

### Adding new unit conversions

If a new ECOTOX file has an unrecognized unit, add it to `UNIT_TO_NG_G`:

```python
UNIT_TO_NG_G = {
    "ng/g": 1,
    "ug/g": 1_000,
    "your_new_unit": conversion_factor,  # ← add here
    ...
}
```

---

## How to Push to GitHub

### First time setup

```bash
# 1. Install git if you don't have it
brew install git

# 2. Configure your identity
git config --global user.name "Your Name"
git config --global user.email "your@email.com"

# 3. Go to github.com → New Repository → name it "pfas-pipeline"
#    Make it Public, don't add README (we have one)

# 4. Initialize git in your project folder
cd /Users/navink.admin/Desktop
git init
git remote add origin https://github.com/YOUR_USERNAME/pfas-pipeline.git
```

### What to commit (and what NOT to)

```bash
# Create a .gitignore to exclude large data files
cat > .gitignore << 'EOF'
ecotox_exports/
comptox_snapshot.csv
pfas_bioaccumulation_dataset.csv
__pycache__/
.DS_Store
EOF
```

### First commit

```bash
git add pfas_pipeline.py
git add README.md
git add pfas_gap_heatmap.png
git add feature_importance.png
git add model_predictions.png
git add cross_species_validation.png
git commit -m "Initial commit: PFAS bioaccumulation pipeline v1"
git branch -M main
git push -u origin main
```

### Every time you make changes

```bash
git add pfas_pipeline.py
git add README.md
git commit -m "describe what you changed here"
git push
```

### Good commit message examples
```
"Add terrestrial ECOTOX data, dataset grows to 1273 rows"
"Add Tier 2 PFAS: PFDA, PFUnDA, PFHpA"
"Fix unit conversion for soil concentration units"
"Add XGBoost model, R2 improves to 0.31"
```

---

## Citation

If using this pipeline in academic work:

- EPA ECOTOX Knowledgebase: https://cfpub.epa.gov/ecotox/
- EPA CompTox Dashboard: https://comptox.epa.gov/dashboard/

---

## Author

PFAS Environmental Informatics Research Project, 2026
