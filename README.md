# PFAS Bioaccumulation Data Pipeline (ECOTOX + CompTox Integration)

A reproducible data pipeline for building a machine-learning-ready dataset of PFAS bioaccumulation across species, tissues, and exposure conditions. Integrates EPA ECOTOX biological data with EPA CompTox chemical properties.

---

## Project Overview

Per- and polyfluoroalkyl substances (PFAS) are persistent environmental contaminants that bioaccumulate across food webs. This pipeline:

1. Ingests raw ECOTOX export files (aquatic + terrestrial)
2. Harmonizes units, species taxonomy, and exposure metadata
3. Merges biological observations with chemical features from CompTox
4. Outputs a clean ML-ready dataset
5. Generates data gap heatmaps and Random Forest model outputs

---

## Outputs

| File | Description |
|---|---|
| `pfas_bioaccumulation_dataset.csv` | Cleaned, ML-ready dataset (697+ rows) |
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
| `Trophic_Level` | int | 1 (plant) to 5 (human) |
| `Is_Aquatic` | int | 1 = aquatic, 0 = terrestrial |
| `Tissue` | string | Tissue or response site |
| `Exposure Route` | string | Exposure type (e.g. Static, Flow-through) |
| `Duration_days` | float | Exposure duration in days |
| `Concentration_ng_g` | float | Measured concentration (ng/g) |
| `log_concentration` | float | Log₁₀ of concentration — ML target variable |

---

## PFAS Chemicals Included

### Tier 1 (core)
| PFAS | CASRN | Class | Chain Length |
|---|---|---|---|
| PFOS | 1763-23-1 | Sulfonate | 8 |
| PFOA | 335-67-1 | Carboxylate | 8 |
| PFHxS | 355-46-4 | Sulfonate | 6 |
| PFNA | 375-95-1 | Carboxylate | 9 |
| PFBS | 375-73-5 | Sulfonate | 4 |

### Tier 2 (extended)
| PFAS | CASRN | Class | Chain Length |
|---|---|---|---|
| PFDA | 335-76-2 | Carboxylate | 10 |
| PFUnDA | 2058-94-8 | Carboxylate | 11 |
| PFDoDA | 307-55-1 | Carboxylate | 12 |
| PFHpA | 375-85-9 | Carboxylate | 7 |
| PFHxA | 307-24-4 | Carboxylate | 6 |

### Tier 3 (emerging)
| PFAS | CASRN | Class |
|---|---|---|
| GenX (HFPO-DA) | 13252-13-6 | Carboxylate |
| ADONA | 958445-44-8 | Carboxylate |
| F53B | 73606-19-6 | Sulfonate |

---

## Data Sources

### ECOTOX Knowledgebase
- URL: https://cfpub.epa.gov/ecotox/
- Search by chemical name → filter by Effect: Accumulation → Export as XLSX
- Place exports in: `ecotox_exports/`

### EPA CompTox Dashboard
- URL: https://comptox.epa.gov/dashboard/batch-search
- Batch search by CASRN list → Export Physical/Chemical Properties
- Save as: `comptox_snapshot.csv`

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

Install all dependencies:

```bash
pip3 install pandas numpy matplotlib seaborn scikit-learn openpyxl requests
```

### Directory Structure

```
Desktop/
├── pfas_pipeline.py
├── comptox_snapshot.csv          ← CompTox bulk download
├── ecotox_exports/
│   ├── ECOTOX_PFOS.xlsx
│   ├── ECOTOX_PFOA.xlsx
│   └── ...
└── outputs/
    ├── pfas_bioaccumulation_dataset.csv
    ├── pfas_gap_heatmap.png
    ├── feature_importance.png
    ├── model_predictions.png
    └── cross_species_validation.png
```

---

## Usage

### 1. Edit input paths in `pfas_pipeline.py`

```python
ECOTOX_EXPORT_DIR = "/Users/yourname/Desktop/ecotox_exports/"
COMPTOX_SNAPSHOT  = "/Users/yourname/Desktop/comptox_snapshot.csv"
OUTPUT_DIR        = "/Users/yourname/Desktop/"
```

### 2. Run the pipeline

```bash
python3 pfas_pipeline.py
```

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
Joined ML-Ready Dataset
        │
   ┌────┴────┐
   ▼         ▼
Gap       Random Forest
Analysis  (BCF / log conc prediction)
```

---

## Model Results (current)

| Metric | Value |
|---|---|
| Training rows | 501 |
| Features | 11 |
| R² | 0.279 |
| RMSE | 1.268 log₁₀ ng/g |

### Top features by importance
1. Exposure Route
2. LogKow
3. Molecular Weight
4. Chain Length
5. Is Fish

---

## Research Questions

1. Can PFAS tissue concentration be predicted from chemical and species traits alone?
2. Which physicochemical properties (LogKow, chain length, MW) matter most?
3. Which PFAS–species combinations have insufficient data for reliable prediction?
4. Can a model trained on fish generalize to mammals or plants?

---

## Known Data Gaps (current dataset)

| PFAS | Fish | Other | Mammal | Plant | Human |
|---|---|---|---|---|---|
| PFOA | 150 | 153 | 0 | 0 | 0 |
| PFOS | 22 | 89 | 0 | 0 | 0 |
| PFNA | 8 | 20 | 0 | 0 | 0 |
| PFHxS | 0 | 23 | 0 | 0 | 0 |
| PFBS | 3 | 33 | 0 | 0 | 0 |

PFHxS has **zero fish observations** — a significant research gap.
Mammal, Plant, and Human data are entirely missing from current exports.

---

## Next Steps

- [ ] Add terrestrial ECOTOX exports (mammal, plant data)
- [ ] Add Tier 2 and Tier 3 PFAS exports
- [ ] Integrate CompTox bulk download for real physicochemical properties
- [ ] Push dataset to 1000+ rows
- [ ] Add uncertainty scoring per observation
- [ ] Add BCF as an alternative target variable
- [ ] Build interactive dashboard (Streamlit or Dash)

---

## Citation

If using this pipeline in academic work, cite:

- EPA ECOTOX Knowledgebase: https://cfpub.epa.gov/ecotox/
- EPA CompTox Dashboard: https://comptox.epa.gov/dashboard/

---

## Author

Built as part of a PFAS environmental informatics research project.
