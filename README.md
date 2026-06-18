# PFAS Bioaccumulation Research Pipeline v3.0

> *"PFAS chemicals don't break down. They move up food chains, accumulate in tissues, and end up in human blood. This pipeline maps where they go — and where our scientific understanding runs out."*

A reproducible, multi-source data pipeline for studying PFAS bioaccumulation across aquatic and terrestrial species. Integrates EPA ECOTOX biological exposure data with EPA CompTox chemical properties to build a machine-learning-ready dataset, identify critical data gaps, and predict bioaccumulation from chemical structure alone.

**Current dataset: 1,524 observations | 8 PFAS | 4 species groups | 2 ML models**

---

## Why This Research Matters

### PFAS Are Everywhere — And They Don't Leave

Per- and polyfluoroalkyl substances (PFAS) are a class of over 12,000 synthetic chemicals used in non-stick cookware, food packaging, firefighting foam, waterproof clothing, and hundreds of industrial applications. They are called "forever chemicals" for a reason: the carbon-fluorine bond is one of the strongest in chemistry. PFAS do not break down in the environment. They do not break down in the human body.

They accumulate.

### The Bioaccumulation Problem

When PFAS enter an ecosystem — through industrial discharge, agricultural runoff, or contaminated groundwater — they are absorbed by plants and small organisms at the base of the food chain. As larger animals eat smaller ones, PFAS concentrations multiply at each trophic level. This process, called biomagnification, means that a fish at the top of an aquatic food chain can carry concentrations thousands of times higher than the water it swims in.

Humans sit at the top of the food chain.

### What We Know — and What We Don't

Studies have detected PFAS in the blood of 97% of Americans. PFAS exposure has been linked to thyroid disease, immune suppression, certain cancers, reproductive harm, and developmental delays in children. The EPA has set maximum contaminant levels for several PFAS in drinking water at 4 parts per trillion — a level so low it required new analytical methods to even measure.

Yet despite this, our scientific understanding of *how* PFAS move through ecosystems remains fragmented. Data is scattered across hundreds of studies, measured in inconsistent units, tested on different species, and reported under different conditions. No single database cleanly maps PFAS bioaccumulation from soil → plant → fish → mammal → human.

**That is the gap this project addresses.**

---

## What This Pipeline Does

1. Ingests raw ECOTOX biological exposure data (aquatic + terrestrial)
2. Harmonizes units, species taxonomy, and exposure metadata across 18+ data files
3. Merges with EPA CompTox chemical property data
4. Outputs a clean, ML-ready dataset with 1,524 observations
5. Trains two Random Forest models — one predicting tissue concentration, one predicting BCF
6. Generates data gap heatmaps showing exactly where scientific knowledge is missing
7. Runs cross-species validation to test whether PFAS behavior in one species can predict another

---

## Version History

### v3.0 (current) — June 2026
**What's new:**
- Added BCF (Bioconcentration Factor) as a second ML target variable
- BCF model achieves R²=0.408 vs concentration model R²=0.174 — a 2.3× improvement
- Added PFAS class (Sulfonate/Carboxylate) as a model feature
- Added terrestrial species data (soil, plant, diet exposure studies)
- Dataset grew from 697 → 1,524 rows (+119%)
- ML training set grew from 501 → 1,123 rows (+124%)
- Clean pipeline rewrite — both models now run inside a single integrated script
- All 7 output files generated in one command

**Key insight from v3.0:**
BCF is a significantly better prediction target than raw tissue concentration because it normalizes for exposure dose, removing a major source of cross-study noise.

---

### v2.1 — June 2026
- Added terrestrial ECOTOX exports (soil, plant, mammal studies)
- Dataset grew from 697 → 1,273 rows
- R² improved from 0.110 → 0.176 after adding species diversity
- Identified that exposure route was dominating model (data leakage) — removed
- Plant species group appeared in data for first time

---

### v2.0 — June 2026
- Added aquatic exports for PFNA, PFHxS, PFBS, PFOA
- First working gap heatmap
- First cross-species validation
- R²=0.279 with Route_encoded (later identified as leakage)

---

### v1.0 — June 2026
- Initial pipeline with PFOS only (411 rows)
- Basic harmonization layer
- First Random Forest model

---

## Key Findings

### Finding 1 — BCF is a better prediction target than raw concentration (R²=0.408 vs 0.174)

Raw tissue concentration varies enormously between studies because it depends on how much chemical the organism was exposed to — which differs by experiment. BCF normalizes this by dividing tissue concentration by water concentration, making it comparable across studies. Our Random Forest predicts BCF more than twice as accurately as raw concentration, confirming that study-level variation in exposure dose is a major source of noise in the literature.

**Implication:** Future PFAS bioaccumulation research should prioritize BCF measurements over raw tissue concentration for cross-study comparability.

---

### Finding 2 — Molecular weight, LogKow, and chain length are the top chemical predictors

Across both models, the same three chemical properties dominate feature importance:

1. **Molecular Weight (MW)** — heavier molecules partition more strongly into lipid-rich tissues
2. **LogKow** — higher octanol-water partition coefficient = greater hydrophobicity = more bioaccumulation
3. **Chain Length** — longer fluorinated carbon chains are more persistent and more bioaccumulative

This is consistent with established bioaccumulation theory and validates that our pipeline is capturing real chemical mechanisms, not statistical artifacts.

---

### Finding 3 — Data gaps are severe, systematic, and chemically biased

| PFAS | Fish | Mammal | Plant | Other |
|---|---|---|---|---|
| PFOA | 150 | 0 | 104 | 326 |
| PFOS | 22 | 1 | 114 | 118 |
| PFBS | 6 | 0 | 0 | 66 |
| PFHxS | 0 | 0 | 0 | 46 |
| PFNA | 16 | 0 | 0 | 40 |
| PFDA | 0 | 0 | 0 | 27 |
| PFHpA | 0 | 0 | 0 | 1 |

**PFHxS has zero fish and zero mammal tissue concentration records** despite being one of the most detected PFAS in human blood globally. PFBS, a common replacement for PFOS, has almost no mammal data. All Tier 2 PFAS (PFDA, PFUnDA, PFHpA) lack fish and mammal data entirely.

This is not a limitation of our pipeline — it is a reflection of the actual state of the scientific literature.

---

### Finding 4 — Mammalian bioaccumulation data is nearly unusable

When we searched ECOTOX for mammalian PFAS data, the vast majority of records were dose-response studies (mg/kg/day administered) rather than tissue residue measurements (ng/g in tissue). These two measurement types answer fundamentally different questions. The former tells you how much PFAS kills a rat. The latter tells you how much PFAS a rat accumulates in its liver from environmental exposure.

Only 1 mammalian tissue residue record exists in our current dataset across all 8 PFAS. This means cross-species prediction from fish to mammals — which is exactly what regulators need to set safe exposure limits for humans — is currently impossible with available data.

---

### Finding 5 — Cross-species prediction fails with current data

Our leave-one-species-out validation shows negative R² for all species groups, meaning a model trained on fish cannot predict plant accumulation, and vice versa. This is itself a finding: PFAS bioaccumulation is so species- and tissue-specific that chemical properties alone cannot bridge the gap. Environmental fate models that assume consistent cross-species behavior may significantly underestimate risk for some organisms.

---

## Model Results

### Concentration Model (log₁₀ ng/g)

| Metric | Value |
|---|---|
| Training rows | 1,123 |
| Features | 11 |
| R² | 0.174 |
| RMSE | ~1.27 log₁₀ ng/g |

### BCF Model (log₁₀ BCF)

| Metric | Value |
|---|---|
| Training rows | 321 |
| Features | 10 |
| R² | 0.408 |
| RMSE | ~0.85 log₁₀ BCF |

### Feature Importance (both models, ranked)
1. Molecular Weight
2. LogKow
3. Chain Length
4. Trophic Level
5. Is Aquatic
6. PFAS Class (Sulfonate vs Carboxylate)

---

## Outputs

| File | Description |
|---|---|
| `pfas_bioaccumulation_dataset.csv` | 1,524 rows, fully cleaned and ML-ready |
| `pfas_gap_heatmap.png` | Observations per PFAS × species group |
| `feature_importance.png` | Concentration model feature importances |
| `model_predictions.png` | Concentration predicted vs actual (R²=0.174) |
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

---

## Setup

```bash
pip3 install pandas numpy matplotlib seaborn scikit-learn openpyxl requests
```

Set paths in `pfas_pipeline_v3.py`:
```python
ECOTOX_EXPORT_DIR = "/path/to/ecotox_exports/"
COMPTOX_SNAPSHOT  = "/path/to/comptox_snapshot.csv"
OUTPUT_DIR        = "/path/to/outputs/"
```

Run:
```bash
python3 pfas_pipeline_v3.py
```

---

## Data Sources

| Source | URL | What it provides |
|---|---|---|
| EPA ECOTOX | https://cfpub.epa.gov/ecotox/ | Species, tissue, concentration, BCF |
| EPA CompTox | https://comptox.epa.gov/dashboard/batch-search | MW, LogKow, chemical properties |
| CDC NHANES | https://www.cdc.gov/nchs/nhanes/ | Human blood serum PFAS (planned) |

---

## How to Add New ECOTOX Data

1. Go to https://cfpub.epa.gov/ecotox/
2. Search a chemical → filter **Effect → Accumulation** if >10,000 results
3. Export as XLSX
4. Move to `ecotox_exports/`:
```bash
mv ~/Downloads/ECOTOX-*.xlsx /path/to/ecotox_exports/
```
5. Run pipeline — all files are picked up automatically

---

## How to Push to GitHub

### First time
```bash
brew install git
git config --global user.name "Your Name"
git config --global user.email "your@email.com"

cd /Users/navink.admin/Desktop
git init
git remote add origin https://github.com/YOUR_USERNAME/pfas-pipeline.git

# Create .gitignore
cat > .gitignore << 'EOF'
ecotox_exports/
comptox_snapshot.csv
pfas_bioaccumulation_dataset.csv
__pycache__/
.DS_Store
EOF

git add pfas_pipeline_v3.py README.md *.png
git commit -m "v3.0: dual model pipeline, BCF R2=0.408, 1524 rows"
git branch -M main
git push -u origin main
```

### Every update
```bash
git add pfas_pipeline_v3.py README.md *.png
git commit -m "describe what changed"
git push
```

---

## Roadmap

### Phase 1 — Data Expansion (June 2026)
- [x] Add PFAS class as ML feature
- [x] Add BCF as second target variable
- [x] Add terrestrial species data
- [x] Add Tier 2 PFAS (PFDA, PFUnDA, PFHpA)
- [ ] Integrate CDC NHANES human blood serum data
- [ ] Add real CompTox physicochemical properties
- [ ] Add actual trophic levels per species

### Phase 2 — Better Models (July 2026)
- [ ] Add XGBoost / Gradient Boosting
- [ ] Add Linear Regression baseline
- [ ] Train per-PFAS models
- [ ] Add prediction confidence intervals
- [ ] BCF cross-species validation

### Phase 3 — Outputs (August 2026)
- [ ] Build Streamlit interactive dashboard
- [ ] Design research poster
- [ ] Write full methods + results report
- [ ] Publish GitHub repository

---

## Citation

EPA ECOTOX Knowledgebase: https://cfpub.epa.gov/ecotox/
EPA CompTox Dashboard: https://comptox.epa.gov/dashboard/

---

## Author

PFAS Environmental Informatics Research Project, 2026
