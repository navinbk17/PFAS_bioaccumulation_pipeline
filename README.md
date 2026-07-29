# PFAS Bioaccumulation Research Pipeline v17.0

Reproducible pipeline integrating EPA ECOTOX tissue residue data, EPA CompTox chemical properties, and CDC NHANES human biomonitoring data to investigate PFAS bioaccumulation mechanisms across fish and human matrices.

**v17.0 dataset: 24,690 rows (40 ECOTOX tissue residue + 23,532 NHANES human serum) | 13 curated PFAS | Primary paper: Findings 14–20 — mechanistic elimination series for sulfonate BCF | Human serum model: R²=0.658 (compound-identity prediction) | Headline evaluation: LOCO R²=−1.419**

---

## Table of Contents
- [What This Project Shows](#what-this-project-shows)
- [Key Findings](#key-findings)
- [Setup](#setup)
- [Data Sources](#data-sources)
- [Pipeline Architecture](#pipeline-architecture)
- [Model Results](#model-results)
- [Data Gaps](#data-gaps)
- [Version History](#version-history)
- [Roadmap](#roadmap)
- [Citation](#citation)

---

## What This Project Shows

**Human serum:** PFAS blood serum concentration is predictable by compound identity (R²=0.658) but not by chemical structure. The model is a lookup table of per-compound means — LOCO R²=−1.419 confirms chemistry cannot predict an unseen compound's accumulation. Between-compound variance is large and structured; within-compound chemistry signal is zero at current sample sizes.

**Fish BCF mechanistic:** The Arnot-Gobas passive-diffusion framework was exhaustively evaluated across 7 parameter dimensions (Findings 14–20). The single-compartment model accurately predicts carboxylate BCF but fails structurally for sulfonates. Every tunable parameter was exhausted — none closes the sulfonate triad simultaneously. The failure is a structural consequence of the 45× Koc range across PFBS/PFHxS/PFOS: no shared scalar parameter can simultaneously satisfy compounds at the Koc extremes. The blocking data gap is compound-specific fish-tissue NLOM partition coefficients.

**One-line thesis (paper):** Still deciding 

---

## Key Findings

### Finding 14 — Arnot-Gobas model systematically under-predicts all sulfonates
The PFAS-adapted Arnot-Gobas mass balance model under-predicts BCF for all sulfonates (PFOS −91%, PFHxS −85%, PFBS −84% at v12.x baseline) while carboxylates fit better (PFOA +5%, PFDA +2%). The sulfonate under-prediction is systematic and class-wide, pointing to protein-binding mechanisms not captured by the standard Koc-based Kprot term.

### Finding 15 — Human serum model predicts compound means, not chemistry
The dedicated human-only model (NHANES rows, n=23,532) achieves R²=0.658, matching the PFAS-mean baseline exactly. Per-PFAS R² is near zero for every individual compound. The model learns which PFAS it is, not within-compound chemistry. Effective sample size is 6 compounds, not 23,532 rows.

### Finding 16 — Sulfonate Kprot scaling is structurally inadequate; tuning cannot fix it
A sweep of the sulfonate Kprot scale factor across 7 values (0.05–0.30) shows PFOS error moves only 3 percentage points across the full range. The Koc-based protein binding term is wrong in structure, not just magnitude. Simple-constant tuning is a dead end.

### Finding 17 — Gill membrane permeability is not the source of the sulfonate gap
A P_mem sensitivity sweep (8 values, 1.0→0.05) falsifies the hypothesis that reduced gill permeability for ionized PFAS explains the error. PFOS error worsens monotonically as P_mem falls. Together with Finding 16, both tunable rate-constant parameters are exhausted. The error is structural and lives in Kfish.

### Finding 18 — Direct albumin binding is mechanistically correct but does not close the gap
Replacing Kprot = scale × Koc with a direct Ka_albumin term (Bischel et al. 2010, Beesoon & Martin 2015) is the correct mechanistic representation. At v13.0: PFOS −25%, PFBS −68%, PFOA −43%. The remaining error is in k1 (gill uptake), not Kfish (tissue partitioning).

### Finding 19 — Single-compartment framework is structurally closed
The augmented k1 formulation (k1 = k1_passive + K_FAC × Ka_albumin × [albumin]_blood × GV, sulfonates only) was swept across 10 K_FAC values. The sweep reveals step-saturation: error is flat across nine orders of K_FAC after saturation at K_FAC≈0.001. No single K_FAC closes PFOS, PFHxS, and PFBS simultaneously. The single-compartment parameter space is exhausted.

### Finding 20 — Two-compartment model does not improve sulfonates vs v14.1; Koc heterogeneity is the structural constraint
The two-compartment TK model (blood/tissue separation) was implemented and the NLOM tissue factor swept from 0.035 to 0. Results at v17.0 with fish-only observed BCF:

| PFAS | v14.1 (1-comp) | v15.0 (2-comp) | n (fish) | Verdict |
|---|---|---|---|---|
| PFOS | −25% | +9421% | 13 | DEGRADED |
| PFBS | −68% | −99% | 2 | DEGRADED |
| PFOA | −43% | +756% | 110 | DEGRADED |

The two-compartment architecture as currently parameterized does not improve on v14.1. The structural constraint is the 45× Koc range across the sulfonate triad (PFBS=47, PFHxS=560, PFOS=2100): optimal nlom_factor spans three orders of magnitude and cannot be reconciled with a shared scalar. Note: PFOS observed BCF (median=10^−0.004, IQR spanning 4 log units, n=13) is itself an unreliable reference. Both the model failure and the data sparsity are publishable findings. Compound-specific fish-tissue NLOM partition coefficients are the enabling measurement for the follow-up paper.

### Supporting findings

**Finding 3/8 (revised post unit-fix):** Prior fish and plant model results (R²=−0.008, −0.018) were artifacts of mixed-unit contamination. After the v17.0 unit fix, ECOTOX tissue rows drop to 40 total; fish and plant models return 0 rows. The environmental data sparsity problem is worse than previously reported.

**Finding 9 — Leakage audit:** Removing all group-identity features drops pooled R² from 0.710 to 0.490. The missing 0.22 R² was the model learning which measurement protocol was used, not PFAS chemistry.

**Finding 12 — NHANES half-life (cross-sectional only):** Apparent population half-lives from two NHANES waves exceed published clinical values (PFOS 42.3y vs 5.4y). These are cross-sectional estimates from two time points with zero residual degrees of freedom and cannot be separated from exposure trends or cohort effects. They are consistent with ongoing population exposure but are not elimination rate estimates.

---

## Setup

### Requirements
```bash
pip install pandas numpy scikit-learn xgboost matplotlib seaborn openpyxl
```

### Directory structure
```
project/
├── pfas_pipeline_v17.py
├── ecotox_exports/           # ECOTOX xlsx files (accumulation/residue endpoints only)
│   └── *.xlsx
├── nhanes_pfas_processed.csv
├── comptox_snapshot.csv      # optional; falls back to curated PFAS_FEATURES table
└── outputs/                  # created automatically
```

### Run
```bash
mkdir -p outputs
python3 pfas_pipeline_v17.py
```

### ECOTOX export guidance
Download from https://cfpub.epa.gov/ecotox/. For fish tissue residue studies, filter by Effect → Accumulation or Residue, and Result Unit → ng/g, ug/g, mg/kg wet wt (mass/mass tissue only). Water concentration and dose-response exports are automatically routed to separate columns and excluded from the ML target.

---

## Data Sources

| Source | Coverage | Rows (v17.0) |
|---|---|---|
| EPA ECOTOX | Fish/plant/mammal tissue residues + BCF values | 40 tissue + 1,118 medium + 445 BCF |
| CDC NHANES 2015–2016 | Human blood serum, 6 PFAS | ~11,766 |
| CDC NHANES 2017–2018 | Human blood serum, 6 PFAS | ~11,766 |
| EPA CompTox | Chemical properties (13 PFAS) | curated table |

**Note on dataset size:** The v17.0 target variable fix separated tissue residues, water concentrations, administered doses, and soil concentrations — previously mixed into a single column. ECOTOX tissue rows dropped from ~1,500 to 40 after enforcing mass/mass tissue units only. Additional ECOTOX exports filtered to tissue residue endpoints are needed for fish/plant modeling.

---

## Pipeline Architecture

```
EPA ECOTOX exports        CDC NHANES             EPA CompTox
(18 xlsx files)     (nhanes_pfas_processed.csv)  (or curated table)
       |                      |                        |
       v                      |                        |
  harmonize()                 |                        |
  Unit routing:               |                        |
  tissue -> Concentration_ng_g|                        |
  medium -> Medium_Conc_ng_L  |                        |
  dose/soil -> excluded       |                        |
       |                      |                        |
       +----------+-----------+                        |
                  |                                    |
           build_dataset() <--------------------------|
                  |
         24,690 row combined dataset
                  |
     +------------+--------------------+
     v            v                   v
Pooled ML    Human-Only Model   Arnot-Gobas BCF
(LOCO CV)    (NHANES only)      mechanistic series
     |            |                   |
Multi-seed   Calibrated         Kprot sweep (F16)
R2+/-SD      intervals          P_mem sweep (F17)
LOCO R2=     per PFAS           Ka-albumin Kfish (F18)
-1.419       R2=0.658           k1_fac sweep (F19)
                                Two-comp + NLOM (F20)
```

---

## Model Results

### Headline — LOCO R² (honest chemistry evaluation)

Leave-One-Compound-Out CV holds out all rows for one PFAS and trains on the remainder. With 6 compounds and chemistry features constant within a compound, this is the only valid test of structure-to-accumulation prediction.

| PFAS held out | n | LOCO R² | RMSE |
|---|---|---|---|
| PFOA | 3,945 | −0.172 | 0.466 |
| PFUnDA | 3,922 | −0.288 | 0.354 |
| PFDA | 3,922 | −1.633 | 0.565 |
| PFOS | 3,936 | −2.350 | 0.818 |
| PFNA | 3,922 | −5.288 | 0.885 |
| PFHxS | 3,922 | −22.971 | 1.918 |
| **Pooled LOCO** | **23,572** | **−1.419** | **0.982** |

Pooled random-split R²=0.637±0.010 (20 seeds) for comparison — inflated because train/test share identical chemistry feature vectors within the same compound. Effective N = 6 chemistry profiles, not 23,572 rows.

### Human-Only Model (NHANES rows, chemistry features)

| Model | R² | RMSE |
|---|---|---|
| Random Forest | 0.658 | 0.360 |
| XGBoost | 0.658 | 0.360 |
| Linear Regression | 0.649 | 0.365 |
| PFAS-mean Baseline | **0.658** | 0.360 |

RF matches the PFAS-mean baseline exactly. Per-PFAS R²≈0 for all 6 compounds. The model is a compound-identity lookup, not a structure-activity predictor. Calibrated intervals: 80% coverage≈79–83%, width=0.812 log10 ng/g; 95% coverage≈92–96%, width=1.374 log10 ng/g. Intervals describe spread of individuals around the per-compound mean — not uncertainty in the structure-accumulation relationship.

### Arnot-Gobas Mechanistic BCF (v13.0, fish-only observed BCF)

| PFAS | Predicted log BCF | Observed log BCF | n (fish) | IQR | % error |
|---|---|---|---|---|---|
| PFOA | 0.912 | 1.155 | 110 | [0.415, 1.444] | −43% |
| PFOS | 1.877 | −0.004 | 13 | [−1.097, 3.192] | +7502%* |
| PFBS | 0.377 | 2.909 | 2 | [2.654, 3.163] | −100% |
| PFUnDA | 2.178 | 3.465 | 2 | [3.413, 3.517] | −95% |

*PFOS observed BCF IQR spans 4 log units — median is not a stable reference point.

---

## Data Gaps

| PFAS | Fish tissue rows | Human rows | BCF records |
|---|---|---|---|
| PFOA | ~5 | 3,945 | 110 |
| PFOS | ~10 | 3,936 | 13 (IQR unreliable) |
| PFNA | 0 | 3,922 | 0 |
| PFHxS | 0 | 3,922 | 0 — **critical gap** |
| PFDA | 0 | 3,922 | 0 |
| PFUnDA | 0 | 3,922 | 2 |
| PFBS | ~3 | 0 | 2 |
| GenX/ADONA/F53B/PFDoDA/PFHpA/PFHxA | 0 | 0 | 0 |

PFHxS is detected in virtually every American's blood and has zero fish records. Fish and plant models currently return 0 rows after the v17.0 unit fix — additional ECOTOX exports filtered to tissue residue endpoints are the highest-priority data need.

---

## Version History

### v17.0 — July 2026
- **Target variable fix:** `UNIT_TO_NG_G` replaced with three typed dictionaries — tissue units → `Concentration_ng_g` (ML target); medium units → `Medium_Concentration_ng_L` (BCF denominator); dose/soil units → excluded with counts reported. 78 dose rows and 340 soil rows now correctly excluded. ECOTOX tissue rows drop to 40; fish/plant models return 0 rows (prior results were mixed-unit artifacts).
- **BCF restricted to fish-only:** `observed_bcf` filter changed from `isin(["Fish","Other"])` to `== "Fish"`. "Other" is a catch-all for unmatched taxa including invertebrates and algae.
- **Per-PFAS leakage fixed:** Species group dummies removed from `run_per_pfas_models()`. Per-compound R²≈0 is the honest finding.
- **Multi-seed R²:** 20-seed mean±SD replaces single-split R²=0.633. Result: R²=0.637±0.010.
- **Verdict logic fix:** Finding 18/19 verdicts now compare absolute errors with 5pp margin; signed-delta bug caused PFOS (−91%→+7502%) to be labelled IMPROVED.
- **NHANES transparency:** Non-detect count reported; survey weight warning; weighted means if WTSB2YR present; ng/mL≈ng/g assumption stated.
- **Half-life framing corrected:** Cross-sectional caveat block in console; zero residual df stated explicitly.
- **Two-comp circularity acknowledged:** Finding 20 prints calibration caveat — NLOM_TISSUE_FACTOR was fit to these same points; residuals are in-sample.

### v16.0 — July 2026
LOCO CV as headline evaluation; relative paths replacing hardcoded Desktop paths.

### v15.0 — July 2026
Two-compartment TK model; NLOM tissue factor sweep; Finding 20.

### v14.x — July 2026
Class-conditional protein-facilitated k1 (sulfonates only); K_FAC sweep; Finding 19. Bug fixes: Kfood formula, BCF coalescing columns.

### v13.0 — July 2026
Direct Ka_albumin tissue partition coefficient replacing Koc proxy; Finding 18.

### v12.x — July 2026
Kprot sweep (Finding 16); P_mem sweep (Finding 17); class-specific KPROT_SCALE.

### v11.x — July 2026
Koc and AlbuminBinding_pKa features; feature ablation ΔR²=0.000 (Finding 13); Arnot-Gobas BCF (Finding 14); human-only model (Finding 15).

### v10.x — June–July 2026
CASRN salt mapping; leakage fix FIX 4; apparent half-life (Finding 12).

### v7.0–v9.0 — June 2026
XGBoost; leakage fixes FIX 1–3; stratified split; NHANES 2015–2016; per-PFAS models; calibrated intervals.

### v1.0–v6.0 — April–June 2026
Initial pipeline through held-out validation framework.

---

## Roadmap

### Complete ✅
- Full mechanistic elimination series: Findings 14–20
- All leakage fixes (group-identity, duration, species dummies)
- Target variable fix (tissue/medium/dose/soil routing)
- LOCO CV as headline metric
- Multi-seed R² with uncertainty
- BCF fish-only restriction with n and IQR

### In Progress
- [ ] **Manuscript: Findings 14–20** — EST Letters submission
- [ ] Additional ECOTOX exports (fish tissue residue endpoints)
- [ ] Repo cleanup: requirements.txt, sample data, clone-and-run test

### Future
- [ ] NHANES survey weights (WTSB2YR) and non-detect imputation (LOD/√2)
- [ ] Earlier NHANES cycles (1999–2014) for half-life fitting
- [ ] Follow-up paper: compound-specific fish-tissue NLOM partition measurements
- [ ] Streamlit dashboard

---

## Citation

- EPA ECOTOX: https://cfpub.epa.gov/ecotox/
- CDC NHANES: https://wwwn.cdc.gov/nchs/nhanes/
- Arnot & Gobas (2004) Environ. Toxicol. Chem. 23:1523
- Kelly et al. (2004) Environ. Sci. Technol.
- Gobas et al. (2003) Environ. Sci. Technol.
- Guelfo & Higgins (2013) Environ. Sci. Technol.
- Bischel et al. (2010) Environ. Sci. Technol. 44:5770
- Beesoon & Martin (2015) Environ. Sci. Technol. 49:5758
- Ng & Hungerbühler (2013) Environ. Sci. Technol. 47:7214
- Conder et al. (2008) Environ. Sci. Technol. 42:9283
- Barber (2003) Chemosphere 53:1099
- Farrell (1991) J. Exp. Biol. 159:213
- ATSDR Toxicological Profile for Perfluoroalkyls (2021)

---

## Author
PFAS Environmental Informatics Research Project (2026)
