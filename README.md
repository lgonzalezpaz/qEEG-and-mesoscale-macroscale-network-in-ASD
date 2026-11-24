# Electroencephalographic criticality and mesoscale–macroscale network decoupling distinguish autism spectrum disorder from neurological disorders: a multi-cohort study

**Open-Science Computational & Statistical Reproducibility Archive**

> **Note on Manuscript Availability:** The accompanying research manuscript is currently undergoing peer review and is intentionally excluded from this public repository and its Zenodo snapshot. Only analytical code, derived statistical summaries, metadata provenance, summary tables, and publication figures are hosted here.

## Overview

This repository provides the refined `R` and `Python` analytical pipelines, statistical audits, empirical Bayes harmonization checks, and publication figure generation scripts supporting our multi-site investigation across **4,834 EEG recordings** (**4,609 valid recordings** spanning **2,569 valid / 2,576 total unique biological participants** and 26,096 ROI observations) harmonized from eleven open-access neurodevelopmental, normative, and clinical neurology repositories (`complete_metadata.xlsx` inside `complete_metadata.zip`).

### Key Scientific & Methodological Findings

1. **Supercritical Hyper-Persistence in ASD:** Frontal cortex (`n = 687`) operates overwhelmingly in a hyper-persistent regime (`H = 0.8805 ± 0.0491`, `95.78% > 0.80`, only `0.44%` within the physiological self-organized criticality [SOC] band `[0.70, 0.80]`), paired with suppressed spectral Shannon entropy (`S = 2.8677 ± 1.2752` vs. TD `4.6409 ± 1.1121`).
2. **Lateral Inferior Frontal Hub Trapping & Modular Collapse:** `88.03%` (`618/702`) of ASD participants trap peak shortest-path betweenness centrality at lateral inferior frontal sensors `F7/F8` ($\chi^2_2 = 1384.16$, $P = 2.71 \times 10^{-301}$), accompanied by collapsed Louvain modularity (`Q = 0.1752 ± 0.1454` vs. TD `0.4142 ± 0.1958`).
3. **Mesoscale–Macroscale Decoupling:** ASD severs the physiological inverse association between frontal Shannon entropy and global efficiency (`r = -0.0409` vs. TD `r = -0.4912`, Fisher `z = 11.71`, $P = 1.18 \times 10^{-31}$; non-ASD `r = -0.6950`, `z = 16.29`), replicating within the balanced single-site Mendeley cohort (`r = +0.0218` vs. TD `-0.4401`, `z = 8.46`, $P = 2.70 \times 10^{-17}$) and in single-condition `restEO` (`n = 31` ASD vs. `31` TD unique participants: `r = +0.2806` vs. `-0.5696`, `z = 3.50`, $P = 4.66 \times 10^{-4}$).
4. **Algorithmic Audits & Leakage-Free Validation:** The refined pipeline eliminates `NA_ratio.at_roi` missing-`O2` artifacts, fixes rank-deficient structural-zero `Group × ROI` interactions (`IFG`, `BG/SN/Th`, `Amygdala` absent in `non-ASD`), validates ComBat batch correction (`ICC_site` reduced `90.4%–99.1%` for bounded features vs. inflation on un-logged power), and confirms `94.58%` hold-out accuracy (macro OVR AUC `0.9945`; ASD vs. non-ASD AUC `0.9986`) using strictly artifact-free universal ROI predictors.
5. **Two-Tier Biological-Participant Sensitivity Verification:** Explicitly mapping recording IDs (`Metadata_ID`) to source-repository subject identifiers (`Subject_ID` in `complete_metadata.xlsx`, `N_bio = 2,569` valid unique participants: `144` ASD, `2,198` TD, `227` non-ASD) preserves all findings under biological-participant-clustered REML LME (`N_obs = 26,096` nested in `N_bio = 2,569`, Group $\chi^2_2 = 136.22$ for Shannon and `1,872.59` for Hurst), participant-averaged frontal decoupling (`N_bio = 2,554`, Fisher `z = 4.09`, $P = 4.32 \times 10^{-5}$ vs. TD), and biological-participant-disjoint Random Forest hold-out evaluation (`93.38%` accuracy across `515` unseen biological participants, ASD vs. non-ASD AUC `0.9995`).

---

## Repository Structure

```text
qEEG-and-mesoscale-macroscale-network-in-ASD/
├── R.zip                                   # Extracted to R/: refined R scripts (01–05) + changelog
│   ├── 01_build_master_database.R          # Deterministic cohort assembly & NA/nan hygiene
│   ├── 02_harmonization_and_roi_aggregation.R # Empirical Bayes ComBat & missingness audit
│   ├── 03_unsupervised_pca_hcpc.R          # 65-variable PCA & 7-cluster HCPC phenotyping
│   ├── 04_multilevel_lme_and_decoupling.R  # Rank-corrected REML LME & Fisher z decoupling
│   ├── 05_supervised_random_forest_validation.R # Leakage-free balanced RF & ablation
│   └── REFINEMENT_CHANGELOG.md             # Line-by-line audit vs. legacy Session001-006b scripts
├── python.zip                              # Extracted to python/: statistical audit & figure scripts
│   ├── 01_statistical_audit_and_modeling.py # Full Python reproduction of R workspace (.RData)
│   └── 02_generate_publication_figures.py  # Colorblind-safe vector PDF & 600-DPI PNG generator
├── data_summary.zip                        # Extracted to data_summary/: derived non-PII summaries & ROC JSON
├── complete_metadata.zip                   # Three-sheet metadata provenance workbook (complete_metadata.xlsx)
├── Tables.zip                              # Extracted to Tables/: publication-ready summary CSV tables
├── Figures.zip                             # Extracted to Figures/: vector PDF & 600-DPI PNG figures (Fig1–Fig6)
├── environment.yml                         # Reproducible Conda environment specification
├── requirements.txt                        # Exact pinned Python package versions
├── CITATION.cff                            # Machine-readable software citation metadata
├── .zenodo.json                            # Zenodo archive & dataset citation metadata
├── LICENSE                                 # MIT License
└── README.md                               # This guide
```

---

## Reproducing the Analysis & Figures

### Python Pipeline (Standalone Verification)

Extract the zipped folders into the repository root (`R/`, `python/`, `data_summary/`, `Tables/`, `Figures/`), then install dependencies via pip (`pip install -r requirements.txt`) or Conda (`conda env create -f environment.yml`):

```bash
python python/01_statistical_audit_and_modeling.py
python python/02_generate_publication_figures.py
```

- `01_statistical_audit_and_modeling.py` reads the original read-only R workspace (`.RData`) alongside the canonical three-sheet metadata registry (`complete_metadata.xlsx` / `complete_metadata.zip`) and regenerates all derived summary tables, `Metadata_Complete_Recording_to_Participant_Map.csv`, `Biological_Participant_Frontal_Summary.csv`, ComBat audits, record- and biological-participant-clustered REML LME Wald $\chi^2$ statistics, Fisher *r*-to-*z* decoupling comparisons (including single-condition `restEO` `1:1` sensitivity), and 8-arm Random Forest ablation & participant-disjoint ROC benchmarks into `data_summary/` and `Tables/`.
- `02_generate_publication_figures.py` produces all six multi-panel figures (`Fig1`–`Fig6`) in vector `.pdf` and 600 DPI `.png` formats inside `Figures/`.

### Refined R Pipeline (`R/01`–`R/05`)

The modernized, artifact-free R scripts inside `R/` replace interactive `.RData` dependencies (`Session001`–`Session006b`) with explicit CSV/RDS checkpoints, automatic structural-zero column pruning in `lme4`/`car::Anova`, log-transform pre-conditioning before `sva::ComBat`, and stratified 6-arm `caret`/`randomForest` cross-validation (documented in `R/REFINEMENT_CHANGELOG.md`).

---

## Data Provenance

No proprietary clinical data or identifiable patient metadata are hosted in this repository. All underlying raw EEG cohorts remain archived at their primary public repositories: OpenNeuro (`ds006777`, `ds006780`, `ds004504`, `ds004584`, `ds003555`), PhysioNet (`CHB-MIT`, `10.13026/C2K01R`), Mendeley Data (`10.17632/spwnt8t25y.1`), University of Sheffield ORDA (`10.15131/shef.data.16840351`), Figshare (`10.6084/m9.figshare.4810705`), BCIAUT-P300 (`10.3389/fnins.2020.568104`), and NMT Scalp EEG Dataset (`10.3389/fnins.2021.755817`).

## License

Code and derived reproducible tables are released under the MIT License (`LICENSE`).
