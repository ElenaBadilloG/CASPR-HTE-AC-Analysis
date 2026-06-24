# CASPR — Heterogeneity of Anticoagulation Treatment Response in ESUS

Code for the analysis in:

> **Evaluation of Recurrent Stroke Risk and Anticoagulation Treatment Response among
> Phenotypically Similar Patients with ESUS: Evidence from the CASPR Registry.**

This repository contains the end-to-end R pipeline used to (1) derive data-driven
patient phenotypes from baseline clinical, imaging, and laboratory features, and
(2) estimate **heterogeneous treatment effects (HTE)** of anticoagulation (AC) vs.
antiplatelet-only therapy on a composite outcome of recurrent ischemic stroke,
major bleeding, or death, using a **causal survival forest**.

---

## Overview

`hte-cfs-code-public.R` runs the complete analysis in sequence:

| Stage | What it does | Key packages |
|-------|--------------|--------------|
| 1. Data prep & cohort | Load Stata data, restrict to the analysis cohort (`include == 1`), define treatment and the composite survival outcome, recode factors. | `haven`, `dplyr` |
| 2. Missingness table | Compare included vs. complete-case-excluded patients (Table S). | `gtsummary`, `flextable` |
| 3. Feature scaling | Min–max normalize continuous covariates. | base R |
| 4. Phenotype clustering | Correlation-weighted **Gower distance** → **PAM** clustering; choose *k* by silhouette; selected **k = 3**. | `cluster`, `DescTools` |
| 5. Cluster visualization | Consensus **t-SNE** over 100 runs (most-stable run selected) + bootstrap **Jaccard** cluster-stability analysis. | `Rtsne`, `fpc` |
| 6. Propensity & IPW Cox | Estimate propensity scores, build IPW weights, fit a PS-weighted Cox model with `tx × cluster` interaction. | `WeightIt`, `survival` |
| 7. Causal survival forest | Fit `grf::causal_survival_forest` (horizon = 800 days); estimate ATE and per-patient **CATEs** with CIs; tree diagnostics. | `grf` |
| 8. HTE interrogation | SHAP-style dependence plots, ordered CATE plots by key covariates, **best-linear-projection** HTE tests (Table). | `fastshap`, `patchwork`, `flextable` |
| 9. Survival curves | **Kaplan–Meier** curves by cluster and by cluster × treatment. | `survminer` |
| 10. Summary tables | Cluster-characteristics table; PFO / HFrEF subgroup table; back-transformed (un-scaled) descriptive stats. | `gtsummary`, `flextable` |

Methodology for the causal survival forest follows Sverdrup & Wager (2024).

---

## Requirements

- **R ≥ 4.2** (developed with the `tidyverse` and `grf` ecosystems).

Install the required packages:

```r
install.packages(c(
  "haven", "tidyverse", "twang", "MatchIt", "cobalt", "bcf",
  "survival", "gtsummary", "WeightIt", "grf", "fpc", "DescTools",
  "cluster", "survminer", "ggthemes", "Rtsne", "fastshap",
  "patchwork", "flextable"
))
```

---

## Data

> **The CASPR registry data are not included in this repository** and are not publicly
> available due to patient-privacy restrictions. Access requires appropriate
> data-use agreements with the registry custodians.

The script expects two Stata (`.dta`) files:

| File | Used for |
|------|----------|
| `data.dta` | Main analysis dataset (loaded at the top of the script). |
| `data/CASPR_lvef_092424_FINAL.dta` | Reference dataset used at the end to back-transform (un-scale) covariates for descriptive tables. |

Expected key columns include: an `include` cohort flag; the treatment indicator
`tx_ap_no_ac` (1 = antiplatelet, no AC); the survival outcome `time_to_event` and
event indicator `stroke_bleed_death2`; and the 21 baseline clustering covariates
(`lv_injury`, `mrs_pta`, `age`, `sex`, `race`, `hispanic`, `insurance`, `hx_*`
history flags, `lae_severe`, `pfo`, `lvef`, `wma`, `mri_dwi_multifocal`,
`dwi_cortical`, `bnihss`).

**Treatment coding:** `tx = 1` → anticoagulation; `tx = 0` → antiplatelet, no anticoagulation.

---

## Running the analysis

1. Place `data.dta` in the repository root and the LVEF reference file under `data/`.
2. Create an output folder for figures:
   ```sh
   mkdir -p figures
   ```
3. Run the script from the repository root:
   ```sh
   Rscript hte-cfs-code-public.R
   ```
   or open it in RStudio and run top to bottom. A fixed seed (`set.seed(427)`)
   is used for reproducibility of clustering, t-SNE, and the forest.

---

## Outputs

Figures and tables are written to the working directory (and `figures/`):

**Tables (`.docx`)**
- `t1_missing_comparison.docx` — included vs. excluded (complete-case) comparison
- `Cluster_Characteristics_ISC.docx` — cluster characteristics
- `Table_BLP_HTE_Analysis.docx` — best-linear-projection HTE tests
- `t1-rr.docx` — cluster-based sample summary
- `t1_pfo_ef_comparison.docx` — PFO / HFrEF subgroup comparison

**Figures (`.pdf`)**
- `fS1.pdf` — silhouette scores by *k*
- `figures/f2.pdf` — consensus t-SNE of patient clusters
- `figures/f4.pdf` — CATE / heterogeneity plots
- `shap_values.pdf` — SHAP dependence plots
- `2_…`–`6_cate_estims_ordered-*.pdf` — CATEs ordered by LV injury, baseline NIHSS, age, Fazekas, LVEF
- `kaplan_meier_by_cluster_isc.pdf`, `kaplan_meier_by_cluster_and_tx_isc.pdf`,
  `kaplan_meier_by_cluster_faceted_isc.pdf` — survival curves
- `fS3.pdf` — composite-outcome component breakdown
- `cluster_stab_jaccard.pdf` — bootstrap cluster-stability (Jaccard)

---

## Notes & caveats

- This is **observational** data, so treatment effects are adjusted via estimated
  propensity scores (IPW for the Cox model; `W.hat` learned internally by the
  forest). For an RCT, set `W.hat = mean(W)` in `causal_survival_forest`.
- The script is written as a sequential analysis pipeline rather than a package;
  some output file paths (e.g. the LVEF reference `.dta`) are hard-coded and may
  need to be adjusted for your environment.

## Citation

If you use this code, please cite the paper above.

## License

See [LICENSE](LICENSE).
