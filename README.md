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
