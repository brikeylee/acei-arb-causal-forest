# Individualized Perioperative ACEI/ARB Management in Cardiac Surgery

A Causal Forest Analysis for Heterogeneous Treatment Effect Estimation

## Overview

This study uses Causal Forest methods to analyze perioperative ACEI/ARB management strategies in cardiac surgery patients, developing individualized treatment decision support tools.

**Target Journal:** Annals of Intensive Care

## Project Structure

```
CF/
├── README.md
├── baseline/                               # Raw data
│   └── ACE_baseline250309.xlsx
├── src/                                    # Source code
│   ├── main.R                              # Main pipeline entry point
│   ├── core/                               # Core modules
│   │   ├── setup.R                         # Package management & themes
│   │   ├── functions.R                     # Helper functions
│   │   ├── data_prep.R                     # Data preprocessing (MI m=5)
│   │   ├── matching.R                      # Propensity score matching
│   │   └── causal_forest.R                 # Causal forest (multi-outcome)
│   └── analysis/                           # Analysis modules
│       ├── validation_analysis.R           # Bootstrap, RATE, hyperparameters
│       ├── sensitivity_analysis.R          # IPTW, E-value, caliper, interactions
│       ├── treatment_heterogeneity.R       # Heterogeneity testing
│       ├── individualized_treatment.R      # Clustering, recommendations
│       └── clinical_decision_support.R     # Score system, calibration, NNT
├── results/                                # Analysis outputs
│   ├── outcome_summary.csv
│   ├── main_analysis/
│   │   ├── matching/                       # Baseline tables, balance plots
│   │   ├── causal_forest/                  # ATE, CATE, subgroups, figures
│   │   ├── models/                         # Saved CF model (.rds)
│   │   └── benefit_risk/                   # NNT/NNH analysis
│   ├── enhanced_analysis/
│   │   ├── individualized_treatment/       # Clusters, recommendations
│   │   ├── heterogeneity_analysis/         # Heterogeneity metrics
│   │   ├── clinical_decision_support/      # Score validation, calibration
│   │   ├── validation_analysis/            # Bootstrap, stability
│   │   └── visualizations/                 # Publication-quality figures
│   └── sensitivity_analysis/               # IPTW, E-value, caliper, interactions
└── docs/                                   # Documentation & manuscript
    ├── STROBE_checklist.md
    ├── METHODS.md
    └── manuscripts/
```

## Quick Start

```r
setwd("/path/to/CF")
source("src/main.R")
all_results <- run_full_analysis()
```

## Key Methodological Features

- **Multiple imputation** (m=5, PMM) with pooled estimates
- **Propensity score matching** (1:1 nearest-neighbor, caliper=0.1)
- **Causal Forest** (grf, 2000 trees, honest splitting, auto-tuning)
- **Multi-outcome analysis** (composite + death, AKI stages)
- **Bootstrap validation** (500 iterations)
- **RATE** (Rank-weighted Average Treatment Effect) for heterogeneity
- **Sensitivity analyses**: IPTW, E-value, caliper sensitivity, subgroup interactions (FDR-corrected)
- **Clinical decision score** with 5-fold cross-validation
- **STROBE-compliant** reporting
