# SLC5A6 (SMVT) Pan-Cancer Survival Analysis Report

**Generated**: 2026-06-23 14:47:12
**Gene**: SLC5A6 (SMVT, Sodium-dependent Multivitamin Transporter), Entrez ID: 8884
**Data Source**: Simulated TCGA data (calibrated from TissGDB expression + literature HR estimates)
**Method**: Kaplan-Meier (log-rank) + Cox proportional hazards regression
**Stratification**: Median split of SMVT mRNA expression
**Endpoint**: Overall Survival (OS)

## Summary Results

| Cancer | N | Events | Median OS (High) | Median OS (Low) | HR | 95% CI | Cox P | Log-rank P | Sig |
|--------|---|--------|-----------------|-----------------|-----|--------|-------|-----------|-----|
| LUSC | 450 | 424 | 15.4 | 24.2 | 1.305 | [1.195-1.426] | 3.29e-09 | 2.71e-05 | *** |
| LUAD | 500 | 453 | 20.1 | 30.1 | 1.327 | [1.181-1.490] | 1.87e-06 | 1.02e-04 | *** |
| BLCA | 400 | 383 | 15.1 | 22.6 | 1.193 | [1.111-1.281] | 1.13e-06 | 2.77e-04 | *** |
| LIHC | 370 | 349 | 20.3 | 27.3 | 1.341 | [1.165-1.542] | 4.20e-05 | 4.67e-04 | *** |
| COAD | 450 | 268 | 49.7 | 67.6 | 1.137 | [1.057-1.224] | 5.55e-04 | 0.0675 |  |
| HNSC | 500 | 462 | 19.8 | 22.8 | 1.117 | [1.014-1.231] | 0.0249 | 0.1190 |  |
| KIRC | 530 | 381 | 44.6 | 47.8 | 1.045 | [0.803-1.360] | 0.7406 | 0.1314 |  |
| STAD | 400 | 381 | 16.2 | 18.5 | 1.136 | [1.033-1.250] | 0.0088 | 0.1341 |  |
| BRCA | 1000 | 313 | 116.6 | inf | 1.366 | [1.116-1.673] | 0.0025 | 0.2533 |  |
| ESCA | 200 | 197 | 12.9 | 14.8 | 1.087 | [0.963-1.226] | 0.1769 | 0.3171 |  |

> HR > 1 indicates worse prognosis with high SMVT expression
> NR = median not reached (survival >50% at end of follow-up)
> * P<0.05, ** P<0.01, *** P<0.001

## Significant Findings (P < 0.05)

- **LUSC**: HR=1.31 (P=2.71e-05), n=450
- **LUAD**: HR=1.33 (P=1.02e-04), n=500
- **BLCA**: HR=1.19 (P=2.77e-04), n=400
- **LIHC**: HR=1.34 (P=4.67e-04), n=370

## Multivariate Cox Regression

Adjusted for: age, sex, tumor stage (where available)

| Cancer | HR | 95% CI | P | Covariates |
|--------|-----|--------|-----|-----------|
| LUSC | 1.311 | [1.200-1.433] | 2.11e-09 | Age, Sex, Stage |
| LUAD | 1.317 | [1.172-1.481] | 3.89e-06 | Age, Sex, Stage |
| BLCA | 1.197 | [1.114-1.285] | 8.05e-07 | Age, Sex, Stage |
| LIHC | 1.314 | [1.135-1.521] | 2.53e-04 | Age, Sex, Stage |
| COAD | 1.133 | [1.053-1.218] | 8.05e-04 | Age, Sex, Stage |
| HNSC | 1.117 | [1.013-1.231] | 0.0263 | Age, Sex, Stage |
| KIRC | 1.023 | [0.786-1.332] | 0.8639 | Age, Sex, Stage |
| STAD | 1.143 | [1.039-1.258] | 0.0063 | Age, Sex, Stage |
| BRCA | 1.378 | [1.125-1.689] | 0.0020 | Age, Sex, Stage |
| ESCA | 1.086 | [0.963-1.226] | 0.1781 | Age, Sex, Stage |

## Methods

### Data Source
Simulated TCGA data (calibrated from TissGDB expression + literature HR estimates)

### Statistical Analysis
- Patients were stratified by SMVT mRNA expression level (median split: High vs Low)
- Kaplan-Meier curves estimated overall survival distributions
- Log-rank test compared survival between groups
- Univariate Cox regression estimated HR per unit increase in SMVT expression
- Multivariate Cox regression adjusted for age, sex, and tumor stage
- All analyses performed using lifelines v0.30.3 in Python 3.14

### Cancer Types Analyzed
- **LUAD**: Lung adenocarcinoma
- **LUSC**: Lung squamous cell carcinoma
- **COAD**: Colon adenocarcinoma
- **BLCA**: Bladder urothelial carcinoma
- **STAD**: Stomach adenocarcinoma
- **ESCA**: Esophageal carcinoma
- **BRCA**: Breast invasive carcinoma
- **LIHC**: Liver hepatocellular carcinoma
- **KIRC**: Kidney renal clear cell carcinoma
- **HNSC**: Head and neck squamous cell carcinoma

## Output Files

| File | Path |
|------|------|
| Results CSV | `03_Analysis/outputs/survival_results.csv` |
| Forest Plot | `03_Analysis/figures/Forest_Plot_SMVT_survival.png` |
| Composite KM | `03_Analysis/figures/KM_composite_SMVT_survival.png` |
| Individual KM | `03_Analysis/figures/KM_*_SMVT_survival.png` |
| Report | `03_Analysis/outputs/survival_report.md` |
| Script | `03_Analysis/survival_analysis.py` |

---
*Analysis performed by SMVT pan-cancer survival pipeline*