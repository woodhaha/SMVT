# SLC5A6 (SMVT) COAD-Specific Survival Analysis

**Generated**: 2026-06-23 16:37:56
**Gene**: SLC5A6 (SMVT, Sodium-dependent Multivitamin Transporter), Entrez ID: 8884
**Cancer**: Colon Adenocarcinoma (COAD / CRC)
**Data Source**: Simulated COAD survival data, calibrated to:
    - SMVT Log2FC = +1.51 in COAD (from existing expression analysis)
    - COAD median OS ~60 months (TCGA)
    - MSI-H frequency ~17% (TCGA COAD Nature 2012)
    - Effect HR ~1.14 (from pan-cancer survival pipeline)
    - Event rate ~40% at 5-year follow-up
**Method**: Kaplan-Meier (log-rank) + Cox proportional hazards
**Stratification**: Median split of SMVT mRNA expression
**Endpoints**: Overall Survival (OS) + Disease-Free Survival (DFS)
**Covariates**: Age, Stage, Gender, MSI Status

## 1. Summary

| Metric | Value |
|--------|-------|
| N patients | 450 |
| OS events | 303 (67.3%) |
| DFS events | 360 (80.0%) |
| Median SMVT exp. | 14.70 |
| Log2FC vs normal | 1.75 |
| MSI-H | 76 / 450 (16.9%) |

## 2. Overall Survival

### Kaplan-Meier

| Metric | High SMVT | Low SMVT |
|--------|-----------|----------|
| Median OS (mo) | 58.7 [nan-nan] | 70.1 [nan-nan] |
| 1Y survival | N/A | N/A |
| 3Y survival | N/A | N/A |
| 5Y survival | N/A | N/A |

0.4227

### Univariate Cox

| Variable | HR | 95% CI | Z | P |
|----------|-----|--------|-----|-----|
| SMVT (cont.) | 1.013 | [1.002-1.024] | 2.23 | 2.57e-02 |

## 3. Disease-Free Survival

| Metric | High SMVT | Low SMVT |
|--------|-----------|----------|
| Median DFS (mo) | 34.8 | 48.9 |
| 1Y DFS | N/A | N/A |
| 3Y DFS | N/A | N/A |
| 5Y DFS | N/A | N/A |

0.0193

## 4. Multivariate Cox (SMVT + Age + Stage + Gender + MSI)

| Variable | HR | 95% CI | P |
|----------|-----|--------|-----|
| SMVT expr | 1.013 | [1.001-1.024] | 0.0304 * |
| AGE | 1.004 | [0.993-1.015] | 0.5008 |
| STAGE | 0.967 | [0.853-1.096] | 0.5968 |
| GENDER (Male) | 1.004 | [0.800-1.260] | 0.9734 |
| MSI-H | 0.767 | [0.557-1.058] | 0.1057 |

**C-index**: 0.547, **N**: 450

## 5. Subgroup Analysis

| Subgroup | N | HR | 95% CI | P |
|----------|-----|-----|--------|-----|
| Stage_I-II | 220 | 1.012 | [0.996-1.028] | 0.1457 |
| Stage_III-IV | 230 | 1.014 | [0.998-1.031] | 0.0806 |
| MSI_MSI-H | 76 | 1.036 | [1.005-1.069] | 0.0234 * |
| MSI_MSS | 374 | 1.009 | [0.997-1.022] | 0.1485 |

## 6. Comparison with Pan-Cancer Results

Pan-cancer COAD (simulated, n=450):
- Univariate HR = 1.137 [1.057-1.224], P = 5.55e-04
- Multivariate HR = 1.133 [1.053-1.218], P = 8.05e-04
- Log-rank P = 0.0675

COAD deep-dive (this analysis, with MSI adjustment):
- Univariate HR = 1.013 [1.002-1.024], P = 2.57e-02

## 7. Methods

### Data Source
Simulated COAD survival data, calibrated to:
    - SMVT Log2FC = +1.51 in COAD (from existing expression analysis)
    - COAD median OS ~60 months (TCGA)
    - MSI-H frequency ~17% (TCGA COAD Nature 2012)
    - Effect HR ~1.14 (from pan-cancer survival pipeline)
    - Event rate ~40% at 5-year follow-up

### Simulation Calibration
- SMVT expr: Log-normal with FC calibrated to Log2FC = +1.51 (COAD vs normal)
- Survival: Weibull rho=1.15, median OS ~60mo
- MSI freq: 17%, MSI-H HR ~0.65
- Age: N(66,11), Stage: I 18%, II 32%, III 35%, IV 15%

### Output Files

| File | Path |
|------|------|
| Results CSV | `outputs/coad_survival_results.csv` |
| KM OS | `figures/KM_COAD_SMVT_survival.png` |
| KM DFS | `figures/KM_COAD_SMVT_disease_free.png` |
| Subgroup | `figures/Forest_COAD_SMVT_subgroup.png` |
| Report | `outputs/coad_survival_report.md` |