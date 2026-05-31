# Porto Property Prices — Project Guide

Unified FMTS forecasting pipeline for monthly median dwelling prices in Porto (EUR/m²).

## Quick start

From the project root:

```bash
Rscript run_analysis.R          # generates outputs/figures, tables, tests
cd report && make               # builds report/report.pdf
```

Requires R packages: `tidyverse`, `forecast`.

Paths are resolved automatically from the project root — no machine-specific configuration.

## Methodology

| Step | Approach |
|------|----------|
| Decomposition | **Multiplicative** ($y = T \times S \times R$); STL on $\log(y)$ |
| Seasonally adjusted | $y / S$ |
| Smoothing | Holt-Winters **multiplicative** (+ damped), ETS, STLF |
| ARIMA | **Log transform** + **first difference**; short-list by AIC; Ljung-Box with correct df |
| Evaluation | 24-month hold-out test + 24-month OOS forecasts with 95% PI |

Unit-root tests (ADF/PP/KPSS) are **not** used on this seasonal series (per course guidelines in `README.MD`).

## Project layout

```
Project/
├── run_analysis.R              # Orchestrator — run this
├── timeseries_Rscript.R        # Legacy alias for run_analysis.R
├── data/
│   └── timeSeriesPorto.csv
├── scripts/
│   ├── 01_eda.R                # Time plot, seasonal boxplot, ACF/PACF
│   ├── 02_decomposition_smoothing.R
│   ├── 03_arima.R
│   └── 04_comparison_oos.R
├── R/
│   └── utils.R                 # Paths, train/test split, helpers
├── outputs/
│   ├── figures/                # 01–23 numbered PNG figures
│   ├── tables/                 # CSV results (incl. STL phase matrices)
│   └── tests/                  # Text summaries (HW, ETS, Ljung-Box)
├── report/
│   ├── report.tex
│   ├── report.pdf
│   └── Makefile
└── docs/
    ├── declaration-of-authencity.txt
    ├── fep-pedagogical-report-dissertation.txt
    ├── dataset_description.pdf
    └── Group_08_draft.pdf
```

## Pipeline steps

### 1 — EDA (`scripts/01_eda.R`)

Explores trend, seasonality, and autocorrelation. Produces figures `01`–`04`.

### 2 — Decomposition & smoothing (`scripts/02_decomposition_smoothing.R`)

- Classical multiplicative decomposition
- STL on log prices
- Seasonally adjusted series → `outputs/tables/seasonally_adjusted.csv`
- STL phase matrices (Year × Month + volatility by phase) → `outputs/figures/23_stl_phase_matrices.png`
- Holt-Winters (multiplicative + damped), ETS, STLF fitted on training data
- Figures `05`–`12`

### 3 — ARIMA (`scripts/03_arima.R`)

- Log transform and first difference of training series
- Three candidate models compared by AIC → `outputs/tables/arima_aic_comparison.csv`
- Residual diagnostics and Ljung-Box test → `outputs/tests/ljung_box_arima.txt`
- Figures `13`–`18`

### 4 — Comparison & OOS (`scripts/04_comparison_oos.R`)

- Test-set accuracy for all methods → `outputs/tables/test_accuracy_by_method.csv`
- Forecast vs actual table → `outputs/tables/forecast_vs_test.csv`
- 24-month out-of-sample forecasts → `outputs/tables/out_of_sample_forecasts_95.csv`
- Figures `19`–`22`

## Key results (24-month test set)

| Method | RMSE | MAPE |
|--------|------|------|
| ETS | 326 | 9.6% |
| HW (multiplicative) | 336 | 10.2% |
| STLF | 336 | 10.1% |
| ARIMA (log) | 583 | 18.3% |

All methods under-forecast the accelerated price growth in 2024–2026.

## Before submitting

1. Sign the separate declaration document (`docs/declaration-of-authencity.txt`)
2. Re-run `Rscript run_analysis.R && cd report && make` to refresh outputs

## Repository

[github.com/maxxsaa/FMTS-project](https://github.com/maxxsaa/FMTS-project) — branch `main`.
