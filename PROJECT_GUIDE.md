# Porto Property Prices — Project Guide

Unified FMTS forecasting pipeline (merged methodology from `main` + `improvedVersion`).

## Quick start

```bash
Rscript run_analysis.R          # full analysis → outputs/
cd report && make               # build report.pdf
```

Works on any machine: paths are resolved relative to the project root (no hardcoded user paths).

## Methodology (theory-aligned)

| Step | Approach |
|------|----------|
| Decomposition | **Multiplicative** ($y = T \times S \times R$); STL on $\log(y)$ |
| Seasonally adjusted | $y / S$ |
| Smoothing | Holt-Winters **multiplicative** (+ damped variant), ETS, STLF |
| ARIMA | **Log transform** + **first difference**; short-list by AIC; Ljung-Box with correct df |
| Evaluation | 24-month hold-out test + 24-month OOS forecasts with 95% PI |

Unit-root tests (ADF/PP/KPSS) are **not** used on this seasonal series (per course guidelines).

## Pipeline structure

```
run_analysis.R              # orchestrator
R/utils.R                     # portable paths, helpers
01_EDA.R                      # exploration
02_decomposition_smoothing.R  # multiplicative decomp + smoothing
03_arima.R                    # log-ARIMA + diagnostics
04_comparison_oos.R           # test comparison + OOS
```

## Outputs

| Folder | Contents |
|--------|----------|
| `outputs/figures/` | 22 PNG figures |
| `outputs/tables/` | accuracy, forecast vs test, OOS CSV |
| `outputs/tests/` | model summaries, Ljung-Box |
| `report/report.pdf` | submission report |

## Repository

[github.com/maxxsaa/FMTS-project](https://github.com/maxxsaa/FMTS-project) — branch `main` (unified); `improvedVersion` (Pedro's original split scripts).

## Before submitting

1. Fill in student numbers in `report/report.tex`
2. Sign the Declaration of Authenticity
3. Run `Rscript run_analysis.R && cd report && make` to refresh all outputs
