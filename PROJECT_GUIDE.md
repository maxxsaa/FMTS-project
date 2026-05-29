# Porto Property Prices — Simple Project Guide

This document explains the FMTS forecasting project in plain language.  
For the formal submission, see `report/report.pdf`.

---

## What is this project?

You must **forecast future monthly property prices in Porto** (euros per square metre) using three types of methods taught in the course:

1. **Smoothing** — Holt-Winters and ETS  
2. **Decomposition** — split the series into trend, season, and noise; forecast with STLF  
3. **Statistical models** — ARIMA / SARIMA with formal diagnostics  

Then you **compare** all methods on a test period the model has never seen, and produce **24-month ahead forecasts** with **95% prediction intervals**.

---

## The data

| Item | Detail |
|------|--------|
| **Variable** | Median dwelling price in Porto |
| **Unit** | EUR/m² |
| **Frequency** | Monthly |
| **Period** | January 2011 – February 2026 (182 months) |
| **File** | `timeSeriesPorto.csv` |
| **Source** | Portuguese housing price statistics (series code 11A1312: Porto) |

**Main patterns:**
- Prices were flat (~1,000 EUR/m²) until about 2015.  
- Strong upward trend since 2016, faster after 2020.  
- Small but regular seasonal pattern (roughly ±15 EUR/m²).  
- A sharp spike around 2022–2023, then a drop and recovery.

---

## How the analysis works (step by step)

### Step 0 — Run the code

Everything is in `timeseries_Rscript.R`. It reads the CSV, cleans it, and saves **15 figures** and **2 tables** in `outputs/`.

```bash
Rscript timeseries_Rscript.R
```

### Step 1 — Explore the series

**Goal:** Understand trend, seasonality, and whether the series is stationary.

**Outputs:**
- `01_time_plot.png` — the raw series over time  
- `02_acf.png` — autocorrelation (slow decay + spikes at 12, 24 → trend + yearly season)  
- `03_pacf.png` — partial autocorrelation (helps choose ARIMA orders)

**Simple reading:** The line goes up a lot over time. You cannot model it in levels without differencing.

### Step 2 — Train / test split

- **Training:** January 2011 – February 2024 (158 months)  
- **Test:** March 2024 – February 2026 (24 months)  

Models are fit **only on training data**, then checked against the test period. This is fair comparison.

### Step 3 — Smoothing methods

#### Holt-Winters (additive seasonality)

Updates three things each month:
- **Level** (current price)  
- **Trend** (direction)  
- **Seasonal indices** (12 monthly adjustments)

**Key parameters (training sample):**
- α = 0.859, β = 0.022, γ = 1.000  

**Output:** `04_forecast_hw_vs_test.png` — blue = forecast, red = actual test, grey = 95% interval.

#### ETS (automatic exponential smoothing)

R picks the best Error / Trend / Seasonal structure by AIC.

**Selected:** ETS(A,A,N) — additive errors, additive trend, **no seasonality** (seasonal effect was too small for ETS to keep on the training window).

**Output:** `05_forecast_ets_vs_test.png`

### Step 4 — Decomposition

**Goal:** Write  
`Observed = Trend + Seasonal + Remainder`

#### STL (robust, preferred)

**Outputs:**
- `06_stl_decomposition.png` — four panels: data, seasonal, trend, remainder  
- `07_seasonally_adjusted.png` — original minus seasonal component  

**Simple reading:** Trend explains almost everything. Seasonal wiggles are small (~±15). Remainder shows big shocks in 2022–2023.

#### Classical additive decomposition

**Output:** `09_classical_decomposition_additive.png` — required by the assignment; same idea as STL but less robust to outliers.

#### STLF forecast

Decompose → forecast trend with ETS → add season back.

**Output:** `08_forecast_stlf_vs_test.png`

### Step 5 — ARIMA / SARIMA

**Box–Cox check:** λ ≈ −0.008 (near log transform). We kept **levels** for comparability with other methods.

**Model search:** `auto.arima()` with seasonal search, plus manual neighbours compared by AIC.

**Selected model:** **ARIMA(0,1,3) with drift**
- One difference removes trend  
- MA(3) on the differenced series  
- Drift ≈ 8.4 EUR/m² per month on the differenced scale  
- AIC = 1576.6  

No seasonal ARIMA terms were kept — seasonality is small relative to trend on the training sample.

**Diagnostics:**
- `11_sarima_residual_diagnostics.png` — residuals should look like random noise  
- **Ljung-Box test:** p = 0.28 (lag 24, df = 20 after subtracting 4 estimated parameters) → residuals look fine  

**Forecast output:** `10_forecast_sarima_vs_test.png`

### Step 6 — Compare methods on the test set

**Table:** `outputs/test_accuracy_by_method.csv`

| Method | RMSE | MAPE | Meaning |
|--------|------|------|---------|
| Holt-Winters | 302 | 8.8% | **Best** (lowest error) |
| ETS | 326 | 9.6% | Second |
| STLF | 336 | 10.1% | Third |
| SARIMA | 350 | 10.5% | Fourth |

**Important:** All methods **under-forecast** (actual prices rose faster than history suggested). Holt-Winters wins, but every method struggled with the 2024–2026 acceleration.

Positive **ME** (~258–306) = forecasts systematically too low.

### Step 7 — Out-of-sample forecasts

Models are **refit on all data**, then forecast the next **24 months** (March 2026 – February 2028).

**Outputs:**
- `12_oos_holtwinters.png` … `15_oos_sarima.png`  
- `outputs/out_of_sample_forecasts_95.csv` — point forecasts and 95% bounds for all methods  

Near-term forecasts cluster around **3,170 EUR/m²** for March 2026, rising toward **~3,350** by early 2028. Intervals widen further ahead.

---

## All output files (quick reference)

| File | What it is |
|------|------------|
| `01_time_plot.png` | Raw series |
| `02_acf.png`, `03_pacf.png` | Exploratory correlograms |
| `04`–`05` | Smoothing vs test |
| `06`–`07`, `09` | Decomposition & seasonally adjusted |
| `08` | STLF vs test |
| `10`–`11` | SARIMA forecast & diagnostics |
| `12`–`15` | Out-of-sample forecasts |
| `test_accuracy_by_method.csv` | Accuracy comparison |
| `out_of_sample_forecasts_95.csv` | Future forecasts table |

---

## How to build the PDF report

```bash
cd report
make
```

This exports LaTeX tables from R, then compiles `report.pdf`.

Before submitting, edit `report/report.tex` and fill in your **full name** and **student number** on the title page and declaration.

---

## What to say in the oral discussion

1. **Why Holt-Winters won:** It explicitly models seasonality; ETS dropped it; all methods missed the post-2024 speed-up.  
2. **Why SARIMA is still useful:** Formal diagnostics (Ljung-Box, residual ACF) even if test RMSE is higher.  
3. **Limitations:** Models assume the past repeats; none capture structural breaks well; intervals assume Gaussian errors.  
4. **Seasonality:** Small in euros but visible; STL shows it clearly.  
5. **No ADF/PP/KPSS on seasonal data** — as per course instructions; we used ACF/PACF and differencing instead.

---

## File map

```
Project/
├── timeSeriesPorto.csv      # Data
├── timeseries_Rscript.R     # Full analysis
├── outputs/                 # Figures and CSV tables
├── PROJECT_GUIDE.md         # This file
└── report/
    ├── report.tex           # LaTeX source
    ├── report.pdf           # Submission PDF (after make)
    └── Makefile             # Build script
```
