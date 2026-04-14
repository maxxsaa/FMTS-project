# =====================================================
# Porto property prices — FMTS full analysis pipeline
# =====================================================

rm(list = ls())

# --- 0. Paths (CSV + figures next to this script) --------------------------------

resolve_project_paths <- function(csv_name = "timeSeriesPorto.csv") {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[1])
    root <- dirname(normalizePath(script_path, mustWork = TRUE))
  } else {
    root <- normalizePath(getwd(), mustWork = TRUE)
  }
  csv_path <- file.path(root, csv_name)
  if (!file.exists(csv_path)) {
    stop(
      "Cannot find ", csv_name, " in ", root, ".\n",
      "Place the CSV next to this script or set working directory to the project folder.",
      call. = FALSE
    )
  }
  out_dir <- file.path(root, "outputs")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  list(root = root, csv_path = normalizePath(csv_path), out_dir = out_dir)
}

paths <- resolve_project_paths()
CSV_PATH <- paths$csv_path
OUT_DIR <- paths$out_dir

# --- 1. Packages -----------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(forecast)
})

# forecast::accuracy() returns Training + Test rows when comparing to a hold-out vector
test_accuracy_row <- function(fc, actual) {
  a <- accuracy(fc, actual)
  if ("Test set" %in% rownames(a)) {
    return(a["Test set", , drop = FALSE])
  }
  if (nrow(a) >= 2) {
    return(a[2, , drop = FALSE])
  }
  a[1, , drop = FALSE]
}

# --- 2. Import & clean -----------------------------------------------------------

raw_data <- read_csv(
  CSV_PATH,
  col_names = "raw_line",
  col_types = "c",
  show_col_types = FALSE
)

porto_prices <- raw_data %>%
  separate(
    raw_line,
    into = c("month_year", "series_type", "region_code", "euros_str", "extra"),
    sep = ";",
    extra = "drop",
    fill = "right"
  ) %>%
  mutate(
    month_name = str_extract(month_year, "^[A-Za-z]+"),
    year = as.numeric(str_extract(month_year, "\\d{4}$")),
    month_number = case_when(
      tolower(month_name) == "january" ~ 1,
      tolower(month_name) == "february" ~ 2,
      tolower(month_name) == "march" ~ 3,
      tolower(month_name) == "april" ~ 4,
      tolower(month_name) == "may" ~ 5,
      tolower(month_name) == "june" ~ 6,
      tolower(month_name) == "july" ~ 7,
      tolower(month_name) == "august" ~ 8,
      tolower(month_name) == "september" ~ 9,
      tolower(month_name) == "october" ~ 10,
      tolower(month_name) == "november" ~ 11,
      tolower(month_name) == "december" ~ 12,
      TRUE ~ NA_real_
    ),
    euros_per_sqm = as.numeric(gsub(",", ".", euros_str)),
    date = as.Date(paste0(year, "-", sprintf("%02d", month_number), "-01"))
  ) %>%
  filter(!is.na(date), !is.na(euros_per_sqm)) %>%
  distinct(date, .keep_all = TRUE) %>%
  select(date, euros_per_sqm) %>%
  arrange(date)

porto_ts <- ts(
  porto_prices$euros_per_sqm,
  start = c(
    as.numeric(format(min(porto_prices$date), "%Y")),
    as.numeric(format(min(porto_prices$date), "%m"))
  ),
  frequency = 12
)

cat("\n=== Data summary ===\n")
print(head(porto_prices, 3))
print(tail(porto_prices, 3))
cat("n =", nrow(porto_prices), "| start =", paste(start(porto_ts), collapse = "/"),
    "| end =", paste(end(porto_ts), collapse = "/"), "\n")

# --- 3. Train / test split (hold-out for model comparison) -----------------------

h_test <- 24L
if (length(porto_ts) <= h_test + 36) {
  h_test <- max(12L, floor(0.15 * length(porto_ts)))
}
n_tot <- length(porto_ts)
t_time <- time(porto_ts)
train <- window(porto_ts, end = t_time[n_tot - h_test])
test <- window(porto_ts, start = t_time[n_tot - h_test + 1])
cat("Train ends:", paste(end(train), collapse = "/"), "| Test h =", h_test, "\n")

# --- 4. Exploratory plots (features of the series) -----------------------------

png(filename = file.path(OUT_DIR, "01_time_plot.png"), width = 900, height = 500, res = 120)
plot(
  porto_ts,
  main = "Porto median dwelling price (EUR/m²) — monthly",
  ylab = "EUR per m²",
  xlab = "Time"
)
dev.off()

png(filename = file.path(OUT_DIR, "02_acf.png"), width = 900, height = 500, res = 120)
acf(porto_ts, main = "ACF — full series")
dev.off()

png(filename = file.path(OUT_DIR, "03_pacf.png"), width = 900, height = 500, res = 120)
pacf(porto_ts, main = "PACF — full series")
dev.off()

# --- 5. Smoothing methods (estimation on training sample) ----------------------

# 5a. Holt-Winters (additive seasonality; state-space equivalent to triple exponential smoothing)
fit_hw <- HoltWinters(train, seasonal = "additive")
cat("\n=== Holt-Winters (additive) — fitted on train ===\n")
print(fit_hw)

fc_hw <- forecast(fit_hw, h = h_test, level = c(95))

# 5b. ETS (automatic error/trend/seasonal structure; reports AIC etc.)
fit_ets <- ets(train)
cat("\n=== ETS (selected on train) ===\n")
print(fit_ets)

fc_ets <- forecast(fit_ets, h = h_test, level = c(95))

png(filename = file.path(OUT_DIR, "04_forecast_hw_vs_test.png"), width = 900, height = 500, res = 120)
plot(fc_hw, main = "Holt-Winters forecast vs test set")
lines(test, col = "red", lwd = 2)
legend(
  "topleft",
  legend = c("Forecast", "95% PI", "Actual (test)"),
  col = c("blue", "lightblue", "red"),
  lty = c(1, 1, 1),
  lwd = c(2, 6, 2),
  bty = "n"
)
dev.off()

png(filename = file.path(OUT_DIR, "05_forecast_ets_vs_test.png"), width = 900, height = 500, res = 120)
plot(fc_ets, main = "ETS forecast vs test set")
lines(test, col = "red", lwd = 2)
dev.off()

# --- 6. Decomposition — trend, seasonal, irregular; seasonally adjusted --------

# STL: robust seasonal-trend decomposition (preferred for monthly data with varying seasonality)
stl_full <- stl(porto_ts, s.window = "periodic", robust = TRUE)
cat("\n=== STL decomposition (full series, s.window = periodic) ===\n")
print(summary(stl_full))

sa_series <- seasadj(stl_full) # y - seasonal component (same index as porto_ts)

png(filename = file.path(OUT_DIR, "06_stl_decomposition.png"), width = 1000, height = 900, res = 120)
plot(stl_full, main = "STL decomposition")
dev.off()

png(filename = file.path(OUT_DIR, "07_seasonally_adjusted.png"), width = 900, height = 500, res = 120)
plot(
  cbind(Original = porto_ts, Seasonally_adjusted = sa_series),
  main = "Original vs STL seasonally adjusted",
  ylab = "EUR per m²"
)
dev.off()

# Forecast from decomposition approach: STLF (STL + ETS on seasonally adjusted + seasonal naive)
fc_stlf <- stlf(train, method = "ets", h = h_test, level = c(95))

png(filename = file.path(OUT_DIR, "08_forecast_stlf_vs_test.png"), width = 900, height = 500, res = 120)
plot(fc_stlf, main = "STLF (decomposition-based) forecast vs test set")
lines(test, col = "red", lwd = 2)
dev.off()

# Classical additive decomposition (for report: show trend / seasonal / residual)
decomp_add <- decompose(porto_ts, type = "additive")
png(filename = file.path(OUT_DIR, "09_classical_decomposition_additive.png"), width = 1000, height = 900, res = 120)
plot(decomp_add)
dev.off()

# --- 7. SARIMA — specification search, diagnostics, test-set performance -------

# Transformations: for this series, check whether a Box–Cox / log stabilizes variance.
# (Course note: seasonal unit-root tests such as ADF/PP/KPSS are not the focus here.)
lambda_guerrero <- BoxCox.lambda(train, method = "guerrero")
cat(
  "\nBox–Cox lambda (train, Guerrero) =", round(lambda_guerrero, 4),
  "— near 1.0 suggests levels or mild transform; far from 1 may warrant BoxCox(y, lambda) in modeling.\n"
)

# Short list: auto.arima + two sensible manual competitors (same d, D from auto)
fit_auto <- auto.arima(
  train,
  seasonal = TRUE,
  stepwise = FALSE,
  approximation = FALSE,
  trace = FALSE
)
ord <- arimaorder(fit_auto)
cat("\n=== auto.arima (train) ===\n")
print(fit_auto)
cat("AIC =", fit_auto$aic, "\n")

# Manual alternatives around the chosen orders (adjust if auto fails)
p0 <- ord[1]; d0 <- ord[2]; q0 <- ord[3]
P0 <- ord[4]; D0 <- ord[5]; Q0 <- ord[6]

fit_alt1 <- tryCatch(
  Arima(train, order = c(max(0, p0 - 1), d0, q0), seasonal = c(P0, D0, Q0)),
  error = function(e) NULL
)
fit_alt2 <- tryCatch(
  Arima(train, order = c(p0, d0, max(0, q0 - 1)), seasonal = c(P0, D0, Q0)),
  error = function(e) NULL
)

sarima_candidates <- list(auto = fit_auto)
if (!is.null(fit_alt1)) sarima_candidates$p_minus_1 <- fit_alt1
if (!is.null(fit_alt2)) sarima_candidates$q_minus_1 <- fit_alt2

cat("\n=== SARIMA candidate AIC (train) ===\n")
aic_tab <- tibble::tibble(
  model = names(sarima_candidates),
  AIC = vapply(sarima_candidates, function(m) m$aic, numeric(1))
) %>% arrange(AIC)
print(aic_tab)

fit_sarima <- sarima_candidates[[aic_tab$model[1]]]
cat("\nSelected SARIMA for test evaluation:", aic_tab$model[1], "\n")

fc_sarima <- forecast(fit_sarima, h = h_test, level = c(95))

png(filename = file.path(OUT_DIR, "10_forecast_sarima_vs_test.png"), width = 900, height = 500, res = 120)
plot(fc_sarima, main = "SARIMA forecast vs test set")
lines(test, col = "red", lwd = 2)
dev.off()

# Residual diagnostics (includes ACF of residuals + histogram)
png(filename = file.path(OUT_DIR, "11_sarima_residual_diagnostics.png"), width = 1000, height = 1000, res = 120)
checkresiduals(fit_sarima)
dev.off()

# Ljung-Box on residuals: fitdf = number of estimated parameters in the conditional mean
# (AR, MA, seasonal AR/MA, intercept/drift, etc.); conservative vs omitting scale.
n_arima_coef <- length(coef(fit_sarima))
lb_lag <- min(2 * frequency(train), length(residuals(fit_sarima)) - 1)
lb <- Box.test(
  residuals(fit_sarima, type = "innovation"),
  lag = lb_lag,
  type = "Ljung-Box",
  fitdf = n_arima_coef
)
cat("\n=== Ljung-Box on SARIMA innovations ===\n")
cat("lag =", lb_lag, "| fitdf =", n_arima_coef, "(estimated mean-structure coefficients)\n")
print(lb)

# --- 8. Compare methods on the test set ----------------------------------------

cmp <- rbind(
  HoltWinters = test_accuracy_row(fc_hw, test),
  ETS = test_accuracy_row(fc_ets, test),
  STLF = test_accuracy_row(fc_stlf, test),
  SARIMA = test_accuracy_row(fc_sarima, test)
)
rownames(cmp) <- c("HoltWinters", "ETS", "STLF", "SARIMA")
cat("\n=== Test-set accuracy (rows sorted by RMSE) ===\n")
cmp_sorted <- cmp[order(cmp[, "RMSE"]), , drop = FALSE]
print(round(cmp_sorted, 3))
best_row <- rownames(cmp_sorted)[1]
cat("\nBest RMSE on test set:", best_row, "\n")

utils::write.csv(
  round(cmp_sorted, 4),
  file.path(OUT_DIR, "test_accuracy_by_method.csv"),
  row.names = TRUE
)

# --- 9. Out-of-sample forecasts (refit on full series) -------------------------

h_oos <- 24L
fit_hw_full <- HoltWinters(porto_ts, seasonal = "additive")
fc_hw_oos <- forecast(fit_hw_full, h = h_oos, level = c(95))

fit_ets_full <- ets(porto_ts)
fc_ets_oos <- forecast(fit_ets_full, h = h_oos, level = c(95))

fc_stlf_oos <- stlf(porto_ts, method = "ets", h = h_oos, level = c(95))

# Refit SARIMA using same structure as selected train model
fit_sarima_full <- Arima(porto_ts, model = fit_sarima)
fc_sarima_oos <- forecast(fit_sarima_full, h = h_oos, level = c(95))

cat("\n=== Out-of-sample point forecasts (next", h_oos, "months) — head ===\n")
print(head(as.data.frame(fc_sarima_oos), 6))

png(filename = file.path(OUT_DIR, "12_oos_holtwinters.png"), width = 900, height = 500, res = 120)
plot(fc_hw_oos, main = "Holt-Winters — out-of-sample (95% interval)")
dev.off()

png(filename = file.path(OUT_DIR, "13_oos_ets.png"), width = 900, height = 500, res = 120)
plot(fc_ets_oos, main = "ETS — out-of-sample (95% interval)")
dev.off()

png(filename = file.path(OUT_DIR, "14_oos_stlf.png"), width = 900, height = 500, res = 120)
plot(fc_stlf_oos, main = "STLF — out-of-sample (95% interval)")
dev.off()

png(filename = file.path(OUT_DIR, "15_oos_sarima.png"), width = 900, height = 500, res = 120)
plot(fc_sarima_oos, main = "SARIMA — out-of-sample (95% interval)")
dev.off()

last_obs_date <- max(porto_prices$date)
oos_month_start <- seq.Date(last_obs_date, by = "1 month", length.out = h_oos + 1L)[-1L]

oos_tbl <- tibble::tibble(
  month = oos_month_start,
  sarima_point = as.numeric(fc_sarima_oos$mean),
  sarima_lo95 = as.numeric(fc_sarima_oos$lower[, 1]),
  sarima_hi95 = as.numeric(fc_sarima_oos$upper[, 1]),
  hw_point = as.numeric(fc_hw_oos$mean),
  hw_lo95 = as.numeric(fc_hw_oos$lower[, 1]),
  hw_hi95 = as.numeric(fc_hw_oos$upper[, 1]),
  ets_point = as.numeric(fc_ets_oos$mean),
  ets_lo95 = as.numeric(fc_ets_oos$lower[, 1]),
  ets_hi95 = as.numeric(fc_ets_oos$upper[, 1]),
  stlf_point = as.numeric(fc_stlf_oos$mean),
  stlf_lo95 = as.numeric(fc_stlf_oos$lower[, 1]),
  stlf_hi95 = as.numeric(fc_stlf_oos$upper[, 1])
)
utils::write.csv(oos_tbl, file.path(OUT_DIR, "out_of_sample_forecasts_95.csv"), row.names = FALSE)

# --- 10. Brief conclusions (for report; expand with interpretation) ---------------

cat("\n=== Summary for report (automated text) ===\n")
cat(
  "- Smoothing: Holt-Winters and ETS fitted on training data; parameter/state estimates printed above.\n",
  "- Decomposition: STL and classical additive decomposition plotted; seasonally adjusted series from STL saved in R object `sa_series`.\n",
  "- SARIMA: auto.arima shortlist vs manual neighbors by AIC; residual ACF/PACF via checkresiduals(); Ljung-Box with fitdf =",
  n_arima_coef, " (all estimated mean parameters).\n",
  "- Comparison: test-set metrics table above; best by RMSE:", best_row, ".\n",
  "- Out-of-sample: 95% intervals for HW, ETS, STLF, SARIMA over h =", h_oos, ".\n",
  sep = ""
)
cat(
  "\nBenefits / limitations (draft):\n",
  "- Smoothing (HW/ETS): fast, strong baselines for trend+seasonality; less explicit uncertainty about ARIMA structure.\n",
  "- STLF: combines STL with robust seasonal structure; depends on stable seasonal pattern.\n",
  "- SARIMA: interpretable dynamics + formal diagnostics; misspecification shows up in residual correlation; intervals assume Gaussian innovations.\n",
  sep = ""
)
cat("\nFigures + tables written to:", OUT_DIR, "\n")
