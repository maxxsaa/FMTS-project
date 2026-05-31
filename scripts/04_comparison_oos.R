# 04 — Test-set comparison and out-of-sample forecasts

exp_arima_114 <- exp_forecast(fc_arima_log)
exp_arima_314 <- exp_forecast(fc_arima_314_log)

cmp_list <- list(
  HW_Multiplicative = extract_test_metrics(fc_hw_mult, test),
  HW_Damped         = extract_test_metrics(fc_hw_damped, test),
  ETS               = extract_test_metrics(fc_ets, test),
  STLF              = extract_test_metrics(fc_stlf, test),
  ARIMA_114_log     = extract_test_metrics(exp_arima_114$mean, test),
  ARIMA_314_log     = extract_test_metrics(exp_arima_314$mean, test)
)
cmp <- do.call(rbind, cmp_list)
cmp_sorted <- cmp[order(cmp[, "RMSE"]), , drop = FALSE]
cat("\n=== Test-set accuracy (sorted by RMSE) ===\n")
print(round(cmp_sorted, 3))
best_method <- rownames(cmp_sorted)[1]
cat("\nBest RMSE:", best_method, "\n")

readr::write_csv(
  tibble::rownames_to_column(as.data.frame(round(cmp_sorted, 4)), "method"),
  file.path(CFG$out_tab, "test_accuracy_by_method.csv")
)

# Forecast vs actual table (all methods)
test_dates <- seq.Date(
  as.Date(paste(end(train)[1], round(end(train)[2]), "01", sep = "-")) + 32,
  by = "month",
  length.out = h_test
)
fc_vs_test <- tibble::tibble(
  month = format(test_dates, "%Y-%m"),
  actual = as.numeric(test),
  hw_mult = round(as.numeric(fc_hw_mult$mean), 1),
  hw_mult_lo95 = round(as.numeric(fc_hw_mult$lower[, 1]), 1),
  hw_mult_hi95 = round(as.numeric(fc_hw_mult$upper[, 1]), 1),
  hw_damped = round(as.numeric(fc_hw_damped$mean), 1),
  hw_damped_lo95 = round(as.numeric(fc_hw_damped$lower[, 1]), 1),
  hw_damped_hi95 = round(as.numeric(fc_hw_damped$upper[, 1]), 1),
  ets = round(as.numeric(fc_ets$mean), 1),
  ets_lo95 = round(as.numeric(fc_ets$lower[, 1]), 1),
  ets_hi95 = round(as.numeric(fc_ets$upper[, 1]), 1),
  stlf = round(as.numeric(fc_stlf$mean), 1),
  stlf_lo95 = round(as.numeric(fc_stlf$lower[, 1]), 1),
  stlf_hi95 = round(as.numeric(fc_stlf$upper[, 1]), 1),
  arima_114 = round(exp_arima_114$mean, 1),
  arima_114_lo95 = round(exp_arima_114$lower, 1),
  arima_114_hi95 = round(exp_arima_114$upper, 1),
  arima_314 = round(exp_arima_314$mean, 1),
  arima_314_lo95 = round(exp_arima_314$lower, 1),
  arima_314_hi95 = round(exp_arima_314$upper, 1)
)
readr::write_csv(fc_vs_test, file.path(CFG$out_tab, "forecast_vs_test.csv"))

# --- Out-of-sample: refit on full series --------------------------------------
fit_hw_mult_full <- HoltWinters(porto_ts, seasonal = "multiplicative")
fc_hw_mult_oos <- forecast::forecast(fit_hw_mult_full, h = H_OOS, level = 95)

fc_hw_damped_oos <- forecast::hw(
  porto_ts, seasonal = "multiplicative", damped = TRUE, h = H_OOS, level = 95
)

fit_ets_full <- ets(porto_ts)
fc_ets_oos <- forecast::forecast(fit_ets_full, h = H_OOS, level = 95)
fc_stlf_oos <- forecast::stlf(porto_ts, h = H_OOS, level = 95)

porto_log_full <- log(porto_ts)
fit_arima_114_full <- forecast::Arima(porto_log_full, order = c(1, 1, 4))
fit_arima_314_full <- forecast::Arima(porto_log_full, order = c(3, 1, 4))
fc_arima_114_oos <- forecast::forecast(fit_arima_114_full, h = H_OOS, level = 95)
fc_arima_314_oos <- forecast::forecast(fit_arima_314_full, h = H_OOS, level = 95)
exp_oos_arima_114 <- exp_forecast(fc_arima_114_oos)
exp_oos_arima_314 <- exp_forecast(fc_arima_314_oos)

plot_forecast_oos(
  fc_hw_mult_oos,
  "Holt-Winters (multiplicative) — out-of-sample",
  file.path(CFG$out_fig, "19_oos_hw_mult.png")
)
plot_forecast_oos(
  fc_hw_damped_oos,
  "Holt-Winters (multiplicative, damped) — out-of-sample",
  file.path(CFG$out_fig, "20_oos_hw_damped.png")
)
plot_forecast_oos(
  fc_ets_oos,
  paste0("ETS (", fit_ets_full$method, ") — out-of-sample"),
  file.path(CFG$out_fig, "21_oos_ets.png")
)
plot_forecast_oos(
  fc_stlf_oos,
  "STLF — out-of-sample",
  file.path(CFG$out_fig, "22_oos_stlf.png")
)
plot_exp_forecast_oos(
  fc_arima_114_oos, porto_ts,
  "ARIMA(1,1,4) — out-of-sample (EUR/m²)",
  file.path(CFG$out_fig, "23_oos_arima_114.png")
)
plot_exp_forecast_oos(
  fc_arima_314_oos, porto_ts,
  "ARIMA(3,1,4) — out-of-sample (EUR/m²)",
  file.path(CFG$out_fig, "24_oos_arima_314.png")
)

last_date <- max(porto_prices$date)
oos_dates <- seq.Date(last_date, by = "month", length.out = H_OOS + 1L)[-1L]
oos_tbl <- tibble::tibble(
  month = oos_dates,
  hw_mult_point = as.numeric(fc_hw_mult_oos$mean),
  hw_mult_lo95 = as.numeric(fc_hw_mult_oos$lower[, 1]),
  hw_mult_hi95 = as.numeric(fc_hw_mult_oos$upper[, 1]),
  hw_damped_point = as.numeric(fc_hw_damped_oos$mean),
  hw_damped_lo95 = as.numeric(fc_hw_damped_oos$lower[, 1]),
  hw_damped_hi95 = as.numeric(fc_hw_damped_oos$upper[, 1]),
  ets_point = as.numeric(fc_ets_oos$mean),
  ets_lo95 = as.numeric(fc_ets_oos$lower[, 1]),
  ets_hi95 = as.numeric(fc_ets_oos$upper[, 1]),
  stlf_point = as.numeric(fc_stlf_oos$mean),
  stlf_lo95 = as.numeric(fc_stlf_oos$lower[, 1]),
  stlf_hi95 = as.numeric(fc_stlf_oos$upper[, 1]),
  arima_114_point = exp_oos_arima_114$mean,
  arima_114_lo95 = exp_oos_arima_114$lower,
  arima_114_hi95 = exp_oos_arima_114$upper,
  arima_314_point = exp_oos_arima_314$mean,
  arima_314_lo95 = exp_oos_arima_314$lower,
  arima_314_hi95 = exp_oos_arima_314$upper
)
readr::write_csv(oos_tbl, file.path(CFG$out_tab, "out_of_sample_forecasts_95.csv"))

cat("\n=== Pipeline complete ===\n")
cat("Figures:", CFG$out_fig, "\n")
cat("Tables:", CFG$out_tab, "\n")
