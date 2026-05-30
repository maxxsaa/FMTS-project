# 04 — Test-set comparison and out-of-sample forecasts

cmp_list <- list(
  HW_Multiplicative = extract_test_metrics(fc_hw_mult, test),
  HW_Damped         = extract_test_metrics(fc_hw_damped, test),
  ETS               = extract_test_metrics(fc_ets, test),
  STLF              = extract_test_metrics(fc_stlf, test),
  ARIMA_log         = extract_test_metrics(exp_forecast(fc_arima_log)$mean, test)
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
exp_arima <- exp_forecast(fc_arima_log)
fc_vs_test <- tibble::tibble(
  month = format(test_dates, "%Y-%m"),
  actual = as.numeric(test),
  hw_mult = round(as.numeric(fc_hw_mult$mean), 1),
  hw_mult_lo95 = round(as.numeric(fc_hw_mult$lower[, 1]), 1),
  hw_mult_hi95 = round(as.numeric(fc_hw_mult$upper[, 1]), 1),
  arima_log = round(exp_arima$mean, 1),
  arima_lo95 = round(exp_arima$lower, 1),
  arima_hi95 = round(exp_arima$upper, 1),
  ets = round(as.numeric(fc_ets$mean), 1),
  stlf = round(as.numeric(fc_stlf$mean), 1)
)
readr::write_csv(fc_vs_test, file.path(CFG$out_tab, "forecast_vs_test.csv"))

# --- Out-of-sample: refit on full series --------------------------------------
fit_hw_mult_full <- HoltWinters(porto_ts, seasonal = "multiplicative")
fc_hw_mult_oos <- forecast::forecast(fit_hw_mult_full, h = H_OOS, level = 95)

fit_hw_damped_oos <- forecast::hw(porto_ts, seasonal = "multiplicative", damped = TRUE, h = H_OOS, level = 95)
fc_hw_damped_oos <- fit_hw_damped_oos

fit_ets_full <- ets(porto_ts)
fc_ets_oos <- forecast::forecast(fit_ets_full, h = H_OOS, level = 95)
fc_stlf_oos <- forecast::stlf(porto_ts, h = H_OOS, level = 95)

fit_arima_full <- forecast::Arima(log(porto_ts), model = fit_arima)
fc_arima_log_oos <- forecast::forecast(fit_arima_full, h = H_OOS, level = 95)
exp_oos_arima <- exp_forecast(fc_arima_log_oos)

plot_forecast_vs_test(fc_hw_mult_oos, window(porto_ts, start = time(porto_ts)[1]),
  "Holt-Winters (multiplicative) — out-of-sample",
  file.path(CFG$out_fig, "19_oos_hw_mult.png"))

plot_forecast_vs_test(fc_ets_oos, window(porto_ts, start = time(porto_ts)[1]),
  paste0("ETS (", fit_ets_full$method, ") — out-of-sample"),
  file.path(CFG$out_fig, "20_oos_ets.png"))

plot_forecast_vs_test(fc_stlf_oos, window(porto_ts, start = time(porto_ts)[1]),
  "STLF — out-of-sample",
  file.path(CFG$out_fig, "21_oos_stlf.png"))

# ARIMA OOS in EUR
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
  arima_point = exp_oos_arima$mean,
  arima_lo95 = exp_oos_arima$lower,
  arima_hi95 = exp_oos_arima$upper
)
readr::write_csv(oos_tbl, file.path(CFG$out_tab, "out_of_sample_forecasts_95.csv"))

grDevices::png(file.path(CFG$out_fig, "22_oos_arima_log.png"), width = 900, height = 500, res = 120)
plot(porto_ts, main = paste0(arima_model_label(fit_arima_full), " — out-of-sample (EUR/m²)"),
     ylab = "EUR per m²", xlim = c(start(porto_ts)[1], max(time(porto_ts)) + H_OOS / 12))
t_oos <- seq(max(time(porto_ts)) + 1 / 12, by = 1 / 12, length.out = H_OOS)
lines(t_oos, exp_oos_arima$mean, col = "blue", lwd = 2)
polygon(c(t_oos, rev(t_oos)), c(exp_oos_arima$lower, rev(exp_oos_arima$upper)),
        col = rgb(0.3, 0.5, 0.9, 0.25), border = NA)
legend("topleft", legend = c("Forecast", "95% PI"), col = c("blue", "lightblue"),
       lty = 1, lwd = c(2, 6), bty = "n")
dev.off()

cat("\n=== Pipeline complete ===\n")
cat("Figures:", CFG$out_fig, "\n")
cat("Tables:", CFG$out_tab, "\n")
