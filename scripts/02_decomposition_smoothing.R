# 02 — Multiplicative decomposition and smoothing methods

# --- Decomposition (multiplicative) -------------------------------------------
# Classical: y_t = T_t × S_t × R_t
decomp_mult <- decompose(porto_ts, type = "multiplicative")
png(file.path(CFG$out_fig, "05_decomposition_multiplicative.png"), width = 1000, height = 900, res = 120)
plot(decomp_mult)
dev.off()

# STL on log(y): additive STL on logs ≈ multiplicative decomposition on levels
porto_log <- log(porto_ts)
stl_log <- stl(porto_log, s.window = "periodic", robust = TRUE)
png(file.path(CFG$out_fig, "06_stl_decomposition_log.png"), width = 1000, height = 900, res = 120)
plot(stl_log, main = "STL decomposition of log(prices)")
dev.off()

# Seasonally adjusted (multiplicative): y_t / S_t
sa_mult <- porto_ts / decomp_mult$seasonal
png(file.path(CFG$out_fig, "07_seasonally_adjusted.png"), width = 900, height = 500, res = 120)
plot(cbind(Original = porto_ts, Seasonally_adjusted = sa_mult),
     main = "Original vs multiplicatively seasonally adjusted",
     ylab = "EUR per m²")
dev.off()

readr::write_csv(
  tibble::tibble(
    date = porto_prices$date,
    original = as.numeric(porto_ts),
    seasonally_adjusted = as.numeric(sa_mult)
  ),
  file.path(CFG$out_tab, "seasonally_adjusted.csv")
)

# --- Smoothing (training sample) ----------------------------------------------
fit_hw_mult <- HoltWinters(train, seasonal = "multiplicative")
cat("\n=== Holt-Winters (multiplicative) ===\n")
print(fit_hw_mult)

fit_ets <- ets(train)
cat("\n=== ETS ===\n")
print(fit_ets)

fc_hw_mult   <- forecast::forecast(fit_hw_mult, h = h_test, level = 95)
fc_hw_damped <- forecast::hw(train, seasonal = "multiplicative", damped = TRUE, h = h_test, level = 95)
fc_ets       <- forecast::forecast(fit_ets, h = h_test, level = 95)
fc_stlf      <- forecast::stlf(train, h = h_test, level = 95)

capture.output(summary(fit_hw_mult), file = file.path(CFG$out_test, "hw_mult_summary.txt"))
capture.output(summary(fit_ets), file = file.path(CFG$out_test, "ets_summary.txt"))

plot_forecast_vs_test(fc_hw_mult, test,
  "Holt-Winters (multiplicative) vs test set",
  file.path(CFG$out_fig, "08_forecast_hw_mult_vs_test.png"))

plot_forecast_vs_test(fc_hw_damped, test,
  "Holt-Winters (multiplicative, damped) vs test set",
  file.path(CFG$out_fig, "09_forecast_hw_damped_vs_test.png"))

plot_forecast_vs_test(fc_ets, test,
  paste0("ETS (", fit_ets$method, ") vs test set"),
  file.path(CFG$out_fig, "10_forecast_ets_vs_test.png"))

plot_forecast_vs_test(fc_stlf, test,
  "STLF vs test set",
  file.path(CFG$out_fig, "11_forecast_stlf_vs_test.png"))

# Combined zoom (2020 onward)
start_zoom <- 2020
end_zoom <- max(time(test))
png(file.path(CFG$out_fig, "12_combined_smoothing_zoomed.png"), width = 1200, height = 1000, res = 120)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (item in list(
  list(fc_hw_mult, "Holt-Winters (multiplicative)"),
  list(fc_hw_damped, "Holt-Winters (damped)"),
  list(fc_ets, paste("ETS (", fit_ets$method, ")", sep = "")),
  list(fc_stlf, "STLF")
)) {
  plot(item[[1]], main = item[[2]], ylab = "EUR per m²", xlim = c(start_zoom, end_zoom))
  lines(test, col = "red", lwd = 2)
}
dev.off()
par(mfrow = c(1, 1))
