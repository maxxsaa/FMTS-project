# 01 — Exploratory data analysis
# Note: unit-root tests (ADF/PP/KPSS) are not applied to this seasonal series
# (see course guidelines); stationarity is assessed via ACF/PACF and differencing.

years <- seq(
  lubridate::year(min(porto_prices$date)),
  lubridate::year(max(porto_prices$date)),
  by = 1
)

png(file.path(CFG$out_fig, "01_time_plot.png"), width = 900, height = 500, res = 120)
plot(porto_ts, main = "Porto median dwelling price (EUR/m²) — monthly",
     ylab = "EUR per m²", xlab = "Time")
dev.off()

# Detrended series: multiplicative STL trend removed, re-scaled to EUR/m²
stl_log_eda <- stats::stl(log(porto_ts), s.window = "periodic", robust = TRUE)
trend_eda <- exp(stl_log_eda$time.series[, "trend"])
seasonal_eda <- exp(stl_log_eda$time.series[, "seasonal"])
detrended_ts <- (porto_ts / trend_eda) * mean(trend_eda)
detrended_ts <- ts(detrended_ts, start = start(porto_ts), frequency = frequency(porto_ts))
png(file.path(CFG$out_fig, "01_detrended_time_plot.png"), width = 900, height = 500, res = 120)
plot(detrended_ts, main = "Porto median dwelling price (EUR/m²) — detrended (monthly)",
     ylab = "EUR per m² (trend removed)", xlab = "Time")
dev.off()

# Deseasonalized series: multiplicative STL seasonal factor removed
deseasonalized_ts <- porto_ts / seasonal_eda
deseasonalized_ts <- ts(deseasonalized_ts, start = start(porto_ts), frequency = frequency(porto_ts))
png(file.path(CFG$out_fig, "01_deseasonalized_time_plot.png"), width = 900, height = 500, res = 120)
plot(deseasonalized_ts, main = "Porto median dwelling price (EUR/m²) — deseasonalized (monthly)",
     ylab = "EUR per m² (seasonality removed)", xlab = "Time")
dev.off()

png(file.path(CFG$out_fig, "02_seasonal_boxplot.png"), width = 900, height = 500, res = 120)
boxplot(porto_ts ~ cycle(porto_ts), names = month.abb,
        main = "Seasonal pattern by month", ylab = "EUR per m²")
dev.off()

png(file.path(CFG$out_fig, "03_acf_pacf_full.png"), width = 900, height = 700, res = 120)
par(mfrow = c(2, 1))
acf(porto_ts, lag.max = 48, main = "ACF — full series (levels)")
pacf(porto_ts, lag.max = 48, main = "PACF — full series (levels)")
dev.off()
par(mfrow = c(1, 1))

png(file.path(CFG$out_fig, "04_lag_plots.png"), width = 900, height = 500, res = 120)
stats::lag.plot(porto_ts, lags = 12, do.lines = FALSE, diag.col = "red",
                main = "Lag plots (k = 1 … 12)")
dev.off()

cat("=== EDA complete ===\n")
cat("n =", length(porto_ts), "| train ends:", paste(end(train), collapse = "/"),
    "| test h =", h_test, "\n")
