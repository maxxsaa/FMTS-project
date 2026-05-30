# 03 — ARIMA / SARIMA on log-transformed, differenced series

train_log <- log(train)
train_diff <- diff(train_log)

png(file.path(CFG$out_fig, "13_log_series_train.png"), width = 900, height = 500, res = 120)
plot(train_log, main = "Log-transformed training series", ylab = "log(EUR per m²)")
dev.off()

png(file.path(CFG$out_fig, "14_diff_levels_train.png"), width = 900, height = 500, res = 120)
plot(diff(train), main = "First difference of levels (training)", ylab = "Δ EUR per m²")
dev.off()

png(file.path(CFG$out_fig, "15_diff_log_train.png"), width = 900, height = 500, res = 120)
plot(train_diff, main = "First difference of log series (training)", ylab = "Δ log(EUR per m²)")
dev.off()

png(file.path(CFG$out_fig, "16_acf_pacf_diff_log.png"), width = 1200, height = 600, res = 120)
par(mfrow = c(1, 2))
acf(train_diff, lag.max = 36, main = "ACF — differenced log series")
pacf(train_diff, lag.max = 36, main = "PACF — differenced log series")
dev.off()
par(mfrow = c(1, 1))

# Short-list of candidate models (fixed d = 1 on log scale)
cand_a <- forecast::auto.arima(
  train_log, d = 1, seasonal = FALSE,
  stepwise = FALSE, approximation = FALSE
)
cand_b <- forecast::auto.arima(
  train_log, d = 1, seasonal = TRUE,
  stepwise = FALSE, approximation = FALSE
)
cand_c <- forecast::Arima(
  train_log,
  order = c(3, 1, 0),
  seasonal = list(order = c(0, 0, 1), period = 12)
)

sarima_candidates <- list(
  auto_arima = cand_a,
  auto_sarima = cand_b,
  manual_sarima = cand_c
)

aic_tab <- tibble::tibble(
  model = names(sarima_candidates),
  label = vapply(sarima_candidates, arima_model_label, character(1)),
  AIC = vapply(sarima_candidates, function(m) m$aic, numeric(1)),
  BIC = vapply(sarima_candidates, function(m) m$bic, numeric(1))
) %>% dplyr::arrange(AIC)

readr::write_csv(aic_tab, file.path(CFG$out_tab, "arima_aic_comparison.csv"))
cat("\n=== ARIMA candidate models (log scale, d = 1) ===\n")
print(aic_tab)

fit_arima <- sarima_candidates[[aic_tab$model[1]]]
cat("\nSelected model:", aic_tab$label[1], "\n")

fc_arima_log <- forecast::forecast(fit_arima, h = h_test, level = 95)

png(file.path(CFG$out_fig, "17_arima_residual_diagnostics.png"), width = 1000, height = 1000, res = 120)
checkresiduals(fit_arima)
dev.off()

n_coef <- length(coef(fit_arima))
lb_lag <- min(2 * frequency(train), length(residuals(fit_arima)) - 1)
lb_test <- Box.test(
  residuals(fit_arima, type = "innovation"),
  lag = lb_lag,
  type = "Ljung-Box",
  fitdf = n_coef
)
cat("\n=== Ljung-Box (innovation residuals) ===\n")
cat("lag =", lb_lag, "| fitdf =", n_coef, "\n")
print(lb_test)
capture.output(lb_test, file = file.path(CFG$out_test, "ljung_box_arima.txt"))

plot_exp_forecast_vs_test(
  fc_arima_log, train, test,
  paste0(arima_model_label(fit_arima), " vs test set (EUR/m²)"),
  file.path(CFG$out_fig, "18_forecast_arima_vs_test.png")
)
