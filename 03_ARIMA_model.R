# ==============================================================================
# Script 3: ARIMA and SARIMA Modeling
# ==============================================================================

library(tidyverse)
library(forecast)

# 1. Paths, Directories, and Helper Function
input_path <- "/Users/pedromcorreia/Desktop/MPST/timeSeriesPorto.csv"
output_dir <- "/Users/pedromcorreia/Desktop/MPST/outputs"
out_fig <- file.path(output_dir, "figures")
out_tab <- file.path(output_dir, "tables")
out_test <- file.path(output_dir, "tests")

test_accuracy_row <- function(fc, actual) {
  a <- accuracy(fc, actual)
  if ("Test set" %in% rownames(a)) return(a["Test set", , drop = FALSE])
  if (nrow(a) >= 2) return(a[2, , drop = FALSE])
  a[1, , drop = FALSE]
}

# 2. Data Processing & Train/Test Split
raw_data <- read_csv(input_path, col_names = "raw_line", col_types = "c", show_col_types = FALSE)
porto_prices <- raw_data %>%
  separate(raw_line, into = c("month_year", "series_type", "region_code", "euros_str", "extra"), sep = ";", extra = "drop", fill = "right") %>%
  mutate(
    year = as.numeric(str_extract(month_year, "\\d{4}$")),
    month_number = match(tolower(str_extract(month_year, "^[A-Za-z]+")), tolower(month.name)),
    euros_per_sqm = as.numeric(gsub(",", ".", euros_str)),
    date = as.Date(paste0(year, "-", sprintf("%02d", month_number), "-01"))
  ) %>% filter(!is.na(date), !is.na(euros_per_sqm)) %>% arrange(date)

porto_ts <- ts(porto_prices$euros_per_sqm, start = c(min(porto_prices$year), min(porto_prices$month_number)), frequency = 12)
years <- seq(from = min(porto_prices$year), to = max(porto_prices$year), by = 1)

h_test <- 12
porto_ts_train <- head(porto_ts, length(porto_ts) - h_test)
porto_ts_test <- tail(porto_ts, h_test)

train_log <- log(porto_ts_train)
train_diff <- diff(train_log)

# 3. Transformations and Stationarity
png(file.path(out_fig, "Fig15_plot_log.png"), width = 800, height = 600)
plot(train_log, main = "Log-Transformed Training Data", ylab = "log(EUR per sqm)", xlab = "Time", xaxt = "n")
axis(1, at = years, labels = years)
dev.off()

png(file.path(out_fig, "Fig16_plot_diff.png"), width = 800, height = 600)
plot(diff(porto_ts_train), main = "Raw First Diff (Training)", ylab = "Change in EUR per sqm", xlab = "Time", xaxt = "n")
axis(1, at = years, labels = years)
dev.off()

png(file.path(out_fig, "Fig17_plot_diff_log.png"), width = 800, height = 600)
plot(train_diff, main = "First Differenced Log-Series (Training)", ylab = "Change in log(EUR/sqm)", xlab = "Time", xaxt = "n")
axis(1, at = years, labels = years)
dev.off()

png(file.path(out_fig, "Fig18_acf_differenced.png"), width = 800, height = 600)
acf(train_diff, lag.max = 36, main = "ACF of Differenced Data")
dev.off()

png(file.path(out_fig, "Fig19_pacf_differenced.png"), width = 800, height = 600)
pacf(train_diff, lag.max = 36, main = "PACF of Differenced Data")
dev.off()

# 4. Model Selection
cand_A_arima <- auto.arima(train_log, d = 1, seasonal = FALSE, stepwise = FALSE, approximation = FALSE)
cand_B_sarima <- auto.arima(train_log, d = 1, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
cand_C_manual <- Arima(train_log, order = c(3, 1, 0), seasonal = list(order = c(0, 0, 1), period = 12))

model_comparison <- data.frame(
  Model = c(
    paste("Cand A: ARIMA", paste(cand_A_arima$arma[c(1,6,2)], collapse=",")),
    paste("Cand B: SARIMA", paste(cand_B_sarima$arma[c(1,6,2)], collapse=",")),
    "Cand C: SARIMA(3,1,0)(0,0,1)[12]"
  ),
  AIC = c(cand_A_arima$aic, cand_B_sarima$aic, cand_C_manual$aic),
  BIC = c(cand_A_arima$bic, cand_B_sarima$bic, cand_C_manual$bic)
)
write.csv(model_comparison, file.path(out_tab, "Tab2_aic_bic_comparison.csv"), row.names = FALSE)

# Proceeding with Candidate A for diagnostics
winning_model <- cand_A_arima

# 5. Diagnostic Checking
png(file.path(out_fig, "Fig20_arima_residuals.png"), width = 800, height = 800)
checkresiduals(winning_model)
dev.off()

df_sarima <- length(winning_model$coef)
lb_test <- Box.test(residuals(winning_model), type = "Ljung-Box", lag = 24, fitdf = df_sarima)
capture.output(lb_test, file = file.path(out_test, "Test7_ljung_box.txt"))

# 6. Predictive Performance on Test Set
forecast_test <- forecast(winning_model, h = h_test)
exp_point_forecasts <- exp(forecast_test$mean)
accuracy_test <- test_accuracy_row(exp_point_forecasts, porto_ts_test)

write.csv(round(accuracy_test, 4), file.path(out_tab, "Tab3_arima_test_accuracy.csv"), row.names = TRUE)

png(file.path(out_fig, "Fig21_forecast_vs_actual.png"), width = 800, height = 600)
plot(forecast_test, main = "Out-of-Sample Performance: ARIMA vs Actuals", ylab = "log(EUR per sqm)", xlab = "Time")
lines(log(porto_ts_test), col = "red", lwd = 2)
legend("topleft", legend = c("Forecast", "Actuals (Test Set)"), col = c("blue", "red"), lty = 1, lwd = 2)
dev.off()