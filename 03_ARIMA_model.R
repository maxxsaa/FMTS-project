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

# ==============================================================================
# 4. Model Selection: Main Benchmarks vs. Annex Sensitivity Analysis
# ==============================================================================

# 1. Main Benchmarks (Automated + Your Initial Manual Choice)
main_models <- list(
  "Benchmark: Auto-ARIMA" = auto.arima(train_log, d = 1, seasonal = FALSE, stepwise = FALSE, approximation = FALSE),
  "Benchmark: Auto-SARIMA" = auto.arima(train_log, d = 1, seasonal = TRUE, stepwise = FALSE, approximation = FALSE),
  "Manual Choice: SARIMA(3,1,0)(0,0,1)[12]" = Arima(train_log, order = c(3, 1, 0), seasonal = list(order = c(0, 0, 1), period = 12)),
  "Manual Choice: ARIMA(3,1,4)" = Arima(train_log, order = c(3, 1, 4))
)

# 2. Annex Sensitivity Analysis (The "Exploration" models)
annex_models <- list(
  "Annex: ARIMA(1,1,3)" = Arima(train_log, order = c(1, 1, 3)),
  "Annex: ARIMA(2,1,4)" = Arima(train_log, order = c(2, 1, 4)),
  "Annex: ARIMA(0,1,4)" = Arima(train_log, order = c(0, 1, 4)),
  "Annex: SARIMA(1,1,4)(0,0,1)[12]" = Arima(train_log, order = c(1, 1, 4), seasonal = list(order = c(0, 0, 1), period = 12)),
  "Annex: ARIMA(3,1,0)" = Arima(train_log, order = c(3, 1, 0)),
  "Annex: SARIMA(2,1,0)(0,0,1)[12]" = Arima(train_log, order = c(2, 1, 0), seasonal = list(order = c(0, 0, 1), period = 12)),
  "Annex: SARIMA(3,1,1)(0,0,1)[12]" = Arima(train_log, order = c(3, 1, 1), seasonal = list(order = c(0, 0, 1), period = 12))
)

# 3. Combine for comprehensive comparison
all_models <- c(main_models, annex_models)

# 4. Create and save the table
model_comparison <- data.frame(
  Model = names(all_models),
  AIC = sapply(all_models, function(m) m$aic),
  BIC = sapply(all_models, function(m) m$bic)
)

# Sort by AIC
model_comparison_sorted <- model_comparison[order(model_comparison$AIC), ]
write.csv(model_comparison_sorted, file.path(out_tab, "Tab2_aic_bic_comparison_full.csv"), row.names = FALSE)

# 5. Selection Logic: Choose model with best BIC for final diagnostics
bic_values <- sapply(all_models, function(m) m$bic)
winning_model_name <- names(which.min(bic_values))
winning_model <- all_models[[winning_model_name]]

cat("Diagnostic tests performed on:", winning_model_name, 
    "(selected for superior BIC/parsimony).\n")

# ==============================================================================
# 5 & 6. Diagnostic Checking & Predictive Performance (BIC vs. AIC)
# ==============================================================================

# 1. Identify the two winners
bic_winner_name <- names(which.min(sapply(all_models, function(m) m$bic)))
aic_winner_name <- names(which.min(sapply(all_models, function(m) m$aic)))

winners <- list(BIC = all_models[[bic_winner_name]], AIC = all_models[[aic_winner_name]])
winner_names <- list(BIC = bic_winner_name, AIC = aic_winner_name)

# 2. Loop through both models to generate diagnostics and accuracy tables
for (w in names(winners)) {
  model_obj <- winners[[w]]
  model_name <- winner_names[[w]]
  
  # Diagnostic Checking
  png(file.path(out_fig, paste0("Fig20_", w, "_residuals.png")), width = 800, height = 800)
  checkresiduals(model_obj, main = paste("Residuals:", w, "-Winner"))
  dev.off()
  
  df_model <- length(model_obj$coef)
  lb_test <- Box.test(residuals(model_obj), type = "Ljung-Box", lag = 24, fitdf = df_model)
  capture.output(lb_test, file = file.path(out_test, paste0("Test7_", w, "_ljung_box.txt")))
  
  # Predictive Performance
  forecast_test <- forecast(model_obj, h = h_test)
  # Remember: we are working with log-transformed data, so we exponentiate
  exp_point_forecasts <- exp(forecast_test$mean)
  accuracy_test <- test_accuracy_row(exp_point_forecasts, porto_ts_test)
  
  write.csv(round(accuracy_test, 4), file.path(out_tab, paste0("Tab3_", w, "_test_accuracy.csv")), row.names = TRUE)
  
  # Plot Forecast vs Actual
  png(file.path(out_fig, paste0("Fig21_forecast_", w, ".png")), width = 800, height = 600)
  plot(forecast_test, main = paste("Out-of-Sample:", w, "-Winner"), ylab = "log(EUR per sqm)", xlab = "Time")
  lines(log(porto_ts_test), col = "red", lwd = 2)
  legend("topleft", legend = c("Forecast", "Actuals (Test Set)"), col = c("blue", "red"), lty = 1, lwd = 2)
  dev.off()
  
  cat("Completed diagnostics and forecast for:", model_name, "(", w, " Winner)\n")
}