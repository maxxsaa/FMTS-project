# ==============================================================================
# Script 1: Exploratory Data Analysis (EDA)
# ==============================================================================

library(tidyverse)
library(forecast)
library(pracma)
library(tseries)

# 1. Paths and Directories
input_path <- "/Users/pedromcorreia/Desktop/MPST/timeSeriesPorto.csv"
output_dir <- "/Users/pedromcorreia/Desktop/MPST/outputs"

out_fig <- file.path(output_dir, "figures")
out_tab <- file.path(output_dir, "tables")
out_test <- file.path(output_dir, "tests")

dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(out_test, recursive = TRUE, showWarnings = FALSE)

# 2. Data Processing
raw_data <- read_csv(input_path, col_names = "raw_line", col_types = "c", show_col_types = FALSE)
porto_prices <- raw_data %>%
  separate(raw_line, into = c("month_year", "series_type", "region_code", "euros_str", "extra"), sep = ";", extra = "drop", fill = "right") %>%
  mutate(
    year = as.numeric(str_extract(month_year, "\\d{4}$")),
    month_number = match(tolower(str_extract(month_year, "^[A-Za-z]+")), tolower(month.name)),
    euros_per_sqm = as.numeric(gsub(",", ".", euros_str)),
    date = as.Date(paste0(year, "-", sprintf("%02d", month_number), "-01"))
  ) %>%
  filter(!is.na(date), !is.na(euros_per_sqm)) %>%
  arrange(date)

porto_ts <- ts(porto_prices$euros_per_sqm, start = c(min(porto_prices$year), min(porto_prices$month_number)), frequency = 12)
years <- seq(from = min(porto_prices$year), to = max(porto_prices$year), by = 1)

# 3. Exploratory Plots
png(file.path(out_fig, "Fig1_time_plot.png"), width = 800, height = 600)
plot(porto_ts, main = "Porto Property Prices", ylab = "EUR per sqm", xlab = "Time", xaxt = "n")
axis(1, at = years, labels = years)
dev.off()

png(file.path(out_fig, "Fig2_seasonal_boxplot.png"), width = 800, height = 600)
boxplot(porto_ts ~ cycle(porto_ts), names = month.abb, main = "Seasonal Boxplot")
dev.off()

png(file.path(out_fig, "Fig3_seasonal_boxplot_subset.png"), width = 800, height = 600)
porto_ts_subset <- window(porto_ts, end = c(2025, 12))
boxplot(porto_ts_subset ~ cycle(porto_ts_subset), names = month.abb, main = "Seasonal Boxplot (Excl. Jan/Feb 2026)")
dev.off()

# 4. Stationarity Tests
adf_test <- adf.test(porto_ts, alternative = "stationary")
capture.output(adf_test, file = file.path(out_test, "Test1_adf_results.txt"))

h_result <- hurstexp(as.numeric(porto_ts), display = FALSE)
capture.output(print(paste("Hurst Exponent:", h_result$Hs)), file = file.path(out_test, "Test2_hurst_results.txt"))

# 5. Autocorrelation & Lags
png(file.path(out_fig, "Fig4_lag_plots.png"), width = 800, height = 600)
lag.plot(porto_ts, lags = 4, do.lines = FALSE, main = "Lag Plots (k=1 to 4)")
dev.off()

png(file.path(out_fig, "Fig5_lag_plots_extended.png"), width = 800, height = 600)
lag.plot(porto_ts, lags = 12, do.lines = FALSE, diag.col = "red", main = "Lag Plots (1 to 12)")
dev.off()

png(file.path(out_fig, "Fig6_acf_pacf_plots.png"), width = 800, height = 600)
par(mfrow = c(2, 1))
acf(porto_ts, main = "Autocorrelation Function")
pacf(porto_ts, main = "Partial Autocorrelation Function")
dev.off()

png(file.path(out_fig, "Fig7_acf_seasonality_check.png"), width = 800, height = 600)
acf(porto_ts, lag.max = 24, main = "ACF: Check for Seasonality at Lag 12")
dev.off()