# ==============================================================================
# Script 2: Decomposition and Smoothing Methods
# ==============================================================================

library(tidyverse)
library(forecast)

# 1. Paths, Directories, and Helper Function
input_path <- "/Users/pedromcorreia/Desktop/MPST/timeSeriesPorto.csv"
output_dir <- "/Users/pedromcorreia/Desktop/MPST/outputs"
out_fig <- file.path(output_dir, "figures")
out_tab <- file.path(output_dir, "tables")
out_test <- file.path(output_dir, "tests")

# Function to extract test set accuracy
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

# 3. Decomposition
# Classical Multiplicative Decomposition (for baseline comparison)
porto_decomp <- decompose(porto_ts, type = "multiplicative")
png(file.path(out_fig, "Fig8_decomposition_multiplicative.png"), width = 800, height = 600)
plot(porto_decomp)
dev.off()

# STL Decomposition (Log-Transformed for Multiplicative Reality)
# We log the data so the additive STL can model the proportional (multiplicative) relationships
porto_log <- log(porto_ts)
porto_stl <- stl(porto_log, s.window = "periodic")

png(file.path(out_fig, "Fig9_stl_decomposition_logged.png"), width = 800, height = 600)
# Note: The Y-axis here will now show log values
plot(porto_stl, main = "STL Decomposition (Log-Transformed Data)")
dev.off()

# Seasonally Adjusted Series (using the Classical Decomposition)
porto_adj <- porto_ts / porto_decomp$seasonal
png(file.path(out_fig, "Fig10_seasonally_adjusted.png"), width = 800, height = 600)
plot(porto_adj, main = "Seasonally Adjusted Prices", ylab = "EUR per sqm", xlab = "Time", xaxt = "n")
axis(1, at = years, labels = years)
dev.off()

# 4. Smoothing Methods (Fitted on Training Set)
years_fc <- seq(from = min(porto_prices$year), to = max(porto_prices$year) + ceiling(h_test/12), by = 1)

# Method 1: Holt-Winters (Multiplicative)
hw_mult <- hw(porto_ts_train, seasonal = "multiplicative", h = h_test)
capture.output(summary(hw_mult), file = file.path(out_test, "Test3_hw_mult_summary.txt"))
png(file.path(out_fig, "Fig11_forecast_hw_mult.png"), width = 800, height = 600)
plot(hw_mult, main = "Holt-Winters (Multiplicative) vs Test Set", ylab = "EUR per sqm", xlab = "Time", xaxt = "n")
lines(porto_ts_test, col = "red", lwd = 2)
axis(1, at = years_fc, labels = years_fc)
dev.off()

# Method 2: Holt-Winters (Damped)
hw_damped <- hw(porto_ts_train, seasonal = "multiplicative", damped = TRUE, h = h_test)
capture.output(summary(hw_damped), file = file.path(out_test, "Test4_hw_damped_summary.txt"))
png(file.path(out_fig, "Fig12_forecast_hw_damped.png"), width = 800, height = 600)
plot(hw_damped, main = "Holt-Winters (Damped) vs Test Set", ylab = "EUR per sqm", xlab = "Time", xaxt = "n")
lines(porto_ts_test, col = "red", lwd = 2)
axis(1, at = years_fc, labels = years_fc)
dev.off()

# Method 3: ETS
ets_model <- ets(porto_ts_train)
ets_fc <- forecast(ets_model, h = h_test)
capture.output(summary(ets_model), file = file.path(out_test, "Test5_ets_summary.txt"))
png(file.path(out_fig, "Fig13_forecast_ets.png"), width = 800, height = 600)
plot(ets_fc, main = paste("ETS (", ets_model$method, ") vs Test Set", sep=""), ylab = "EUR per sqm", xlab = "Time", xaxt = "n")
lines(porto_ts_test, col = "red", lwd = 2)
axis(1, at = years_fc, labels = years_fc)
dev.off()

# Method 4: STLF
stlf_fc <- stlf(porto_ts_train, h = h_test)
capture.output(summary(stlf_fc$model), file = file.path(out_test, "Test6_stlf_summary.txt"))
png(file.path(out_fig, "Fig14_forecast_stlf.png"), width = 800, height = 600)
plot(stlf_fc, main = "STLF vs Test Set", ylab = "EUR per sqm", xlab = "Time", xaxt = "n")
lines(porto_ts_test, col = "red", lwd = 2)
axis(1, at = years_fc, labels = years_fc)
dev.off()

# 5. Accuracy Comparison Table
cmp_smoothing <- rbind(
  HW_Multiplicative = test_accuracy_row(hw_mult, porto_ts_test),
  HW_Damped         = test_accuracy_row(hw_damped, porto_ts_test),
  ETS_Automated     = test_accuracy_row(ets_fc, porto_ts_test),
  STLF_Loess        = test_accuracy_row(stlf_fc, porto_ts_test)
)
cmp_smoothing_sorted <- cmp_smoothing[order(cmp_smoothing[, "RMSE"]), , drop = FALSE]

write.csv(round(cmp_smoothing_sorted, 4), file.path(out_tab, "Tab1_smoothing_accuracy.csv"), row.names = TRUE)