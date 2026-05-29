#!/usr/bin/env Rscript
# Export CSV/LaTeX tables used by report/report.tex

suppressPackageStartupMessages({
  library(tidyverse)
  library(forecast)
})

root <- normalizePath("..", mustWork = TRUE)
csv_path <- file.path(root, "timeSeriesPorto.csv")
out_dir <- file.path(root, "report", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw_data <- read_csv(csv_path, col_names = "raw_line", col_types = "c", show_col_types = FALSE)

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

h_test <- 24L
n_tot <- length(porto_ts)
t_time <- time(porto_ts)
train <- window(porto_ts, end = t_time[n_tot - h_test])
test <- window(porto_ts, start = t_time[n_tot - h_test + 1])

fit_hw <- HoltWinters(train, seasonal = "additive")
fit_ets <- ets(train)
fc_hw <- forecast(fit_hw, h = h_test, level = 95)
fc_ets <- forecast(fit_ets, h = h_test, level = 95)
fc_stlf <- stlf(train, method = "ets", h = h_test, level = 95)

fit_auto <- auto.arima(train, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
fc_sarima <- forecast(fit_auto, h = h_test, level = 95)

test_dates <- seq(max(porto_prices$date) - (h_test - 1) * 31, by = "month", length.out = h_test)
test_dates <- seq.Date(
  as.Date(paste(end(train)[1], round(end(train)[2]), "01", sep = "-")) + 32,
  by = "month",
  length.out = h_test
)

fc_vs_test <- tibble(
  month = format(test_dates, "%Y-%m"),
  actual = as.numeric(test),
  HW = round(as.numeric(fc_hw$mean), 1),
  HW_lo95 = round(as.numeric(fc_hw$lower[, 1]), 1),
  HW_hi95 = round(as.numeric(fc_hw$upper[, 1]), 1),
  ETS = round(as.numeric(fc_ets$mean), 1),
  SARIMA = round(as.numeric(fc_sarima$mean), 1),
  SARIMA_lo95 = round(as.numeric(fc_sarima$lower[, 1]), 1),
  SARIMA_hi95 = round(as.numeric(fc_sarima$upper[, 1]), 1)
)

write_csv(fc_vs_test, file.path(out_dir, "forecast_vs_test.csv"))
write_csv(read_csv(file.path(root, "outputs", "test_accuracy_by_method.csv"), show_col_types = FALSE),
          file.path(out_dir, "test_accuracy_by_method.csv"))

stl_full <- stl(porto_ts, s.window = "periodic", robust = TRUE)
sa_series <- seasadj(stl_full)
sa_tbl <- tibble(
  date = porto_prices$date,
  original = as.numeric(porto_ts),
  seasonally_adjusted = as.numeric(sa_series)
)
write_csv(sa_tbl, file.path(out_dir, "seasonally_adjusted.csv"))

cat("Exported tables to", out_dir, "\n")
