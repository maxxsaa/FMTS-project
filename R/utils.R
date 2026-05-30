# Shared utilities for the FMTS Porto analysis pipeline.
# Portable paths: works from any machine when run via run_analysis.R or Rscript.

H_TEST <- 24L
H_OOS  <- 24L

resolve_project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

init_project <- function(csv_name = "timeSeriesPorto.csv") {
  root <- resolve_project_root()
  csv_path <- file.path(root, csv_name)
  if (!file.exists(csv_path)) {
    stop(
      "Cannot find ", csv_name, " in ", root, ".\n",
      "Place the CSV in the project root or set the working directory there.",
      call. = FALSE
    )
  }
  out_dir <- file.path(root, "outputs")
  out_fig <- file.path(out_dir, "figures")
  out_tab <- file.path(out_dir, "tables")
  out_test <- file.path(out_dir, "tests")
  for (d in c(out_dir, out_fig, out_tab, out_test)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
  list(
    root = root,
    csv_path = normalizePath(csv_path),
    out_dir = out_dir,
    out_fig = out_fig,
    out_tab = out_tab,
    out_test = out_test
  )
}

load_porto_prices <- function(csv_path) {
  raw_data <- readr::read_csv(
    csv_path,
    col_names = "raw_line",
    col_types = "c",
    show_col_types = FALSE
  )
  raw_data %>%
    tidyr::separate(
      raw_line,
      into = c("month_year", "series_type", "region_code", "euros_str", "extra"),
      sep = ";",
      extra = "drop",
      fill = "right"
    ) %>%
    dplyr::mutate(
      month_name = stringr::str_extract(month_year, "^[A-Za-z]+"),
      year = as.numeric(stringr::str_extract(month_year, "\\d{4}$")),
      month_number = match(tolower(month_name), tolower(month.name)),
      euros_per_sqm = as.numeric(gsub(",", ".", euros_str)),
      date = as.Date(paste0(year, "-", sprintf("%02d", month_number), "-01"))
    ) %>%
    dplyr::filter(!is.na(date), !is.na(euros_per_sqm)) %>%
    dplyr::distinct(date, .keep_all = TRUE) %>%
    dplyr::arrange(date) %>%
    dplyr::select(date, euros_per_sqm)
}

make_porto_ts <- function(porto_prices) {
  ts(
    porto_prices$euros_per_sqm,
    start = c(
      as.numeric(format(min(porto_prices$date), "%Y")),
      as.numeric(format(min(porto_prices$date), "%m"))
    ),
    frequency = 12
  )
}

split_train_test <- function(y, h = H_TEST) {
  n <- length(y)
  if (n <= h + 36) h <- max(12L, floor(0.15 * n))
  t_time <- time(y)
  train <- window(y, end = t_time[n - h])
  test  <- window(y, start = t_time[n - h + 1])
  list(train = train, test = test, h = h)
}

test_accuracy_row <- function(fc, actual) {
  a <- forecast::accuracy(fc, actual)
  if ("Test set" %in% rownames(a)) return(a["Test set", , drop = FALSE])
  if (nrow(a) >= 2) return(a[2, , drop = FALSE])
  a[1, , drop = FALSE]
}

# Back-transform log-scale forecast to EUR/m² (point + 95% PI)
exp_forecast <- function(fc_log) {
  structure(
    list(
      mean  = exp(as.numeric(fc_log$mean)),
      lower = exp(as.numeric(fc_log$lower[, 1])),
      upper = exp(as.numeric(fc_log$upper[, 1])),
      level = fc_log$level
    ),
    class = "exp_forecast"
  )
}

extract_test_metrics <- function(fc, actual) {
  if (inherits(fc, "forecast")) {
    a <- test_accuracy_row(fc, actual)
  } else {
    a <- forecast::accuracy(as.numeric(fc), actual)
    if ("Test set" %in% rownames(a)) a <- a["Test set", , drop = FALSE]
    else if (nrow(a) >= 2) a <- a[2, , drop = FALSE]
    else a <- a[1, , drop = FALSE]
  }
  c(
    ME = as.numeric(a[, "ME"]),
    RMSE = as.numeric(a[, "RMSE"]),
    MAE = as.numeric(a[, "MAE"]),
    MAPE = as.numeric(a[, "MAPE"])
  )
}

plot_forecast_vs_test <- function(fc, test, main, file, ylab = "EUR per m²") {
  grDevices::png(file, width = 900, height = 500, res = 120)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(fc, main = main, ylab = ylab, xlab = "Time")
  lines(test, col = "red", lwd = 2)
  legend(
    "topleft",
    legend = c("Forecast", "95% PI", "Actual (test)"),
    col = c("blue", "lightblue", "red"),
    lty = 1,
    lwd = c(2, 6, 2),
    bty = "n"
  )
}

forecast_log_to_levels <- function(fc_log, history_ts) {
  fc <- fc_log
  fc$mean <- ts(exp(fc_log$mean), start = time(fc_log)[1], frequency = frequency(fc_log))
  fc$lower <- ts(exp(fc_log$lower), start = time(fc_log)[1], frequency = frequency(fc_log))
  fc$upper <- ts(exp(fc_log$upper), start = time(fc_log)[1], frequency = frequency(fc_log))
  fc$x <- history_ts
  fc
}

plot_exp_forecast_vs_test <- function(fc_log, history_ts, test, main, file) {
  fc_eur <- forecast_log_to_levels(fc_log, history_ts)
  plot_forecast_vs_test(fc_eur, test, main, file)
}

arima_model_label <- function(m) {
  ord <- forecast::arimaorder(m)
  if (length(ord) < 6L) ord <- c(ord, rep(0L, 6L - length(ord)))
  p <- ord[1]; d <- ord[2]; q <- ord[3]
  P <- ord[4]; D <- ord[5]; Q <- ord[6]
  if (any(c(P, D, Q) > 0, na.rm = TRUE)) {
    sprintf("SARIMA(%d,%d,%d)(%d,%d,%d)[12]", p, d, q, P, D, Q)
  } else {
    sprintf("ARIMA(%d,%d,%d)", p, d, q)
  }
}
