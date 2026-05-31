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
  csv_path <- file.path(root, "data", csv_name)
  if (!file.exists(csv_path)) {
    csv_path <- file.path(root, csv_name)
  }
  if (!file.exists(csv_path)) {
    stop(
      "Cannot find ", csv_name, " in data/ or project root.\n",
      "Place the CSV in data/ or set the working directory to the project folder.",
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

plot_forecast_oos <- function(fc, main, file, ylab = "EUR per m²") {
  grDevices::png(file, width = 900, height = 500, res = 120)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(fc, main = main, ylab = ylab, xlab = "Time", shaded = TRUE)
  legend(
    "topleft",
    legend = c("Forecast", "95% PI"),
    col = c("blue", "lightblue"),
    lty = 1,
    lwd = c(2, 6),
    bty = "n"
  )
}

plot_exp_forecast_oos <- function(fc_log, history_ts, main, file) {
  fc_eur <- forecast_log_to_levels(fc_log, history_ts)
  plot_forecast_oos(fc_eur, main, file)
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

# --- STL phase matrices (multiplicative LOESS-style decomposition) ------------

calculate_decomp_metrics <- function(actual, fitted) {
  actual <- as.numeric(actual)
  fitted <- as.numeric(fitted)
  errors <- actual - fitted
  mape <- mean(abs(errors / actual), na.rm = TRUE) * 100
  ss_res <- sum(errors^2, na.rm = TRUE)
  ss_tot <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  cv <- sd(actual, na.rm = TRUE) / mean(actual, na.rm = TRUE) * 100
  c(
    Mean_Error = mean(errors, na.rm = TRUE),
    MAE = mean(abs(errors), na.rm = TRUE),
    MSE = mean(errors^2, na.rm = TRUE),
    `MAPE (%)` = mape,
    Correlation = stats::cor(actual, fitted, use = "complete.obs"),
    r_square = r2,
    STD_Error = stats::sd(errors, na.rm = TRUE),
    `CV(%)` = cv
  )
}

build_year_month_matrix <- function(df, value_col, month_col = "Month_label") {
  wide <- df %>%
    dplyr::select(Year, Month = dplyr::all_of(month_col), value = dplyr::all_of(value_col)) %>%
    tidyr::pivot_wider(names_from = Month, values_from = value)
  month_order <- intersect(month.abb, names(wide))
  wide <- wide %>% dplyr::select(Year, dplyr::all_of(month_order))
  mat <- as.data.frame(wide)
  rownames(mat) <- mat$Year
  mat$Year <- NULL
  mat
}

plot_decomp_heatmap <- function(mat, title, fill_label, mid = NULL) {
  df <- mat %>%
    tibble::rownames_to_column("Year") %>%
    tidyr::pivot_longer(-Year, names_to = "Month", values_to = "value")
  df$Month <- factor(df$Month, levels = month.abb)
  df$Year <- factor(df$Year, levels = sort(unique(as.character(df$Year))))

  fill_scale <- if (is.null(mid)) {
    ggplot2::scale_fill_viridis_c(option = "C", name = fill_label)
  } else {
    ggplot2::scale_fill_gradient2(
      midpoint = mid,
      low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
      name = fill_label
    )
  }

  ggplot2::ggplot(df, ggplot2::aes(x = Month, y = Year, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    fill_scale +
    ggplot2::labs(title = title, x = "Month", y = "Year") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
}

plot_matrix_table_grob <- function(mat, title, digits = 1) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required for STL phase matrix tables.", call. = FALSE)
  }
  mat_num <- as.matrix(mat)
  mat_num[is.na(mat_num)] <- 0
  tab <- format(round(mat_num, digits), nsmall = digits)
  rownames(tab) <- rownames(mat)
  gridExtra::tableGrob(
    tab,
    rows = rownames(tab),
    cols = colnames(tab),
    theme = gridExtra::ttheme_minimal(
      core = list(fg_params = list(cex = 0.55)),
      rowhead = list(fg_params = list(cex = 0.55, fontface = "bold")),
      colhead = list(fg_params = list(cex = 0.55, fontface = "bold"))
    )
  )
}

plot_volatility_table_grob <- function(vol_report) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required for STL phase matrix tables.", call. = FALSE)
  }
  tab <- format(round(as.matrix(vol_report), 6), nsmall = 6, scientific = FALSE)
  gridExtra::tableGrob(
    tab,
    rows = rownames(tab),
    cols = colnames(tab),
    theme = gridExtra::ttheme_minimal(
      core = list(fg_params = list(cex = 0.65)),
      rowhead = list(fg_params = list(cex = 0.65, fontface = "bold")),
      colhead = list(fg_params = list(cex = 0.55, fontface = "bold"))
    )
  )
}

format_volatility_report <- function(vol_df) {
  metrics <- setdiff(colnames(vol_df), "Method")
  mat <- t(as.matrix(vol_df[, metrics, drop = FALSE]))
  colnames(mat) <- vol_df$Method
  rownames(mat) <- metrics
  as.data.frame(mat)
}

print_volatility_report <- function(vol_report) {
  cat("\n--- RELATÓRIO DE VOLATILIDADE POR FASE (LOESS) ---\n")
  idx <- seq_len(ncol(vol_report)) - 1L
  cat(sprintf("%-18s", ""))
  for (i in idx) cat(sprintf("%32d", i))
  cat("\n")
  cat(sprintf("%-18s", "Método"))
  for (nm in colnames(vol_report)) cat(sprintf("%32s", nm))
  cat("\n")
  for (i in seq_len(nrow(vol_report))) {
    cat(sprintf("%-18s", rownames(vol_report)[i]))
    for (j in seq_len(ncol(vol_report))) {
      cat(sprintf("%32.6f", vol_report[i, j]))
    }
    cat("\n")
  }
}

write_volatility_report_txt <- function(vol_report, path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines("--- RELATÓRIO DE VOLATILIDADE POR FASE (LOESS) ---", con)
  idx_line <- paste(c("", as.character(seq_len(ncol(vol_report)) - 1L)), collapse = "\t")
  writeLines(idx_line, con)
  metodo_line <- paste(c("Método", colnames(vol_report)), collapse = "\t")
  writeLines(metodo_line, con)
  for (i in seq_len(nrow(vol_report))) {
    line <- paste(
      c(rownames(vol_report)[i], sprintf("%.6f", unlist(vol_report[i, , drop = TRUE]))),
      collapse = "\t"
    )
    writeLines(line, con)
  }
}

porto_stl_phase_matrices <- function(porto_prices, porto_ts, cfg) {
  serie_ts <- porto_ts
  serie_ts[serie_ts <= 0] <- 0.001

  stl_fit <- stats::stl(log(serie_ts), s.window = "periodic", robust = TRUE)
  trend_real <- exp(stl_fit$time.series[, "trend"])
  seasonal_real <- exp(stl_fit$time.series[, "seasonal"])
  remainder_real <- exp(stl_fit$time.series[, "remainder"])
  fitted_mult <- trend_real * seasonal_real
  resid_real <- as.numeric(serie_ts) - fitted_mult

  dates <- porto_prices$date
  df_base <- tibble::tibble(
    date = dates,
    Raw_Data = as.numeric(serie_ts),
    Trend_LOESS = trend_real,
    Seasonal_LOESS = seasonal_real,
    Remainder_LOESS = remainder_real,
    Resid_LOESS = resid_real,
    Year = lubridate::year(dates),
    Month = lubridate::month(dates),
    Month_label = month.abb[Month]
  )

  phases <- list(
    Phase_1_Stable = list(
      start = "2011-01-01", end = "2015-12-31",
      title = "Phase 1 Stable (2011-2015)",
      label_fig = "Stable"
    ),
    Phase_2_Growth = list(
      start = "2016-01-01", end = "2019-12-31",
      title = "Phase 2 Growth (2016-2019)",
      label_fig = "Growth"
    ),
    Phase_3_Pandemic = list(
      start = "2020-01-01", end = "2021-12-31",
      title = "Phase 3 Pandemic (2020-2021)",
      label_fig = "Pandemic"
    ),
    Phase_4_Acceleration = list(
      start = "2022-01-01", end = "2026-02-01",
      title = "Phase 4 Acceleration (2022-2026)",
      label_fig = "Acceleration"
    )
  )

  phase_details <- list()
  metrics_list <- list()

  for (nm in names(phases)) {
    ph <- phases[[nm]]
    idx <- df_base$date >= as.Date(ph$start) & df_base$date <= as.Date(ph$end)
    sub <- df_base[idx, ]
    sub$Phase <- ph$title
    phase_details[[nm]] <- sub

    act <- sub$Raw_Data
    fit <- sub$Trend_LOESS * sub$Seasonal_LOESS
    metrics_list[[nm]] <- c(Method = ph$title, calculate_decomp_metrics(act, fit))
  }

  df_detail <- dplyr::bind_rows(phase_details)
  vol_df <- as.data.frame(do.call(rbind, metrics_list), stringsAsFactors = FALSE)
  rownames(vol_df) <- names(phases)
  vol_df$Method <- as.character(vol_df$Method)
  num_cols <- setdiff(colnames(vol_df), "Method")
  vol_df[num_cols] <- lapply(vol_df[num_cols], as.numeric)

  mat_trend <- build_year_month_matrix(df_base, "Trend_LOESS")
  mat_seasonal <- build_year_month_matrix(df_base, "Seasonal_LOESS")
  mat_remainder <- build_year_month_matrix(df_base, "Remainder_LOESS")
  mat_resid <- build_year_month_matrix(df_base, "Resid_LOESS")
  vol_report <- format_volatility_report(vol_df)
  fig_labels <- vapply(phases, function(p) p$label_fig, character(1))
  vol_report_fig <- vol_report
  colnames(vol_report_fig) <- fig_labels[match(colnames(vol_report_fig), vol_df$Method)]

  decomp_full <- df_base %>%
    dplyr::transmute(
      date,
      Observed = Raw_Data,
      Trend = Trend_LOESS,
      Seasonal = Seasonal_LOESS,
      Remainder = Remainder_LOESS
    )

  readr::write_csv(decomp_full, file.path(cfg$out_tab, "stl_full_decomposition.csv"))
  readr::write_csv(
    tibble::rownames_to_column(as.data.frame(vol_report), "Metric"),
    file.path(cfg$out_tab, "stl_phase_volatility.csv")
  )
  readr::write_csv(
    tibble::rownames_to_column(as.data.frame(mat_trend), "Year"),
    file.path(cfg$out_tab, "stl_trend_matrix.csv")
  )
  readr::write_csv(
    tibble::rownames_to_column(as.data.frame(mat_seasonal), "Year"),
    file.path(cfg$out_tab, "stl_seasonal_matrix.csv")
  )
  readr::write_csv(
    tibble::rownames_to_column(as.data.frame(mat_remainder), "Year"),
    file.path(cfg$out_tab, "stl_remainder_matrix.csv")
  )
  readr::write_csv(
    tibble::rownames_to_column(as.data.frame(mat_resid), "Year"),
    file.path(cfg$out_tab, "stl_resid_matrix.csv")
  )
  readr::write_csv(df_detail, file.path(cfg$out_tab, "stl_phase_detail.csv"))
  write_volatility_report_txt(
    vol_report,
    file.path(cfg$out_test, "stl_phase_volatility_report.txt")
  )

  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    warning("Install 'gridExtra' to render stl_phase_matrices.png")
    print_volatility_report(vol_report)
    return(invisible(list(
      volatility = vol_report,
      trend = mat_trend,
      seasonal = mat_seasonal,
      resid = mat_resid
    )))
  }

  grDevices::png(
    file.path(cfg$out_fig, "23_stl_volatility_by_phase.png"),
    width = 1400, height = 500, res = 120
  )
  grid::grid.newpage()
  gridExtra::grid.arrange(
    grobs = list(plot_volatility_table_grob(vol_report_fig)),
    top = grid::textGrob(
      "Volatility by phase (STL decomposition)",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  )
  grDevices::dev.off()

  grDevices::png(
    file.path(cfg$out_fig, "23_stl_phase_matrices.png"),
    width = 1800, height = 1400, res = 120
  )
  grid::grid.newpage()
  gridExtra::grid.arrange(
    grobs = list(
      ggplot2::ggplotGrob(plot_decomp_heatmap(mat_trend, "Trend (EUR/m²)", "EUR/m²")),
      ggplot2::ggplotGrob(plot_decomp_heatmap(mat_seasonal, "Seasonal factor", "Factor", mid = 1)),
      ggplot2::ggplotGrob(plot_decomp_heatmap(mat_remainder, "Remainder factor", "Factor", mid = 1)),
      plot_volatility_table_grob(vol_report_fig)
    ),
    ncol = 2,
    top = grid::textGrob(
      "STL multiplicative decomposition — year × month heatmaps (Porto prices)",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  )
  grDevices::dev.off()

  for (item in list(
    list(mat_trend, "Trend (EUR/m²)", "EUR/m²", NULL, "24_stl_trend_heatmap.png"),
    list(mat_seasonal, "Seasonal factor", "Factor", 1, "24_stl_seasonal_heatmap.png"),
    list(mat_remainder, "Remainder factor", "Factor", 1, "24_stl_remainder_heatmap.png")
  )) {
    grDevices::png(file.path(cfg$out_fig, item[[5]]), width = 900, height = 600, res = 120)
    print(plot_decomp_heatmap(item[[1]], item[[2]], item[[3]], mid = item[[4]]))
    grDevices::dev.off()
  }

  print_volatility_report(vol_report)

  invisible(list(
    volatility = vol_report,
    trend = mat_trend,
    seasonal = mat_seasonal,
    remainder = mat_remainder,
    resid = mat_resid,
    detail = df_detail,
    decomposition = decomp_full
  ))
}
