#!/usr/bin/env Rscript
# =====================================================
# FMTS Porto property prices — unified analysis pipeline
# Run: Rscript run_analysis.R
# =====================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(forecast)
})

source("R/utils.R")

CFG <- init_project()
assign("CFG", CFG, envir = .GlobalEnv)

porto_prices <- load_porto_prices(CFG$csv_path)
porto_ts <- make_porto_ts(porto_prices)
split <- split_train_test(porto_ts, H_TEST)

train   <- split$train
test    <- split$test
h_test  <- split$h

assign("porto_prices", porto_prices, envir = .GlobalEnv)
assign("porto_ts", porto_ts, envir = .GlobalEnv)
assign("train", train, envir = .GlobalEnv)
assign("test", test, envir = .GlobalEnv)
assign("h_test", h_test, envir = .GlobalEnv)

cat("Project root:", CFG$root, "\n")
cat("Data:", CFG$csv_path, "\n\n")

source(file.path(CFG$root, "scripts", "01_eda.R"))
source(file.path(CFG$root, "scripts", "02_decomposition_smoothing.R"))
source(file.path(CFG$root, "scripts", "03_arima.R"))
source(file.path(CFG$root, "scripts", "04_comparison_oos.R"))
