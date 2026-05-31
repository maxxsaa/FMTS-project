# ==============================================================================
# Script 4: Final Comparison of All Methods (Corrected Labeling)
# ==============================================================================

library(tidyverse)

output_dir <- "/Users/pedromcorreia/Desktop/MPST/outputs"
out_tab <- file.path(output_dir, "tables")

# Function to read a file
load_file <- function(filename) {
  df <- read.csv(file.path(out_tab, filename))
  colnames(df)[1] <- "Method" # Ensure first col is Method
  return(df)
}

# 1. Load data
df_smoothing <- load_file("Tab1_smoothing_accuracy.csv")
df_bic       <- load_file("Tab3_BIC_test_accuracy.csv")
df_aic       <- load_file("Tab3_AIC_test_accuracy.csv")

# 2. Apply labels BEFORE binding
df_smoothing$Method <- c("HW Multiplicative", "HW Damped", "ETS Automated", "STLF Loess")
df_bic$Method       <- "ARIMA BIC Winner"
df_aic$Method       <- "ARIMA AIC Winner"

# 3. Combine: bind_rows() now sees the updated labels
master_comparison <- bind_rows(df_smoothing, df_bic, df_aic)

# 4. Sort by RMSE
master_comparison_sorted <- master_comparison[order(master_comparison$RMSE), ]

# 5. Round ONLY numeric columns
numeric_cols <- sapply(master_comparison_sorted, is.numeric)
master_comparison_sorted[, numeric_cols] <- round(master_comparison_sorted[, numeric_cols], 4)

# 6. Save
write.csv(master_comparison_sorted, file.path(out_tab, "Tab4_Final_Master_Comparison.csv"), row.names = FALSE)

# Print for your review
print(master_comparison_sorted)