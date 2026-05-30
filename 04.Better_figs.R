# ==============================================================================
# Script 4: Combined Zoomed Forecast Figures
# ==============================================================================

# Note: Ensure you have run Script 02 first so the models are in your environment!

# 1. Paths and Figure Settings
output_dir <- "/Users/pedromcorreia/Desktop/MPST/outputs"
out_fig <- file.path(output_dir, "figures")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

png(filename = file.path(out_fig, "Fig_Combined_Smoothing_Zoomed.png"), width = 1200, height = 1000, res = 120)

# Set up the 2x2 grid
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1) + 0.1)

# Define zoom parameters (from 2020 to slightly past the test set)
start_zoom <- 2020.0
end_zoom <- max(time(porto_ts_test))

# --- Panel 1: Holt-Winters (Multiplicative) ---
plot(
  hw_mult, # Corrected variable name
  main = "Holt-Winters (Multiplicative) Zoomed",
  ylab = "EUR per sqm",
  xlab = "Time",
  xaxt = "n",
  xlim = c(start_zoom, end_zoom)
)
lines(porto_ts_test, col = "red", lwd = 2)
axis(1, at = years_fc, labels = years_fc)
legend("topleft", legend = c("Forecast", "95% PI", "Actual (test)"), col = c("blue", "lightblue", "red"), lty = 1, lwd = c(2, 6, 2), bty = "n")

# --- Panel 2: Holt-Winters (Damped) ---
plot(
  hw_damped, # Corrected variable name
  main = "Holt-Winters (Damped) Zoomed",
  ylab = "EUR per sqm",
  xlab = "Time",
  xaxt = "n",
  xlim = c(start_zoom, end_zoom)
)
lines(porto_ts_test, col = "red", lwd = 2)
axis(1, at = years_fc, labels = years_fc)
legend("topleft", legend = c("Forecast", "95% PI", "Actual (test)"), col = c("blue", "lightblue", "red"), lty = 1, lwd = c(2, 6, 2), bty = "n")

# --- Panel 3: ETS ---
plot(
  ets_fc,
  main = paste("ETS (", ets_model$method, ") Zoomed", sep=""),
  ylab = "EUR per sqm",
  xlab = "Time",
  xaxt = "n",
  xlim = c(start_zoom, end_zoom)
)
lines(porto_ts_test, col = "red", lwd = 2)
axis(1, at = years_fc, labels = years_fc)
legend("topleft", legend = c("Forecast", "95% PI", "Actual (test)"), col = c("blue", "lightblue", "red"), lty = 1, lwd = c(2, 6, 2), bty = "n")

# --- Panel 4: STLF ---
plot(
  stlf_fc,
  main = "STLF Zoomed",
  ylab = "EUR per sqm",
  xlab = "Time",
  xaxt = "n",
  xlim = c(start_zoom, end_zoom)
)
lines(porto_ts_test, col = "red", lwd = 2)
axis(1, at = years_fc, labels = years_fc)
legend("topleft", legend = c("Forecast", "95% PI", "Actual (test)"), col = c("blue", "lightblue", "red"), lty = 1, lwd = c(2, 6, 2), bty = "n")

# Close the device and reset layout
dev.off()
par(mfrow = c(1, 1))

# Instead of two separate images (Fig18 and Fig19), we create one combined image:
png(file.path(out_fig, "Fig_Combined_acf_pacf_differenced_combined.png"), width = 1200, height = 600) # Increased width to fit two plots

# Set up a 1x2 grid (1 row, 2 columns)
par(mfrow = c(1, 2))

# Plot 1 (Left Side)
acf(train_diff, lag.max = 36, main = "ACF of Differenced Data")

# Plot 2 (Right Side)
pacf(train_diff, lag.max = 36, main = "PACF of Differenced Data")

# Close the device and reset the layout parameters
dev.off()
par(mfrow = c(1, 1))