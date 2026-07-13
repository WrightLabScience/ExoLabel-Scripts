# This script calculates and plots the absolute and percentage deviations
# between original edge weights and post-compression edge weights in ExoLabel.

# Detect directory and define paths
SOURCE_DIR <- "."

pdf_path <- file.path(SOURCE_DIR, "DeviationsPlot.pdf")

# ExoLabel Compression Constants
BITS_FOR_WEIGHT <- 16
BITS_FOR_EXP <- 4
MAX_POSSIBLE_WEIGHT <- 2^BITS_FOR_WEIGHT - 1     # 65535
MAX_POSSIBLE_EXPONENT <- 2^BITS_FOR_EXP - 1       # 15

# Helper: Index of most significant bit (0-indexed)
get_msb32 <- function(v) {
  msb <- integer(length(v))
  valid <- v > 0
  msb[valid] <- as.integer(floor(log2(v[valid])))
  return(msb)
}

# R version of compressEdgeValues (weight part only)
compress_weight <- function(weight) {
  w_prime <- log2(weight + 1)
  w_prime <- pmin(w_prime, MAX_POSSIBLE_WEIGHT)

  floor_w <- floor(w_prime)
  msb <- get_msb32(floor_w)

  to_shift <- MAX_POSSIBLE_EXPONENT - msb
  to_shift <- pmax(to_shift, 0)

  scaled_w <- w_prime * (2^to_shift)
  w_comp <- as.numeric(floor(scaled_w))
  w_comp <- bitwAnd(as.integer(w_comp), MAX_POSSIBLE_WEIGHT)

  to_shift_comp <- bitwShiftL(as.integer(to_shift), BITS_FOR_WEIGHT)
  compressed <- bitwOr(to_shift_comp, w_comp)
  return(compressed)
}

# R version of decompressEdgeValue (weight part only)
decompress_weight <- function(compressed) {
  w_comp <- bitwAnd(compressed, MAX_POSSIBLE_WEIGHT)
  bits_shifted <- bitwAnd(bitwShiftR(compressed, BITS_FOR_WEIGHT), MAX_POSSIBLE_EXPONENT)

  decompressed_log_w <- w_comp / (2^bits_shifted)
  decompressed_w <- 2^decompressed_log_w - 1
  return(decompressed_w)
}

# 1. Generate log-spaced weights from 1e-5 to 10000 (2000 points to keep the plot clean)
cat("Generating weights and simulating compression...\n")
n_points <- 2000
weights <- 10^seq(-5, 4, length.out = n_points)

# 2. Run compression and decompression
comp <- compress_weight(weights)
decomp <- decompress_weight(comp)

# 3. Calculate absolute and percentage deviations
abs_deviation <- abs(decomp - weights)
pct_deviation <- (abs_deviation / weights) * 100

# 4. Generate the plots
cat("Creating plots...\n")
pdf(pdf_path, width=8, height=8, onefile=TRUE)

# Set grid parameters for a 2-panel vertical layout
par(mfrow=c(2, 1), mar=c(0.5, 4.5, 2.5, 1.5), las=1, mgp=c(2.5, 0.7, 0), tck=-0.015)

# --- PANEL A: ABSOLUTE DEVIATION (LOG-LOG SCALE) ---
# Y-limit expanded down to 1e-9 so no positive error points are clipped offscreen
plot(NULL, xlim=c(1e-5, 10000), ylim=c(1e-9, 10), log="xy",
     xlab="", ylab="Raw Deviation from Original",
     main="ExoLabel Weight Compression Error", axes=FALSE)

# Custom log-spaced axes
x_at <- 10^(-5:4)
axis(1, at=x_at, labels=FALSE, lwd=0, lwd.ticks=0.5, col.ticks="gray60")
y_at_abs <- 10^(-9:1)
axis(2, at=y_at_abs, labels=parse(text=paste0("10^", -9:1)), lwd=0, lwd.ticks=0.5, col.ticks="gray60")

# Draw background grid
abline(h=y_at_abs, col="gray90", lty="solid", lwd=0.5)
abline(v=x_at, col="gray90", lty="solid", lwd=0.5)

# Add line plot (only plot positive values to avoid log(0) warnings/issues)
valid_abs <- abs_deviation > 0
lines(weights[valid_abs], abs_deviation[valid_abs], col="#0E68C5", lwd=1.2)
box(col="gray80", lwd=0.8)
# Add horizontal dashed line at 0.025%
abline(h=0.000025, col="red", lty="dashed", lwd=1.2)
text(x=1e3, y=0.000025, labels=expression(2.5 %*% 10^-5 ~ "error"), col="red", font=1, cex=0.8,
     adj=c(0.1,-0.5))

#legend("bottomright", legend="Absolute Deviation", col="#0E68C5", lty=1, bty="n", cex=0.85, text.col="gray30")


# --- PANEL B: PERCENTAGE DEVIATION (LOG-LOG SCALE) ---
# Y-limit expanded down to 1e-6 so no positive error points are clipped offscreen
par(mar=c(4.5, 4.5, 0.5, 1.5))
plot(NULL, xlim=c(1e-5, 10000), ylim=c(1e-6, 150), log="xy",
     xlab="Original Weight", ylab="Percentage Deviation from Original",
     main="", axes=FALSE)

# Custom axes
axis(1, at=x_at, labels=parse(text=paste0("10^", -5:4)), lwd=0, lwd.ticks=0.5, col.ticks="gray60")
y_at_pct <- 10^(-6:2)
axis(2, at=y_at_pct, labels=parse(text=paste0("10^", -6:2)), lwd=0, lwd.ticks=0.5, col.ticks="gray60")

# Draw background grid
abline(h=y_at_pct, col="gray90", lty="solid", lwd=0.5)
abline(v=x_at, col="gray90", lty="solid", lwd=0.5)

# Add horizontal dashed line at 0.025%
abline(h=0.02, col="red", lty="dashed", lwd=1.2)
text(x=1e3, y=0.02, labels="0.02% error", col="red", font=1, cex=0.8,
     adj=c(0,-0.5))

# Add line plot
valid_pct <- pct_deviation > 0
lines(weights[valid_pct], pct_deviation[valid_pct], col="#824484", lwd=1.2)
box(col="gray80", lwd=0.8)

dev.off()

cat(paste("Successfully saved plot to:", pdf_path, "\n"))
