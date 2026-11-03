# Source PDF: Debarshi_day_4_011.pdf (superset of day_2_Debarshi_Adhikary.pdf and day_3_Debarshi.pdf)
# Stage 1: read CEL files, raw PM/MM inspection, log transform, RMA normalisation, QC plots

# Define the directory path containing our CEL files
celpath <- "/home/ibab/bio-env"

# Load packages for microarray processing
library(affy)

library(limma)

library(tidyverse)

# Read raw microarray data from the specified directory
yeast_microchip <- ReadAffy(celfile.path = celpath, compress = TRUE)

# Extract gene expression matrix into a standard data frame
df_yeast <- data.frame(exprs(yeast_microchip))

# Preview the top rows of the expression matrix
head(df_yeast)

# Extract Perfect Match (PM) probe signals
pm_yeast <- data.frame(affy::pm(yeast_microchip))

# Extract Mismatch (MM) probe signals
mm_yeast <- data.frame(affy::mm(yeast_microchip))

# Scatter plot of PM intensities for the first Wild Type sample
plot(pm_yeast$EWE11465_0001_WT1.CEL, col = "green", main = "PM vs MM Signals")
# Overlay MM intensities on the same plot
points(mm_yeast$EWE11465_0001_WT1.CEL, col = "pink")

## Performed log-transformation on PM and MM intensity data

pmlog <- log(pm_yeast)
mmlog <- log(mm_yeast)

# Set up layout for two side-by-side plots
par(mfrow = c(1, 2))

# Histogram of raw PM values
hist(pm_yeast$EWE11465_0001_WT1.CEL, main = "Raw PM Data", xlab = "Intensity", col = "red")

# Histogram of log-transformed PM values
hist(pmlog$EWE11465_0001_WT1.CEL, main = "Log-transformed PM Data", xlab = "Log Intensity", col = "blue")

# Extract individual sample log-PM vectors
wt1_log <- pmlog$EWE11465_0001_WT1.CEL
wt2_log <- pmlog$EWE11465_0002_WT2.CEL
wt3_log <- pmlog$EWE11465_0003_WT3.CEL
drna1_log <- pmlog$EWE11465_0004_drna1.CEL
drna2_log <- pmlog$EWE11465_0005_drna2.CEL
drna3_log <- pmlog$EWE11465_0006_drna3.CEL

# Plot boxplots across samples
boxplot(wt1_log, wt2_log, wt3_log, drna1_log, drna2_log, drna3_log,
        names = c("WT1", "WT2", "WT3", "dRNA1", "dRNA2", "dRNA3"),main = "Log Transformed Dataset")

##Applied RMA

yeast_microchip_rma <- rma(yeast_microchip)

##Extraced gene expression matrix into a standard data frame

df_rma_yeast <- data.frame(exprs(yeast_microchip_rma))

#the top rows of the expression matrix
head(df_rma_yeast)

##Stored individual sample vectors from the RMA-normalized data frame

wt1_norm <- df_rma_yeast$EWE11465_0001_WT1.CEL
wt2_norm <- df_rma_yeast$EWE11465_0002_WT2.CEL
wt3_norm <- df_rma_yeast$EWE11465_0003_WT3.CEL
drna1_norm <- df_rma_yeast$EWE11465_0004_drna1.CEL
drna2_norm <- df_rma_yeast$EWE11465_0005_drna2.CEL
drna3_norm <- df_rma_yeast$EWE11465_0006_drna3.CEL

##Boxplot of RMA-normalized expression values

boxplot(wt1_norm, wt2_norm, wt3_norm, drna1_norm, drna2_norm, drna3_norm,
        names = c("WT1", "WT2", "WT3", "dRNA1", "dRNA2", "dRNA3"),main = "Normalized Dataset")

##Prepared raw vectors for comparison

wt1_og <- pm_yeast$EWE11465_0001_WT1.CEL
wt2_og <- pm_yeast$EWE11465_0002_WT2.CEL
wt3_og <- pm_yeast$EWE11465_0003_WT3.CEL
drna1_og <- pm_yeast$EWE11465_0004_drna1.CEL
drna2_og <- pm_yeast$EWE11465_0005_drna2.CEL
drna3_og <- pm_yeast$EWE11465_0006_drna3.CEL

# Display three side-by-side boxplots
par(mfrow = c(1, 3))

# 1. Original raw dataset
boxplot(wt1_og, wt2_og, wt3_og, drna1_og, drna2_og, drna3_og,
        names = c("WT1", "WT2", "WT3", "dRNA1", "dRNA2", "dRNA3"),
        ylim = c(0, 500), main = "Original Dataset", las = 2)

# 2. Log-transformed dataset
boxplot(wt1_log, wt2_log, wt3_log, drna1_log, drna2_log, drna3_log,
        names = c("WT1", "WT2", "WT3", "dRNA1", "dRNA2", "dRNA3"),main = "Log Transformed Dataset", las = 2)

# 3. RMA-normalized dataset
boxplot(wt1_norm, wt2_norm, wt3_norm, drna1_norm, drna2_norm, drna3_norm,
        names = c("WT1", "WT2", "WT3", "dRNA1", "dRNA2", "dRNA3"),main = "Normalized Dataset", las = 2)

##compared scatter plots

par(mfrow = c(1, 3))

# Raw scatter plot
plot(wt1_og, drna1_og, main = "Original Dataset", xlim = c(0, 30000))

# Log-transformed scatter plot
plot(wt1_log, drna1_log, main = "Log-transformed Dataset")

# Normalized scatter plot
plot(wt1_norm, drna1_norm, main = "Normalized Dataset")

##Generated MA plots

par(mfrow = c(1, 3))

# MA Calculation: Raw Data
A_og <- (log(wt1_og) + log(drna1_og)) / 2
M_og <- log(wt1_og / drna1_og)
plot(A_og, M_og, main = "Original Dataset", col='brown')
abline(h = 0, col = "red")

# MA Calculation: Log-Transformed Data
A_log <- (wt1_log + drna1_log) / 2
M_log <- wt1_log - drna1_log
plot(A_log, M_log, main = "Log Transformed Dataset", col='green')
abline(h = 0, col = "red")

# MA Calculation: Normalized Data
A_norm <- (wt1_norm + drna1_norm) / 2
M_norm <- wt1_norm - drna1_norm
plot(A_norm, M_norm, main = "Normalized Dataset", col='blue')
abline(h = 0, col = "red")
