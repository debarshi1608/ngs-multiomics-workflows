# Source: diffbind_analysis_Debarshi_011.pdf  (run inside `pixi run R`; DiffBind differential binding WT vs Mutant)
library(DiffBind)

tam = dba(sampleSheet="sample_sheet.csv")
tam <- dba.count(tam)
info <- dba.show(tam)
libsizes <- cbind(
  LibReads = info$Reads,
  FRiP = info$FRiP,
  PeakReads = round(info$Reads * info$FRiP)
)
rownames(libsizes) <- info$ID
libsizes
plot(tam)

# 1. Normalize the count data (default: sequencing depth)
tam <- dba.normalize(tam)
# 2. Retrieve normalization details
norm <- dba.normalize(tam, bRetrieve = TRUE)
norm
# 3. View table of raw library sizes vs. normalized library sizes
info <- dba.show(tam)
normlibs <- cbind(
  FullLibSize = norm$lib.sizes,
  NormFacs = round(norm$norm.factors, 4),
  NormLibSize = round(norm$lib.sizes / norm$norm.factors)
)
rownames(normlibs) <- info$ID
normlibs

# Establishing contrast
tam <- dba.contrast(
  tam,
  categories = DBA_CONDITION,
  minMembers = 2
)
tam

# ------------------------------------------------------------
# 5. Differential binding analysis
# ------------------------------------------------------------
tam_analyzed <- dba.analyze(
  tam,
  method = DBA_DESEQ2,
  bParallel = FALSE
)
# Check the results summary
dba.show(tam_analyzed, bContrasts = TRUE)
tam_analyzed$contrast

# ------------------------------------------------------------
# 6. Export all regions (significant + non-significant)
# ------------------------------------------------------------
if (!dir.exists("DiffBind_results")) {
  dir.create("DiffBind_results")
}
res <- dba.report(
  tam_analyzed,
  th = 1,
  bCounts = TRUE
)
res_df <- as.data.frame(res)
write.csv(
  res_df,
  "DiffBind_results/DiffBind_all_regions.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Significant regions (FDR < 0.05)
# ------------------------------------------------------------
res_sig <- res_df[
  !is.na(res_df$FDR) &
  res_df$FDR < 0.05,
]
write.csv(
  res_sig,
  "DiffBind_results/DiffBind_significant_regions_FDR0.05.csv",
  row.names = FALSE
)
cat("Total analyzed regions: ", nrow(res_df), "\n")
cat("Significant regions (FDR < 0.05): ", nrow(res_sig), "\n")
cat(" - Enriched in Mutant (Fold > 0): ", sum(res_sig$Fold > 0), "\n")
cat(" - Enriched in WT (Fold < 0): ", sum(res_sig$Fold < 0), "\n")

# 8. Classify significant regions
res_sig$Direction <- ifelse(
  res_sig$Fold > 0,
  "dGcn4_enriched",
  "WT_enriched"
)
# NOTE: the file name in the next write.csv() is cut off in the source PDF; reconstructed from the other output names
write.csv(
  res_sig,
  "DiffBind_results/DiffBind_significant_regions_with_direction.csv",
  row.names=FALSE
)

# ------------------------------------------------------------
# 9. Summary & Save Object
# ------------------------------------------------------------
cat("\n====================================\n")
cat("DIFFERENTIAL BINDING SUMMARY\n")
cat("====================================\n")
cat("Total regions:", nrow(res_df), "\n")
cat("Significant regions FDR < 0.05:", nrow(res_sig), "\n")
cat("dGcn4-enriched:", sum(res_sig$Direction == "dGcn4_enriched"), "\n")
cat("WT-enriched:", sum(res_sig$Direction == "WT_enriched"), "\n")
saveRDS(
  tam_analyzed,
  "DiffBind_results/DiffBind_final.rds"
)

### Peak length distributions
library(ggplot2)
library(dplyr)
peak_files <- list(
  "WT_Rep1" = "WT_Rep1_narrow_peaks.narrowPeak",
  "WT_Rep2" = "WT_Rep2_narrow_peaks.narrowPeak",
  "Mutant_Rep1" = "Mutant_Rep1_narrow_peaks.narrowPeak",
  "Mutant_Rep2" = "Mutant_Rep2_narrow_peaks.narrowPeak"
)
# Function to read a narrowPeak file and calculate length
read_peak_length <- function(filepath, sample_name) {
  df <- read.table(filepath, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  # Col 2: chromStart, Col 3: chromEnd
  data.frame(
    Sample = sample_name,
    Condition = ifelse(grepl("WT", sample_name), "WT", "Mutant"),
    Length = df$V3 - df$V2
  )
}
all_peak_lengths <- bind_rows(
  lapply(names(peak_files), function(nm) {
    read_peak_length(peak_files[[nm]], nm)
  })
)
all_peak_lengths %>%
  group_by(Sample) %>%
  summarise(
    Peak_Count = n(),
    Min_Length = min(Length),
    Median_Length = median(Length),
    Mean_Length = round(mean(Length), 1),
    Max_Length = max(Length)
  )

# ------------------------------------------------------------
# Convert significant regions from CSV to BED (input for HOMER)
# ------------------------------------------------------------
res_sig <- read.csv("DiffBind_results/DiffBind_significant_regions_FDR0.05.csv")
write.table(
  res_sig[, c("seqnames", "start", "end")],
  file = "DiffBind_significant_regions_FDR0.05.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
cat("BED file successfully generated!\n")

# ------------------------------------------------------------
# Memory-safe re-run used later (export_results.R), recoverable part only.
# The remainder of this heredoc is garbled in the source PDF.
# ------------------------------------------------------------
library(BiocParallel)
register(SerialParam())  # single-core processing to prevent RAM spikes
output_dir <- "DiffBind_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
cat("\n[1/5] Loading sample sheet...\n")
tam <- dba(sampleSheet = "sample_sheet.csv")
cat("[2/5] Counting reads safely...\n")
tam <- dba.count(tam, bParallel = FALSE)
cat("[3/5] Normalizing data and establishing contrast...\n")
tam <- dba.normalize(tam)
tam <- dba.contrast(tam, categories = DBA_CONDITION, minMembers = 2)
cat("[4/5] Running differential binding analysis...\n")
tam_analyzed <- dba.analyze(tam, method = DBA_DESEQ2, bParallel = FALSE)
