# Source: day_5-1_Debarshi.pdf
#RNA-Seq Gene Expression Normalization Using RPKM and TPM and HeatmapVisualization
counts <- read.table("raw_count.xlsx", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
anno <- read.table("Scer_annotation.txt", header = FALSE, sep = "\t",stringsAsFactors = FALSE)
head(anno)
colnames(anno) <- c("Chr", "Feature", "Start", "End", "Strand", "Gene")
anno$Length <- anno$End - anno$Start + 1
gene_length <- anno[, c("Gene", "Length")]
head(anno)
data <- merge(counts, gene_length, by.x = "Genes", by.y = "Gene")
head(data)
library_size <- colSums(data[,2:5])
library_size
#Calculated RPKM
data$RPKM_ConA_R1 <- (data$ConA_R1 * 1e9) / (data$Length *library_size["ConA_R1"])
data$RPKM_ConA_R2 <- (data$ConA_R2 * 1e9) / (data$Length * library_size["ConA_R2"])
data$RPKM_ConB_R1 <- (data$ConB_R1 * 1e9) / (data$Length * library_size["ConB_R1"])
data$RPKM_ConB_R2 <- (data$ConB_R2 * 1e9) / (data$Length * library_size["ConB_R2"])
data$Length_kb <- data$Length/1000
head(data$RPKM_ConA_R1)
88
# Calculated RPK
data$RPK_ConA_R1 <- data$ConA_R1 / data$Length_kb
data$RPK_ConA_R2 <- data$ConA_R2 / data$Length_kb
data$RPK_ConB_R1 <- data$ConB_R1 / data$Length_kb
data$RPK_ConB_R2 <- data$ConB_R2 / data$Length_kb
# Calculated Sums
sumA1 <- sum(data$RPK_ConA_R1)
sumA2 <- sum(data$RPK_ConA_R2)
sumB1 <- sum(data$RPK_ConB_R1)
sumB2 <- sum(data$RPK_ConB_R2)
#Calculated TPM
data$TPM_ConA_R1 <- (data$RPK_ConA_R1 / sumA1) * 1e6
data$TPM_ConA_R2 <- (data$RPK_ConA_R2 / sumA2) * 1e6
data$TPM_ConB_R1 <- (data$RPK_ConB_R1 / sumB1) * 1e6
data$TPM_ConB_R2 <- (data$RPK_ConB_R2 / sumB2) * 1e6
#Saved the results
write.csv(data, "RPKM_TPM_results.csv", row.names=FALSE)
colSums(data[, c("RPKM_ConA_R1", "RPKM_ConA_R2", "RPKM_ConB_R1", "RPKM_ConB_R2")])
colSums(data[, c("TPM_ConA_R1", "TPM_ConA_R2", "TPM_ConB_R1", "TPM_ConB_R2")])
# 1. Extract TPM columns
expr <- data[, c("TPM_ConA_R1", "TPM_ConA_R2", "TPM_ConB_R1", "TPM_ConB_R2")]
# 2. Set gene names as matrix row names
rownames(expr) <- data$Genes
# 3. Convert to numeric matrix
expr <- as.matrix(expr)
mode(expr) <- "numeric"
# 4. Log2 transformation with pseudo-count
expr_log <- log2(expr + 1)
head(expr_log)
#the heatmap
heatmap(expr_log, Rowv = NA, Colv = NA, scale = "none", col = heat.colors(100))
