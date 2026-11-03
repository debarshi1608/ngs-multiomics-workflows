# Source: Debarshi_15_8.pdf (superset of Debarshi_13_08_part2.pdf; adds GO enrichment)
library(Rsubread)
setwd("/home/ibab/bio-env")
gtf_file <- "GCF_000146045.2_R64_genomic.gtf"
out_FeatCount_conA_rep1 <- featureCounts(
files = "conA_rep1.sam",
annot.ext = gtf_file,
isGTFAnnotationFile = TRUE,
GTF.featureType = "exon",
GTF.attrType = "gene_id"
)
out_FeatCount_conA_rep2 <- featureCounts(
files = "conA_rep2.sam",
annot.ext = gtf_file,
isGTFAnnotationFile = TRUE,
GTF.featureType = "exon",
GTF.attrType = "gene_id"
)
out_FeatCount_conB_rep1 <- featureCounts(
files = "conB_rep1.sam",
annot.ext = gtf_file,
isGTFAnnotationFile = TRUE,
GTF.featureType = "exon",
GTF.attrType = "gene_id"
)
out_FeatCount_conB_rep2 <- featureCounts(
files = "conB_rep2.sam",
annot.ext = gtf_file,
isGTFAnnotationFile = TRUE,
GTF.featureType = "exon",
GTF.attrType = "gene_id"
)
out_FeatCount_conA_rep1_count = data.frame(out_FeatCount_conA_rep1$counts)
out_FeatCount_conA_rep2_count = data.frame(out_FeatCount_conA_rep2$counts)
out_FeatCount_conB_rep1_count = data.frame(out_FeatCount_conB_rep1$counts)
out_FeatCount_conB_rep2_count = data.frame(out_FeatCount_conB_rep2$counts)
count_matrix <- cbind(
out_FeatCount_conA_rep1_count,
out_FeatCount_conA_rep2_count,
out_FeatCount_conB_rep1_count,
out_FeatCount_conB_rep2_count
)
write.table(count_matrix, file = "count_matrix.txt",
sep="\t")
library(edgeR)
group=c(rep("conA",2),rep("conB",2))
y=DGEList(counts=count_matrix,group=group)
y
# Filter lowly expressed genes #filterByExpr() identifies genes with sufficient expression to be reliably tested. Applies the cpm and gives cutoff
keep=filterByExpr(y)
head(keep)
# how many genes survived
sum(keep)
# to see the total number of genes before filtering.
length(keep)
library(edgeR)
# 1. Ensure design matrix is defined first
design1 <- model.matrix(~0 + group)
# 2. Recalculate 'keep' using filterByExpr on object 'y'
keep <- filterByExpr(y, design = design1)
# 3. Filter the DGEList (lengths will now match)
y <- y[keep, , keep.lib.sizes = FALSE]
dim(y)
# 4. Calculate normalization factors (TMM)
y <- calcNormFactors(y)
# 5. Estimate dispersion
y <- estimateDisp(y, design1)
# 6. Test for differential expression
# Note: exactTest requires treatment group factors; ensure group is set on y$samples$group
DEet <- exactTest(y, pair = c("conA", "conB"))
# View top results
head(DEet$table)
#Identifying the up and down regulated gene at least 2-fold higher expression and strong statistical significance
et_table = DEet$table
up_et = subset(et_table,et_table$logFC >= 1 & et_table$PValue <= 0.001)
down_et = subset(et_table,et_table$logFC <= -1 & et_table$PValue <= 0.001)
# up regulated gene
head(up_et)
# Downregulated genes
head(down_et)
#No of gene upregulated and down regulated
nrow(up_et)
nrow(down_et)
library(ggplot2)
# Create a color column
et_table$color <- "grey"
# Upregulated genes
et_table$color[et_table$logFC >= 1 & et_table$PValue <= 0.001] <- "red"
# Downregulated genes
et_table$color[et_table$logFC <= -1 & et_table$PValue <= 0.001] <- "blue"
# Calculate -log10(P-value)
et_table$negLogP <- -log10(et_table$PValue)
# Volcano plot
plot(
et_table$logFC,
et_table$negLogP,
pch = 16,
cex = 0.6,
col = et_table$color,
xlab = "Log2 Fold Change",
ylab = "-Log10 P-value",
main = "Volcano Plot: conB vs conA"
)
# Vertical cutoff lines
abline(
v = c(-1, 1),
lty = 2,
col = "black"
)
# Horizontal cutoff line
abline(
h = -log10(0.001),
lty = 2,
col = "black"
)
# RPKM_Normalization for raw read counts
gene_length=out_FeatCount_conA_rep1$annotation$Length
head(gene_length)
gene_id <- out_FeatCount_conA_rep1$annotation$GeneID
head(gene_id)
gene_info <- data.frame(Geneid = gene_id,Length = gene_length)
head(gene_info)
gene_info <- gene_info[match(rownames(count_matrix), gene_info$Geneid),]
count_data <- cbind(
Geneid = rownames(count_matrix),
Length = gene_info$Length,
count_matrix,
row.names = NULL
)
write.table(count_data, file = "count_data.txt",
sep="\t")
dim(count_data)
library_size <- colSums(count_data[, 3:6])
library_size
scaling_factor <- library_size / 1000000
scaling_factor
gene_length_kb <- count_data$Length / 1000
#Calculate RPKM
RPKM_conA_rep1<-count_data$conA_rep1.sam/(gene_length_kb*scaling_factor[1])
RPKM_conA_rep2<-count_data$conA_rep2.sam/(gene_length_kb*scaling_factor[2])
RPKM_conB_rep1<-count_data$conB_rep1.sam/(gene_length_kb*scaling_factor[3])
RPKM_conB_rep2<-count_data$conB_rep2.sam/(gene_length_kb*scaling_factor[4])
#Create the RPKM matrix
rpkm <- data.frame(
Geneid = count_data$Geneid,
RPKM_conA_rep1 = RPKM_conA_rep1,
RPKM_conA_rep2 = RPKM_conA_rep2,
RPKM_conB_rep1 = RPKM_conB_rep1,
RPKM_conB_rep2 = RPKM_conB_rep2
)
rownames(rpkm) <- rpkm$Geneid
rpkm$Geneid <- NULL
#Log transform for heatmap
log_rpkm <- log2(rpkm + 1)
#For DEG list perform RPKM Normalization
up_genes <- rownames(up_et)
down_genes <- rownames(down_et)
deg_genes <- c(up_genes, down_genes)
length(deg_genes)
rpkm_deg <- rpkm[deg_genes, ]
dim(rpkm_deg)
log_rpkm_deg <- log2(rpkm_deg + 1)
library(pheatmap)
pheatmap(
log_rpkm_deg,
scale = "row",
cluster_rows = TRUE,
cluster_cols = TRUE,
main = "RPKM Heatmap of Differentially Expressed Genes"
)
# Quasi-likelihood method
group=c(rep("conA",2),rep("conB",2))
group
#DGEList is an edgeR object that stores your RNA-seq count data along with informationabout the experim
y=DGEList(counts=count_matrix,group=group)
keep=filterByExpr(y)
y=y[keep,,keep.lib.sizes=FALSE]
dim(y)
y <- calcNormFactors(y)
design=model.matrix(~0+group)
design
# This fits a quasi-likelihood negative binomial generalized linear model to your RNA-seq count data.
y_glm=glmQLFit(y,design)
BvsA=makeContrasts(groupconB - groupconA, levels=design)
#This performs the actual quasi-likelihood F-test for the contrast:
res = glmQLFTest(y_glm, contrast = BvsA)
#topTags(res)
Ftest_table <- res$table
up_Ftest <- subset(Ftest_table,logFC >= 1 & PValue <= 0.01)
down_Ftest <- subset(Ftest_table,logFC <= -1 & PValue <= 0.01)
nrow(up_Ftest)
nrow(down_Ftest)
Ftest_table$color <- "grey"
# Upregulated genes
Ftest_table$color[Ftest_table$logFC >= 1 & Ftest_table$PValue <= 0.01] <- "red"
# Downregulated genes
Ftest_table$color[Ftest_table$logFC <= -1 & Ftest_table$PValue <= 0.01] <- "blue"
# Calculate -log10(P-value)
Ftest_table$negLogP <- -log10(Ftest_table$PValue)
# Volcano plot
plot(
Ftest_table$logFC,
Ftest_table$negLogP,
pch = 16,
cex = 0.6,
col = Ftest_table$color,
xlab = "Log2 Fold Change",
ylab = "-Log10 P-value",
main = "Volcano Plot: conB vs conA"
)
# Vertical cutoff lines
abline(
v = c(-1, 1),
lty = 2,
col = "black"
)
# Horizontal cutoff line
abline(
h = -log10(0.001),
lty = 2,
col = "black"
)
#Load the gprofiler2 package for Gene Ontology (GO) enrichment analysis
library(gprofiler2)
#Convert the upregulated gene table into a data frame
up_classic=data.frame(up_et)
#Add the gene IDs from the row names as a separate column
up_classic$Genes=rownames(up_classic)
head(up_classic$Genes)
# Remove the row names from the data frame
row.names(up_classic)=NULL
# Convert the downregulated gene table into a data frame
down_classic=data.frame(down_et)
# Add the gene IDs from the row names as a separate column
down_classic$Genes=rownames(down_classic)
head(down_classic$Genes)
#Remove the row names from the data frame
row.names(down_classic)=NULL
write.table(up_classic, file="UP_regulated.csv",row.names = FALSE)
# We have to perform GO enrichment analysis for the upregulated genes organism = “scerevisiae” specifies Saccharomyces cerevisiae (yeast) Bonferroni correction is used to correct for multiple testing
# Run functional profiling on upregulated genes
res_up <- gost(
query = up_classic$Genes,
organism = "scerevisiae",
exclude_iea = TRUE,
correction_method = "bonferroni"
)
# View the top enrichment results
head(res_up$result)
# We have to extract the GO enrichment results into a data frame
df_GO=res_up$result
head(df_GO)
# Select the first 10 GO terms from the enrichment results
go_up_top10 <- df_GO[1:10, ]
head(go_up_top10)
# We have to convert GO P-values to -log10(P-value) Smaller P-values will have larger -log10(P-value) values
go_up_top10$negLogP <- -log10(go_up_top10$p_value)
#Sort the GO terms according to their P-value .The most statistically significant termscome first
go_up_top10 <- go_up_top10[order(go_up_top10$p_value),]
#Select top 10 after sorting
top_10 <- go_up_top10[1:10, ]
# Rename the long GO term for better visualization
top_10$term_name[top_10$term_name ==
"maturation of SSU-rRNA from tricistronic rRNA transcript (SSU-rRNA, 5.8S rRNA, LSU-rRNA)"] <-
"SSU-rRNA transcript maturation"
# We have to plot showing the top 10 enriched GO terms
ggplot(top_10, aes(
x = negLogP,
y = reorder(term_name, negLogP)
)) +
geom_col(
fill = "pink",
width = 0.75
) +
geom_point(
aes(size = intersection_size),
shape = 21,
fill = "blue",
color = "blue"
) +
labs(
title = "Top 10 GO Enriched Biological Processes",
x = "-log10(P-value)",
y = "GO Biological Process",
size = "Gene Count"
) +
theme_bw()
# Figure Description & Interpretation
# The x-axis displays statistical significance using −log10​
(P-value), meaning points further to the right carry higher statistical confidence. Pink bars indicate the overallenrichment level for each GO term, while the blue circles show the number of upregulated genes mapped to that term, where larger circles represent a higher gene count.
# Functional enrichment analysis reveals that genes upregulated in conB relative to conA are heavily focused on ribosome assembly and RNA processing. Ribosome biogenesis stands out as the single most enriched biological process, along with ribonucleoprotein complex assembly, rRNA metabolism, and pre-rRNA processing. Both small and large ribosomalsubunit maturation terms are also strongly represented. Overall, these results indicatethat switching from conA to conB primarily drives cellular investment into building ribosomes, processing rRNA, and scaling up the core translational machinery.
