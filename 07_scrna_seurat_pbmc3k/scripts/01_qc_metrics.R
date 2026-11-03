# Source PDF: Seurat_Debarshi-1.pdf

# Mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Show QC metrics for the first 5 cells (was missing)
head(pbmc@meta.data, 5)

# Violin plot of QC metrics
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
