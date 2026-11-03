# Source PDF: Seurat_Debarshi-1.pdf

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
