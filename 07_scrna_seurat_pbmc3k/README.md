# Seurat Workflow on PBMC3k: Single-cell RNA-seq from Raw Counts to Labelled Cell Types

**Aim** To analyse a single-cell RNA-seq dataset of about 3,000 human peripheral blood mononuclear cells (PBMC3k) with the Seurat package in R, moving step by step from the raw count matrix to biologically labelled cell populations: load the 10x matrix into a Seurat object, perform QC and remove low-quality cells, normalise / find highly variable genes / scale, reduce dimensionality with PCA and choose the number of PCs, cluster and visualise with UMAP, find cluster marker genes, and assign cell-type names.

**Data / organism** PBMC3k, 10x Genomics, human (hg19); archive `pbmc3k_filtered_gene_bc_matrices.tar.gz` (barcodes.tsv, genes.tsv, matrix.mtx).

**Tools** Terminal (tar, R started from the `bio-env` directory), R, Seurat (v5.x; "v5.6" banner shown in Output 1), Louvain clustering, UMAP.

## Procedure

The source PDF is a written report with console screenshots; it does not contain the full code. The only code recoverable is what is visible in the screenshots (see note under Files).

**Part A (terminal)** - the archive `pbmc3k_filtered_gene_bc_matrices.tar.gz` was checked and extracted into `data/pbmc3k`; the three expected files (barcodes.tsv, genes.tsv, matrix.mtx) were confirmed inside `filtered_gene_bc_matrices/hg19/`; a results folder was created and R was started from the bio-env directory. (No commands are shown in the PDF, so there is no script for this part.)

**Part B (R, Seurat)**
1. Read the data with `Read10X` and create the Seurat object (`CreateSeuratObject`, min.cells = 3, min.features = 200).
2. Quality control: calculate percent.mt, view violin and scatter plots, filter cells (nFeature_RNA between 200 and 2500, percent.mt below 5) - `scripts/01_qc_metrics.R` (visible lines only).
3. Normalisation with `NormalizeData` (LogNormalize, scale factor 10,000).
4. 2,000 highly variable genes with `FindVariableFeatures` (vst).
5. Scaling of all genes with `ScaleData`.
6. PCA with `RunPCA`, then `VizDimLoadings`, `DimPlot`, `DimHeatmap`.
7. Dimensionality choice with `ElbowPlot`.
8. Graph-based clustering with `FindNeighbors` (dims 1:10) and `FindClusters` (resolution 0.5).
9. UMAP with `RunUMAP` (dims 1:10).
10. Cluster biomarkers with `FindMarkers` and `FindAllMarkers`, checked with violin plots, feature plots and a heatmap - `scripts/02_clusters_markers.R` (visible lines only).
11. Renaming clusters with cell-type labels (`RenameIdents`) and final labelled UMAP.

## Results

### 1. Seurat object
13,714 genes (features) across 2,700 cells in one assay (RNA), only the raw counts layer present. Warning: feature names cannot contain underscores, replaced with dashes.

![Output 1: Seurat object after loading](results/out01_seurat_object_loaded.png)

Raw counts of CD3D, TCL1A and MS4A1 in the first 30 cells; dense matrix 709,591,472 bytes (~710 MB) vs sparse 29,905,192 bytes (~30 MB), ratio 23.7.

![Output 2: raw counts, dense vs sparse size](results/out02_raw_counts_dense_vs_sparse.png)

### 2. Quality control
QC metrics of the first five cells (e.g. cell 1: 2,419 counts, 779 genes, 3.02% mt; cell 3: 3,147 counts, 1,129 genes, 0.89% mt).

![Output 3: QC metrics, first five cells](results/out03_qc_metrics_first5_cells.png)

![Figure 1: violin plots of nFeature_RNA, nCount_RNA, percent.mt](results/fig01_violin_qc_metrics.png)

Correlations: percent.mt vs nCount_RNA = -0.05; nFeature_RNA vs nCount_RNA = 0.95.

![Figure 2: FeatureScatter plots](results/fig02_featurescatter_qc.png)

After filtering: 13,714 genes and 2,638 cells (62 of 2,700 cells, about 2.3%, removed).

### 3. Highly variable features
2,000 variable genes (red) and 11,714 non-variable genes (black); labelled top genes include PPBP, PF4, S100A8, GNLY, IGLL5, FTL, FTH1.

![Figure 3: variable feature plot](results/fig03_variable_features.png)

### 4. PCA
Top genes per PC:

| PC | Positive side | Negative side | Interpretation |
|---|---|---|---|
| PC1 | CST3, TYROBP, LST1, AIF1, FTL | MALAT1, LTB, IL32, IL7R, CD2 | Monocytes/myeloid cells vs T cells |
| PC2 | NKG7, PRF1, CST7, GZMB, GZMA | CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1 | Cytotoxic NK/T cells vs B cells |
| PC3 | PPBP, PF4, SDPR, SPARC, GNG11 | HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPB1 | Platelets vs B cells |
| PC4 | HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1 | VIM, IL7R, S100A6, IL32, S100A8 | B cells vs T cells/monocytes |
| PC5 | GZMB, NKG7, S100A8, FGFBP2, GNLY | LTB, IL7R, CKB, VIM, MS4A7 | NK cells (and some monocyte genes) vs T cells |

![Figure 4: VizDimLoadings PC1 and PC2](results/fig04_vizdimloadings_pc1_pc2.png)

![Figure 5: PCA plot (PC1 vs PC2)](results/fig05_pca_pc1_pc2.png)

![Figure 6: DimHeatmap of PC1 (500 cells)](results/fig06_dimheatmap_pc1.png)

![Figure 7: DimHeatmap of PC1 to PC15](results/fig07_dimheatmap_pc1_15.png)

### 5. Dimensionality
Elbow plot standard deviations: PC1 about 7.1, PC2 about 4.5, PC3 about 3.9, PC4 about 3.75, PC5 about 3.2, PC6 about 2.5, PC7 about 2.1, flat at about 1.8 after PC9. 10 PCs (dims = 1:10) were used.

![Figure 8: elbow plot](results/fig08_elbow_plot.png)

### 6. Clustering
Graph with 2,638 nodes and 95,927 edges; Louvain, resolution 0.5; 9 communities; maximum modularity in 10 random starts 0.8728.

![Output 4: FindNeighbors / FindClusters](results/out04_findneighbors_findclusters_louvain.png)

First five cells fall in clusters 2, 3, 2, 1, 6; levels 0 to 8.

![Output 5: cluster IDs of first five cells](results/out05_cluster_ids_first5_cells.png)

![Figure 9: UMAP coloured by cluster](results/fig09_umap_clusters.png)

### 7. Cluster biomarkers
Cluster 2 (FindMarkers, top 5): IL32 (avg_log2FC 1.31), LTB (1.33), CD3D (1.06), IL7R (1.44), LDHB (0.99), adjusted p-values about 1e-63 to 1e-86. IL7R: pct.1 75% vs pct.2 32.6%; IL32: 94.7% vs 46.5%.

![Output 6: FindMarkers cluster 2](results/out06_findmarkers_cluster2.png)

Cluster 5 vs clusters 0 and 3: FCGR3A (6.79), IFITM3 (6.19), CFD (6.02), CD68 (5.53), RP11-290F20.3 (6.30); adjusted p about 1e-182 to 1e-200; expressed in 84 to 97.5% of cluster 5 cells vs 1.7 to 4.9% in clusters 0 and 3.

![Output 7: FindMarkers cluster 5 vs 0 and 3](results/out07_findmarkers_cluster5_vs_0_3.png)

FindAllMarkers (positive, avg_log2FC > 1), cluster 0 top 10: LDHB (1.21), CCR7 (2.40), CD3D (1.06), CD3E (1.04), LEF1 (2.14), NOSIP, PIK3IP1, PRKCQ-AS1, FHIT (2.70), MAL (1.96); 7,009 more rows for other clusters.

![Output 8: FindAllMarkers cluster 0](results/out08_findallmarkers_cluster0.png)

![Figure 10: violin MS4A1 and CD79A](results/fig10_violin_ms4a1_cd79a.png)

![Figure 11: violin NKG7 and PF4 (raw counts, log scale)](results/fig11_violin_nkg7_pf4_counts.png)

![Figure 12: feature plots of nine marker genes](results/fig12_featureplots_markers.png)

![Figure 13: heatmap of top 10 markers per cluster](results/fig13_heatmap_top10_markers.png)

### 8. Cluster identities

| Cluster | Cell type | Key supporting markers |
|---|---|---|
| 0 | Naive CD4 T | CCR7, LEF1, MAL, CD3D, CD3E |
| 1 | CD14+ Mono | CD14, LYZ, S100A8, S100A9, FCN1 |
| 2 | Memory CD4 T | IL7R, IL32, LTB, CD3D |
| 3 | B | MS4A1, CD79A, HLA-DQA1 |
| 4 | CD8 T | CD8A, NKG7, CD3E |
| 5 | FCGR3A+ Mono | FCGR3A, CD68, IFITM3 |
| 6 | NK | GNLY, NKG7, GZMB, PRF1 |
| 7 | DC | FCER1A, CLEC10A |
| 8 | Platelet | PPBP, PF4 |

![Figure 14: final UMAP with cell-type labels](results/fig14_umap_celltype_labels.png)

## Interpretation

**Output 1.** The object was created with 13,714 genes across 2,700 cells in one assay (RNA), with only the raw counts layer present. The underscore warning is only a naming fix and does not change the data.

**Output 2.** Most entries in the count table are dots; a dot is a zero, meaning the gene was not detected in that cell. CD3D (a T cell gene) has non-zero counts in many of the 30 cells, while TCL1A and MS4A1 (B cell genes) are mostly zero but show a few cells with counts, and MS4A1 reaches high values (such as 36) in what are probably B cells. This is typical of single-cell data, where each cell expresses only a fraction of all genes and dropouts are common. The sparse format is about 24 times smaller than the dense one; Seurat stores only the non-zero values and their positions, which is why large single-cell datasets can be analysed on an ordinary computer.

**Output 3.** The metadata table gives nCount_RNA (total RNA molecules), nFeature_RNA (number of different genes detected) and percent.mt (percentage of reads from mitochondrial genes). All five cells have low mitochondrial percentages (about 0.9 to 3.8%), which indicates healthy, intact cells.

**Figure 1.** nFeature_RNA: most cells have between about 500 and 1,200 detected genes, widest part near 800; a few reach about 2,000 to 2,400. Cells with very few genes are usually empty droplets or damaged cells, and cells with very many genes may be doublets, hence the cut-offs 200 and 2,500. nCount_RNA: most cells have about 1,500 to 3,500 counts, widest part near 2,200, with a long tail up to about 8,000; the shape follows nFeature_RNA. percent.mt: values between 0 and 5%, most between 1 and 3%; a high mitochondrial percentage indicates dying cells whose cytoplasmic RNA has leaked out, so the 5% cut-off removes them. The upper limits visible on the plot match the filtering thresholds, so these violins represent the cells retained after filtering.

**Figure 2.** percent.mt vs nCount_RNA (correlation -0.05): mitochondrial percentage does not depend on sequencing depth, so filtering on percent.mt does not bias the data towards low or high count cells. nFeature_RNA vs nCount_RNA (correlation 0.95): a very strong positive relationship, with points along a tight, almost linear band, as expected for good quality data; cells far from the band would be suspicious. After filtering 62 of the original 2,700 cells (about 2.3%) were removed, so most cells were of good quality.

**Section 3 (Normalisation).** LogNormalize divides each gene's counts in a cell by the total counts of that cell, multiplies by 10,000 and applies log(1 + x). This removes the effect of different sequencing depths so expression can be compared across cells, and the log transformation reduces the influence of very highly expressed genes. No plot is produced at this step.

**Figure 3.** Average expression is on the x-axis (log scale) and standardised variance on the y-axis. The 2,000 red dots are the highly variable genes and the 11,714 black dots are non-variable. Most genes sit near a standardised variance of about 1; the red genes lie above this band, differ more between cells than expected, and carry information about cell types. The labelled genes include PPBP and PF4 (platelet), S100A8 (monocyte), GNLY (NK cells), IGLL5 (B cell lineage), and FTL and FTH1 (ferritin genes, high in monocytes), showing that the main immune cell types drive the variability. Overlapping tick labels and gene labels are only a display problem. Using only these 2,000 genes reduces noise and computing time.

**Section 5 (Scaling).** ScaleData was run on all genes; each gene is centred to mean zero and scaled to variance one across cells, giving all genes equal weight in PCA. Scaling all genes (not only the variable ones) is needed for the heatmaps later, although PCA itself used only the variable genes.

**PCA table.** Each PC is associated with a real biological difference between cell types and not with noise, confirming that the variable genes contain the information needed to separate the cell populations.

**Figure 4.** For PC1, the positive genes (CST3, TYROBP, LST1, AIF1, FTL, LYZ, FCN1, S100A9) are myeloid or monocyte genes, while the single strong negative gene is MALAT1, relatively higher in lymphocytes. For PC2, the positive genes (NKG7, PRF1, CST7, GZMB, GZMA, GNLY, CCL5) are cytotoxic genes of NK and CD8 T cells, and the negative genes (CD79A, CD79B, MS4A1, TCL1A, HLA-DRA, HLA-DQA1, HLA-DQB1) are B cell and antigen-presentation genes. Loadings of the leading genes are of similar size (about 0.1 to 0.15), so the PCs are driven by groups of co-expressed genes and not by a single gene.

**Figure 5.** The cells separate into three arms: a large dense group on the left (negative PC1, PC2 about 0 to 15), the lymphocyte group (T cells, with cells at high PC2 being cytotoxic NK/CD8 cells); a smaller tight group at the bottom left (negative PC1, PC2 near -8), B cells; and a long spread on the right (positive PC1), monocytes and myeloid cells. The first two PCs recover the main lineages (T/NK, B, monocytes), but smaller subtypes are not yet clearly resolved, which is why clustering and UMAP follow.

**Figure 6.** Cells are ordered by PC1 score and genes with the strongest loadings are on the y-axis (yellow high, purple low). In cells on one side, CST3, TYROBP, LST1, AIF1, FTL, LYZ, S100A9, FCN1 are high and MALAT1, LTB, IL32, IL7R, CD2 are low; on the other side the pattern is reversed. PC1 separates myeloid cells from T lymphocytes, and the sharp boundary between the blocks indicates a strong, clear source of variation.

**Figure 7.** PC1 to about PC5 show clear and distinct blocks; PC6 to about PC9 show weaker, patchier structure; from about PC10 onwards (especially PC13 to PC15) the heatmaps become noisy and diffuse, so these PCs mostly carry noise or very small effects. This supports using only the first 10 or so PCs.

**Figure 8.** PC1 is highest (about 7.1), then PC2 (about 4.5), PC3 (about 3.9), PC4 (about 3.75), PC5 (about 3.2); the values drop to about 2.5 at PC6 and about 2.1 at PC7, and after PC9 the curve is flat at about 1.8. The curve bends in the region of PC5 to PC9, and PCs after about PC9 or PC10 add very little information. Together with Figure 7 this shows that the first 10 PCs capture almost all the real biological signal, so dims = 1:10 was used. Using a few more PCs than the strict elbow is safer than using too few, because it does not discard weaker but real populations such as the small DC and platelet groups.

**Output 4.** A nearest-neighbour graph was built with 2,638 nodes and 95,927 edges; Louvain grouped the cells into 9 communities at resolution 0.5. The maximum modularity in 10 random starts was 0.8728. Values near 1 indicate very clear community structure, so 0.87 means the clusters are well separated.

**Output 5.** The first five cells were assigned to clusters 2, 3, 2, 1 and 6. The levels are 0 to 8, so there are nine clusters in total.

**Figure 9.** Clusters 0, 2, 4 and 6 form one large connected region in the upper right; clusters 0 and 2 lie side by side and partly touch (very similar cell types); cluster 4 continues from this region and leads to cluster 6 at the lower right tip, a continuum from T cells to cytotoxic cells. Clusters 1 and 5 form a second group on the left, cluster 5 directly next to cluster 1. Cluster 3 forms an isolated island at the bottom. Clusters 7 and 8 are very small: cluster 7 lies just above the cluster 1 group, cluster 8 is a tiny isolated group at the far left. The broad separation into three regions (T/NK cells, monocytes, B cells) agrees with the PCA plot. Distances between separate islands in UMAP should not be interpreted quantitatively; only grouping and neighbourhood are meaningful.

**Output 6.** The top markers of cluster 2 are IL32, LTB, CD3D, IL7R and LDHB with extremely small adjusted p-values. CD3D is a core T cell gene, and IL7R together with LTB and IL32 is characteristic of memory CD4 T cells. Cluster 2 is therefore identified as Memory CD4 T cells. The log fold changes are modest (about 1 to 1.4) because other clusters also contain T cells that share these genes.

**Output 7.** FCGR3A (CD16), CD68 and IFITM3 are markers of non-classical (FCGR3A+) monocytes. The fold changes are very large because cluster 5 was compared with clusters 0 (T cells) and 3 (B cells), which are very different cell types; the values would be lower against other monocytes. Cluster 5 is therefore identified as FCGR3A+ monocytes.

**Output 8.** CD3D and CD3E show that these are T cells, while CCR7, LEF1 and MAL are classical markers of naive T cells, which have not yet met an antigen. Cluster 0 is therefore identified as Naive CD4 T cells. This also explains why clusters 0 and 2 lie side by side in the UMAP: both are CD4 T cells but in different states (naive versus memory).

**Figure 10.** Both genes are strongly expressed only in cluster 3: MS4A1 has a median of about 2 to 2.5 and CD79A about 3, while in all other clusters the violins are flat at zero with only a few scattered cells. MS4A1 (CD20) and CD79A (a B cell receptor component) are classical B cell markers, so cluster 3 is identified as B cells.

**Figure 11.** NKG7 is highest in cluster 6 (up to about 40 counts), also high in cluster 4, and low or absent in the others; this is typical of cytotoxic cells, with cluster 6 being NK cells and cluster 4 being CD8 T cells. PF4 is expressed only in cluster 8 (counts of about 7 to 45). PF4 (and PPBP) are platelet-specific genes, so cluster 8 is identified as platelets. The horizontal rows of dots at 1, 2 and 3 appear because raw counts are whole numbers.

**Figure 12.** MS4A1: only the isolated bottom island, the B cell cluster. GNLY: the lower right tip of the T/NK region, the NK cluster (6). CD3E: throughout the large upper right region, confirming clusters 0, 2, 4 and 6 are T lineage cells (lower in the NK tip). CD14: the left group, mainly cluster 1, classical monocytes. FCER1A: a very small group near the monocytes, the dendritic cell cluster (7). FCGR3A: the smaller part of the monocyte group (cluster 5) and the NK tip, since NK cells also express CD16. LYZ: throughout the monocyte group and the small DC group. PPBP: a single tiny group, the platelet cluster (8). CD8A: the lower part of the T cell region (cluster 4) and partly the NK tip. Each marker gene is restricted to the expected cluster, so the cell-type identification is supported by the expression data and not only by the clustering.

**Figure 13.** A clear yellow block appears along the diagonal for most clusters, meaning each cluster has its own set of strongly expressed genes. Clusters 0 and 2 share T cell genes (CD3D, IL7R, LTB, LDHB), and cluster 0 shows CCR7, LEF1 and MAL; the difference between them is smaller than between the others. Cluster 1 shows S100A8, S100A9, LGALS1, FCN1, TYROBP, CST3 and LYZ (classical monocytes). Cluster 3 shows CD79A, MS4A1, HLA-DQA1, HLA-DQB1 and HLA-DRA (B cells). Clusters 4 and 6 share cytotoxic genes (NKG7, GZMA, CST7, CTSW); cluster 6 additionally shows GZMB, FGFBP2, PRF1 and GNLY (NK cells). Cluster 5 shows FCGR3A, MS4A7, IFITM3 and related genes (FCGR3A+ monocytes), with some monocyte genes shared with cluster 1. Clusters 7 and 8 are very narrow because they contain few cells; their specific genes (FCER1A and CLEC10A for cluster 7, PPBP and PF4 for cluster 8) appear at the bottom as small bright blocks. The heatmap confirms that the clusters are real biological groups with distinct expression programmes.

**Figure 14.** The final plot shows the same UMAP as Figure 9, with each cluster carrying its cell-type name. The T cell region (upper right) contains Naive CD4 T, Memory CD4 T, CD8 T and NK cells in a connected arrangement, expected for closely related lymphocytes. The myeloid cells (CD14+ Mono and FCGR3A+ Mono) form the left group, with the small DC cluster next to them. B cells form a separate island at the bottom, and platelets a tiny isolated group on the left. The cluster sizes look reasonable for blood: T cells are the most numerous, followed by monocytes, then B cells and NK cells, while DCs and platelets are rare.

**Overall conclusion.** Starting from 2,700 cells and 13,714 genes, quality control retained 2,638 good quality cells. Normalisation, selection of 2,000 variable genes and PCA captured the main variation, and the elbow plot and PC heatmaps showed that about the first 10 PCs contain the real signal. Louvain clustering (modularity 0.87) gave 9 clusters, and the marker genes, violin plots, feature plots and heatmap identified them as Naive CD4 T cells, CD14+ monocytes, Memory CD4 T cells, B cells, CD8 T cells, FCGR3A+ monocytes, NK cells, dendritic cells and platelets. These are the expected cell populations of human PBMC, which shows that the workflow was carried out correctly and that the dataset can be resolved into its major immune cell types.

## Files

```
07_scrna_seurat_pbmc3k/
  README.md
  reports/Seurat_Debarshi-1.pdf
  results/   fig01..fig14 (plots), out01..out08 (console screenshots)
  scripts/01_qc_metrics.R        (only lines visible in the PDF screenshots)
  scripts/02_clusters_markers.R  (only lines visible in the PDF screenshots)
```

Note: the PDF does not include a full code listing; scripts contain only code that is legible in its console screenshots.
