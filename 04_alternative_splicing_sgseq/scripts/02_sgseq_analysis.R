# Source: Alternative_Splicing_Analysis_Debarshi_Das_Adhikary-1.pdf (Part B: RStudio, SGSeq splice graph analysis)

# Step 8: Install and load SGSeq (install only the first time)
setwd("~/bio-env")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("SGSeq")
library(SGSeq)

# Step 9: Prepare the input data
# The column names must be exactly sample_name and file_bam.
file_bam <- c("A_1_sorted.bam", "A_2_sorted.bam", "B_1_sorted.bam", "B_2_sorted.bam")
sample_name <- c("WT1", "WT2", "MT1", "MT2")
v <- data.frame(sample_name, file_bam)
v

# Step 10: Extract BAM file information
bam_info <- getBamInfo(v)
bam_info

# Step 11: Import the annotation
transcripts <- importTranscripts(file = "GCF_000146045.2_R64_genomic.gtf")
tx_features <- convertToTxFeatures(transcripts)

# Step 12: Predict features and quantify expression
AnaFeat <- analyzeFeatures(bam_info, features = tx_features)
FPKM(AnaFeat)
rowRanges(AnaFeat)

# Step 13: Visualise splice features
plotFeatures(AnaFeat, geneName = "YAL030W")

# For a better plot, a modified version of the function was used. It matches features by featureID first, then
# by name, and drops unmatched (NA) entries, which avoids the indexing problem of the default function.
plotFeatures_fixed <- function(x, geneID = NULL, geneName = NULL, which = NULL,
    tx_view = FALSE, cex = 1, assay = "FPKM",
    include = c("junctions", "exons", "both"),
    transform = function(x) log2(x + 1),
    Rowv = NULL, distfun = dist, hclustfun = hclust, margin = 0.2,
    RowSideColors = NULL, square = FALSE, cexRow = 1, cexCol = 1,
    labRow = colnames(x),
    col = colorRampPalette(c("black", "gold"))(256),
    zlim = NULL, heightPanels = c(1, 2), ...)
{
    include <- match.arg(include)
    x <- restrictFeatures(x, geneID = geneID, which = which, geneName = geneName)
    if (nrow(x) == 0) return()
    features <- rowRanges(x)
    if (!is.null(RowSideColors) && !is.list(RowSideColors)) {
        RowSideColors <- list(RowSideColors)
    }
    n_sample <- ncol(x)
    include_type <- switch(include, junctions = "J", exons = "E", both = c("E", "J"))
    n_feature <- length(which(type(features) %in% include_type))
    pars <- getLayoutParameters(n_sample, n_feature, margin, heightPanels,
                                RowSideColors, square, tx_view)
    layout(pars$mat, heights = pars$hei, widths = pars$wid)
    par(cex = cex, mai = c(pars$mai[1], 0, pars$mai[3], 0))
    df <- plotSpliceGraph(x = features, label = "id", tx_view = tx_view,
                          short_output = FALSE, ...)
    i_df <- switch(include,
                   exons = which(df$type == "E"),
                   junctions = which(df$type == "J"),
                   both = seq_len(nrow(df)))
    i_df <- i_df[!duplicated(df$name[i_df])]
    # FIX: match on featureID first, then fall back to name, then drop NAs
    i <- match(df$featureID[i_df], featureID(features))
    miss <- is.na(i)
    if (any(miss)) {
        i[miss] <- match(df$name[i_df][miss], feature2name(features))
    }
    keep <- !is.na(i)
    i_df <- i_df[keep]
    i <- i[keep]
    X <- transform(assay(x[i, ], assay))
    labCol <- df$id[i_df]
    colLabCol <- df$color[i_df]
    plotImage(X, Rowv = Rowv, distfun = distfun, hclustfun = hclustfun,
              RowSideColors = RowSideColors, cexRow = cexRow, cexCol = cexCol,
              labRow = labRow, labCol = labCol, colLabCol = colLabCol,
              col = col, zlim = zlim)
    frame()
    df <- df[c("id", "name", "type", "featureID")]
    invisible(df)
}

# let it see SGSeq's internal functions (getLayoutParameters, plotImage, etc.)
environment(plotFeatures_fixed) <- asNamespace("SGSeq")

while (dev.cur() > 1) dev.off()
plotFeatures_fixed(AnaFeat, geneID = 21, include = "both")

# Step 14: Find the gene ID, then plot the splice graph and coverage
rr_df <- as.data.frame(rowRanges(AnaFeat))
View(rr_df)    # find the geneID of your gene of interest (e.g. YAL030W)

par(mfrow = c(5, 1), mar = c(1, 3, 1, 1))
plotSpliceGraph(rowRanges(AnaFeat), geneID = 505, toscale = "none", color_novel = "red")
for (j in 1:4) {
  plotCoverage(AnaFeat[, j], geneID = 505, toscale = "none")
}

# save the figure as a PNG file
png("splicegraph_coverage.png", width = 1800, height = 2000, res = 200)
par(mfrow = c(5, 1), mar = c(1, 3, 1, 1))
plotSpliceGraph(rowRanges(AnaFeat), geneID = 505, toscale = "none", color_novel = "red")
for (j in 1:4) plotCoverage(AnaFeat[, j], geneID = 505, toscale = "none")
dev.off()
