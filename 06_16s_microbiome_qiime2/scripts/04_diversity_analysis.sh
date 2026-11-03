# Source: Debarshi_Das_Adhkary_255hsbb011_16S_rRNA_QIIME2-1.pdf  --  5(D) Diversity analysis
qiime diversity alpha-rarefaction --i-table table-dada2.qza --i-phylogeny rooted-tree.qza --p-max-depth 9820 --m-metadata-file sample-metadata.tsv --o-visualization alpha-rarefaction.qzv
qiime tools view alpha-rarefaction.qzv
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table-dada2.qza --p-sampling-depth 897 --m-metadata-file sample-metadata.tsv --output-dir core-metrics-results
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results/faith-pd-body-site.qzv
qiime tools view core-metrics-results/faith-pd-body-site.qzv
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/shannon_vector.qza --m-metadata-file sample-metadata.tsv --o-visualization core-metrics-results/shannon-body-site.qzv
qiime tools view core-metrics-results/shannon-body-site.qzv
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza --m-metadata-file sample-metadata.tsv --m-metadata-column body-site --o-visualization core-metrics-results/bray-curtis-body-site.qzv --p-pairwise
qiime tools view core-metrics-results/bray-curtis-body-site.qzv
qiime tools view core-metrics-results/bray_curtis_emperor.qzv
