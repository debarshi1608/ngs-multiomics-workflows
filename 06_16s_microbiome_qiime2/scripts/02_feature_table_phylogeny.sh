# Source: Debarshi_Das_Adhkary_255hsbb011_16S_rRNA_QIIME2-1.pdf  --  5(B) Feature table and phylogeny
qiime feature-table summarize --i-table table-dada2.qza --o-visualization table.qzv
qiime feature-table tabulate-seqs --i-data rep-seqs-dada2.qza --o-visualization rep-seqs.qzv
qiime tools view table.qzv
qiime tools view rep-seqs.qzv
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs-dada2.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza
