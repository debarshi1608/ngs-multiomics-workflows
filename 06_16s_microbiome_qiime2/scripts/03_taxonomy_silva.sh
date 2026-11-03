# Source: Debarshi_Das_Adhkary_255hsbb011_16S_rRNA_QIIME2-1.pdf  --  5(C) Taxonomic analysis (SILVA 138, 99%)
mkdir -p ~/Downloads/SILVA2026
cd ~/Downloads/SILVA2026
wget https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza
cd ~/bio-env
qiime feature-classifier classify-sklearn --i-classifier ~/Downloads/SILVA2026/silva-138-99-nb-classifier.qza --i-reads rep-seqs-dada2.qza --o-classification taxonomy.qza --p-reads-per-batch 2000
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
qiime tools view taxonomy.qzv
qiime taxa barplot --i-table table-dada2.qza --i-taxonomy taxonomy.qza --m-metadata-file sample-metadata.tsv --o-visualization taxa-bar-plots.qzv
qiime tools view taxa-bar-plots.qzv
