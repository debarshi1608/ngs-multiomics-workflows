# Source: Debarshi_Das_Adhkary_255hsbb011_16S_rRNA_QIIME2-1.pdf  --  5(E) Gut-specific differential abundance (ANCOM-BC)
qiime feature-table filter-samples --i-table table-dada2.qza --m-metadata-file sample-metadata.tsv --p-where "[body-site]='gut'" --o-filtered-table gut-table.qza
qiime composition ancombc --i-table gut-table.qza --m-metadata-file sample-metadata.tsv --p-formula "subject" --o-differentials ancombc-subject.qza
qiime composition da-barplot --i-data ancombc-subject.qza --o-visualization da-barplot-subject.qzv
qiime tools view da-barplot-subject.qzv
qiime taxa collapse --i-table gut-table.qza --i-taxonomy taxonomy.qza --p-level 7 --o-collapsed-table gut-table-l7.qza
qiime composition ancombc --i-table gut-table-l7.qza --m-metadata-file sample-metadata.tsv --p-formula "subject" --o-differentials l7-ancombc-subject.qza --verbose
qiime composition da-barplot --i-data l7-ancombc-subject.qza --o-visualization l7-da-barplot-subject.qzv
qiime tools view l7-da-barplot-subject.qzv
qiime tools export --input-path l7-ancombc-subject.qza --output-path ancombc-genus-export
