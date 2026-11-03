# Source: Debarshi_Das_Adhkary_255hsbb011_16S_rRNA_QIIME2-1.pdf  --  5(A) Data preparation: import, demultiplex, DADA2
mkdir emp-single-end-sequences
cp sequences.fastq.gz barcodes.fastq.gz emp-single-end-sequences/
ls
qiime tools import --type EMPSingleEndSequences --input-path emp-single-end-sequences --output-path emp-single-end-sequences.qza
qiime demux emp-single --i-seqs emp-single-end-sequences.qza --m-barcodes-file sample-metadata.tsv --m-barcodes-column barcode-sequence --o-per-sample-sequences demux.qza --o-error-correction-details demux-details.qza
qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
qiime tools view demux.qzv
qiime tools export --input-path demux.qza --output-path demux_export
ls demux_export
qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza --p-trim-left 0 --p-trunc-len 120 --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza --o-denoising-stats stats-dada2.qza
qiime metadata tabulate --m-input-file stats-dada2.qza --o-visualization stats-dada2.qzv
qiime tools view stats-dada2.qzv
