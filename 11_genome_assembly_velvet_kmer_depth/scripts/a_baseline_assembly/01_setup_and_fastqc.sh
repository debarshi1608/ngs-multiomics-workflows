# Source PDF: Debarshi_day_20_8.pdf
pixi add --channel conda-forge --channel bioconda velvet quast augustus gffread ragtag busco seqtk

pixi run fastqc SR1.fastq SR2.fastq
