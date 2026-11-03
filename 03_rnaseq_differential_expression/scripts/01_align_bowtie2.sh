# Source: Debarshi_13_8_part1.pdf  (RNA-seq: QC + alignment; run in a folder with the 4 reads and the S. cerevisiae R64 FASTA/GTF)
mv conA_rep1.fq.gz conA_rep1.fq
mv conA_rep2.fq.gz conA_rep2.fq
mv conB_rep1.fq.gz conB_rep1.fq
mv conB_rep2.fq.gz conB_rep2.fq
file * .fq

# FastQC
fastqc -o fastqc_results \
conA_rep1.fq \
conA_rep2.fq \
conB_rep1.fq \
conB_rep2.fq

# Index reference genome
mkdir -p bowtie2_index
bowtie2-build GCF_000146045.2_R64_genomic.fna S288C

# Alignment - to be done for all four files
bowtie2 \
-x bowtie2_index/S288C \
-U conA_rep1.fq \
-S conA_rep1.sam
