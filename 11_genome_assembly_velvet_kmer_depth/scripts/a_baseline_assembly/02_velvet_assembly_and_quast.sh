# Source PDF: Debarshi_day_20_8.pdf
velveth Assembly_step1 31 -shortPaired -fastq -separate SR1.fastq SR2.fastq

velvetg Assembly_step1 -cov_cutoff 5 -unused_reads yes

grep -c ">" Assembly_step1/contigs.fa

quast.py -o quast_contigs Assembly_step1/contigs.fa
