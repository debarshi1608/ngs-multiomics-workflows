# Source: Alternative_Splicing_Analysis_Debarshi_Das_Adhikary-1.pdf (Part A: terminal, HISAT2 + SAMtools)

# Step 1: Create and activate the conda environment
conda create -n altsplice -c conda-forge -c bioconda hisat2 samtools -y
conda activate altsplice

# Step 2: Check the installed versions
hisat2 --version
samtools --version

# Step 3: Index the reference genome
hisat2-build GCF_000146045.2_R64_genomic.fna Hisat_index
ls Hisat_index*

# Step 4: Align the single-end reads
hisat2 -x Hisat_index -U conA_rep1.fq -S A_1.sam
hisat2 -x Hisat_index -U conA_rep2.fq -S A_2.sam
hisat2 -x Hisat_index -U conB_rep1.fq -S B_1.sam
hisat2 -x Hisat_index -U conB_rep2.fq -S B_2.sam

# Step 5: Convert SAM to BAM
samtools view -bS A_1.sam -o A_1.bam
samtools view -bS A_2.sam -o A_2.bam
samtools view -bS B_1.sam -o B_1.bam
samtools view -bS B_2.sam -o B_2.bam

# Step 6: Sort the BAM files
samtools sort -o A_1_sorted.bam A_1.bam
samtools sort -o A_2_sorted.bam A_2.bam
samtools sort -o B_1_sorted.bam B_1.bam
samtools sort -o B_2_sorted.bam B_2.bam

# Step 7: Index the sorted BAM files
samtools index A_1_sorted.bam
samtools index A_2_sorted.bam
samtools index B_1_sorted.bam
samtools index B_2_sorted.bam
ls *_sorted.bam*
