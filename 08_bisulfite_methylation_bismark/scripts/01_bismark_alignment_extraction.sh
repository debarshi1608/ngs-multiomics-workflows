# Source: Debarshi_011_bismark.pdf  (whole-genome bisulfite sequencing, Bismark pipeline)

# Step 1: reference genome preparation (in-silico bisulfite conversion + indexes)
mkdir BS_seq_WGS
cp ref.fa BS_seq_WGS/
bismark_genome_preparation BS_seq_WGS/

# Step 2: bisulfite alignment (Bowtie2); R1 and R2 aligned separately
bismark --bowtie2 BS_seq_WGS/ BS_seq_WGS/R1.fastq.gz
bismark --bowtie2 BS_seq_WGS/ BS_seq_WGS/R2.fastq.gz

# Step 3: sort the BAM files
samtools sort R1_bismark_bt2.bam R1_bismark_bt2_sorted
samtools sort R2_bismark_bt2.bam R2_bismark_bt2_sorted

# Step 4: index the sorted BAM files
samtools index R1_bismark_bt2_sorted.bam
samtools index R2_bismark_bt2_sorted.bam

# Step 5: methylation extraction + bedGraph
bismark_methylation_extractor --bedGraph --counts R1_bismark_bt2_sorted.bam
bismark_methylation_extractor --bedGraph --counts R2_bismark_bt2_sorted.bam

# Step 6: IGV (manual): load ref.fa, the sorted BAMs and the bedGraph files

# Step 7: M-bias inspection
head -20 R1_bismark_bt2_sorted.M-bias.txt
