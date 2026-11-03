# Source: diffbind_analysis_Debarshi_011.pdf  (paired-end WT vs Mutant ChIP-seq: alignment, MACS peak calling, BigWig tracks)

# ---- Environment (Pixi) ----
pixi project channel add bioconda
pixi add bwa samtools macs2 deeptools

# ---- Step 1: index reference ----
bwa index reference.fna
samtools faidx reference.fna

# ---- Step 2: align, sort, index ----
# Mutant replicate 1 (ChIP and input)
bwa mem -t 10 reference.fna MC1_R1.fastq.gz MC1_R2.fastq.gz | samtools sort -@ 4 -m 2G -o MC1.sorted.bam && samtools index -@ 4 MC1.sorted.bam
bwa mem -t 10 reference.fna MI1_R1.fastq.gz MI1_R2.fastq.gz | samtools sort -@ 4 -m 2G -o MI1.sorted.bam && samtools index -@ 4 MI1.sorted.bam
# Mutant replicate 2
bwa mem -t 10 reference.fna MC2_R1.fastq.gz MC2_R2.fastq.gz | samtools sort -@ 4 -m 2G -o MC2.sorted.bam && samtools index -@ 4 MC2.sorted.bam
bwa mem -t 10 reference.fna MI2_R1.fastq.gz MI2_R2.fastq.gz | samtools sort -@ 4 -m 2G -o MI2.sorted.bam && samtools index -@ 4 MI2.sorted.bam
# WT replicate 1
bwa mem -t 10 reference.fna C1_R1.fastq.gz C1_R2.fastq.gz | samtools sort -@ 4 -m 2G -o C1.sorted.bam && samtools index -@ 4 C1.sorted.bam
bwa mem -t 10 reference.fna I1_R1.fastq.gz I1_R2.fastq.gz | samtools sort -@ 4 -m 2G -o I1.sorted.bam && samtools index -@ 4 I1.sorted.bam
# WT replicate 2
bwa mem -t 10 reference.fna C2_R1.fastq.gz C2_R2.fastq.gz | samtools sort -@ 4 -m 2G -o C2.sorted.bam && samtools index -@ 4 C2.sorted.bam
bwa mem -t 10 reference.fna I2_R1.fastq.gz I2_R2.fastq.gz | samtools sort -@ 4 -m 2G -o I2.sorted.bam && samtools index -@ 4 I2.sorted.bam

# ---- Step 3.1: narrow peak calling ----
# (the PDF notes: replace -g hs with the right effective genome size for non-human organisms, e.g. -g 1.2e7 for yeast)
macs3 callpeak -t MC1.sorted.bam -c MI1.sorted.bam -f BAMPE -g hs -n Mutant_Rep1_narrow -B
macs3 callpeak -t MC2.sorted.bam -c MI2.sorted.bam -f BAMPE -g hs -n Mutant_Rep2_narrow -B
macs3 callpeak -t C1.sorted.bam -c I1.sorted.bam -f BAMPE -g hs -n WT_Rep1_narrow -B
macs3 callpeak -t C2.sorted.bam -c I2.sorted.bam -f BAMPE -g hs -n WT_Rep2_narrow -B

# ---- Step 3.2: broad peak calling ----
macs3 callpeak -t MC1.sorted.bam -c MI1.sorted.bam -f BAMPE -g hs --broad -n Mutant_Rep1_broad -B
macs3 callpeak -t MC2.sorted.bam -c MI2.sorted.bam -f BAMPE -g hs --broad -n Mutant_Rep2_broad -B
macs3 callpeak -t C1.sorted.bam -c I1.sorted.bam -f BAMPE -g hs --broad -n WT_Rep1_broad -B
macs3 callpeak -t C2.sorted.bam -c I2.sorted.bam -f BAMPE -g hs --broad -n WT_Rep2_broad -B

# ---- Step 4: RPKM-normalised BigWig tracks ----
bamCoverage -b MC1.sorted.bam -o MC1.rpkm.bw --normalizeUsing RPKM --binSize 10 -p 12
bamCoverage -b MC2.sorted.bam -o MC2.rpkm.bw --normalizeUsing RPKM --binSize 10 -p 12
bamCoverage -b C1.sorted.bam -o C1.rpkm.bw --normalizeUsing RPKM --binSize 10 -p 12
bamCoverage -b C2.sorted.bam -o C2.rpkm.bw --normalizeUsing RPKM --binSize 10 -p 12

# ---- Step 5: IGV (manual) ----
# Load reference.fna, the four *.rpkm.bw tracks and the *_peaks.narrowPeak / *_peaks.broadPeak files;
# set an identical Data Range (e.g. 0-50) on all coverage tracks.
