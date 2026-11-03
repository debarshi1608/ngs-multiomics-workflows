# Source: Debarshi_chip_seq_analysis_24_09.pdf  (E. coli BW25113 transcription-factor ChIP-seq, samples C and D)

# ---- Align to reference and sort ----
# 1. ChIP sample C
bowtie2 -p 8 -x E_coli_BW25113_index -U SO_4933_C_CHR1_R1.fastq.gz | samtools sort -@ 4 -o C_CHR1_sorted.bam -
samtools index C_CHR1_sorted.bam
# 2. Input control C (INR)
bowtie2 -p 8 -x E_coli_BW25113_index -U SO_4933_C_INR1_R1.fastq.gz | samtools sort -@ 4 -o C_INR1_sorted.bam -
samtools index C_INR1_sorted.bam
# 3. ChIP sample D
bowtie2 -p 8 -x E_coli_BW25113_index -U SO_4933_D_CHR1_R1.fastq.gz | samtools sort -@ 4 -o D_CHR1_sorted.bam -
samtools index D_CHR1_sorted.bam
# 4. Input control D (INR)
bowtie2 -p 8 -x E_coli_BW25113_index -U SO_4933_D_INR1_R1.fastq.gz | samtools sort -@ 4 -o D_INR1_sorted.bam -
samtools index D_INR1_sorted.bam

# ---- MACS peak calling ----
# Sample C
pixi run macs3 callpeak \
 -t C_CHR1_sorted.bam \
 -c C_INR1_sorted.bam \
 -f BAM \
 -g 4.6e6 \
 --outdir SampleC \
 -n sorted_C \
 -B \
 --nomodel \
 --extsize 150 \
 -p 0.001
# Sample D
pixi run macs3 callpeak \
 -t D_CHR1_sorted.bam \
 -c D_INR1_sorted.bam \
 -f BAM \
 -g 4.6e6 \
 --outdir SampleD \
 -n sorted_D \
 -B \
 --nomodel \
 --extsize 150 \
 -p 0.001

# ---- IGV ----
cd ~/bio-env
pixi run igv

# ---- Merge peaks ----
pixi run mergePeaks ./SampleC/sorted_C_peaks.narrowPeak -d 300 > merge300.txt
pixi run mergePeaks ./SampleD/sorted_D_peaks.narrowPeak -d 300 > merge300_D.txt

# ---- Annotation (HOMER) ----
pixi run annotatePeaks.pl ./SampleC/sorted_C_peaks.narrowPeak E.coli_BW25113.fasta -gff E.coli_BW25113_annotation.gff > annotated_peaks.tsv
pixi run annotatePeaks.pl ./SampleD/sorted_D_peaks.narrowPeak E.coli_BW25113.fasta -gff E.coli_BW25113_annotation.gff > annotated_peaks_D.tsv

# ---- Closest genes (bedtools) ----
pixi run bedtools closest -a ./SampleC/sorted_C_peaks.narrowPeak -b E.coli_BW25113_annotation.bed | awk -F "\t" '{if($18 == "gene"){print}}' > C_closest_genes.bed
pixi run bedtools closest -a ./SampleD/sorted_D_peaks.narrowPeak -b E.coli_BW25113_annotation.bed | awk -F "\t" '{if($18 == "gene"){print}}' > D_closest_genes.bed

# ---- Extract peak DNA sequences ----
pixi run samtools faidx E.coli_BW25113.fasta
pixi run bedtools getfasta -fi E.coli_BW25113.fasta -bed ./SampleC/sorted_C_peaks.narrowPeak -fo C_seq.fasta
pixi run bedtools getfasta -fi E.coli_BW25113.fasta -bed ./SampleD/sorted_D_peaks.narrowPeak -fo D_seq.fasta

# ---- Motif discovery (MEME) ----
pixi run meme C_seq.fasta -dna -nmotifs 3 -oc meme_out_C
pixi run meme D_seq.fasta -dna -nmotifs 3 -oc meme_out_D
