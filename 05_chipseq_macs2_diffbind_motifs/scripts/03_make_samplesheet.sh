# Source: diffbind_analysis_Debarshi_011.pdf  (DiffBind sample sheet; paths as used in the lab)
cat << 'EOF2' > sample_sheet.csv
SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
S1,yeast,ChIP,WT,1,/home/ibab/bio-env/mutant fastq files/C1.sorted.bam,I1,/home/ibab/bio-env/mutant fastq files/I1.sorted.bam,/home/ibab/bio-env/mutant fastq files/WT_Rep1_narrow_peaks.narrowPeak,narrow
S2,yeast,ChIP,WT,2,/home/ibab/bio-env/mutant fastq files/C2.sorted.bam,I2,/home/ibab/bio-env/mutant fastq files/I2.sorted.bam,/home/ibab/bio-env/mutant fastq files/WT_Rep2_narrow_peaks.narrowPeak,narrow
S3,yeast,ChIP,Mutant,1,/home/ibab/bio-env/mutant fastq files/MC1.sorted.bam,MI1,/home/ibab/bio-env/mutant fastq files/MI1.sorted.bam,/home/ibab/bio-env/mutant fastq files/Mutant_Rep1_narrow_peaks.narrowPeak,narrow
S4,yeast,ChIP,Mutant,2,/home/ibab/bio-env/mutant fastq files/MC2.sorted.bam,MI2,/home/ibab/bio-env/mutant fastq files/MI2.sorted.bam,/home/ibab/bio-env/mutant fastq files/Mutant_Rep2_narrow_peaks.narrowPeak,narrow
EOF2
pixi run R
