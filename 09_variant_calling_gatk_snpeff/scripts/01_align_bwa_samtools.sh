# Source: Debarshi_011_Variant_call.pdf  (Phase 1: QC and alignment; single-end yeast reads conditionA_Rep2)

# 1. Quality control
fastqc conditionA_Rep2.fastq

# 2. Index the reference genome for BWA
bwa index GCF_000146045.2_R64_genomic.fa

# 3. Align with BWA-MEM, including read groups (GATK needs @RG: ID, SM, PL, LB)
bwa mem -R '@RG\tID:Flowcell.Lane\tSM:conditionA_Rep2\tPL:ILLUMINA\tLB:LibraryName' \
  GCF_000146045.2_R64_genomic.fa conditionA_Rep2.fastq > conditionA_Rep2.sam

# 4. SAM -> BAM
samtools view -bS conditionA_Rep2.sam > conditionA_Rep2.bam

# 5. Coordinate-sort
samtools sort conditionA_Rep2.bam conditionA_Rep2.sorted

# 6. Index the sorted BAM
samtools index conditionA_Rep2.sorted.bam
