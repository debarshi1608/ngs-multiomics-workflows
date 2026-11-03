# Source: Debarshi_011_Variant_call.pdf  (Phases 2-4: GATK preparation, duplicates, BQSR, variant calling, filtering)

# 7A: sequence dictionary
gatk CreateSequenceDictionary -R GCF_000146045.2_R64_genomic.fa
# 7B: FASTA index
samtools faidx GCF_000146045.2_R64_genomic.fa

# 8. Mark PCR duplicates
picard MarkDuplicates \
   I=conditionA_Rep2.sorted.bam \
   O=conditionA_Rep2_markdup.bam \
   M=marked_dup_metrics.txt
# Index the markdup BAM
samtools index conditionA_Rep2_markdup.bam

# 9. Initial variant calling (draft set for BQSR)
gatk HaplotypeCaller \
  -R GCF_000146045.2_R64_genomic.fa \
  -I conditionA_Rep2_markdup.bam \
  -O conA_variants.vcf

# 10. Base quality score recalibration
gatk BaseRecalibrator \
  -R GCF_000146045.2_R64_genomic.fa \
  -I conditionA_Rep2_markdup.bam \
  --known-sites conA_variants.vcf \
  -O BQSR_out_conA.table

# 11. Apply BQSR
gatk ApplyBQSR \
  -R GCF_000146045.2_R64_genomic.fa \
  -I conditionA_Rep2_markdup.bam \
  --bqsr-recal-file BQSR_out_conA.table \
  -O conA_BQSR_applied_reads.bam

# 12. Final variant calling
gatk HaplotypeCaller \
  -R GCF_000146045.2_R64_genomic.fa \
  -I conA_BQSR_applied_reads.bam \
  -O conA_variants_final.vcf

# 13. Hard-filter low-quality variants
gatk VariantFiltration \
  -R GCF_000146045.2_R64_genomic.fa \
  -V conA_variants_final.vcf \
  --filter-expression "QD < 10.0" --filter-name "LowQD" \
  --filter-expression "FS > 60.0" --filter-name "HighFS" \
  --filter-expression "MQ < 40.0" --filter-name "LowMQ" \
  --filter-expression "MQRankSum < -12.5" --filter-name "LowMQRankSum" \
  -O conA_variants_final.filtered.vcf
