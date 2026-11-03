# Source: Debarshi_011_Functional Annotation.pdf  (section 2: heuristic pileup calling with VarScan)

# Step 1: reference FASTA index
samtools faidx GCF_000146045.2_R64_genomic.fa
# Step 2: base pileup from the sorted BAM
samtools mpileup -f GCF_000146045.2_R64_genomic.fa conditionA_Rep2.sorted.bam > aln.raw.pileup
# Step 3: raw SNPs
varscan pileup2snp aln.raw.pileup > snp_list
# Step 4: raw indels
varscan pileup2indel aln.raw.pileup > indel_list
# Step 5: allele-frequency filter (>= 20 %)
varscan pileup2snp aln.raw.pileup --min-var-freq 0.20 > freq_filtered_snps
# Step 6: Fisher's exact test p-value filter
varscan pileup2snp aln.raw.pileup --p-value 0.1 > pvalue_filtered_snps
