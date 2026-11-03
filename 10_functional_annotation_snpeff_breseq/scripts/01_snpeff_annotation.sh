# Source: Debarshi_011_Functional Annotation.pdf  (section 1: SnpEff annotation of the GATK VCF)

# Step 1: inspect the input VCF
bcftools view -H conA_variants_final.vcf | head -3
bcftools view -h conA_variants_final.vcf | grep "^##reference"
bcftools view -h conA_variants_final.vcf | grep "^##contig"
bcftools view -h conA_variants_final.vcf | grep "^##source"

# Step 2: install / verify SnpEff
conda install -c bioconda -c conda-forge snpeff=5.2
snpEff -version

# Steps 3-4: fetch the reference database (yeast: R64-1-1.99)
wget https://snpeff-public.s3.amazonaws.com/databases/v5_0/snpEff_v5_0_R64-1-1.99.zip
unzip snpEff_v5_0_R64-1-1.99.zip

# Step 5: inspect database contig naming (Roman numerals, not NCBI accessions)
snpEff dump -dataDir "$PWD/data" R64-1-1.99 | head -n 50

# Steps 6-7: harmonise chromosome identifiers (NCBI accession -> SnpEff name)
cat << 'EOF2' > chr_map.txt
NC_001133.9	I
NC_001134.8	II
NC_001135.5	III
NC_001136.10	IV
NC_001137.3	V
NC_001138.5	VI
NC_001139.9	VII
NC_001140.6	VIII
NC_001141.2	IX
NC_001142.9	X
NC_001143.9	XI
NC_001144.5	XII
NC_001145.3	XIII
NC_001146.8	XIV
NC_001147.6	XV
NC_001148.4	XVI
NC_001224.1	Mito
EOF2
bcftools annotate --rename-chrs chr_map.txt -O v -o conA_renamed.vcf conA_variants_final.vcf

# Step 8: verify the renamed VCF
bcftools view -H conA_renamed.vcf | head -5

# Step 9: functional annotation
snpEff -dataDir "$PWD/data" R64-1-1.99 conA_renamed.vcf > conA_annotated.vcf

# Steps 10-11: confirm ANN subfields were added
bcftools view -h conA_annotated.vcf | grep "ID=ANN"
bcftools view -H conA_annotated.vcf | head -5

# Step 12: check for unmatched chromosomes (expected output: 0)
bcftools view -H conA_annotated.vcf | grep -c "ERROR_CHROMOSOME_NOT_FOUND"

# Steps 13-14: extract a readable table
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/ANN\n' conA_annotated.vcf > conA_annotation.txt
column -t conA_annotation.txt | less -S

# Steps 15-16: high-priority candidates
# Isolate severe loss-of-function variants
grep "HIGH" conA_annotation.txt > conA_HIGH_impact.txt
# Isolate amino-acid substitutions
grep "missense_variant" conA_annotation.txt > conA_missense.txt
