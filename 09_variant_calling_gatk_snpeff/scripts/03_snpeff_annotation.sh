# Source: Debarshi_011_Variant_call.pdf  (steps 14-15: harmonise chromosome names, SnpEff annotation)

# 14.1 Verify SnpEff database names
snpEff dump R64-1-1.99 | head -n 50
# 14.2 Check VCF chromosome names
grep -v "^#" conA_variants_final.filtered.vcf | cut -f1 | sort -u
# 14.3 Rename with bcftools (chr_map.txt = tab-delimited map, e.g. NC_001133.9 -> SnpEff name)
bcftools annotate --rename-chrs chr_map.txt \
   -o conA_variants_snpEff_input.vcf -Ov \
   conA_variants_final.filtered.vcf

# 15. Functional annotation with SnpEff (+ HTML stats)
snpEff -stats conA_snpEff_summary.html \
  R64-1-1.99 \
  conA_variants_snpEff_input.vcf > conA_variants_annotated.vcf
