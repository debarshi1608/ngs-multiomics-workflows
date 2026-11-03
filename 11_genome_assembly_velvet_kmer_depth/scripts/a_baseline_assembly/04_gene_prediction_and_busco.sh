# Source PDF: Debarshi_day_20_8.pdf
augustus --species=E_coli_K12 ragtag_output/ragtag.scaffold.fasta > scaffold_1.gff

gffread scaffold_1.gff -g ragtag_output/ragtag.scaffold.fasta -w multi.fasta

busco -i ragtag_output/ragtag.scaffold.fasta -m geno -o busco_output_bact -l bacteria_odb12
