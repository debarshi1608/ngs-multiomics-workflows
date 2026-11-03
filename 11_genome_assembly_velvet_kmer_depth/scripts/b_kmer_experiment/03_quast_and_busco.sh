# Source PDF: k_mer_change_Debarshi.pdf
pixi run quast.py Assembly_k21/contigs.fa Assembly_k31/contigs.fa Assembly_k41/contigs.fa -o quast_kmer_comparison

pixi run busco -i Assembly_k21/contigs.fa -m geno -o busco_k21 -l bacteria_odb12
pixi run busco -i Assembly_k31/contigs.fa -m geno -o busco_k31 -l bacteria_odb12
pixi run busco -i Assembly_k41/contigs.fa -m geno -o busco_k41 -l bacteria_odb12
