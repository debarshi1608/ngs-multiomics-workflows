# Source PDF: depth_change_Debarshi.pdf
pixi run quast.py Assembly_3M/contigs.fa Assembly_5M/contigs.fa Assembly_10M/contigs.fa -o quast_depth_comparison

pixi run busco -i Assembly_3M/contigs.fa -m geno -o busco_3M -l bacteria_odb12
pixi run busco -i Assembly_5M/contigs.fa -m geno -o busco_5M -l bacteria_odb12
pixi run busco -i Assembly_10M/contigs.fa -m geno -o busco_10M -l bacteria_odb12
