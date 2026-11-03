# Source PDF: Debarshi_day_20_8.pdf
pixi run ragtag.py correct GCA_014131755.1_ASM1413175v1_genomic.fna Assembly_step1/contigs.fa

ragtag.py scaffold GCA_014131755.1_ASM1413175v1_genomic.fna ragtag_output/ragtag.correct.fasta

quast.py -o quast_scaffolds ragtag_output/ragtag.scaffold.fasta
