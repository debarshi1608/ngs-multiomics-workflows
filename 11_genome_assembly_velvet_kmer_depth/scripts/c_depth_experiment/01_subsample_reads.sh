# Source PDF: depth_change_Debarshi.pdf
# Subsampled primary FASTA files into 3M, 5M, and 10M read pair datasets
pixi run seqtk sample -s100 SR1.fastq 3000000 > SR1_3M.fastq
pixi run seqtk sample -s100 SR2.fastq 3000000 > SR2_3M.fastq


pixi run seqtk sample -s100 SR1.fastq 5000000 > SR1_5M.fastq
pixi run seqtk sample -s100 SR2.fastq 5000000 > SR2_5M.fastq


pixi run seqtk sample -s100 SR1.fastq 10000000 > SR1_10M.fastq
pixi run seqtk sample -s100 SR2.fastq 10000000 > SR2_10M.fastq
