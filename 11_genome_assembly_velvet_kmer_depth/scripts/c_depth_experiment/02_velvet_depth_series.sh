# Source PDF: depth_change_Debarshi.pdf
# De Bruijn Graph Construction and Assembly for 3M Reads (k=31)
pixi run velveth Assembly_3M 31 -shortPaired -fastq -separate SR1_3M.fastq SR2_3M.fastq
pixi run velvetg Assembly_3M -cov_cutoff 5 -unused_reads yes


# De Bruijn Graph Construction and Assembly for 5M Reads (k=31)
pixi run velveth Assembly_5M 31 -shortPaired -fastq -separate SR1_5M.fastq SR2_5M.fastq
pixi run velvetg Assembly_5M -cov_cutoff 5 -unused_reads yes


# De Bruijn Graph Construction and Assembly for 10M Reads (k=31)
pixi run velveth Assembly_10M 31 -shortPaired -fastq -separate SR1_10M.fastq SR2_10M.fastq
pixi run velvetg Assembly_10M -cov_cutoff 5 -unused_reads yes
