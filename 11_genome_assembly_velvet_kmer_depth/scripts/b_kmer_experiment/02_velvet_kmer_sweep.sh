# Source PDF: k_mer_change_Debarshi.pdf
# De Bruijn Graph Construction and Assembly for k=21
pixi run velveth Assembly_k21 21 -shortPaired -fastq -separate SR1.fastq SR2.fastq
pixi run velvetg Assembly_k21 -cov_cutoff 5 -unused_reads yes

# De Bruijn Graph Construction and Assembly for k=31
pixi run velveth Assembly_k31 31 -shortPaired -fastq -separate SR1.fastq SR2.fastq
pixi run velvetg Assembly_k31 -cov_cutoff 5 -unused_reads yes

# De Bruijn Graph Construction and Assembly for k=41
pixi run velveth Assembly_k41 41 -shortPaired -fastq -separate SR1.fastq SR2.fastq
pixi run velvetg Assembly_k41 -cov_cutoff 5 -unused_reads yes
