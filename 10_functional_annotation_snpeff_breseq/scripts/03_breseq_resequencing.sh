# Source: Debarshi_011_Functional Annotation.pdf  (section 3: breseq on E. coli, SRR12340679; reference U00096 GenBank file)

# first attempt with a local GenBank reference
breseq -r U00096.gbff SRR12340679_1.fastq.gz SRR12340679_2.fastq.gz

# download the reads from SRA and run with 4 threads
prefetch SRR12340679
fasterq-dump --split-files SRR12340679/SRR12340679.sra
breseq -j 4 -r genomic.gbff SRR12340679_1.fastq SRR12340679_2.fastq
# outputs: output.gd (GenomeDiff) and index.html
