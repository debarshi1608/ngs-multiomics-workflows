# Source: diffbind_analysis_Debarshi_011.pdf  (annotation + motif analysis of the significant differential regions)

# Download the GFF annotation
wget -O saccharomyces_cerevisiae_S288C.gff.gz "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz"
# 1. Check file size
ls -lh saccharomyces_cerevisiae_S288C.gff.gz
# 2. Decompress GFF file
gunzip saccharomyces_cerevisiae_S288C.gff.gz
# 3. Extract unique chromosome/sequence identifiers
grep -v '^#' saccharomyces_cerevisiae_S288C.gff | cut -f1 | sort -u | head -15

# (R step: convert significant regions CSV to BED - see 04_diffbind_analysis.R, last section)

# HOMER genomic region annotation
annotatePeaks.pl DiffBind_significant_regions_FDR0.05.bed reference.fna -gff3 saccharomyces_cerevisiae_S288C.gff > DiffBind_significant_regions_HOMER_annotated.txt

# HOMER motif analysis
findMotifsGenome.pl DiffBind_significant_regions_FDR0.05.bed reference.fna motif_analysis -size 200
