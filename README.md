# NGS & Multi-omics Workflows

Eleven computational biology analyses (microarray, RNA-seq, splicing, ChIP-seq, 16S microbiome, single-cell, methylation, variant calling, genome assembly). Every project has the same layout:

```
<project>/
├── README.md    aim, data, tools, procedure, results (with figures), interpretation
├── scripts/     all commands and code, one file per stage
├── results/     figures and screenshots extracted from the report
└── reports/     the original PDF report(s)
```

| # | Project | Main tools | Interpretation written in source? |
|---|---|---|---|
| 01 | [R distributions basics](01_r_distributions_basics) | R | no (raw notebook) |
| 02 | [Microarray differential expression](02_microarray_yeast_limma) | affy, RMA, limma | method notes only |
| 03 | [RNA-seq: alignment → edgeR → GO](03_rnaseq_differential_expression) | Bowtie2, featureCounts, edgeR, gprofiler2 | GO figure only |
| 04 | [Alternative splicing](04_alternative_splicing_sgseq) | HISAT2, SAMtools, SGSeq | yes |
| 05 | [ChIP-seq peaks, DiffBind, motifs](05_chipseq_macs2_diffbind_motifs) | Bowtie2/BWA, MACS3, DiffBind, HOMER, MEME | yes |
| 06 | [16S microbiome](06_16s_microbiome_qiime2) | QIIME 2, DADA2, SILVA, ANCOM-BC | yes |
| 07 | [Single-cell RNA-seq (PBMC 3k)](07_scrna_seurat_pbmc3k) | Seurat | yes |
| 08 | [Bisulfite methylation](08_bisulfite_methylation_bismark) | Bismark, IGV | yes |
| 09 | [Variant calling and annotation](09_variant_calling_gatk_snpeff) | BWA, GATK, SnpEff | yes |
| 10 | [Functional annotation, VarScan, breseq](10_functional_annotation_snpeff_breseq) | SnpEff, VarScan, breseq | yes |
| 11 | [Genome assembly: baseline, k-mer, depth](11_genome_assembly_velvet_kmer_depth) | Velvet, QUAST, BUSCO, RagTag | yes |


