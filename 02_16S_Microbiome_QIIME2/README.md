# 16S rRNA Microbiome Profiling Pipeline

## Overview
End-to-end amplicon sequencing analysis pipeline for characterization of microbial community structures.

## Pipeline & Methods
* **Quality Control & Denoising:** Processed raw paired-end FASTQ reads using DADA2 for error correction and Amplicon Sequence Variant (ASV) generation.
* **Taxonomic Classification:** Trained naive Bayes classifiers against reference databases (SILVA / Greengenes) for taxonomic assignment.
* **Diversity Metrics:** Evaluated alpha-diversity (Shannon, Faith's PD) and beta-diversity (Bray-Curtis, UniFrac distances) visualized via ordination plots (PCoA).

## Deliverables
* Workflow Report: [`report/16S_rRNA_QIIME2_Analysis.pdf`](./report/16S_rRNA_QIIME2_Analysis.pdf)
