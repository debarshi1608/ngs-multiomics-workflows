# ChIP-Seq Differential Binding Analysis

## Overview
Identification and statistical modeling of differentially bound genomic regions using R and Bioconductor toolsets.

## Pipeline & Methods
* **Peak Calling & Consensus:** Processed aligned BAM files to establish consensus peak sets.
* **Affinity Binding Analysis:** Utilized `DiffBind` to model binding affinity across experimental conditions.
* **Differential Binding Statistics:** Implemented negative binomial generalized linear models via `DESeq2` to isolate statistically significant differential binding sites.

## Deliverables
* Full Analysis Report: [`report/ChIPseq_DiffBind_Analysis.pdf`](./report/ChIPseq_DiffBind_Analysis.pdf)
* Analysis Scripts: `scripts/`
* Processed Outputs: `results/`
