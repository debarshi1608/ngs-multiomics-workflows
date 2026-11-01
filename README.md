# Computational Multi-Omics & NGS Workflows

[![Status](https://img.shields.io/badge/status-active-success.svg)]()
[![Focus](https://img.shields.io/badge/focus-bioinformatics_%7C_genomics-blue.svg)]()

Production-grade bioinformatics pipelines, analytical scripts, and technical reports spanning epigenomics, metagenomics, transcriptomics, and functional enrichment.

---

## 📂 Repository Architecture

| Module | Description | Key Tools / Packages | Report |
| :--- | :--- | :--- | :--- |
| **[01_ChIPseq_DiffBind](./01_ChIPseq_DiffBind)** | ChIP-seq consensus peaks & affinity binding | `DiffBind`, `DESeq2`, R | [PDF Report](./01_ChIPseq_DiffBind/report/ChIPseq_DiffBind_Analysis.pdf) |
| **[02_16S_Microbiome_QIIME2](./02_16S_Microbiome_QIIME2)** | Amplicon denoising & microbial diversity | `QIIME 2`, DADA2, Python | [PDF Report](./02_16S_Microbiome_QIIME2/report/16S_rRNA_QIIME2_Analysis.pdf) |
| **[03_Functional_Annotation](./03_Functional_Annotation)** | GO enrichment & biological network mapping | hypergeometric tests, KEGG | [PDF Report](./03_Functional_Annotation/report/Functional_Enrichment_Analysis.pdf) |
| **[04_RNAseq_and_NGS_Labs](./04_RNAseq_and_NGS_Labs)** | scRNA-seq, Bisulfite-seq & Variant Calling | `Seurat`, `Bismark`, Samtools | [View Reports](./04_RNAseq_and_NGS_Labs/report/) |

---

## 🛠️ Technical Stack & Environment
* **Languages:** R, Python, Bash / Shell Scripting
* **Core Libraries:** Bioconductor ecosystem (`DESeq2`, `DiffBind`, `Seurat`), `QIIME 2`, `scikit-learn`
* **Version Control & Reproducibility:** Git, Conda, Pixi
