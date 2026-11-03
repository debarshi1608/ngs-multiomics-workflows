# Microarray differential expression: WT vs Mutant (affy RMA + limma)

**Aim.** Process Affymetrix CEL files, check raw/log/RMA-normalised intensities, then test WT vs Mutant for differential expression with limma (moderated t-test), plot a volcano plot, attach probe annotation and export up- and down-regulated probe lists.
**Data / organism.** Six CEL files `EWE11465_0001_WT1`, `_0002_WT2`, `_0003_WT3`, `_0004_drna1`, `_0005_drna2`, `_0006_drna3` (3 WT, 3 Mutant), read from `/home/ibab/bio-env`; the object is named `yeast_microchip`. Note: the package warnings mention `ecoli2cdf` and the annotation file (`EWE11465_18380_Results_021517.MAS5.TXT`) lists IDs such as `MG1655_mrr_b4351` and `EDL933_Z0312`.
**Tools.** R, affy, limma, tidyverse (ggplot2, dplyr).


## Procedure
1. Read CEL files with `ReadAffy(celfile.path = celpath, compress = TRUE)`, take `exprs()` as a data frame and preview it (`scripts/01_affy_qc_rma.R`).
2. Extract Perfect Match (`affy::pm`) and Mismatch (`affy::mm`) probe signals; plot PM vs MM for WT1 (`scripts/01_affy_qc_rma.R`).
3. Natural-log transform PM/MM (`log()`); compare histograms of raw vs log PM for WT1; boxplot per sample (`scripts/01_affy_qc_rma.R`).
4. Normalise with `rma(yeast_microchip)` (background correction, normalisation, expression calculation); boxplot of RMA values (`scripts/01_affy_qc_rma.R`).
5. Compare original / log / RMA data side by side: boxplots, WT1-vs-dRNA1 scatter plots, and MA plots (A = mean, M = difference) (`scripts/01_affy_qc_rma.R`).
6. Define groups `c("WT","WT","WT","Mutant","Mutant","Mutant")` as a factor with WT as reference; build the no-intercept design `model.matrix(~0+Groups)`; `lmFit`, `makeContrasts(GroupsMutant-GroupsWT)`, `contrasts.fit`, `eBayes` (`scripts/02_limma_de_annotation.R`).
7. `topTable(fit12, number=Inf, adjust.method="BH")` for all probes; subset up (logFC >= 0.6 and adj.P.Val <= 0.05) and down (logFC <= -0.6 and adj.P.Val <= 0.05) (`scripts/02_limma_de_annotation.R`).
8. Volcano plot (logFC vs -log10 adj.P.Val; Up: logFC >= 1 and P.Value < 0.05, Down: logFC <= -1 and P.Value < 0.05; guide lines at logFC +/-1 and adj.P = 0.05) (`scripts/02_limma_de_annotation.R`).
9. Read the annotation file, keep `Probe Set ID`, `Representative Public ID`, `Gene Symbol`, `Gene Title`, `left_join` onto the results by probe ID, write `diffexp_annotated.csv`, and write `Upregulated_Genes.csv` (logFC >= 1) and `Downregulated_Genes.csv` (logFC <= -1) (`scripts/02_limma_de_annotation.R`).

## Results

**Raw intensities (first rows of `exprs()`, WT1 / WT2 / WT3 / drna1 columns).** Probe 1: 75 / 95 / 83 / 71; probe 2: 4530 / 4222 / 4416 / 4636; probe 3: 74 / 91 / 78 / 85; probe 4: 4918 / 4362 / 4714 / 4557.

PM vs MM signal for WT1 (PM green, MM pink):

![PM vs MM signal, WT1](results/fig01_pm_vs_mm_signal.png)

Raw vs log-transformed PM histogram, WT1:

![Raw vs log PM histogram](results/fig02_pm_raw_vs_log_histogram.png)

Boxplot of log-transformed PM for the six samples:

![Boxplot log transformed](results/fig03_boxplot_log_transformed.png)

**RMA values (first rows, WT1 / WT2 / WT3):** 1759068_at 3.413791 / 3.015579 / 3.437268; 1759069_at 9.020832 / 9.126874 / 9.048895; 1759070_s_at 7.771956 / 7.623827 / 7.432955; 1759071_s_at 6.148551 / 5.695259 / 5.836691; 1759072_s_at 3.670370 / 3.580104 / 3.662050; 1759073_at 2.925485 / 3.253104 / 2.651791.

![Boxplot RMA normalized](results/fig04_boxplot_rma_normalized.png)

Original vs log vs RMA boxplots:

![Original, log, RMA boxplots](results/fig05_boxplots_original_log_rma.png)

WT1 vs dRNA1 scatter plots (original, log, normalised):

![Scatter WT1 vs dRNA1](results/fig06_scatter_wt1_vs_drna1.png)

MA plots (original, log, normalised):

![MA plots](results/fig07_ma_plots.png)

**limma model.** Design columns `GroupsWT`, `GroupsMutant` (samples 1-3 = WT, 4-6 = Mutant); `df.residual` 4 per probe; eBayes `df.prior` 2.791656, `s2.prior` 0.03315848, `var.prior` 112.7717, `proportion` 0.01; `df.total` 6.791656; 10,203 further probes beyond the five printed. Example contrast coefficients (Mutant - WT): 1759068_at -0.24249026, 1759069_at -0.52167654, 1759070_s_at 1.29906777, 1759071_s_at -0.27378326, 1759072_s_at -0.04937133.

**Top 10 probes by `topTable` (BH-adjusted).** The table viewer reports 10,000 rows.

| Probe | logFC | AveExpr | t | P.Value | adj.P.Val | B |
|---|---|---|---|---|---|---|
| 1766647_s_at | -5.933308 | 5.554242 | -54.193405 | 3.280059e-10 | 3.348284e-06 | 11.564632 |
| AFFX-Ec-U57609-2_s_at | 7.155444 | 6.466314 | 39.603087 | 2.744066e-09 | 1.400571e-05 | 10.714404 |
| 1769071_s_at | -2.820358 | 10.795486 | -25.627512 | 5.180489e-08 | 1.241138e-04 | 8.933794 |
| 1765239_s_at | -2.983624 | 9.268889 | -25.254639 | 5.717466e-08 | 1.241138e-04 | 8.861778 |
| 1766702_s_at | 4.506791 | 9.075793 | 24.208565 | 7.599151e-08 | 1.241138e-04 | 8.649837 |
| 1764012_s_at | 2.838458 | 8.984151 | 23.819899 | 8.472505e-08 | 1.241138e-04 | 8.567176 |
| 1762388_s_at | 3.329537 | 5.449938 | 23.408790 | 9.523913e-08 | 1.241138e-04 | 8.477307 |
| 1763982_s_at | 5.848421 | 7.352742 | 23.187517 | 1.015131e-07 | 1.241138e-04 | 8.427870 |
| 1766364_s_at | 2.570210 | 8.176837 | 22.673929 | 1.179888e-07 | 1.241138e-04 | 8.310147 |
| 1762983_s_at | 3.976983 | 10.713461 | 22.572749 | 1.215848e-07 | 1.241138e-04 | 8.286452 |

**Gene counts.** logFC >= 0.6 and adj.P.Val <= 0.05: 636 probes up; logFC <= -0.6 and adj.P.Val <= 0.05: 554 probes down. After annotation, `up_annotated` (logFC >= 1) has 438 rows.

![Volcano plot](results/fig08_volcano_plot.png)

**Annotation.** The annotation file has 14 columns (`Probe Set ID`, six MAS5 signal columns, `Representative Public ID`, `Gene Title`, `Gene Symbol`, Gene Ontology BP/CC/MF, `Pathway`). Examples: 1759069_at = MG1655_mrr_b4351, gene symbol `mrr`; 1759071_s_at = MG1655_b2503, `ECs3365 /// yfgF`; 1759073_at = MG1655_ycfK_b1154, `stfP`. The joined table has 11 columns.



## Files
```
02_microarray_yeast_limma/
  README.md
  scripts/  01_affy_qc_rma.R  02_limma_de_annotation.R
  results/  fig01 ... fig08 (.png)
  reports/  day_2_Debarshi_Adhikary.pdf  day_3_Debarshi.pdf  Debarshi_day_4_011.pdf
```
