# Source PDF: Debarshi_day_4_011.pdf (superset of day_2_Debarshi_Adhikary.pdf and day_3_Debarshi.pdf)
# Stage 2: limma differential expression (WT vs Mutant), volcano plot, probe annotation, up/down CSVs

##This will act as an input for the rest of the part

yeast_microchip_rma <- rma(yeast_microchip)

library(limma)

##read the RMA normalized values , printing that we can get to know the exact order

pData_yest = pData(yeast_microchip_rma)
pData_yest

##We are grouping this,it groups first three as WT group and next three as Mutant, if it is in alternate order , we can write as WT, mutant, WT,
#mutant like this… ##Groups =pData_yeast$Treatnment, order matters here, we can see two levels,if we write WT1 by mistake it will take that as
#another level, so instead of two here it will become three and WT1 wll not have replicates.

Groups = c("WT","WT","WT","Mutant","Mutant","Mutant")
Groups

##Purpose: Convert the character vector into a factor ##factors take the parameters in a alphanumeric way,which automatically detects what to
#compare with whom,,,order does not matter here ##relevel tells what we are going to use as reference..it prevents use of that alphanumeric
#way,,,and specifically tells what to use as reference ##but use of factor is important to tell how many levels are here.

Groups = factor(Groups, levels = c("WT","Mutant"))
Groups = relevel(Groups, ref= "WT")
Groups

##We have to create a design matrix, here first three columns, it will take as wt because it has values of 1 for gene 1 and so on and will take the
#avg..

design = model.matrix(~0+Groups)

##design = model.matrix(~0+Groups)..the first 0 ...it tells to take intercept as 0. No intercept model is also called as means model..without 0,,it will keep the intercept and that will be the mean reference model..there I will get only one column that is mutant.
##IN the means model,,if i get 9 in WT and * in Mutant ,,that states the avg intensity,,but in mean reference model the values are in ratio..

design

##Here the first three columns is 1 for WT,,,so while calculating the avg intensities from gene 1 for wild type it will take the first three column’s
#intensities and put in the formula of B0WT+B1MT,,where it will multiply 0 with the values of next three columns which is effectively mutant
#types,,,so we will get avg intensity values for WT for gene1. but the avg will be calculated using lmfit…

fit_y = lmFit(yeast_microchip_rma, design)
contrast.matrix = makeContrasts(GroupsMutant-GroupsWT, levels = design)
##Considering WT as reference, and we are comparing against mutants. This is important, if the values are positive it will say that the genes are upregulated in mutants.
data.fit.con = contrasts.fit(fit_y,contrast.matrix)
##contrast matrix and design matrix does not have the avg intensity values, so we are using fit_y here,,fit_y have the avg intensity values, variances , SE etc of all the columns , defined by the design matrix. And if thre are multiple groups like WT1, MT1, MT2--to compare and specify which two we are going to compare we can use contrast matrix.
fit12 = eBayes(data.fit.con)
##Here by eBayes we are using the moderated t-test. It borrows variance info from all the other genes and then calulate the avg and adjust the variances according to that.
contrast.matrix

fit_y

##makeContrasts is for deifing contrasts, where it is telling to take the differences between mutant and wild type. So it is like specifically defining
#condition of interest for us.

fit12

##from thr eBayes we got log of Odds, pvalue etc…top table is used

diffexp= topTable(fit12, number=Inf,adjust.method="BH")
# adjusted p value calulated basd on Benjamin hotchback correction and we are writing Inf to calulate fold diff and adj p value for all the genes/probes. It prints in the order of higest expression difference to lowest. If I have given 10 in the "number"..it will calculate for only 10 highest,,,
diffexp

diffexp = topTable(fit12,
number = Inf,
adjust.method = "BH")
diffexp

#the fold change is in log2 scale, 1 means 2 fold upreguation and -1 means 2 fold downregualtion ,,,,0.6 means 1.5 fold upregualtion.

# now plot a volcano plot between log10(Adjpvalue) and
# logFC… to get to now that pvalue cut off suits the most for
# our data.
# Upregulated genes
up = subset(diffexp,
logFC >= 0.6 &
adj.P.Val <= 0.05)
up

# we dont have gene id we have prob id,,there is a annotation
# file we have downloded, there it is present.
# Downregulated genes
down = subset(diffexp,
logFC <= -0.6 &
adj.P.Val <= 0.05)
down

diffexp$Status <- "NS"
diffexp$Status[
  diffexp$logFC >= 1 & diffexp$P.Value < 0.05
] <- "Up"
diffexp$Status[
diffexp$logFC <= -1 & diffexp$P.Value < 0.05
] <- "Down"
ggplot(diffexp,
aes(logFC, -log10(adj.P.Val), color = Status)) +
geom_point(size = 2) +
scale_color_manual(values = c(
"Up" = "red",
"Down" = "blue",
"NS" = "grey70"
)) +
geom_vline(xintercept = c(-1, 1), linetype = 2) +
geom_hline(yintercept = -log10(0.05), linetype = 2) +
theme_bw()

annotation <- read.delim(
"~/bio-env/EWE11465_18380_Results_021517.MAS5.TXT",
header = TRUE,
stringsAsFactors = FALSE,
check.names = FALSE
)
colnames(annotation)

annotation_small <- annotation[, c("Probe Set ID",
"Representative Public ID",
"Gene Symbol",
"Gene Title")]
diffexp$`Probe Set ID` <- rownames(diffexp)
annotation_small

library(dplyr)
diffexp_annotated <- left_join(diffexp,
annotation_small,
by = "Probe Set ID")
diffexp_annotated

write.csv(diffexp_annotated,
"diffexp_annotated.csv",
row.names = FALSE)

up_annotated <- subset(diffexp_annotated, logFC >= 1)
up_annotated

down_annotated <- subset(diffexp_annotated, logFC <= -1)

write.csv(up_annotated, "Upregulated_Genes.csv", row.names = FALSE)

write.csv(down_annotated, "Downregulated_Genes.csv", row.names = FALSE)
