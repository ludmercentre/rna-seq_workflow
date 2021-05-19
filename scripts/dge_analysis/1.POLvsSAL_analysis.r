library(limma)
library(Glimma)
library(edgeR)
library(statmod)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(knitr)
# library(kableExtra)
# library(biomaRt)

#### Script Results ####
## POLvsSAL
# 14618    72
# |         | Down| NotSig|  Up|
# |:--------|----:|------:|---:|
# |POLvsSAL |  567|  13730| 321|

## POLvsSAL with voomWithQualityWeights
# |         | Down| NotSig|  Up|
# |:--------|----:|------:|---:|
# |POLvsSAL | 1072|  12640| 906|

## POLvsSAL with roi
# 14971    72
# |         | Down| NotSig|  Up|
# |:--------|----:|------:|---:|
# |POLvsSAL |  394|  14310| 267|
  
# POLvsSAL with sex
# 14857    72
# |         | Down| NotSig|  Up|
# |:--------|----:|------:|---:|
# |POLvsSAL |  626|  13864| 367|

## POLvsSAL with roi and sex
# 15025    72
# |         | Down| NotSig|  Up|
# |:--------|----:|------:|---:|
# |POLvsSAL |  560|  14225| 240|

## POLvsSAL with roi and sex : voomWQW
# 15025    72
# |         | Down| NotSig|  Up|
# |:--------|----:|------:|---:|
# |POLvsSAL |  971|  13350| 704|

## POLvsSAL with roi and sex - 1 outlier (X01_09_ACC)
# 15048    71
# |         | Down| NotSig|  Up|
# |:--------|----:|------:|---:|
# |POLvsSAL |  482|  14399| 167|

## POLvsSAL with roi and sex - 1 outlier (X01_09_ACC) : voomWQW
# 15048    71
# |         | Down| NotSig|  Up|
# |:--------|----:|------:|---:|
# |POLvsSAL |  962|  13418| 668|

## POLvsSAL with roi and sex - 2 outliers (X01_09_ACC & X08_09_dHIP)
# 15046    70
# |         | Down| NotSig|  Up|
# |:--------|----:|------:|---:|
# |POLvsSAL |  462|  14426| 158|

## POLvsSAL with roi and sex - 2 outliers (X01_09_ACC & X08_09_dHIP) : voomWQW
# 15046    70
# |         | Down| NotSig|  Up|
# |:--------|----:|------:|---:|
# |POLvsSAL |  929|  13483| 634|


## Conclusion:
########################           


setwd("C:/Users/pedro/OneDrive - McGill University/MIA_analysis/limma_voom/scripts/")

#### Import required files
# Experimental variables:
targets <- read.csv('../data/conditions_df.csv')

# STAR gene count output results dataframe:
star_results_df = read.csv('../data/star_results.csv', row.names=1)
# star_results_df = read.csv('../data/star_results_no_mm.csv', row.names=1)

# Biomart gene annotation file http://useast.ensembl.org/biomart/martview/:
annotations_df = read.csv("../data/mm10_biomart_annotated_genes.csv")

# # Import Biomart annotation with Entrez Gene ids
# annotations_df = read.table(file = "../data/mm10_biomart_annotated_genes_16-07-2020_with_entrez_ids.tsv", sep = '\t', quote="", header = TRUE)
# # ids_df <- annotations_df[c("ensembl_gene_id", "entrez_gene_id")]

# dataframe with factors:
targets = data.frame(sample_id = targets$sample_id, subject = substr(targets$sample_id, 1, 5), treat = targets$treatment, sex = targets$sex, roi = targets$region)


#### Group
group = factor(targets$treat)


#### DGEList
# remove.zeros=TRUE gets rid of genes with 0 counts
dge <- DGEList(counts=star_results_df, group=group, remove.zeros=TRUE)

# add factors to DGEList object
dge$samples$subject <- targets$subject
dge$samples$treat <- targets$treat
dge$samples$sex <- targets$sex
dge$samples$roi <- targets$roi


#### Annotations
# Same order as in DGE object gene ids:
annotations_df <- annotations_df[match(rownames(dge), annotations_df$ensembl_gene_id), ]
dge$genes <- annotations_df


#### Filtering and Design Matrix
# ## Filter out 1 outliers samples:
# dge <- dge[, which(colnames(dge$counts)!="X01_09_ACC")]
# 
# ## Filter out 2 outliers samples:
# dge <- dge[, which(colnames(dge$counts)!="X01_09_ACC" & colnames(dge$counts)!="X08_09_dHIP")]

#### Filtering
# ## Keep protein coding genes only:
keep <- dge$genes$gene_biotype == 'protein_coding'
dge <- dge[keep, ]

#### Design Matrix
group = dge$samples$group

# POLvsSAL
# design = model.matrix(~ 0 + group, data = dge$samples)
# # POLvsSAL with roi
# design = model.matrix(~ 0 + group + roi, data = dge$samples)
# # POLvsSAL with sex
# design = model.matrix(~ 0 + group + sex, data = dge$samples)
# POLvsSAL with roi and sex
design = model.matrix(~ 0 + group + roi + sex, data = dge$samples)
# design = model.matrix(~ 0 + group, data = dge$samples)

# Remove "group" from design column names:
colnames(design) = sub("group", "", colnames(design))


# ## Remove Lowly Expressed Genes:
keep <- filterByExpr(dge, design=design, min.count=10)
dge <- dge[keep,, keep.lib.sizes=FALSE]

# At least 6 samples with a count of 10 or higher
# dge <- dge[rowSums(dge$counts >= 10) >= 6,]

print(dim(dge)) # [1] 15048    71


#### Normalization
# Compute composition normalization factors to scale the raw library sizes with TMM method:
dge <- calcNormFactors(dge, method = "TMM")


#### MDS Plots
par(mfrow=c(1,2))
# ROI (Brain Regions) on dimensions 1 vs. 2:col.group <- dge$samples$roi # gen color group vector
col.group <- dge$samples$roi # gen color group vector
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1") # assign colors as factors
col.group <- as.character(col.group)
mds <- plotMDS(dge, labels=dge$samples$roi, col=col.group, dim.plot=c(1,2), main="A. ROI dimensions 1 vs. 2")

# Sex on dimensions 4 vs. 5:
levels(dge$samples$sex) <- c(levels(dge$samples$sex), "E") # Add level to sex to color differently the excluded samples
dge$samples['X01_09_ACC',]$sex <- "E" # Replace sex with 'E' for excluded sample
col.group <- dge$samples$sex # gen color group vector
dge$samples['X01_09_ACC',]$sex <- "F" # Assign it back to 'F'
dge$samples$sex <- factor(dge$samples$sex) # Reset factor

# levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1") # assign colors as factors
levels(col.group) <-  c("#E41A1C", "#377EB8", "#984EA3") # purple color for excluded
# levels(col.group) <-  c("#E41A1C", "#377EB8", "#A65628") # brown color for excluded
col.group <- as.character(col.group)
mds <- plotMDS(dge, labels=dge$samples$sex, col=col.group, dim.plot=c(4,5), main="B. Sex dimensions 4 vs. 5")


#### Voom and Linear Modeling
## Use voom() to convert the read counts to log2-cpm, with associated weights, ready for linear modelling:
# v <- voom(dge, design, plot=TRUE)

## If using vWQW:
v <- voomWithQualityWeights(dge, design, plot=TRUE)

# Estimate the correlation between the brain region cell lines that were extracted from the same mouse subject.
corfit = duplicateCorrelation(v, design, block = dge$samples$subject)

## The intra cell line correlation will change the voom weights slightly, so we run voom a second time:
# v <- voom(dge, design, plot=TRUE, block = dge$samples$subject, correlation = corfit$consensus)

## If using vWQW:
v <- voomWithQualityWeights(dge, design, plot=TRUE, block = dge$samples$subject, correlation = corfit$consensus)

# Similary, we run update the correlation for the new voom weights:
corfit = duplicateCorrelation(v, design, block = dge$samples$subject)

fit = lmFit(v, design, block = dge$samples$subject, correlation = corfit$consensus)
# fit = lmFit(v, design)

#### Contrasts
cm = makeContrasts(
  POLvsSAL = POL-SAL,
  levels = design
)


#### Fit Model
fit2 = contrasts.fit(fit, cm)
fit2 = eBayes(fit2, robust=TRUE)


#### Examine DE Genes
dt <- decideTests(fit2, method="separate")
# dt <- decideTests(fit2, method="global")
kable(t(summary(dt)))

# degs[which(degs$logFC > 0), ] [1:50,] # up reg top 50
# degs[which(degs$logFC < 0), ] [1:50,] # down reg top 50

#### Saving gene lists to file
for (c in colnames(dt)) {
  
  print(c)
  
  # BH is the default fdr method.
  degs = topTable(fit2, coef = c, n=Inf, sort.by = 'P', adjust.method='BH')
  
  # Add Annotations:
  # NB: degs is already annotated
  # annotations = dge$genes
  # annotations = annotations[match(rownames(degs), annotations$ensembl_gene_id), ]
  
  # deg_data <- cbind(degs, annotations)
  
  write.csv(degs, file=paste0("../DE/1.POLvsSAL_with_roi_sex-1_outliers/", c, ".csv"), row.names=F)
  
  # write.csv(deg_data, file=paste0("../results/limma_results/POLvsSAL/", c, ".csv"))
  # write.csv(deg_data, file=paste0("../results/limma_results/POLvsSAL_with_roi/", c, ".csv"))
  # write.csv(deg_data, file=paste0("../results/limma_results/POLvsSAL_with_sex/", c, ".csv"))

  # write.csv(deg_data, file=paste0("../results/final_results/voom/1.POLvsSAL_with_roi_sex/", c, ".csv"))
  # write.csv(deg_data, file=paste0("../results/final_results/voomWithQualityWeights/1.POLvsSAL_with_roi_sex/", c, ".csv"))
  # 
  # write.csv(deg_data, file=paste0("../results/final_results/voom/1.POLvsSAL_with_roi_sex-1_outliers/", c, ".csv"))
  # write.csv(deg_data, file=paste0("../results/final_results/voomWithQualityWeights/1.POLvsSAL_with_roi_sex-1_outliers/", c, ".csv"))
  # 
  # write.csv(deg_data, file=paste0("../results/final_results/voom/1.POLvsSAL_with_roi_sex-2_outliers/", c, ".csv"))
  # write.csv(deg_data, file=paste0("../results/final_results/voomWithQualityWeights/1.POLvsSAL_with_roi_sex-2_outliers/", c, ".csv"))
}

# CAMERA:
# load(system.file("extdata", "mouse_c2_v5p1.rda", package = "RNAseq123"))
Mm.c2 <- load(file="C:/Users/pedro/OneDrive - McGill University/MIA_GO_enrichment_analysis/CAMERA/data/mouse_c2_v5p2.rdata")

idx <- ids2indices(Mm.c2,id=v$genes$entrez_gene_id)

cam.POLvsSAL <- camera(v,idx,design,contrast=cm[,1])
head(cam.POLvsSAL,5)

## If more than on comparison,look at different columns of cm:
cam.BasalvsML <- camera(v,idx,design,contrast=contr.matrix[,2])
## etc.

barcodeplot(efit$t[,3], index=idx$LIM_MAMMARY_LUMINAL_MATURE_UP, 
            index2=idx$LIM_MAMMARY_LUMINAL_MATURE_DN, main="LPvsML")
cam.BasalvsML <- camera(v,idx,design,contrast=contr.matrix[,2])