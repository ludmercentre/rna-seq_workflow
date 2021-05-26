library(limma)
library(Glimma)
library(edgeR)
library(RColorBrewer)
library(kableExtra)

setwd("C:/Users/path_to_analysis_directory/rna-seq_workflow/scripts/dge_analysis")

#### Import required files
# Experimental variables:
colData_df = read.csv('../../data_files/limma_voom/conditions_df.csv')

# STAR gene count output results dataframe:
assays_df = read.csv('../../data_files/limma_voom/star_results.csv', row.names=1)

# Biomart gene annotation file http://useast.ensembl.org/biomart/martview/:
rowData_df = read.csv("../../data_files/limma_voom/mm10_biomart_annotated_genes.csv")

# dataframe with factors:
colData_df <- data.frame(sample_id = colData_df$sample_id, 
                         subject = substr(colData_df$sample_id, 1, 5), 
                         treat = colData_df$treatment, 
                         sex = colData_df$sex, 
                         roi = colData_df$region
)


# Group
group = factor(colData_df$treat)

# remove.zeros=TRUE gets rid of genes with 0 counts
dge <- DGEList(counts=assays_df, group=group, remove.zeros=TRUE)

# add factors to DGEList object
dge$samples$subject <- factor(colData_df$subject)
dge$samples$treat <- factor(colData_df$treat)
dge$samples$sex <- factor(colData_df$sex)
dge$samples$roi <- factor(colData_df$roi)


# add annotations, same order as in DGE object gene ids:
annotations_df <- annotations_df[match(rownames(dge), annotations_df$ensembl_gene_id), ]
dge$genes <- annotations_df


#### Filtering and Design Matrix
# Keep protein coding genes only:
keep <- dge$genes$gene_biotype == 'protein_coding'
dge <- dge[keep, ]

# Only keep female samples:
# dge <- dge[, which(dge$samples$sex=="F")]

design = model.matrix(~ 0 + group + roi + sex, data = dge$samples)
# Remove "group" from design column names:
colnames(design) = sub("group", "", colnames(design))

# Remove Lowly Expressed Genes:
keep <- filterByExpr(dge, design=design, min.count=10)
dge <- dge[keep,, keep.lib.sizes=FALSE]


#### Background Gene Lists Generation:
# Save background genes list to file:
write.csv(dge$genes, "../../data_files/limma_voom/bg_list.csv", row.names=FALSE)


#### Normalization
# Compute composition normalization factors to scale the raw library sizes with TMM method:
dge <- calcNormFactors(dge, method = "TMM")


#### MDS Plots
par(mfrow=c(1,2))

# ROI (Brain Regions) on dimensions 1 vs. 2:
col.group <- dge$samples$roi # gen color group vector
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1") # assign colors as factors
col.group <- as.character(col.group)
mds <- plotMDS(dge, labels=dge$samples$roi, col=col.group, dim.plot=c(1,2), main="A. ROI dimensions 1 vs. 2")

# Sex on dimensions 4 vs. 5:
col.group <- dge$samples$sex # gen color group vector
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1") # assign colors as factors
col.group <- as.character(col.group)
mds <- plotMDS(dge, labels=dge$samples$sex, col=col.group, dim.plot=c(4,5), main="B. Sex dimensions 4 vs. 5")

# Glimma:
html_filename = "all_samples"
glMDSPlot(dge, labels=dge$samples$group,
          groups=dge$samples[, 4:7], launch=TRUE, html=html_filename)


# Voom and Linear Modeling
## Use voom() to convert the read counts to log2-cpm, with associated weights, ready for linear modelling:
v <- voom(dge, design, plot=TRUE)

## If using vWQW:
# v <- voomWithQualityWeights(dge, design, plot=TRUE)

# Estimate the correlation between the brain region cell lines that were extracted from the same mouse subject.
corfit = duplicateCorrelation(v, design, block = dge$samples$subject)

## The intra cell line correlation will change the voom weights slightly, so we run voom a second time:
v <- voom(dge, design, plot=TRUE, block = dge$samples$subject, correlation = corfit$consensus)

## If using vWQW:
# v <- voomWithQualityWeights(dge, design, plot=TRUE, block = dge$samples$subject, correlation = corfit$consensus)

# Similary, we run update the correlation for the new voom weights:
corfit = duplicateCorrelation(v, design, block = dge$samples$subject)

fit = lmFit(v, design, block = dge$samples$subject, correlation = corfit$consensus)

#### Contrasts
cm = makeContrasts(
  POLvsSAL = POL-SAL,
  levels = design
)

fit2 = contrasts.fit(fit, cm)
fit3 = eBayes(fit2, robust=TRUE)


#### Examine DE Genes
dt <- decideTests(fit3, method="separate")
# To make the output nice:
kable(t(summary(dt)))

top.table <- topTable(fit3, sort.by="P", n=Inf, adjust.method="BH")

# Saving to file:
write.csv(top.table, file="../../data_files/limma_voom/POLvsSAL.csv", row.names=F)

# Partial lists:
write.csv(top.table[which(top.table$adj.P.Val < 0.05 & top.table$logFC > 0),], file="../../data_files/limma_voom/POLvsSAL_UP_DEGs.csv")

write.csv(top.table[which(top.table$adj.P.Val < 0.05 & top.table$logFC < 0),], file="../../data_files/limma_voom/POLvsSAL_DOWN_DEGs.csv")

# In case of multiple comparisons:
for (c in colnames(dt)) {
  
  degs = topTable(fit2, coef=c, n=Inf, sort.by='P', adjust.method='BH')
  
  write.csv(degs, file=paste0("../DE_results_folder/", c, ".csv"), row.names=F)
}