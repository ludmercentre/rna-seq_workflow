# Pathway Enrichment Analysis

Results from modern omics analyses are often made up of long lists of genes, requiring an impractically large amount of manual literature research to interpret. Pathway enrichment analysis helps researchers gain mechanistic insight into these lists by helping identify biological pathways that are more enriched in a gene list than what would be expected by chance.

Using the 2019 nature protocols paper by Reimand et al. as inspiration, the tutorial will be split in 3 main sections: definition of a gene list from omics data, determination of statistically enriched pathways, and visualization and interpretation of the results.[^1] Note that these principles can be applied to diverse types of omics data.

In this tutorial we will be presenting 4 different tools. 3 for conducting alternative methods of pathway enrichment analysis and on for visualization. These tools, at the time of this being written, would be the go-to, recommended ones to use. [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) and [Enrichr](https://maayanlab.cloud/Enrichr/) to perform gene enrichment analysis, [GSEA](http://www.gsea-msigdb.org/gsea/index.jsp) (Gene-set enrichment analysis) and [Cytoscape](https://cytoscape.org/) for visualizing the results.

## 1. Definition of a gene list
Gene pathway analysis can be applicable to the analysis of lists of genes or biomolecules from any organism derived from large-scale data, including proteomics, genomics, epigenomics and gene-regulation studies, they can come from gene expression microarrays, quantitative proteomics, germline and somatic genome sequencing and global DNA methylation assays amongst others. They usually contain the features (genes) as rows and metadata specific to the experiment that generated it. For DEG analysis there are usually two kinds of gene lists, ranked and unranked. The difference is trivial, an unranked list can be any output list of biological features. The ranked list is ordered based on a specific parameter, for DEG results this is often the FDR-adjusted p-value or the FC (fold change) sign (+ or -)*-log10(p-value), it could also be by the fold change sign only amongst many others. For an example of this ouput see the ___ result file from the STAR comparison described in the [Differential Gene Expression Analysis](https://ludmercentre.github.io/rna-seq_workflow/markdown_files/DGE_analysis.html) tutorial page.

Outputs a ranked list of genes. 
By p-value of the significance of differential expression, Q-value (a.k.a. adjusted p-value, corrected for multiple testing across all genes), effect size and direction of expression change (upregulated genes are positive and at the top of the list and downregulated genes are negative and at the bottom of the list, often expressed as log-transformed fold-change (logFC). 
(e.g., −log10 p-value multiplied by the sign of logFC)

List:
E.g.: All somatically mutated genes in a tumor identified by exome sequencing, or all proteins that interact with a bait in a proteomics experiment.
Ranked List:
Can be partial or filtered by a particular threshold (e.g., FDR-adjusted P value <0.05 and fold-change >2)

## 2. Determination of statistically enriched pathways
Usually performed using Gene Ontology (GO) (biological processes, cellular components and molecular functions) terms or annotations from other databases (e.g., KEGG (everything), Reactome (pathways), TRANSFAC (transcription factors), etc.)

Gene enrichment analysis.
Takes list or partial ranked list as input.
The p-value of the enrichment of a pathway is computed using a Fisher’s exact test and multiple-test correction is applied.
Multiple software: g:Profiler, DAVID, GOrilla, etc.


Gene-set enrichment analysis (GSEA).
In some experiments comparing two conditions, there might be no or only a few genes that are significantly DE, but group of genes slightly DE.
Takes full ranked list as input.
Kolmogorov-Smirnov (KS) test to assign enrichment score (ES) to a group of genes with multiple testing corrections applied.
Performed by the GSEA software from the Broad Institute (2005).

Statistical Tests:
Fisher's exact test compares the expected number of significant genes at random to the observed number of significant genes to arrive at a probability.
The KS test compares the distribution of gene p-values expected at random to the observed distribution of the gene p-values to arrive at a probability.


<br />

---

[^1]: Reimand, J., Isserlin, R., Voisin, V., Kucera, M., Tannus-Lopes, C., Rostamianfar, A., ... & Bader, G. D. (2019). Pathway enrichment analysis and visualization of omics data using g: Profiler, GSEA, Cytoscape and EnrichmentMap. Nature protocols, 14(2), 482-517.
