# Pathway Enrichment Analysis

Results from modern omics analyses are often made up of long lists of genes, requiring an impractically large amount of manual literature research to interpret. Pathway enrichment analysis helps researchers gain mechanistic insight into these lists by helping identify biological pathways that are more enriched in a gene list than what would be expected by chance.

Using the 2019 nature protocols paper by Reimand et al. as inspiration, the tutorial will be split in 3 main sections: definition of a gene list from omics data, determination of statistically enriched pathways, and visualization and interpretation of the results.[^1] Note that these principles can be applied to diverse types of omics data.

In this tutorial we will be presenting 4 different tools. 3 for conducting alternative methods of pathway enrichment analysis and on for visualization. These tools, at the time of this being written, would be the go-to, recommended ones to use. [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) and [Enrichr](https://maayanlab.cloud/Enrichr/) to perform gene enrichment analysis, [GSEA](http://www.gsea-msigdb.org/gsea/index.jsp) (Gene-set enrichment analysis) and [Cytoscape](https://cytoscape.org/) for visualizing the results.

## 1. Definition of a gene list
Gene pathway analysis can be applicable to the analysis of lists of genes or biomolecules from any organism derived from large-scale data, including proteomics, genomics, epigenomics and gene-regulation studies, they can come from gene expression microarrays, quantitative proteomics, germline and somatic genome sequencing and global DNA methylation assays amongst others. They usually contain the features (genes) as rows and metadata specific to the experiment that generated it. For DEG analysis there are usually two kinds of gene lists, ranked and unranked. The difference is trivial, an unranked list can be any output list of biological features. The ranked list is ordered based on a specific parameter, for DEG results this is often the FDR-adjusted p-value or the FC (fold change) sign (+ or -)*-log10(p-value), it could also be by the fold change sign only amongst many others. For an example of this ouput see the ___ result file from the STAR comparison described in the [Differential Gene Expression Analysis](https://ludmercentre.github.io/rna-seq_workflow/markdown_files/DGE_analysis.html) tutorial page.

Examples of unranked lists:
All somatically mutated genes in a tumor identified by exome sequencing, or all proteins that interact with a bait in a proteomics experiment.

Examples of ranked list:
Ouput list of genes ranked y p-value of the significance of differential expression, Q-value (a.k.a. adjusted p-value, corrected for multiple testing across all genes), effect size and direction of expression change (upregulated genes are positive and at the top of the list and downregulated genes are negative and at the bottom of the list, often expressed as log-transformed fold-change (logFC). 
(e.g., −log10 p-value multiplied by the sign of logFC)

**N.B.:** Lists can also be **partial**, i.e., containing only outputs filtered by a particular threshold (e.g., FDR-adjusted P value <0.05 and fold-change >2).

## 2. Determination of statistically enriched pathways
The statistical detection of pathways or other groups of genes showing an over-representation in the gene list of interest in contrast to what would be expected by chance. Usually performed using Gene Ontology (GO) (biological processes, cellular components and molecular functions) terms or annotations from other databases (e.g., KEGG (everything), Reactome (pathways), TRANSFAC (transcription factors), etc.)

There are two presently used method to conduct enrichment analysis, gene enrichment analysis and gene-set enrichment analysis (GSEA).

### 2.1. Gene enrichment analysis.
Takes list or partial ranked list as input. The p-value of the enrichment of a pathway is computed using a Fisher’s exact test and multiple-test correction is applied (usually Benjamini‐Hochberg's method or a custom method). The results are usually computed by referencing to the many databases mentioned above. Multiple software: g:Profiler, DAVID, Enrichr, GOrilla, etc.

In the scope of this tutorial we recommend using two tools for the gene enrichment analysis:

#### 2.1.1. Enrichr
Enrichr is an interesting tool because it hosts and very comprehensive online tool. It also hosts API capabilities for access programmatically. It's very simple to use, accepts multiple types of gene ids and queries multiple databases. It currently contains a collection of ∼400,000 annotated gene sets organized into ∼300 gene‐set libraries. [^2] It also immediately provides the user with tables and different plots that can be used for publications or even easily customized in software such as Adobe Illustrator or the free an open-source Inkscape.

Check the online tool [here](https://maayanlab.cloud/Enrichr/) and the latest publication [here](https://doi.org/10.1002/cpz1.90).[^2]

#### 2.1.2. g:Profiler
g:Profiler uses very similar calculations as Enrichir, it also has a very user friendly [online tool](https://biit.cs.ut.ee/gprofiler/gost) and is also available as an R package with very comprehensive [documentation](https://biit.cs.ut.ee/gprofiler/page/r) and offers [API access](https://biit.cs.ut.ee/gprofiler/page/apis). Check out the [gProfiler_script.R](https://github.com/ludmercentre/rna-seq_workflow/blob/master/scripts/pathway_enrichment_analysis/gProfiler/gProfiler_script.R) for an example of how to use it in R. g:Profiler and Enrichr complete each other well because while Enrichr provides quick results and ready made graphs, g:Profiler is much more customizable, accepts all kinds of list (this include background lists, see why this often recommended and necessary here ___). g:Profiler also supports almost all genomes however, whereas Enrichr only works for a selected few. g:Profiler also supports custom gene set files ([GMT](http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) (Gene Matrix Transposed file format)). You can find their latest publication [here](https://doi.org/10.1093/nar/gkz369).[^3]


### 2.2. Gene-set enrichment analysis (GSEA).
In some experiments comparing two conditions, there might not be any genes or only a few genes that are significantly over-represented in pathways or gene sets, but this doesn't mean that groups of genes aren't enriched. See the figure below:
![diab2](https://ludmercentre.github.io/rna-seq_workflow/markdown_files/images/diab2.png)


 GSEA Takes full ranked list as input. It uses the Kolmogorov-Smirnov (KS) test to assign enrichment score (ES) to a group of genes with multiple testing corrections applied.
Performed by the GSEA software from the Broad Institute (2005).

Statistical Tests:
Fisher's exact test compares the expected number of significant genes at random to the observed number of significant genes to arrive at a probability.
The KS test compares the distribution of gene p-values expected at random to the observed distribution of the gene p-values to arrive at a probability.


<br />

---

[^1]: Reimand, J., Isserlin, R., Voisin, V., Kucera, M., Tannus-Lopes, C., Rostamianfar, A., ... & Bader, G. D. (2019). Pathway enrichment analysis and visualization of omics data using g: Profiler, GSEA, Cytoscape and EnrichmentMap. Nature protocols, 14(2), 482-517.
[^2]: Xie, Z., Bailey, A., Kuleshov, M. V., Clarke, D. J., Evangelista, J. E., Jenkins, S. L., ... & Ma'ayan, A. (2021). Gene Set Knowledge Discovery with Enrichr. Current Protocols, 1(3), e90.
[^3]: Raudvere, U., Kolberg, L., Kuzmin, I., Arak, T., Adler, P., Peterson, H., & Vilo, J. (2019). g: Profiler: a web server for functional enrichment analysis and conversions of gene lists (2019 update). Nucleic acids research, 47(W1), W191-W198.
