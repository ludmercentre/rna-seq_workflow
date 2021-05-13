# Differential Gene Expression Analysis
Multiple tools across multiple programming languages can be used to perform Differential Gene Expression (DGE) analysis. When using R and the [Bioconductor](https://bioconductor.org/) packages, the main three are [limma](http://www.bioconductor.org/packages/release/bioc/html/limma.html), [edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html) and [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

In this guide we will be looking at the limma package. The advantage of using limma over the other tools is its higher flexibility, specially when it comes to more complex experimental designs. Since here we will only be looking at a trivial control vs. condition comparison any of the tools can easily be used, it will still be interesting to see how to use limma in this scenario.

limma stands for **li**near **m**odels and differential expression for **m**icro**a**rray data, and as the name entails it was first designed for the analysis of microarray expression data. It has been updated since to also take into allow for the analysis of RNA-Seq data. limma was created by the same team behind edgeR, therefore some of the functions are common to both tools and it is sometimes required to also import edgeR into the limma pipeline.

Two main resources were used to inspire this pipeline, the main documentation for the [limma package]() and the **RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR** paper by Law et al.[^1]

## 1. Preparing Input Data

Before starting DGE analysis it is necessary assemble the STAR output into data that can easily be use to create a [DGEList object](https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/DGEList), i.e., the data format used by the edgeR and limma software to organize and perform DE analysis. As an example of how this can be achieved, please refer to the [gen_star_results.py script](https://github.com/ludmercentre/rna-seq_workflow/blob/master/scripts/dge_analysis/gen_star_results.py) on this tutorial's github page.

Albeit with the possibility of being even more complex, the DGEList object has to contain at least two data frames, **assays** and **colData**. As picture onf this figure from the Huber et al., 2015 Nature methods paper.[^2]
![Figure 2: The integrative data container SummarizedExperiment](https://ludmercentre.github.io/rna-seq_workflow/markdown_files/images/41592_2015_Article_BFnmeth3252_Fig2_HTML.webp)

<br />

---

[^1]: Law, C. W., Alhamdoosh, M., Su, S., Dong, X., Tian, L., Smyth, G. K., & Ritchie, M. E. (2016). RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR. F1000Research, 5.
[^2]: Huber, W., Carey, V. J., Gentleman, R., Anders, S., Carlson, M., Carvalho, B. S., ... & Morgan, M. (2015). Orchestrating high-throughput genomic analysis with Bioconductor. Nature methods, 12(2), 115.
