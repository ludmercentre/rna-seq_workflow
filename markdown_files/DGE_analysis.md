# Differential Gene Expression Analysis
Multiple tools across multiple programming languages can be used to perform Differential Gene Expression (DGE) analysis. When using R and the [Bioconductor](https://bioconductor.org/) packages, the main three are [limma](http://www.bioconductor.org/packages/release/bioc/html/limma.html), [edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html) and [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

In this guide we will be looking at the limma package. The advantage of using limma over the other tools is its higher flexibility, specially when it comes to more complex experimental designs. Since here we will only be looking at a trivial control vs. condition comparison any of the tools can easily be used, it will still be interesting to see how to use limma in this scenario.

limma stands for **li**near **m**odels and differential expression for **m**icro**a**rray data, as the name entails it was first designed for the analysis of microarray expression data. It has been updated since to also take into allow for the analysis of RNA-Seq data. limma was created by the same team behind edgeR, therefore some of the functions are common to both tools and it is sometimes required to also import edgeR into the limma pipeline.

Two main resources were used to inspire this pipeline, the main documentation for the [limma package]() and the **RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR** paper by Law et al.[^1]



<br />

---

[^1]:Law, C. W., Alhamdoosh, M., Su, S., Dong, X., Tian, L., Smyth, G. K., & Ritchie, M. E. (2016). RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR. F1000Research, 5.
