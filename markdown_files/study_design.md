# Study Design
It is important to determine a study design upfront.

Here are some guidelines to assist in doing so:

**N.B.** This page is inspired in part by the University of Melbourne's bioinformatics tutorials [RNA-Seq Experimental Design](https://www.melbournebioinformatics.org.au/tutorials/tutorials/rna_seq_exp_design/rna_seq_experimental_design/) page and the [Canadian Center for Computational Genomics (C3G)](https://www.computationalgenomics.ca/) online tutorials.

_A typical RNA-seq experiment aims to find differentially expressed genes between two conditions (e.g. up and down-regulated genes in knock-out mice compared to wild-type mice). RNA-seq can also be used to discover new transcripts, splice variants, and fusion genes._

<!-- <br/>

## Table of Contents
1. [Replicates](#1replicates)
2. [Example2](#example2)
3. [Third Example](#third-example)
4. [Fourth Example](#fourth-examplehttpwwwfourthexamplecom)

<br/> -->

## 1. Replicates

**Technical replicates:** Sequences derived from the same library/sample (lanes, flow cells, etc.)

**Biological replicates:** Sequences derived from different samples with the same phenotype/genotype or experimental condition

As a general rule, the number of biological replicates should never be below 3. Biological replicates are important. Technical replicates are often unnecessary.[^1] In fact, 6 biological replicates is the recommend amount according to a 2016 study by Schurch et al.[^2] 

*The results of this study suggest the following should be considered when designing an RNA-seq experiment for DGE:*

* *At least six replicates per condition for all experiments.*

* *At least 12 replicates per condition for experiments where identifying the majority of all DE genes is important.*

* *For experiments with <12 replicates per condition; use edgeR (exact) or DESeq2.*

* *For experiments with >12 replicates per condition; use DESeq.*

* *Apply a fold-change threshold appropriate to the number of replicates per condition between 0.1 ≤ T ≤ 0.5 (see Fig. 2 and the discussion of tool performance as a function of replication).*

From the [Canadian Center for Computational Genomics](https://www.computationalgenomics.ca/):
* *More biological replicates per group are recommended if
samples are expected to have high variation (live specimens) or the effect is expected to be subtle.*
* *More technical replicates are recommended if higher coverage per sample is required*

## 2. Sequencing Depth (Coverage)

Sequencing depth or coverage is the average quantity of reads mapped to a genome, usually defined theoretically as 'LN/G', where L is the read length, N is the number of reads and G is the haploid genome length.[^3] Unlike with replicates, increasing read depth does not impact statistical power in a differential expression analysis. Differential expression analysis generally requires pretty low coverage (unless very low expressed genes are of interest) whereas any kind of de novo build or splicing analysis requires higher read depth. These variables, if possible, can be determined by the researcher the goal of the study needs to be considered. Hence the importance of conducting a pilot study first with a few samples to determine what depth and how many replicates (biological or technical) would be needed to gather enough power in the main analysis. There are published studies providing guideline on doing so:

*With regard to testing of various experimental designs, this work strongly suggests that greater power is gained through the use of biological replicates relative to library (technical) replicates and sequencing depth. Strikingly, sequencing depth could be reduced as low as 15% without substantial impacts on false positive or true positive rates.*[^4]

Illumina provides an online [Sequencing Coverage Calculator](https://support.illumina.com/downloads/sequencing_coverage_calculator.html) as well as [documentation](https://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_coverage_calculation.pdf) to better understand and estimate sequencing coverage.

Some general guidelines from C3G:

| Type of experiment | No. of mapped reads (per sample) | Length of reads
| :---: | :---: | :---: |
| Gene expression profiling | 10-25 million | 50-75 bp |
| Differential analysis and alternative splicing | 40-60 million | 75 bp |
| Transcriptome assembly | > 100 million | > 75 bp |
| miRNA and sRNA analysis | 1-5 million (targeted) | 50 bp (single-end) |

## 3. Read Length

Defining read length is also an important decision. Usually sequencers have a technical limit to the amount of reads they can sequence. Short read sequencing provides high accuracy but only small fragments of data, which gives an incomplete picture.

## 4. Read Types

Part of a study design is also the choice of read types, they can be either single or paired-end. *For basic differential expression analysis RNA-seq experiments, single-end sequencing is recommended to obtain gene transcript counts. In more advanced experiments, paired-ends are useful for determining transcript structure and discovering splice variants.*[^1]

## 5. Strandedness
Difference in alignment and read distribution.
*With a non-directional (unstranded) protocol, there is no way to identify whether a read originated from the coding strand or its reverse complement. Non-directional protocols allow mapping of a read to a genomic location, but not the direction in which the RNA was transcribed. They are therefore used to count transcripts for known genes, and are recommended for basic RNA-seq experiments. Directional protocols (stranded) preserve strand information and are useful for novel transcript discovery.*

## 6. Total RNA Needed and RIN:
*Abundant RNA’s dominate (rRNA), high amounts of unprocessed RNA and genomic DNA. Many sequencing centres such as AGRF recommend at least 250ng of total RNA for RNA sequencing. It is possible to go as low as 100ng of total RNA, but results are not guaranteed. The quality of RNA is also important when making libraries. A RNA Integrity Number (RIN) is a number from 1 (poor) to 10 (good) and can indicate how much degradation there is in the sample. A poor score can lead to over representation at the 3’ end of the transcript and low yield. Samples with low RIN scores (below 8) are not recommended for sequencing. Care should also be taken to ensure RIN is consistent between conditions to avoid confounding this technical effect with the biological question.*[^1]

## 7. Enrichment Method:

*Ribosomal RNA makes up >95% of total cellular RNA, so a preparation for RNA-seq must either enrich for mRNA using poly-A enrichment, or deplete rRNA. The amount of RNA needed for each method differs.*

### 7.1 PolyA selection: 
* Limited transcript representation, low unprocessed RNA and genomic DNA.
* Poly-A enrichment is recommended for most standard RNA-seq experiments, but will not provide information about microRNAs and other non-coding RNA species.
* Minimum of RNA 100ng needed.

### 7.2 rRNA reduction: 
* *Abundant RNAs de-emphasized, still high amounts of unprocessed RNA and genomic DNA.*
* In general, ribo-depleted RNA-seq data will contain more noise, however, the protocol is recommended if you have poor or variable quality of RNA as the 3’ bias of poly-A enrichment will be more pronounced with increased RNA degradation.
* Minimum of 200ng of RNA recommended.

### 7.3 cDNA capture: 
* Targeted transcript representation (using cDNA), all other RNA molecules de-emphasized

## 8. Multiplexing and Batch Effect:
* Multiplexing is an approach to sequence multiple samples in the same sequencing lane. By sequencing all samples in the same lane, multiplexing can also minimise bias from lane effects.

<br />

---

<br />

[^1]: Melbourne Bioinformatics. (2017). RNA-Seq Experimental Design. https://www.melbournebioinformatics.org.au/tutorials/tutorials/rna_seq_exp_design/rna_seq_experimental_design/
[^2]: Schurch, N. J., Schofield, P., Gierliński, M., Cole, C., Sherstnev, A., Singh, V., ... & Barton, G. J. (2016). How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?. Rna, 22(6), 839-851.
[^3]: Sims, D., Sudbery, I., Ilott, N. E., Heger, A., & Ponting, C. P. (2014). Sequencing depth and coverage: key considerations in genomic analyses. Nature Reviews Genetics, 15(2), 121-132.
[^4]: Robles, J. A., Qureshi, S. E., Stephen, S. J., Wilson, S. R., Burden, C. J., & Taylor, J. M. (2012). Efficient experimental design and analysis strategies for the detection of differential expression using RNA-Sequencing. BMC genomics, 13(1), 1-14.
[^5]: Robles, J. A., Qureshi, S. E., Stephen, S. J., Wilson, S. R., Burden, C. J., & Taylor, J. M. (2012). Efficient experimental design and analysis strategies for the detection of differential expression using RNA-Sequencing. BMC genomics, 13(1), 1-14.
