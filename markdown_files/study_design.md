# Study Design
It is important to determine a study design upfront.

Here are some guidelines to do so:

**N.B.** This page is heavily inspired by the University of Melbourne's bioinformatics tutorials [RNA-Seq Experimental Design](https://www.melbournebioinformatics.org.au/tutorials/tutorials/rna_seq_exp_design/rna_seq_experimental_design/) page. 

_A typical RNA-seq experiment aims to find differentially expressed genes between two conditions (e.g. up and down-regulated genes in knock-out mice compared to wild-type mice). RNA-seq can also be used to discover new transcripts, splice variants, and fusion genes._

## 1. Replicates.

As a general rule, the number of biological replicates should never be below 3. In fact, 6 biological replicates is the recommend amount according to a 2016 study by Schurch et al.[^1]

The results of this study suggest the following should be considered when designing an RNA-seq experiment for DGE:

At least six replicates per condition for all experiments.

At least 12 replicates per condition for experiments where identifying the majority of all DE genes is important.

For experiments with <12 replicates per condition; use edgeR (exact) or DESeq2.

For experiments with >12 replicates per condition; use DESeq.

Apply a fold-change threshold appropriate to the number of replicates per condition between 0.1 ≤ T ≤ 0.5 (see Fig. 2 and the discussion of tool performance as a function of replication).



[^1]: Schurch, N. J., Schofield, P., Gierliński, M., Cole, C., Sherstnev, A., Singh, V., ... & Barton, G. J. (2016). How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?. Rna, 22(6), 839-851.
