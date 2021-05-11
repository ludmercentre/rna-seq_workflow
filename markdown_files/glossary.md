# Glossary

**Library Size:** 

The RNA that was sequenced is called the RNA library.  With longer read lengths and more accurate sequencing, these days in most organisms, most of the reads are mapped. 

Library size could mean one of two things: the total number of reads that were sequenced in the run or the total number of mapped reads. We will use the total number of mapped reads as the library size in our analyses.  Normalization of RNA-seq data proceeds by computing an "effective" library size, which is computed from the actual library size and the distribution of the counts. ([source](https://online.stat.psu.edu/stat555/node/13/#:~:text=Library%20Size,total%20number%20of%20mapped%20reads.))

**Sequencing Quality:**

A Phred Quality Score above 30 signifies that the probability of an incorrect base call is 1 in 1000. Refer to the table below ([source](https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf)):

| Phred Quality Score | Probability of Incorrect Base Call | Base Call Accuracy |
| :---: | :---:        | :---: |
| 10    | 1 in 10      | 90%   |
| 20    | 1 in 100     | 99%   |
| 30    | 1 in 1,000   | 99.9% |
| 40    | 1 in 10,000  | 99.99%|
| 50    | 1 in 100,000 | 99.999% |
| 60    | 1 in 1,000,000 | 99.9999% |