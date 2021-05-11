The pre-processing steps starts with raw sequencing data (usually in fastq.gz format)

## 1. Pre-Alignment Quality Control (QC):

Multiple programs can be used to perform quality control. Here we will use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://multiqc.info/docs/).

### 1.1. FastQC:

Run the program on the raw read files to return a QC report. 

Depending on your system, you can use a bash command similar to the following:

`fastqc -t 12 raw_data/* -o fastqc_output/`

Here the raw data is located in the **raw_data/** folder and the output (one for each data file will be in the **fastqc_output/** folder). You can find a example of a FastQC output [here](https://ludmercentre.github.io/rna-seq_workflow/data_files/fastqc_output/raw_data/NS.1223.004.NEBNext_dual_i7_A1---NEBNext_dual_i5_A1.03_32_vHIP_R1_fastqc.html).

As stated in the name, the read file here is of 

## 4. Sequencing Quality

A Phred Quality Score above 30 signifies that the probability of an incorrect base call is 1 in 1000. Refer to the table below[^5]:

| Phred Quality Score | Probability of Incorrect Base Call | Base Call Accuracy |
| :---: | :---:        | :---: |
| 10    | 1 in 10      | 90%   |
| 20    | 1 in 100     | 99%   |
| 30    | 1 in 1,000   | 99.9% |
| 40    | 1 in 10,000  | 99.99%|
| 50    | 1 in 100,000 | 99.999% |
| 60    | 1 in 1,000,000 | 99.9999% |

