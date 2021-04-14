The pre-processing steps starts with raw sequencing data (usually in fastq.gz format)

## 1. Pre-Alignment Quality Control (QC):

Multiple programs can be used to perform quality control. Here we will use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://multiqc.info/docs/).

### 1.1. FastQC:

Run the program on the raw read files to return a QC report. 

Depending on your system, you can use a bash command similar to the following:

`fastqc -t 12 raw_data/* -o fastqc_output/`

Here the raw data in the in the **raw_data/** folder and the output (one for each data file will be in the **fastqc_output/** folder). You can find a example of a FastQC output [here]().





