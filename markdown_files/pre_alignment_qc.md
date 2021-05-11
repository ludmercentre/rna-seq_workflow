The pre-processing steps starts with raw sequencing data (usually in fastq.gz format)

## 1. Pre-Alignment Quality Control (QC):

Multiple programs can be used to perform quality control. Here we will use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://multiqc.info/docs/).

### 1.1. FastQC:

Run the program on the raw read files to return a QC report. 

Depending on your system, you can use a bash command similar to the following:

`fastqc -t 12 raw_data/* -o fastqc_output/`

Here the raw data is located in the **raw_data/** folder and the output (one for each data file will be in the **fastqc_output/** folder). You can find a example of a FastQC output [here](https://ludmercentre.github.io/rna-seq_workflow/data_files/fastqc_output/raw_data/NS.1223.004.NEBNext_dual_i7_A1---NEBNext_dual_i5_A1.03_32_vHIP_R1_fastqc.html).

Sometimes there are information contained in the read file name, in this case it contains the company NEB (New England Biolabs) and the fact that the read are paired-end. R1 means this is the first of the two pairs. See the [Study Design](https://ludmercentre.github.io/rna-seq_workflow/markdown_files/study_design.html) page.

### 1.1. MultiQC:

Run the program from anywhere within your project folder. MultiQC will automatically recognize file extensions and organize its report with the right files in each section.

Depending on your system, you can use a bash command similar to the following:

`multiqc . -o multiqc_output/`

The output will consist of only one html file and will be located in the **multiqc_output/** folder). You can find a example of a MultiQC output [here](https://ludmercentre.github.io/rna-seq_workflow/data_files/multiqc_output/multiqc_report.html).

Sometimes there are information contained in the read file name, in this case it contains the company NEB (New England Biolabs) and the fact that the read are paired-end. R1 means this is the first of the two pairs. See the [Study Design](https://ludmercentre.github.io/rna-seq_workflow/markdown_files/study_design.html) page.

