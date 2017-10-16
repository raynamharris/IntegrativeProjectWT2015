# RNAseq workflow on TACC

## The workflow
* **00_rawdata:** Download the data to scratch on Stampede with `00_gsaf_download.sh`. 
* **01_fastqc:** Evaluate the quality of the reads using the program FastQC.
* **02_filtrimreads:** Filter low quality reads and trim adapters using the program cutadapt.
* **03_fastqc:** Evaluate the quality of the processed reads
* **04_kallisto:** Quantify transcript-level expression using Kallisto