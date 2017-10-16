# IntegrativeProjectWT2015

This repo contains the experiment that I like to call "IntegrativeProjectWT2015" because it is an **integrative** analysis of behavior, electrophysiology, and RNA-seq data collected from **wild type** mice in **2015**. 

## Organization

Three are four main directories for the R scripts, UNIX workflow, figures, and data. 

The file names were chosen to indicated their order of execution in the workflow and provide a brief description of the contents.


### Workflow
 
- Overview
	- **00_methods:** a description of the behavior, electrophysiology, and RNA-seq methods 
- Part 1: Behavior	
	- **01a_behavior_create_dfs:** wrangling the behavior data
	- **01b_behavior_figures:** statistical analysis and behavior data visualization
- Part 2.1: UNIX for RNA-seq
	-  UNIX workflow 
Initial processing of RNA sequencing was performed on the Stampede cluster at the Texas Advanced Computing Facility. My workflow for this is described in the 'UNIXworkflow' directory. 
	- **../UNIXworkflow/00_rawdata:** Download the data to scratch on Stampede with `00_gsaf_download.sh`. 
	- **../UNIXworkflow/01_fastqc:** Evaluate the quality of the reads using the program FastQC.
	- **../UNIXworkflow/02_filtrimreads:** Filter low quality reads and trim adapters using the program cutadapt.
	- **../UNIXworkflow/03_fastqc:** Evaluate the quality of the processed reads
	- **../UNIXworkflow/04_kallisto:** Quantify transcript-level expression using Kallisto
- Part 2.2: R for RNA-seq
	- **02a_rnaseq_makedfs:** converting the Kallisto transcript counts to gene counts and wrangling the categorical data about the samples
	- **02b_rnaseqALL:** analyzing all the RNA-seq data together
	- **02c_rnaseqSubfield:** analyzing the data for each hippocampal subfield separately
	- **02d_rnaseqAvoidance:** combining the two yoked group and the two training groups before analyzing and then analyzing each subfield separately
	- **02e_GO_MWU** an analysis of gene ontology. Note: this directory contains scripts and data
- Part 3: Electrophysiology
	- **03_ephys:** analysis of electrophysiology data	
- Part 4: Integrative analysis
	- **04_integration:** correlations across levels, mostly using PCA data
- Extras
	- **99_colorblindr:** a simulation of how some of my figures will look to colorblind readers
	- **99_forDeep:** a targeted PC analysis of a colleague	

### Results

Quick view of the overall numbers and patterns of differentially expressed beahvioras and genes in response to beahvioral manipulatlions.


#### Figure 1: Experimental overview 

![](./figures/figures-05.png)

#### Figure 2: Cognitive training induces avoidance behavior in conflict and consistently trained animals. See **01b_behavior_figures:Rmd** for the corresponding analysis


<img src="./figures/figures-01.png" />


#### Figure 3: RNA sequencing confirms large differences in DG, CA3, and CA1 hippocampal subfields  
See **02b_rnaseqALL:Rmd** for the corresponding analysis

<img src="./figures/figures-02.png" />


#### Table 1: Numbers of differentially expressed genes by two-way contrast when all data is included in the normalization and quantificaiton. See **02b_rnaseqALL:Rmd** for the corresponding analysis

| contrast | up | down | total |
| --- | --- | --- | --- |
CA3	vs. DG | 1625 | 1361 | 2986
CA1	vs. DG | 968 | 956 | 1924
CA1	vs. CA3 | 651 | 791 | 1442
consistent vs.	yoked consistent | 113 | 7 | 120
yoked conflict vs. yoked consistent | 30 | 1 | 31
conflict vs. yoked conflict | 15 | 24 | 39
conflict vs. consistent | 0 | 0 | 0

#### Table 2: Sample sizes for tissue specific analysis of the resopnse to behavior manipulations in DG, CA3, and CA1 hippocampal subfields. See **02c_rnaseqSubfield.Rmd** for the corresponding analysis


| Training Group | DG | CA3 | CA1 |
|----|:----:|:-----:|:-----:|
| consistent yoked | 4 | 3 | 2 |
| consistent trained | 3 | 2 | 4 |
| conflict yoked | 4 | 3 | 5 |
| conflict trained | 5| 5 | 4 |


#### Figure 4: Cognitive training alters gene expression in DG and CA1 but not CA3. See **02c_rnaseqSubfield.Rmd** for the corresponding analysis

<img src="./figures/figures2-01.png" />


#### Table 1: Numbers of differentially expressed genes by two-way contrast the data is subset by hippocampal tissue before normalization and quantification. See **02c_rnaseqSubfield.Rmd** for the corresponding analysis


| Two-way contrasts between groups | DG | CA3 | CA1 |
|----|:----:|:-----:|:-----:|
| consistent vs. yoked-consisten | 95 | 0 | 196 |
| conflict vs. consistent | 0 | 1 | 0 |
| yoked-conflict vs. yoked-consistent | 2 | 0 | 201 |
| conflict vs. yoked-conflict | 5| 0 | 2 |



#### Figure 5: Training has little effect on electrophysiology. See **003_ephys.Rmd** for the corresponding analysis

<img src="./figures/figures-04.png" />

#### Figure 6: Correlations across levels of biological organization. See **04_integration.Rmd** for the corresponding analysis

<img src="./figures/figures-03.png" />

### Data

This directory contains both raw and intermediate data files. 
- Intermediate data files have alphanumeric prefixes that correspond to the R script that created them. 
- Raw files have only numeric prefix (aka lack an alphabetical character) that indicates whether it is for behavior (01), RNA-seq (02), or ephys (03). 
- Files with more general names were created for public repositories
