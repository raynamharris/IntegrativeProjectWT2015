[![Binder](http://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/raynamharris/IntegrativeProjectWT2015/master?urlpath=rstudio)
*Click the button to launch a Binder R session. Navigate to the
`scripts` directory and open any `.Rmd` file.*

Overview
--------

This repo contains the experiment that I like to call
“IntegrativeProjectWT2015” because it is an **integrative** analysis of
behavior, electrophysiology, and RNA-seq data collected from **wild
type** mice in **2015**.

This research is undergoing peer-review at Hippocampus. A preprint is
available at
<https://www.biorxiv.org/content/10.1101/2020.02.05.935759v1>.

Please cite this GitHub repository as Rayna M Harris. (2017, November
29). raynamharris/IntegrativeProjectWT2015: GitHub Repository for the
Integrative WT 2015 Project (Version v1.3). Zenodo.
[http://doi.org/10.5281/zenodo.1068356](Rayna%20M%20Harris.%20(2017,%20November%2029).%20raynamharris/IntegrativeProjectWT2015:%20GitHub%20Repository%20for%20the%20Integrative%20WT%202015%20Project%20(Version%20v1.3).%20Zenodo.%20http://doi.org/10.5281/zenodo.1068356)

Bioinformatics Workflow
-----------------------

This project is organized into five main subdirectories: -
[`scripts`](scripts): contains all the `knitr`-based R workflows for
statistical analyses and data visualization (source code is in the
`.Rmd` files, and `.md` files are used to visualize the code and the
results - [`UNIXworkflow`](UNIXworkflow): contains and explanation and
all the UNIX commands used to process the raw sequencing data on the
Stampede cluster at the Texas Advanced Computing Facility -
[`figures`](figures): all the figures created from the scripts -
[`data`](data): all the input data and the results - [`docs`](docs):
presentations generated using R

-   **Part 1: Behavioral analysis**
    -   [scripts/00\_behavior\_wrangle](scripts/00_behavior_wrangle.md):
        behavior data wrangling
    -   [scripts/01\_behavior\_analysis](scripts/01_behavior_analysis.md):
        behavior statistics and data visualization
-   **Part 2: RNA sequencing analysis**
    -   **UNIX-based**
        -   [UNIXworkflow/00\_rawdata](UNIXworkflow/00_rawdata.md):
            Download the data to scratch on Stampede with
            `00_gsaf_download.sh`.
        -   [UNIXworkflow/01\_fastqc](UNIXworkflow/01_fastqc.md):
            Evaluate the quality of the reads using the program FastQC.
        -   [UNIXworkflow/02\_filtrimreads](UNIXworkflow/02_filtrimreads.md):
            Filter low quality reads and trim adapters using the program
            cutadapt.
        -   [UNIXworkflow/03\_fastqc](UNIXworkflow/03_fastqc.md):
            Evaluate the quality of the processed reads
        -   [UNIXworkflow04\_kallisto](UNIXworkflow04_kallisto.md):
            Quantify transcript-level expression using Kallisto
    -   **R-based**
        -   [scripts/00\_rnaseq\_wrangle](scripts/00_rnaseq_wrangle.md):
            converting the kallisto transcript counts to gene counts
        -   [scripts/02\_rnaseqQC](scripts/02_rnaseqQC.md): analyzing
            all the RNA-seq data together
        -   [scripts/03\_rnaseqSubfield](scripts/03_rnaseqSubfield.md):
            analyzing the data for each hippocampal subfield separately
        -   [scripts/04\_correlations.Rmd](scripts/04_correlations.md):
            correlations between genes and beahvior
        -   [scripts/05\_GO](scripts/05_GO.md) an analysis of gene
            ontology.
        -   [scripts/06\_candidates](scripts/05_candidates.md) an
            analysis of candidate memory genes.

![](./figures/00_schematics/figure_workflow.png)
