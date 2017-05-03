This repo contains the experiment that I like to call "Integrative Wild Type 2015" because it reflect that behavior, electrophysiology, and RNAseq data were collected from WT mice in 2015 and analyzed in an integrative fashion. The title of the paper is: 


## Workflow

The R markdown files to reproducible run the code are saved in the  RmdFiles subdirectory. Rather than have a single Rmd file for the entire project, the workflow is broken down into pieces. The workflows to recreated the figures are specified by the figure name. 

Click "summary" to view a markdown file with documention and figure. 
Click "Rmd file" to view the bioinformatic workflow.

- Overview - summer the experimental paradigm ([summary](./RmdFiles/01_schematics.md), made with Adobe Illustrator)
- Behavior - analysis of the leanring and memory behavioral data ([analysis](./RmdFiles/01a_beahvior_create_dfs.md), [figures](./RmdFiles/01b_behavior_figures.md)
- Gene expression - anlaysis of all hippocampal transcriptomic data ([analysis](./RmdFiles/02a_rnaseq_makedfs.md), [figures](./RmdFiles/02b_rnaseq_makefigures.md))
- Gene expression - anlaysis of DG, CA1, and CA3 separately ([analysis & figures](./RmdFiles/02c_brainregionspecific.md)

## Working Abstract

The collaborate research project began with hypothesis that memory is supported by the modulation of specific molecular pathways that regulate long-term changes in synaptic function within hippocampal circuits. To address this hypothesis, it is necessary to measure behavior, molecular pathwasy, and synaptic function in the same animals. Here, I present an analysis of animal behavior during the active place avoidance task, synaptic plasticity using ex-vivo slice physiology (in progress, not yet included), and molecular activty using RNA-sequencing. Ongoing analysis will explore how variation in learning and memory behavior is related to variation in synaptic plasticity and molecular activity. 
