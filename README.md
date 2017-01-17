This repo contains the experiment that I like to call "Integrative Wild Type 2015" because it reflect that behavior, electrophysiology, and RNAseq data were collected from WT mice in 2015 and analyzed in an integrative fashion. The title of the paper is: 

The R markdown files to reproducible run the code are saved in [the  RmdFiles subdirectory](./RmdFiles). Rather than have a single Rmd file for the entire project, the workflow is broken down into pieces. The workflows to recreated the figures are specified by the figure name. 

- Figure 1
- Figure 2
- Figure 3

Individual figures were combined into one panel using Adobe Illustrator. The workflows to reproduce the tidy data are outlined in separate files for behavior, physiology, and gene expression. 

- Behavior Data Tidying
- Transcript quantification with Kallisto
- Differential gene expression quantification with DESeq2 


### Talk figures

In a talk for the eisen lab, I presented some slide with the following figures:

![paradigm](.figures/methods_behavior.png "paradigm")

![paradigm](.figures/allenslicephoto-01.tif "allen brain and slice photo")

![spatialtraining1](./RmdFiles/01_behavior_files/figure-markdown_strict/unnamed-chunk-3-3.png "Spatially trained mice (orange & browns boxes) avoid the shock zone")

![spatialtraining2](./RmdFiles/01_behavior_files/figure-markdown_strict/unnamed-chunk-3-2.png "But there is a lot of variability")

![spatialtraining3](./RmdFiles/01_behavior_files/figure-markdown_strict/unnamed-chunk-5-3.png "Heatmap")

![spatialtraining1](./RmdFiles/02_rnaseq_files/heatmap.png "Spatially trained mice (orange & browns boxes) avoid the shock zone")
