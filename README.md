This repo contains the experiment that I like to call "Integrative Wild Type 2015" because it reflect that behavior, electrophysiology, and RNAseq data were collected from WT mice in 2015 and analyzed in an integrative fashion. The R markdown files to reproducible run the code are saved in the  RmdFiles subdirectory. Rather than have a single Rmd file for the entire project, the workflow is broken down into pieces. The workflows to recreated the figures are specified by the figure name. All the files are ordered from 00 to 04 to indicate the order of execution.

## Figure 1: Cognitive training alters spatial approach and avoidance behavior

<img src="./figures/figures-01.png" />

A1) Mice were assigned to one of four groups: consistently-trained (red), yoked-consistent (dark grey), conflict-trained (peach), or yoked-conflict (light grey). Mice were placed on the rotating arena (1 rpm) for training sessions that lasted 10 min and was separated by 2-hour intersession interval or overnight (~17 hrs). Both WT (filled circles) and FMR1 (open circles triangles) min in the consistent training and the conflict training paradigms learn to associate the shock zone with extra-maze cues. A2) Consistently trained (red lines) mice make fewer entrances into the shock zone than yoked-mice (dark grey lines) on all training (T1-T6), restest, and retention (Reten.) sessions but not during the habituation session (Hab.). Conflict-trained mice (peach) and their yoked controls (light grey) show a similar pattern except for that mean number of differences between T1 and T4 do not differ between conflict-trained mice. A3) Maximum time spending avoiding the shock zone shows the reciprocal pattern to the mean number of entrances. This opposing pattern of approach and avoidance behavior is reflected in the hierarchical clustering analysis and the principal component analysis. B1) Hierarchical clustering of 40 behaviors show that approach, latency to approach, and speed are primary behavioral variables captured by our video-tracking software. Clustering easily distinguishes trained and yoked animals but does not provide precise temporal resolution. The color scale shows centered z-s for high (yellow) and low (deep purple) values for each group average behavioral observation. B2) A principal component analysis estimates that 36% of the behavioral variation is explained by cognitive training (red and peach versus dark grey and light grey). Among the top five contributing variables are the number of entrances and the max avoidance time.


## Illustration 1:	Visualization of the hippocampal samples collected for this study.   

<img src="./figures/figures-05.png" />

## Figure 4.2:	New description of transcription in the context of learning and confirmation of known subfields-specific expression patterns. 

<img src="./figures/figures-02.png" />

A) Hierarchical clustering of 3,382 differentially expressed genes by any two-way contrast show that this variation is driven by subfield-specific expression (DG: orange, CA3: green, CA1: purple, yoked-consistent: dark grey, consistent: red, yoked-conflict: light grey, conflict: peach). B1) A principal component analysis estimates that 71% of the variation is capture in PC1 and P2, which visually separate the three hippocampal subfields. B2) PC6 differs significantly with treatment group an explains 1% of the variation in gene expression.


## Figure 3: Integration
<img src="./figures/figures-03.png" />

## Figure 2 Supplement: Brain slice


## Figure 3 Supplement: Electrophysiology

<img src="./figures/figures-04.png" />

