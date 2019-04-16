Setup
-----

    library(plyr)
    library(dplyr)
    library(reshape2)
    library(Hmisc)
    library(corrplot)
    library(viridis)
    library(cowplot)
    library(pheatmap)

    source("figureoptions.R")

    knitr::opts_chunk$set(fig.path = '../figures/04_integration/')

Import Data
-----------

    # behavior
    behavior <- read.csv("../data/01a_behavior.csv", header = T)

    # physiology
    ephys3 <- read.csv("../data/03_ephys.csv", header = T)
    ephys3 <- ephys3[,(c(1:6))]
    names(ephys3)[3] <- "Pre_Potentiation"
    names(ephys3)[4] <- "Early_Potentiation"
    names(ephys3)[5] <- "Late_Potentiation"
    names(ephys3)[6] <- "Max_fEPSP"

    # clearnup the rosetts data and filter extraneous samples
    rossetta <- read.csv("../data/00_rossettastone.csv", header = F)
    names(rossetta)[1] <- "organism"
    names(rossetta)[2] <- "ID"
    names(rossetta)[3] <- "Region"
    names(rossetta)[4] <- "RNAseqID"
    names(rossetta)[5] <- "R1filename"
    rossetta$R2filename <- rossetta$R1filename
    rossetta$R2filename <- gsub("R1", "R2", rossetta$R2filename)
    head(rossetta)

    ##   organism     ID Region   RNAseqID                          R1filename
    ## 1  15-142C 15142C    CA1 142C-CA1-S 142C_CA1_S_S19_L003_R1_001.fastq.gz
    ## 2  15-142C 15142C     DG  142C-DG-S  142C_DG_S_S21_L003_R1_001.fastq.gz
    ## 3  15-143A 15143A    CA3 143A-CA3-1 143A_CA3_1_S35_L002_R1_001.fastq.gz
    ## 4  15-143A 15143A     DG  143A-DG-1  143A_DG_1_S36_L002_R1_001.fastq.gz
    ## 5  15-143B 15143B    CA1 143B-CA1-1 143B_CA1_1_S37_L002_R1_001.fastq.gz
    ## 6  15-143B 15143B     DG  143B-DG-1  143B_DG_1_S38_L002_R1_001.fastq.gz
    ##                            R2filename
    ## 1 142C_CA1_S_S19_L003_R2_001.fastq.gz
    ## 2  142C_DG_S_S21_L003_R2_001.fastq.gz
    ## 3 143A_CA3_1_S35_L002_R2_001.fastq.gz
    ## 4  143A_DG_1_S36_L002_R2_001.fastq.gz
    ## 5 143B_CA1_1_S37_L002_R2_001.fastq.gz
    ## 6  143B_DG_1_S38_L002_R2_001.fastq.gz

    colData <- read.csv("../data/02a_colData.csv", header = T) 
    metadata <- full_join(colData, rossetta)

    ## Joining, by = c("RNAseqID", "ID")

    ## Warning: Column `RNAseqID` joining factors with different levels, coercing
    ## to character vector

    ## Warning: Column `ID` joining factors with different levels, coercing to
    ## character vector

    names(metadata)[1] <- "samplename"
    metadata$title <- as.factor(paste(metadata$ID,metadata$Region, metadata$Group, sep=" "))
    names(metadata)[6] <- "sourcename"
    metadata$char1 <- "Mus musculus"
    metadata$char2 <- "C57BL/6"
    metadata$mol <- "RNA"
    metadata$des <- " "
    metadata$processes <- as.factor(paste(metadata$samplename,"/abundance.txt"))
    metadata$process2 <- "IntegrativeWT2015ColData.csv"
    str(metadata)

    ## 'data.frame':    54 obs. of  19 variables:
    ##  $ samplename: chr  "143A-CA3-1" "143A-DG-1" "143B-CA1-1" "143B-DG-1" ...
    ##  $ Mouse     : Factor w/ 18 levels "15-143A","15-143B",..: 1 1 2 2 3 4 4 5 5 5 ...
    ##  $ Punch     : Factor w/ 3 levels "CA1","CA3","DG": 2 3 1 3 1 1 3 1 2 3 ...
    ##  $ Group     : Factor w/ 3 levels "conflict","consistent",..: 1 1 3 3 2 3 3 1 1 1 ...
    ##  $ Conflict  : Factor w/ 2 levels "Conflict","NoConflict": 1 1 1 1 2 2 2 1 1 1 ...
    ##  $ sourcename: Factor w/ 4 levels "conflict","shocked",..: 1 1 2 2 3 4 4 1 1 1 ...
    ##  $ ID        : chr  "15143A" "15143A" "15143B" "15143B" ...
    ##  $ APA2      : Factor w/ 4 levels "conflict","consistent",..: 1 1 3 3 2 4 4 1 1 1 ...
    ##  $ organism  : Factor w/ 21 levels "15-142C","15-143A",..: 2 2 3 3 4 5 5 6 6 6 ...
    ##  $ Region    : Factor w/ 3 levels "CA1","CA3","DG": 2 3 1 3 1 1 3 1 2 3 ...
    ##  $ R1filename: Factor w/ 54 levels "142C_CA1_S_S19_L003_R1_001.fastq.gz",..: 3 4 5 6 7 10 11 12 13 14 ...
    ##  $ R2filename: chr  "143A_CA3_1_S35_L002_R2_001.fastq.gz" "143A_DG_1_S36_L002_R2_001.fastq.gz" "143B_CA1_1_S37_L002_R2_001.fastq.gz" "143B_DG_1_S38_L002_R2_001.fastq.gz" ...
    ##  $ title     : Factor w/ 54 levels "15142C CA1 NA",..: 3 4 5 6 7 10 11 12 13 14 ...
    ##  $ char1     : chr  "Mus musculus" "Mus musculus" "Mus musculus" "Mus musculus" ...
    ##  $ char2     : chr  "C57BL/6" "C57BL/6" "C57BL/6" "C57BL/6" ...
    ##  $ mol       : chr  "RNA" "RNA" "RNA" "RNA" ...
    ##  $ des       : chr  " " " " " " " " ...
    ##  $ processes : Factor w/ 54 levels "142C-CA1-S /abundance.txt",..: 3 4 5 6 7 10 11 12 13 14 ...
    ##  $ process2  : chr  "IntegrativeWT2015ColData.csv" "IntegrativeWT2015ColData.csv" "IntegrativeWT2015ColData.csv" "IntegrativeWT2015ColData.csv" ...

    metadata <- dplyr::select(metadata, samplename, title, sourcename, organism, char1, char2, mol, des, processes, R1filename, R2filename, process2)

    #write.csv(metadata, "../data/00_metadata.csv")

    # clearnup the rosetts data and filter extraneous samples
    names(rossetta)[1] <- "Mouse"
    names(rossetta)[2] <- "ID"
    names(rossetta)[3] <- "Region"
    names(rossetta)[4] <- "RNAseqID"
    names(rossetta)[5] <- "R1filename"
    rossetta$R1filename <- NULL
    rossetta <- rossetta %>% dplyr::filter(Mouse != "15-100", Mouse != "15-101", Mouse != "15-147")
    rossetta <- unique(rossetta[ , c(1:2,4) ]) # for joining keop only

    ## slim behavior ephy to top 5 pcs and rename the columsn
    behaviorpca <- read.csv("../data/01a_scoresdf.csv", header = T)
    behaviorpca <- behaviorpca[(c(1:4,6,35:36))]
    names(behaviorpca)[names(behaviorpca)=="PC1"] <- "Behavior_PC1"
    names(behaviorpca)[names(behaviorpca)=="PC2"] <- "Behavior_PC2"
    names(behaviorpca)[names(behaviorpca)=="PC3"] <- "Behavior_PC3"
    names(behaviorpca)[names(behaviorpca)=="PC4"] <- "Behavior_PC4"
    names(behaviorpca)[names(behaviorpca)=="PC6"] <- "Behavior_PC6"

    behaviorpca <- behaviorpca %>% dplyr::filter(ID != "15148", ID !=  "15140A", ID !=  "15140B", ID !=  "15140C", ID !=  "15140D", ID !=  "15141C", ID !=  "15141D", ID !=  "15142C", ID !=  "15142D", ID !=  "15142A", ID !=  "15142B", ID !=  "15145C", ID !=  "15145C", ID !=  "15145D", ID !=  "15147A", ID !=  "15147B", ID !=  "15148C", ID !=  "15148D")


    pcadata <- read.csv("../data/02a_pcadata.csv", header = T)
    pcadata$Punch <- ifelse(grepl("DG", pcadata$group), "DG", 
                                            ifelse(grepl("CA3", pcadata$group), "CA3","CA1"))
    names(pcadata)[names(pcadata)=="name"] <- "RNAseqID"
    pcadata <- pcadata[c(11:17,2:10)]

    head(rossetta)

    ##     Mouse     ID   RNAseqID
    ## 1 15-142C 15142C 142C-CA1-S
    ## 2 15-142C 15142C  142C-DG-S
    ## 3 15-143A 15143A 143A-CA3-1
    ## 4 15-143A 15143A  143A-DG-1
    ## 5 15-143B 15143B 143B-CA1-1
    ## 6 15-143B 15143B  143B-DG-1

    #widen then length RNAseq data so each row is an animals
    pcadatabyregion <- dplyr::left_join(pcadata, rossetta, by = "RNAseqID")

    ## Warning: Column `RNAseqID` joining factors with different levels, coercing
    ## to character vector

    names(pcadatabyregion)

    ##  [1] "group"     "Punch"     "avoidance" "RNAseqID"  "wrap"     
    ##  [6] "PunchAPA"  "APA2"      "PC1"       "PC2"       "PC3"      
    ## [11] "PC4"       "PC5"       "PC6"       "PC7"       "PC8"      
    ## [16] "PC9"       "Mouse"     "ID"

    pcadatabyregion <- pcadatabyregion[(c(1:7,17:18,8:16))]

    pcadatabyregion <- melt(pcadatabyregion, id = c(1:9))
    pcadatabyregion$RegionPC <- as.factor(paste(pcadatabyregion$Punch, pcadatabyregion$variable, sep="_"))
    pcadatabyregion <- dcast(pcadatabyregion, Mouse ~ RegionPC)

    alldata <- left_join(pcadatabyregion, rossetta , by="Mouse")
    alldata <- left_join(alldata, behaviorpca, by="ID")

    ## Warning: Column `ID` joining factors with different levels, coercing to
    ## character vector

    alldata <- left_join(alldata, ephys3, by="Mouse")

    ## Warning: Column `Mouse` joining factors with different levels, coercing to
    ## character vector

    alldataslim <- dplyr::filter(alldata, APA2 %in% c("conflict","consistent"))

    names(alldataslim)

    ##  [1] "Mouse"              "CA1_PC1"            "CA1_PC2"           
    ##  [4] "CA1_PC3"            "CA1_PC4"            "CA1_PC5"           
    ##  [7] "CA1_PC6"            "CA1_PC7"            "CA1_PC8"           
    ## [10] "CA1_PC9"            "CA3_PC1"            "CA3_PC2"           
    ## [13] "CA3_PC3"            "CA3_PC4"            "CA3_PC5"           
    ## [16] "CA3_PC6"            "CA3_PC7"            "CA3_PC8"           
    ## [19] "CA3_PC9"            "DG_PC1"             "DG_PC2"            
    ## [22] "DG_PC3"             "DG_PC4"             "DG_PC5"            
    ## [25] "DG_PC6"             "DG_PC7"             "DG_PC8"            
    ## [28] "DG_PC9"             "ID"                 "RNAseqID"          
    ## [31] "Behavior_PC1"       "Behavior_PC2"       "Behavior_PC3"      
    ## [34] "Behavior_PC4"       "Behavior_PC6"       "APA2"              
    ## [37] "Group"              "Pre_Potentiation"   "Early_Potentiation"
    ## [40] "Late_Potentiation"  "Max_fEPSP"

    alldataslim <- alldataslim[,-c(1,5:10,14:19,23:28,29:30,36:37)]
    alldataslim <- sapply( alldataslim, as.numeric )

    cormat <- rcorr(as.matrix(alldataslim))
    res2 <- rcorr(as.matrix(alldataslim))

    corrplot(res2$r, type = "lower", order = "hclust", 
            tl.col = "black", tl.srt = 45)

![](../figures/04_integration/correlationsAvoidance-1.png)

    corrplot(res2$r, type="lower", order="hclust",  tl.col = "black", 
            p.mat = res2$P, sig.level = 0.05, insig = "blank")

![](../figures/04_integration/correlationsAvoidance-2.png)

    DEGes <- as.data.frame(cormat$r)
    DEGes$group <- c("CA1","CA1","CA1","CA3","CA3","CA3","DG","DG","DG","Behavior","Behavior","Behavior","Behavior","Behavior","Physiology","Physiology","Physiology","Physiology")
    df <- as.data.frame(DEGes$group)
    names(df)[1] <- "Level"
    row.names(df) <- row.names(DEGes)
    DEGes$group <- NULL

    paletteLength <- 30
    myBreaks <- c(seq(min(cormat$r), 0, length.out=ceiling(paletteLength/2) + 1),
                  seq(max(cormat$r)/paletteLength, max(cormat$r), length.out=floor(paletteLength/2)))

    pheatmap(DEGes, show_colnames=T, show_rownames = T,
             annotation_col = df, annotation_colors = ann_colorsLevel,
             annotation_row = df, 
             #annotation_legend = FALSE,
             annotation_names_row = FALSE, annotation_names_col = FALSE,
             treeheight_row = 10, treeheight_col = 10,
             fontsize = 8, 
             width=2.5, height=2,
             border_color = "grey60" ,
             color = viridis(30),
             #cellwidth = 10, 
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation",
             clustering_distance_rows="correlation" ,
             main="Cognitively Trainied"
             )

![](../figures/04_integration/correlationsAvoidance-3.png)

    pheatmap(cormat$r, show_colnames=F, show_rownames = T,
             annotation_col = df, annotation_colors = ann_colorsLevel,
             annotation_row = df, 
             #annotation_legend = FALSE,
             annotation_names_row = FALSE, annotation_names_col = FALSE,
             treeheight_row = 10, treeheight_col = 10,
             fontsize = 8, 
             width=4, height=2.25,
             border_color = "grey60" ,
             color = viridis(30),
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation",
             clustering_distance_rows="correlation",
             filename = "../figures/04_integration/corrTrained.pdf"
             )

    head(alldata)

    ##     Mouse   CA1_PC1  CA1_PC2   CA1_PC3   CA1_PC4    CA1_PC5    CA1_PC6
    ## 1 15-143A        NA       NA        NA        NA         NA         NA
    ## 2 15-143A        NA       NA        NA        NA         NA         NA
    ## 3 15-143B -15.81153 15.25303 -3.606462 0.9520639  3.8894245  0.4329035
    ## 4 15-143B -15.81153 15.25303 -3.606462 0.9520639  3.8894245  0.4329035
    ## 5 15-143C -16.30623 17.24287 -7.792934 2.1072940 -0.1733359 -1.1584867
    ## 6 15-143C -16.30623 17.24287 -7.792934 2.1072940 -0.1733359 -1.1584867
    ##        CA1_PC7   CA1_PC8      CA1_PC9   CA3_PC1   CA3_PC2   CA3_PC3
    ## 1           NA        NA           NA -18.11347 -19.38979 -5.184887
    ## 2           NA        NA           NA -18.11347 -19.38979 -5.184887
    ## 3 -0.795806464 -3.027021 -0.006689372        NA        NA        NA
    ## 4 -0.795806464 -3.027021 -0.006689372        NA        NA        NA
    ## 5  0.008446822 -2.806268 -1.884597363        NA        NA        NA
    ## 6  0.008446822 -2.806268 -1.884597363        NA        NA        NA
    ##     CA3_PC4   CA3_PC5   CA3_PC6    CA3_PC7   CA3_PC8    CA3_PC9   DG_PC1
    ## 1 -7.600432 -6.916928 0.2442081 -0.1748939 0.2335554 -0.9225033 30.62595
    ## 2 -7.600432 -6.916928 0.2442081 -0.1748939 0.2335554 -0.9225033 30.62595
    ## 3        NA        NA        NA         NA        NA         NA 26.91100
    ## 4        NA        NA        NA         NA        NA         NA 26.91100
    ## 5        NA        NA        NA         NA        NA         NA       NA
    ## 6        NA        NA        NA         NA        NA         NA       NA
    ##       DG_PC2    DG_PC3    DG_PC4      DG_PC5    DG_PC6    DG_PC7
    ## 1 -0.9656786 -5.455227  4.087336 -0.08806555 2.7988142 -2.090717
    ## 2 -0.9656786 -5.455227  4.087336 -0.08806555 2.7988142 -2.090717
    ## 3 -1.8943162 -2.192852 -0.767995  4.05246079 0.2182273 -2.291432
    ## 4 -1.8943162 -2.192852 -0.767995  4.05246079 0.2182273 -2.291432
    ## 5         NA        NA        NA          NA        NA        NA
    ## 6         NA        NA        NA          NA        NA        NA
    ##       DG_PC8    DG_PC9     ID   RNAseqID Behavior_PC1 Behavior_PC2
    ## 1 -0.6829819 0.7056409 15143A 143A-CA3-1     8.979881    -3.205040
    ## 2 -0.6829819 0.7056409 15143A  143A-DG-1     8.979881    -3.205040
    ## 3  1.6244657 3.8663534 15143B 143B-CA1-1    -7.986703    -4.867634
    ## 4  1.6244657 3.8663534 15143B  143B-DG-1    -7.986703    -4.867634
    ## 5         NA        NA 15143C 143C-CA1-1    12.252855     3.106723
    ## 6         NA        NA 15143C 143C-CA1-S    12.252855     3.106723
    ##   Behavior_PC3 Behavior_PC4 Behavior_PC6           APA2    Group
    ## 1   2.67546000    0.8503198    -6.280920       conflict     <NA>
    ## 2   2.67546000    0.8503198    -6.280920       conflict     <NA>
    ## 3  -0.03632236    2.3240220     8.613456 yoked-conflict  control
    ## 4  -0.03632236    2.3240220     8.613456 yoked-conflict  control
    ## 5   6.88945052   -2.8341330     2.997247     consistent conflict
    ## 6   6.88945052   -2.8341330     2.997247     consistent conflict
    ##   Pre_Potentiation Early_Potentiation Late_Potentiation  Max_fEPSP
    ## 1               NA                 NA                NA         NA
    ## 2               NA                 NA                NA         NA
    ## 3           101.56              95.34             76.90 -0.0061714
    ## 4           101.56              95.34             76.90 -0.0061714
    ## 5           130.98             199.54            195.64 -0.0077112
    ## 6           130.98             199.54            195.64 -0.0077112

    alldata$APA2 <- factor(alldata$APA2, levels = c("yoked-consistent", "consistent",  "yoked-conflict", "conflict"))

    alldata$avoidance <-  ifelse(grepl("yoked", alldata$APA), "no", "yes")

    alldata$wrap <- "All Cognitively Trained" 

    behavDG <- alldata %>%
      filter(avoidance == "yes") %>%
      ggplot(aes(Behavior_PC1,as.numeric(DG_PC2))) + 
      geom_point(size = 2, alpha = 0.75, aes(color=APA2)) +
      scale_color_manual(values = colorvalAPA6) +
        theme_cowplot(font_size = 8, line_size = 0.25)  +
        theme(legend.position="none") +
        stat_smooth(method = "lm", color="red") + 
        ylab(paste0("DG PC1")) +
        xlab(paste0("Behavior PC1"))  +
      facet_wrap(~wrap)
    behavDG

    pdf(file="../figures/04_integration/behavDG.pdf", width=1.75, height=2)
    plot(behavDG)
    dev.off()

    ## pdf 
    ##   3

    alldata <- left_join(pcadatabyregion, rossetta , by="Mouse")
    alldata <- left_join(alldata, behaviorpca, by="ID")
    alldata <- left_join(alldata, ephys3, by="Mouse")

    alldataslim <- dplyr::filter(alldata, APA2 %in% c("yoked-conflict","yoked-consistent"))

    alldataslim <- alldataslim[,-c(1,5:10,14:19,23:28,29:30,36:37)]
    alldataslim <- sapply( alldataslim, as.numeric )

    cormat <- rcorr(as.matrix(alldataslim))
    res2 <- rcorr(as.matrix(alldataslim))

    corrplot(res2$r, type = "lower", order = "hclust", 
            tl.col = "black", tl.srt = 45)

![](../figures/04_integration/correlationsYoked-1.png)

    corrplot(res2$r, type="lower", order="hclust",  tl.col = "black", 
            p.mat = res2$P, sig.level = 0.05, insig = "blank")

![](../figures/04_integration/correlationsYoked-2.png)

    DEGes <- as.data.frame(cormat$r)
    DEGes$group <- c("CA1","CA1","CA1","CA3","CA3","CA3","DG","DG","DG","Behavior","Behavior","Behavior","Behavior","Behavior","Physiology","Physiology","Physiology","Physiology")
    df <- as.data.frame(DEGes$group)
    names(df)[1] <- "Level"
    row.names(df) <- row.names(DEGes)
    DEGes$group <- NULL

    paletteLength <- 30
    myBreaks <- c(seq(min(cormat$r), 0, length.out=ceiling(paletteLength/2) + 1),
                  seq(max(cormat$r)/paletteLength, max(cormat$r), length.out=floor(paletteLength/2)))
    pheatmap(DEGes, show_colnames=T, show_rownames = T,
             annotation_col = df, annotation_colors = ann_colorsLevel,
             annotation_row = df, 
             #annotation_legend = FALSE,
             annotation_names_row = FALSE, annotation_names_col = FALSE,
             treeheight_row = 10, treeheight_col = 10,
             fontsize = 8, 
             width=2.5, height=2,
             border_color = "grey60" ,
             color = viridis(30),
             #cellwidth = 10, 
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation",
             clustering_distance_rows="correlation" ,
             main="All Yoked"
             )

![](../figures/04_integration/correlationsYoked-3.png)

    pheatmap(cormat$r, show_colnames=F, show_rownames = T,
             annotation_col = df, annotation_colors = ann_colorsLevel,
             annotation_row = df, 
             #annotation_legend = FALSE,
             annotation_names_row = FALSE, annotation_names_col = FALSE,
             treeheight_row = 10, treeheight_col = 10,
             fontsize = 8, 
             width=4, height=2.25,
             border_color = "grey60" ,
             color = viridis(30),
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation",
             clustering_distance_rows="correlation",
             filename = "../figures/04_integration/corrYoked.pdf",
             main="All Yoked"
             )

Significant genes

    # read files with pvalues
    CA1trainpvals <- read.csv("../data/02c_CA1_consyokcons.csv", stringsAsFactors = F, row.names = NULL) 
    CA1stresspvals <- read.csv("../data/02c_CA1_yokeconfyokcons.csv", stringsAsFactors = F, row.names = NULL)
    DGtrainpvals <- read.csv("../data/02c_DG_consyokcons.csv", stringsAsFactors = F)

    # remove non-significant genes and genes affected by learning and stresss
    CA1learnstress <- read.csv("../data/02_CA1learningstressgenes.csv", stringsAsFactors = F)
    CA1learnstress <- as.vector(CA1learnstress$gene)

    CA1trainpvals <- CA1trainpvals %>% filter(!gene %in% CA1learnstress,
                                              direction != "NS")
    CA1stresspvals <- CA1stresspvals %>% filter(!gene %in% CA1learnstress,
                                              direction != "NS")
    DGtrainpvals <- DGtrainpvals %>% filter(!gene %in% CA1learnstress,
                                              direction != "NS")

    # read variance stabilized count data
    DGvsd <- read.csv("../data/02c_DGvsd.csv", stringsAsFactors = F, check.names = F) 
    CA1vsd <- read.csv("../data/02c_CA1vsd.csv", stringsAsFactors = F, check.names = F)

    # function to calculate and plot correlations between gene expression and behaivor
    correlationdf <- function(pvalues, vsds){
      
      # add gene to df for joining
      colnames(vsds)[1] <- "gene"
      
      mydf <- left_join(pvalues, vsds)
      row.names(mydf) <- mydf$gene
      mydf <- mydf[-c(1:6)]
      mydf <- as.data.frame(t(mydf))
      
      # create an ID column to join lists of DEGs with their experssion values and the behavior of the animal 
      mydf$sample <-  row.names(mydf)
      mydf$mouse <- sapply(strsplit(as.character(mydf$sample),'-'), "[", 1)
      mydf$ID <- paste("15", mydf$mouse, sep = "")
      mydf$sample <-  NULL
      mydf$mouse <- NULL
      row.names(mydf) <- mydf$ID
      
      # randomely select 15 genes and add ID column back for joining
      mydfrandom <- mydf
      #mydfrandom <- mydfrandom[sample(1:ncol(mydfrandom), 15,replace=FALSE)]
      
      mydfrandom$ID <- row.names(mydfrandom)
      
      # subset the beahvior data to the last time point and only values shows in figure 2
      behavior <- read.csv("../data/01a_behavior.csv", header = T)
      favbehav <- behavior %>% 
        filter(TrainSession == "Retention") %>%
        select(ID, NumEntrances, Time1stEntr)
        
      favbehavgenes <- left_join(mydfrandom, favbehav)
      favbehavgenes$ID <- NULL
      
      mycordf = cor(favbehavgenes)
      #corrplot(mycordf, method = "circle", title = mytitle, tl.col = "black", tl.cex = 0.8)
      mycordf <- as.data.frame(mycordf)
      mycordf$gene <- row.names(mycordf)
      return(mycordf)
    }  

    corCA1stress <- correlationdf(CA1stresspvals, CA1vsd)

    ## Joining, by = "gene"

    ## Joining, by = "ID"

    ## Warning: Column `ID` joining character vector and factor, coercing into
    ## character vector

    corCA1train <- correlationdf(CA1trainpvals, CA1vsd)

    ## Joining, by = "gene"
    ## Joining, by = "ID"

    ## Warning: Column `ID` joining character vector and factor, coercing into
    ## character vector

    corDGtrain <- correlationdf(DGtrainpvals, DGvsd) 

    ## Joining, by = "gene"
    ## Joining, by = "ID"

    ## Warning: Column `ID` joining character vector and factor, coercing into
    ## character vector

    corpointplot <- function(mycordf, mytitle){

      # remove self corelations
      mycordf <- mycordf %>%
        filter(!gene %in% c("NumEntrances", "Time1stEntr"))
      
      p1 <- mycordf %>%
        select(-Time1stEntr) %>% # remove corelations with other behavior
        ggplot(aes(x = reorder(gene, -NumEntrances), y = NumEntrances, color = NumEntrances)) +
        geom_point(stat = "identity") +
        theme_minimal(base_size = 7) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              legend.position = "none")   +
        labs(subtitle = mytitle, x = "signficant genes", y = "Correlation to\nNumber of Entrances")  +
        ylim(-1, 1) +
        scale_color_viridis(limits=c(-1,1))
      
      p2 <- mycordf %>%
        select(-NumEntrances) %>% # remove corelations with other behavior
        ggplot(aes(x =  reorder(gene, Time1stEntr), y = Time1stEntr, color = Time1stEntr)) +
        geom_point(stat = "identity") +
        theme_minimal(base_size = 7) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              legend.position = "none",
              legend.title = element_blank())   +
        labs(subtitle = mytitle, x = "signficant genes", y = "Correlation to\nTime to 1st Entrances")  +
        ylim(-1, 1) +
        scale_color_viridis(limits=c(-1,1)) 
        
      
      p12 <- plot_grid(p2,p1, nrow = 2)
        
      print(p12)
    }

    a <- corpointplot(corCA1train, "CA1 training-responsive genes") 

![](../figures/04_integration/significantgenes-1.png)

    b <- corpointplot(corDGtrain, "DG training-responsive genes")

![](../figures/04_integration/significantgenes-2.png)

    ab <- plot_grid(b,a, nrow = 1)
    ab 

![](../figures/04_integration/significantgenes-3.png)

    pdf(file="../figures/04_integration/correlationpointplots.pdf", width=5, height=4)
    plot(ab)    
    dev.off() 

    ## quartz_off_screen 
    ##                 2

    corDGtrain %>%
      filter(Time1stEntr > 0.75 | NumEntrances < -.75 |  Time1stEntr < -.75 | NumEntrances > 0.75,
             !gene %in% c("Time1stEntr", "NumEntrances")) %>%
      select(gene, Time1stEntr, NumEntrances)

    ##       gene Time1stEntr NumEntrances
    ## 1     Acan   0.7270032   -0.7757527
    ## 2   Amigo2   0.7301423   -0.7621433
    ## 3      Arc   0.7550063   -0.7613195
    ## 4   Armcx5   0.7372492   -0.8373150
    ## 5     Fzd5   0.7954120   -0.6875919
    ## 6    Npas4   0.7686581   -0.7123179
    ## 7    Ptgs2   0.7818914   -0.7566271
    ## 8     Rgs2   0.7664609   -0.7638017
    ## 9  Slc16a1   0.7613746   -0.6036242
    ## 10    Syt4   0.6248840   -0.7540978

    corCA1train %>%
      filter(Time1stEntr > 0.5 | NumEntrances < -.5 |  Time1stEntr < -.5 | NumEntrances > 0.5,
             !gene %in% c("Time1stEntr", "NumEntrances")) %>%
      select(gene, Time1stEntr, NumEntrances)

    ##        gene Time1stEntr NumEntrances
    ## 1     Acad8   0.2212761   -0.5189565
    ## 2     Ahdc1  -0.6498533    0.4807087
    ## 3     Alms1   0.3317366   -0.5903613
    ## 4      Apln   0.6912757   -0.6768894
    ## 5      Aptx   0.5011229   -0.5020300
    ## 6    Arrdc3   0.3842883   -0.6024845
    ## 7    Atpaf2   0.6336303   -0.6909063
    ## 8   Bloc1s5   0.6123718   -0.6593602
    ## 9    Capn10   0.3990502   -0.6733059
    ## 10  Ccdc190   0.5263474   -0.5887538
    ## 11     Cdh7   0.4069231   -0.5462141
    ## 12   Cfap20  -0.5295420    0.3721191
    ## 13     Chuk  -0.5261919    0.5031008
    ## 14 Colgalt2  -0.3676015    0.5855128
    ## 15    Ctcfl  -0.6362308    0.4263709
    ## 16      Dbt   0.3746872   -0.5351618
    ## 17     Dner   0.7111875   -0.5810722
    ## 18    Esyt1   0.4632121   -0.5166430
    ## 19  Fam110b   0.5085691   -0.3963534
    ## 20  Fam19a2   0.6055048   -0.6704304
    ## 21   Fam63b   0.4389843   -0.7527619
    ## 22    Fanci   0.2495205   -0.5278318
    ## 23     Fbn1   0.6268982   -0.5418373
    ## 24     Fgd5   0.3047857   -0.5324103
    ## 25   Filip1   0.6122806   -0.5460707
    ## 26     Flt1   0.3163125   -0.5325390
    ## 27   Glcci1   0.5773992   -0.6595865
    ## 28  Gm21887  -0.3775421    0.5255329
    ## 29   Gm4631   0.5543307   -0.6905110
    ## 30   Gm9821  -0.5657761    0.6178466
    ## 31    Grem2   0.6740555   -0.6294663
    ## 32    Hecw2   0.5761573   -0.6882464
    ## 33    Igbp1   0.6660945   -0.5249595
    ## 34  Igf2bp2  -0.7300574    0.5766138
    ## 35    Inhbb  -0.4958727    0.5084745
    ## 36   Inpp4b   0.3541333   -0.5070030
    ## 37    Kcna4   0.5178216   -0.5925663
    ## 38    Kcnd3   0.5167999   -0.4040270
    ## 39    Kdm5d  -0.5099970    0.5253640
    ## 40    Klkb1  -0.5990080    0.4479896
    ## 41     Kri1  -0.2566445    0.5235254
    ## 42 Mapkapk2   0.4868547   -0.7024112
    ## 43    Marf1  -0.6246670    0.3996785
    ## 44      Mcc   0.4840093   -0.5455392
    ## 45     Med8  -0.3243889    0.5654744
    ## 46     Nat9   0.3625147   -0.6160515
    ## 47    Neto2   0.4211430   -0.6731647
    ## 48    Ntng1   0.2848668   -0.5579322
    ## 49    Ovca2  -0.5373988    0.5294857
    ## 50    Pde6a  -0.7391394    0.5982401
    ## 51    Pias4  -0.5077882    0.5889152
    ## 52  Ppp1r3f   0.6201177   -0.5472079
    ## 53    Ptprm   0.3825078   -0.5491678
    ## 54    Rcan2   0.4660109   -0.5137531
    ## 55     Rdh1   0.3588774   -0.5204037
    ## 56    Rfesd   0.5333210   -0.5376000
    ## 57    Rfwd3   0.4724059   -0.6524796
    ## 58     Rgs2   0.4751647   -0.5122033
    ## 59   Rsc1a1  -0.5193059    0.6224766
    ## 60    Shoc2   0.2692746   -0.5699276
    ## 61 Slc25a38  -0.1698548    0.5691872
    ## 62  Slc2a13   0.3511989   -0.5725884
    ## 63   Slc4a3   0.5608161   -0.7344899
    ## 64    Slx1b  -0.5457649    0.3301103
    ## 65   Snap29  -0.3559440    0.5278961
    ## 66    Snx24   0.6703670   -0.5677200
    ## 67    Sox10   0.5813646   -0.6210244
    ## 68     Sox5   0.2115937   -0.5163476
    ## 69  St8sia4   0.1249904   -0.5032398
    ## 70  Stard13   0.2474957   -0.5596036
    ## 71     Stk3   0.5448689   -0.4969145
    ## 72    Stox2  -0.6174573    0.4853552
    ## 73  Tbc1d30   0.6462628   -0.6726156
    ## 74   Tbc1d9   0.5358903   -0.4989014
    ## 75     Tet3   0.2322443   -0.5865108
    ## 76     Tns2   0.6988295   -0.4519199
    ## 77      Ttn   0.4298739   -0.5784336
    ## 78   Tvp23a   0.4377234   -0.5311110
    ## 79    Uckl1  -0.3703443    0.7326685
    ## 80   Zc3h13   0.4759109   -0.5401196
    ## 81   Zfp446   0.4422391   -0.5452128
    ## 82   Zfp617  -0.2295565    0.5454083
    ## 83   Zfp738  -0.5145745    0.3693506
    ## 84   Zfp831   0.3255999   -0.5092313
    ## 85    Znrd1   0.3535553   -0.6145544

    CA1learn <- read.csv( "../data/02_CA1learninggenes.csv", stringsAsFactors = F)
    CA1learnGO <- left_join(CA1trainpvals,CA1learn)

    ## Joining, by = "gene"

    CA1learnGO

    ##              gene        padj  logpadj       lfc         direction
    ## 1   1810022K09Rik 0.088227456 1.054396 -2.103747 yoked\nconsistent
    ## 2   9430015G10Rik 0.054803146 1.261195 -1.942464 yoked\nconsistent
    ## 3           Aagab 0.087372943 1.058623 -1.129414 yoked\nconsistent
    ## 4            Abl1 0.047560292 1.322755  1.540018        consistent
    ## 5           Acad8 0.099236138 1.003330  1.777297        consistent
    ## 6           Acads 0.019790733 1.703538 -4.336036 yoked\nconsistent
    ## 7           Acot9 0.096286699 1.016434 -1.782720 yoked\nconsistent
    ## 8           Acta1 0.061321606 1.212386  4.198198        consistent
    ## 9           Adam8 0.075039640 1.124709 -2.813243 yoked\nconsistent
    ## 10          Adat1 0.079398573 1.100187  4.687065        consistent
    ## 11         Afg3l2 0.082914792 1.081368  1.156923        consistent
    ## 12          Ahdc1 0.095844557 1.018433 -1.890135 yoked\nconsistent
    ## 13            Aip 0.082783987 1.082054 -1.115433 yoked\nconsistent
    ## 14          Akap9 0.035396074 1.451045  1.252930        consistent
    ## 15          Alg11 0.060863185 1.215645 -1.162234 yoked\nconsistent
    ## 16          Alms1 0.050806429 1.294081  5.364739        consistent
    ## 17         Amigo1 0.077753552 1.109280  1.670131        consistent
    ## 18         Amigo3 0.082305651 1.084570 -2.444226 yoked\nconsistent
    ## 19         Ankfy1 0.093565087 1.028886  1.127166        consistent
    ## 20        Ankrd12 0.087372943 1.058623  1.107307        consistent
    ## 21           Apln 0.059369676 1.226435  5.070182        consistent
    ## 22          Apol8 0.081959305 1.086402 -1.928941 yoked\nconsistent
    ## 23           Aptx 0.089956243 1.045969  1.878595        consistent
    ## 24        Arfgef3 0.089917206 1.046157  1.509322        consistent
    ## 25         Arl13b 0.097862636 1.009383  3.405766        consistent
    ## 26          Arl5a 0.067391414 1.171395  1.144549        consistent
    ## 27         Arntl2 0.090232749 1.044636 -2.648957 yoked\nconsistent
    ## 28         Arrdc3 0.090232749 1.044636  2.406740        consistent
    ## 29          Asxl3 0.052873146 1.276765  3.865796        consistent
    ## 30           Atf5 0.096482969 1.015549 -1.642500 yoked\nconsistent
    ## 31          Atf6b 0.040411491 1.393495 -1.225260 yoked\nconsistent
    ## 32          Atp5e 0.066358213 1.178105 -1.489385 yoked\nconsistent
    ## 33         Atpaf2 0.089917206 1.046157  3.233725        consistent
    ## 34          Baz2a 0.096795675 1.014144  1.089441        consistent
    ## 35          Baz2b 0.090852311 1.041664  2.163200        consistent
    ## 36       BC030500 0.090298624 1.044319 -1.187449 yoked\nconsistent
    ## 37        Bloc1s5 0.082783987 1.082054  4.632711        consistent
    ## 38           Bop1 0.077048887 1.113234 -1.234956 yoked\nconsistent
    ## 39         Btbd10 0.088000145 1.055517  1.252840        consistent
    ## 40          C2cd3 0.067391414 1.171395  1.508586        consistent
    ## 41         Capn10 0.047560292 1.322755  1.992656        consistent
    ## 42        Ccdc190 0.081574331 1.088446  4.719891        consistent
    ## 43         Ccdc43 0.089832225 1.046568 -1.312444 yoked\nconsistent
    ## 44           Ccr5 0.087372943 1.058623 -2.856566 yoked\nconsistent
    ## 45          Cdh20 0.051433838 1.288751  5.358163        consistent
    ## 46           Cdh7 0.038825245 1.410886  5.755961        consistent
    ## 47          Cdk15 0.088083664 1.055105  5.403692        consistent
    ## 48         Cep131 0.067391414 1.171395 -1.608290 yoked\nconsistent
    ## 49         Cfap20 0.075888087 1.119826 -1.244613 yoked\nconsistent
    ## 50           Chuk 0.077048887 1.113234 -1.161201 yoked\nconsistent
    ## 51          Clock 0.075917692 1.119657  1.301507        consistent
    ## 52          Cntrl 0.088083664 1.055105  2.152049        consistent
    ## 53       Colgalt2 0.035396074 1.451045 -2.859570 yoked\nconsistent
    ## 54            Cpq 0.066691534 1.175929 -3.218289 yoked\nconsistent
    ## 55           Crem 0.087372943 1.058623 -2.143971 yoked\nconsistent
    ## 56          Crim1 0.090848135 1.041684  1.086726        consistent
    ## 57          Cstf1 0.081560857 1.088518  1.858374        consistent
    ## 58          Ctcfl 0.088083664 1.055105 -2.887153 yoked\nconsistent
    ## 59     D1Ertd622e 0.050619598 1.295681  1.649572        consistent
    ## 60            Dbt 0.074745860 1.126413  2.516856        consistent
    ## 61         Ddx19b 0.071099580 1.148133  1.678328        consistent
    ## 62          Ddx31 0.075888087 1.119826  4.790523        consistent
    ## 63        Dennd1b 0.095844557 1.018433  2.583597        consistent
    ## 64           Dner 0.019268528 1.715151  1.117948        consistent
    ## 65           Dok6 0.089956243 1.045969  2.097135        consistent
    ## 66           Dpf1 0.050806429 1.294081 -1.152038 yoked\nconsistent
    ## 67         Dusp19 0.094831372 1.023048 -2.064674 yoked\nconsistent
    ## 68         Eif2s1 0.086294718 1.064016 -1.575040 yoked\nconsistent
    ## 69           Emc4 0.023924907 1.621150 -1.419306 yoked\nconsistent
    ## 70          Enkd1 0.045945672 1.337755 -2.992738 yoked\nconsistent
    ## 71          Esrrg 0.089613002 1.047629  2.691423        consistent
    ## 72          Esyt1 0.093565087 1.028886  4.963608        consistent
    ## 73           Ext1 0.079834230 1.097811  1.132101        consistent
    ## 74        Fam110b 0.066775644 1.175382  2.004312        consistent
    ## 75        Fam19a2 0.038825245 1.410886  5.671638        consistent
    ## 76         Fam63a 0.098399309 1.007008  2.575974        consistent
    ## 77         Fam63b 0.091274774 1.039649  1.116486        consistent
    ## 78          Fanci 0.056885301 1.245000  5.231147        consistent
    ## 79           Fbn1 0.071124207 1.147983  2.457364        consistent
    ## 80          Fbxl6 0.099772651 1.000988 -1.239693 yoked\nconsistent
    ## 81           Fgd4 0.061475547 1.211298  1.980883        consistent
    ## 82           Fgd5 0.072409094 1.140207  4.953602        consistent
    ## 83         Filip1 0.064290908 1.191850  4.947799        consistent
    ## 84           Flt1 0.048061879 1.318199  1.570229        consistent
    ## 85           Fmn1 0.053664604 1.270312  2.269599        consistent
    ## 86         Fn3krp 0.099236138 1.003330 -1.919983 yoked\nconsistent
    ## 87         Fndc3a 0.043730707 1.359214  1.499059        consistent
    ## 88          Foxj3 0.047560292 1.322755  1.364845        consistent
    ## 89         Frmpd3 0.088887213 1.051161  1.638085        consistent
    ## 90           Fzd3 0.085048392 1.070334  1.928208        consistent
    ## 91           Get4 0.060926364 1.215195 -1.333426 yoked\nconsistent
    ## 92         Glcci1 0.081560857 1.088518  2.032681        consistent
    ## 93        Gm10053 0.036827325 1.433830 -1.177339 yoked\nconsistent
    ## 94        Gm13889 0.061393493 1.211878  2.432572        consistent
    ## 95        Gm21887 0.043955846 1.356983 -3.630310 yoked\nconsistent
    ## 96         Gm4631 0.030758294 1.512038  1.070944        consistent
    ## 97          Gm527 0.082305651 1.084570 -3.498094 yoked\nconsistent
    ## 98         Gm9821 0.035396074 1.451045 -2.229978 yoked\nconsistent
    ## 99         Gpr161 0.078677396 1.104150  1.069548        consistent
    ## 100        Gpr180 0.051433838 1.288751  2.761881        consistent
    ## 101         Grem2 0.089308034 1.049109  4.761767        consistent
    ## 102         Grik3 0.045066602 1.346145  3.337036        consistent
    ## 103          Grm1 0.052012170 1.283895  1.566786        consistent
    ## 104         Gstt3 0.092552012 1.033614  4.965727        consistent
    ## 105       Gucy1a2 0.075888087 1.119826  2.408653        consistent
    ## 106        Gucy2e 0.089917206 1.046157 -2.966315 yoked\nconsistent
    ## 107          Guk1 0.087297939 1.058996 -1.300820 yoked\nconsistent
    ## 108         Hdac6 0.082783987 1.082054 -1.160096 yoked\nconsistent
    ## 109         Hecw2 0.060863185 1.215645  1.950303        consistent
    ## 110          Helz 0.055649748 1.254537  1.366751        consistent
    ## 111       Herpud2 0.077221657 1.112261  1.230543        consistent
    ## 112          Hexb 0.070044490 1.154626 -1.153432 yoked\nconsistent
    ## 113         Hmox2 0.093506776 1.029157 -1.063294 yoked\nconsistent
    ## 114        Hspa1b 0.035677918 1.447600  1.568311        consistent
    ## 115         Hspa2 0.096927716 1.013552  1.099304        consistent
    ## 116         Igbp1 0.066775644 1.175382  4.435456        consistent
    ## 117       Igf2bp2 0.064309040 1.191728 -4.033622 yoked\nconsistent
    ## 118         Inhbb 0.003540332 2.450956 -3.822015 yoked\nconsistent
    ## 119         Ino80 0.074841461 1.125858 -1.368529 yoked\nconsistent
    ## 120        Inpp4b 0.098321929 1.007350  4.672146        consistent
    ## 121          Irs1 0.089917206 1.046157  2.197218        consistent
    ## 122          Jag2 0.081599684 1.088312  1.445742        consistent
    ## 123          Jak2 0.077800374 1.109018  1.456565        consistent
    ## 124         Josd1 0.077753552 1.109280 -1.290114 yoked\nconsistent
    ## 125         Kat6a 0.053844145 1.268862  1.408021        consistent
    ## 126        Kbtbd3 0.085703922 1.066999  5.112558        consistent
    ## 127         Kcna1 0.075090449 1.124415  1.700055        consistent
    ## 128         Kcna4 0.035677918 1.447600  2.117451        consistent
    ## 129         Kcnc1 0.091949134 1.036452  1.149074        consistent
    ## 130         Kcnc2 0.075888087 1.119826  1.094556        consistent
    ## 131         Kcnd3 0.074745860 1.126413  3.723883        consistent
    ## 132         Kctd5 0.057439505 1.240789 -1.555371 yoked\nconsistent
    ## 133         Kdm5d 0.022055009 1.656493 -1.586114 yoked\nconsistent
    ## 134         Khnyn 0.026827664 1.571417 -2.695515 yoked\nconsistent
    ## 135        Kif26b 0.082305651 1.084570  4.948583        consistent
    ## 136         Klkb1 0.018683187 1.728549 -3.519986 yoked\nconsistent
    ## 137         Kmt5a 0.032866115 1.483252 -1.477646 yoked\nconsistent
    ## 138          Kri1 0.009269714 2.032934 -1.700187 yoked\nconsistent
    ## 139         Lemd3 0.050806429 1.294081  1.598017        consistent
    ## 140         Letm2 0.072367199 1.140458  3.355435        consistent
    ## 141         Lnpep 0.066691534 1.175929  1.773609        consistent
    ## 142         Lrfn4 0.081574331 1.088446 -1.243356 yoked\nconsistent
    ## 143        Lrrc58 0.088083664 1.055105  1.115617        consistent
    ## 144          Ltv1 0.065813600 1.181684 -1.828609 yoked\nconsistent
    ## 145        Man1a2 0.095844557 1.018433  1.049128        consistent
    ## 146         Manba 0.099236138 1.003330  4.491903        consistent
    ## 147      Mapkapk2 0.098399309 1.007008  1.445756        consistent
    ## 148         Marf1 0.045945672 1.337755 -1.401831 yoked\nconsistent
    ## 149          Mbd5 0.090852311 1.041664  1.067662        consistent
    ## 150           Mcc 0.047560292 1.322755  3.327820        consistent
    ## 151          Med8 0.062421003 1.204669 -1.461357 yoked\nconsistent
    ## 152       Mettl16 0.056885301 1.245000  2.210106        consistent
    ## 153      Mettl21e 0.043911797 1.357419 -4.031179 yoked\nconsistent
    ## 154       Mfsd13a 0.071099580 1.148133 -1.286223 yoked\nconsistent
    ## 155        Mmadhc 0.082914792 1.081368 -1.410850 yoked\nconsistent
    ## 156      Mphosph9 0.066039594 1.180196 -1.177832 yoked\nconsistent
    ## 157        Mrgpre 0.089917206 1.046157  4.799626        consistent
    ## 158        Mrpl28 0.071033673 1.148536 -1.277171 yoked\nconsistent
    ## 159        Mrpl48 0.061475547 1.211298 -1.650438 yoked\nconsistent
    ## 160        Msl3l2 0.061475547 1.211298 -1.737967 yoked\nconsistent
    ## 161         Mtfmt 0.064195393 1.192496 -1.990499 yoked\nconsistent
    ## 162          Nat9 0.075888087 1.119826  4.830314        consistent
    ## 163        Ndufs7 0.083011642 1.080861 -1.267415 yoked\nconsistent
    ## 164          Nefm 0.061475547 1.211298  1.074088        consistent
    ## 165         Neto2 0.057596674 1.239603  1.257996        consistent
    ## 166          Nkrf 0.099688753 1.001354  1.756114        consistent
    ## 167          Nle1 0.090852311 1.041664  2.900664        consistent
    ## 168          Nme1 0.077800374 1.109018 -1.386883 yoked\nconsistent
    ## 169           Nov 0.070911743 1.149282 -1.016098 yoked\nconsistent
    ## 170         Npas4 0.072260561 1.141099  2.342853        consistent
    ## 171         Ntng1 0.097012447 1.013173  4.843551        consistent
    ## 172         Ntpcr 0.030758294 1.512038 -2.103366 yoked\nconsistent
    ## 173        Nudt19 0.084452734 1.073386 -1.203811 yoked\nconsistent
    ## 174         Nudt6 0.070354119 1.152710 -4.264652 yoked\nconsistent
    ## 175          Nxt2 0.036827325 1.433830 -1.910268 yoked\nconsistent
    ## 176         Olfm3 0.044213137 1.354449  1.974534        consistent
    ## 177         Ovca2 0.025320503 1.596528 -1.888529 yoked\nconsistent
    ## 178         Patz1 0.081959305 1.086402  1.398037        consistent
    ## 179        Pcdh17 0.049991290 1.301106  1.568507        consistent
    ## 180         Pde6a 0.049487332 1.305506 -3.718742 yoked\nconsistent
    ## 181         Phka1 0.077048887 1.113234  2.975156        consistent
    ## 182         Pias4 0.088154509 1.054755 -1.392803 yoked\nconsistent
    ## 183        Plagl2 0.079596747 1.099105  4.392761        consistent
    ## 184          Plau 0.066755149 1.175515 -2.132120 yoked\nconsistent
    ## 185        Polr2h 0.087372943 1.058623 -1.871598 yoked\nconsistent
    ## 186        Polrmt 0.088083664 1.055105 -1.060275 yoked\nconsistent
    ## 187       Ppp1r3f 0.061475547 1.211298  1.776499        consistent
    ## 188         Psg28 0.081959305 1.086402 -4.158935 yoked\nconsistent
    ## 189         Psmg2 0.057596674 1.239603 -2.399810 yoked\nconsistent
    ## 190         Ptpn4 0.092944106 1.031778  1.284176        consistent
    ## 191         Ptprm 0.086294718 1.064016  2.701829        consistent
    ## 192         Pygo1 0.082783987 1.082054  1.016283        consistent
    ## 193          Rbak 0.057439505 1.240789  4.453020        consistent
    ## 194         Rbbp4 0.027845581 1.555244 -1.325689 yoked\nconsistent
    ## 195         Rcan2 0.079752393 1.098256  1.334951        consistent
    ## 196          Rdh1 0.054803146 1.261195  4.353353        consistent
    ## 197          Rest 0.086294718 1.064016  4.923687        consistent
    ## 198         Rfesd 0.092995319 1.031539  4.377184        consistent
    ## 199         Rfwd3 0.055834943 1.253094  2.668277        consistent
    ## 200          Rgs2 0.050618806 1.295688  1.741528        consistent
    ## 201        Rhbdl3 0.079714765 1.098461 -1.382995 yoked\nconsistent
    ## 202       Rnaseh1 0.025320503 1.596528 -3.119018 yoked\nconsistent
    ## 203        Rnf165 0.079408520 1.100133  1.608227        consistent
    ## 204        Rnf180 0.098533612 1.006416  4.646787        consistent
    ## 205        Rnf216 0.098399309 1.007008 -1.073959 yoked\nconsistent
    ## 206         Rnf25 0.028332466 1.547716 -1.872594 yoked\nconsistent
    ## 207         Rpl10 0.025963596 1.585635 -2.250149 yoked\nconsistent
    ## 208         Rpl36 0.050806429 1.294081 -1.359802 yoked\nconsistent
    ## 209        Rsc1a1 0.049991290 1.301106 -6.696627 yoked\nconsistent
    ## 210        Rsph3a 0.075888087 1.119826 -1.836296 yoked\nconsistent
    ## 211        S100a1 0.090506687 1.043319  4.176658        consistent
    ## 212        Samhd1 0.052247414 1.281935  2.555701        consistent
    ## 213        Sema4a 0.055697684 1.254163  2.290312        consistent
    ## 214         Sgms2 0.090232749 1.044636  5.086672        consistent
    ## 215         Shoc2 0.040649863 1.390941  1.326581        consistent
    ## 216       Slc24a2 0.040888161 1.388402  1.140005        consistent
    ## 217      Slc25a38 0.079049015 1.102104 -1.674412 yoked\nconsistent
    ## 218      Slc25a46 0.035396074 1.451045  1.158403        consistent
    ## 219       Slc26a2 0.081959305 1.086402  3.718625        consistent
    ## 220       Slc2a13 0.095844557 1.018433  1.118219        consistent
    ## 221       Slc35b4 0.049487332 1.305506 -1.323510 yoked\nconsistent
    ## 222       Slc44a1 0.053664604 1.270312  1.364285        consistent
    ## 223        Slc4a3 0.047296750 1.325169  1.015137        consistent
    ## 224        Slc8a1 0.077847562 1.108755  1.535193        consistent
    ## 225         Slx1b 0.064297173 1.191808 -2.154254 yoked\nconsistent
    ## 226         Smek2 0.094182143 1.026031  1.475515        consistent
    ## 227         Smim3 0.073184046 1.135584 -2.173538 yoked\nconsistent
    ## 228        Snap29 0.032289995 1.490932 -1.667227 yoked\nconsistent
    ## 229          Snx2 0.071033673 1.148536 -1.157842 yoked\nconsistent
    ## 230         Snx24 0.089917206 1.046157  4.639860        consistent
    ## 231        Sorcs1 0.087372943 1.058623 -2.851229 yoked\nconsistent
    ## 232         Sox10 0.056487180 1.248050  3.604220        consistent
    ## 233          Sox5 0.062452817 1.204448  4.161665        consistent
    ## 234         Spast 0.090232749 1.044636  1.120841        consistent
    ## 235       Specc1l 0.070201017 1.153657  1.353570        consistent
    ## 236          Srrd 0.047560292 1.322755 -2.740428 yoked\nconsistent
    ## 237    St6galnac4 0.095844557 1.018433  2.278950        consistent
    ## 238       St8sia4 0.070004366 1.154875  4.554064        consistent
    ## 239       Stard13 0.061475547 1.211298  4.597528        consistent
    ## 240          Stk3 0.097871277 1.009345  2.907311        consistent
    ## 241         Stox2 0.069073318 1.160690 -4.228689 yoked\nconsistent
    ## 242          Stx3 0.063703756 1.195835 -1.331315 yoked\nconsistent
    ## 243       Tbc1d30 0.030758294 1.512038  1.604249        consistent
    ## 244        Tbc1d9 0.053664604 1.270312  1.589010        consistent
    ## 245         Tdrd7 0.088000145 1.055517  1.703265        consistent
    ## 246          Tet3 0.053844145 1.268862  1.094635        consistent
    ## 247        Tmem57 0.053844145 1.268862  1.390690        consistent
    ## 248        Tmem65 0.092944106 1.031778 -1.024638 yoked\nconsistent
    ## 249       Tmem88b 0.058480905 1.232986  2.631945        consistent
    ## 250          Tns2 0.077753552 1.109280  5.250512        consistent
    ## 251          Tox2 0.071099580 1.148133  2.103597        consistent
    ## 252         Traf6 0.070622910 1.151054  2.544356        consistent
    ## 253           Ttn 0.094831372 1.023048  5.042780        consistent
    ## 254        Tvp23a 0.053588804 1.270926  1.914888        consistent
    ## 255         Uckl1 0.057965671 1.236829 -1.385668 yoked\nconsistent
    ## 256         Uqcrh 0.023924907 1.621150 -1.273929 yoked\nconsistent
    ## 257         Usmg5 0.023924907 1.621150 -1.350159 yoked\nconsistent
    ## 258         Usp45 0.035677918 1.447600  1.581830        consistent
    ## 259        Usp6nl 0.065511720 1.183681  4.295653        consistent
    ## 260         Xrcc3 0.092944106 1.031778 -2.619414 yoked\nconsistent
    ## 261         Xrcc6 0.057140335 1.243057 -1.818006 yoked\nconsistent
    ## 262         Ylpm1 0.095844557 1.018433  1.188359        consistent
    ## 263           Zak 0.062452817 1.204448  2.743572        consistent
    ## 264       Zc3h12a 0.073328276 1.134729 -2.583596 yoked\nconsistent
    ## 265        Zc3h13 0.035396074 1.451045  1.362548        consistent
    ## 266        Zfp114 0.089451104 1.048414  4.685814        consistent
    ## 267        Zfp395 0.040411491 1.393495 -3.819722 yoked\nconsistent
    ## 268        Zfp414 0.082795261 1.081995 -1.645229 yoked\nconsistent
    ## 269        Zfp446 0.097473286 1.011114  3.317079        consistent
    ## 270        Zfp580 0.071099580 1.148133 -1.492773 yoked\nconsistent
    ## 271        Zfp617 0.061475547 1.211298 -1.440423 yoked\nconsistent
    ## 272        Zfp711 0.056885301 1.245000 -1.956344 yoked\nconsistent
    ## 273        Zfp738 0.071124207 1.147983 -1.918231 yoked\nconsistent
    ## 274        Zfp821 0.056086978 1.251138 -1.211277 yoked\nconsistent
    ## 275        Zfp831 0.025320503 1.596528  1.931827        consistent
    ## 276       Zfyve16 0.089699290 1.047211  2.680013        consistent
    ## 277         Zmym4 0.078256308 1.106481  1.702054        consistent
    ## 278         Znrd1 0.089832225 1.046568  5.159798        consistent

    CA1stress <- read.csv( "../data/02_CA1lstressgenes.csv")
    CA1stressGO <- left_join(CA1stresspvals,CA1stress)

    ## Joining, by = "gene"

    ## Warning: Column `gene` joining character vector and factor, coercing into
    ## character vector

    CA1stressGO

    ##              gene        padj  logpadj       lfc         direction
    ## 1   1110032F04Rik 0.075375751 1.122768  3.700008 yoked\nconsistent
    ## 2   1600002K03Rik 0.057058567 1.243679 -4.505653   yoked\nconflict
    ## 3   2010107G23Rik 0.017584207 1.754877 -2.138287   yoked\nconflict
    ## 4   2210013O21Rik 0.051653076 1.286904 -2.380294   yoked\nconflict
    ## 5   2310009B15Rik 0.076150968 1.118325  4.971213 yoked\nconsistent
    ## 6   5730409E04Rik 0.038796683 1.411205 -1.078415   yoked\nconflict
    ## 7           Abcc4 0.078624137 1.104444  4.083222 yoked\nconsistent
    ## 8           Abhd4 0.094828933 1.023059  1.645698 yoked\nconsistent
    ## 9            Ache 0.099660863 1.001475  1.452986 yoked\nconsistent
    ## 10          Adcy5 0.047614237 1.322263  1.029938 yoked\nconsistent
    ## 11          Adcy6 0.088617852 1.052479  1.772035 yoked\nconsistent
    ## 12         Adgrb3 0.068496484 1.164332  1.071288 yoked\nconsistent
    ## 13          Adpgk 0.069303994 1.159242 -1.721475   yoked\nconflict
    ## 14         Afg3l1 0.079785818 1.098074 -1.440899   yoked\nconflict
    ## 15           Alg9 0.048923039 1.310487  2.236373 yoked\nconsistent
    ## 16         Apcdd1 0.057174055 1.242801  5.127605 yoked\nconsistent
    ## 17        Arfgap3 0.049566295 1.304814 -1.407024   yoked\nconflict
    ## 18        Arhgap4 0.056730474 1.246184  5.212728 yoked\nconsistent
    ## 19      Arhgef10l 0.097742444 1.009917  1.538117 yoked\nconsistent
    ## 20         Arid3b 0.078907879 1.102880 -2.903522   yoked\nconflict
    ## 21          Armt1 0.072768219 1.138058  1.613630 yoked\nconsistent
    ## 22          Arvcf 0.056417749 1.248584 -1.299618   yoked\nconflict
    ## 23           Asph 0.073196550 1.135509  1.234793 yoked\nconsistent
    ## 24          Atp5d 0.049335382 1.306842 -1.177138   yoked\nconflict
    ## 25            Axl 0.032624969 1.486450  2.151214 yoked\nconsistent
    ## 26        B4galt7 0.091946669 1.036464 -1.660627   yoked\nconflict
    ## 27          Bbof1 0.059271847 1.227152 -2.528912   yoked\nconflict
    ## 28           Bbs9 0.042060720 1.376123 -2.488483   yoked\nconflict
    ## 29           Bcl6 0.088333423 1.053875  1.074952 yoked\nconsistent
    ## 30         Bhlhb9 0.051653076 1.286904  1.181709 yoked\nconsistent
    ## 31           Bin3 0.087550965 1.057739 -2.675543   yoked\nconflict
    ## 32          Bnip2 0.094828933 1.023059 -1.392493   yoked\nconflict
    ## 33           Bod1 0.044642665 1.350250 -1.296866   yoked\nconflict
    ## 34            Bok 0.060699106 1.216818  2.881833 yoked\nconsistent
    ## 35           Braf 0.067484160 1.170798  1.059782 yoked\nconsistent
    ## 36  C130074G19Rik 0.077715664 1.109491 -2.312159   yoked\nconflict
    ## 37          C1ql3 0.095675335 1.019200 -1.043961   yoked\nconflict
    ## 38        C1qtnf6 0.027188129 1.565621 -5.378147   yoked\nconflict
    ## 39         C77370 0.079767068 1.098176  1.974578 yoked\nconsistent
    ## 40         Cacfd1 0.088913633 1.051032 -1.332728   yoked\nconflict
    ## 41        Cacna1c 0.058409341 1.233518 -1.189623   yoked\nconflict
    ## 42         Camk1g 0.058259297 1.234635 -1.255724   yoked\nconflict
    ## 43          Carm1 0.033860616 1.470305  1.177495 yoked\nconsistent
    ## 44        Carnmt1 0.072768219 1.138058  5.263765 yoked\nconsistent
    ## 45         Ccdc53 0.096665124 1.014730 -1.931319   yoked\nconflict
    ## 46         Ccdc59 0.056237885 1.249971 -2.030286   yoked\nconflict
    ## 47          Ccng1 0.084901312 1.071086  1.133409 yoked\nconsistent
    ## 48           Ccr2 0.042060720 1.376123 -5.696360   yoked\nconflict
    ## 49           Cd63 0.037079932 1.430861  2.552303 yoked\nconsistent
    ## 50         Cdadc1 0.093160115 1.030770 -1.301332   yoked\nconflict
    ## 51         Cdc123 0.095675335 1.019200 -1.055325   yoked\nconflict
    ## 52       Cdc42bpg 0.092123416 1.035630  4.783444 yoked\nconsistent
    ## 53         Cdkal1 0.077099653 1.112948  4.153773 yoked\nconsistent
    ## 54         Cep135 0.089506596 1.048145  2.545687 yoked\nconsistent
    ## 55          Ces2b 0.059511212 1.225401 -2.400447   yoked\nconflict
    ## 56            Cfp 0.092123416 1.035630 -1.762976   yoked\nconflict
    ## 57           Chd6 0.089163966 1.049811  1.426777 yoked\nconsistent
    ## 58          Ciao1 0.064736293 1.188852 -1.122153   yoked\nconflict
    ## 59         Cirh1a 0.089589175 1.047744  1.638313 yoked\nconsistent
    ## 60        Clec11a 0.092493448 1.033889 -2.709516   yoked\nconflict
    ## 61          Cnbd2 0.069303994 1.159242  1.512554 yoked\nconsistent
    ## 62         Cntrob 0.073259303 1.135137  3.142593 yoked\nconsistent
    ## 63         Commd9 0.092493448 1.033889  2.445479 yoked\nconsistent
    ## 64          Crocc 0.037079932 1.430861 -1.496366   yoked\nconflict
    ## 65          Csmd3 0.063070454 1.200174  1.924787 yoked\nconsistent
    ## 66           Ctsh 0.085297067 1.069066 -2.905902   yoked\nconflict
    ## 67         Cuedc1 0.023273579 1.633137  3.928244 yoked\nconsistent
    ## 68          Cxcl9 0.029104100 1.536046 -6.374635   yoked\nconflict
    ## 69         Cyb561 0.095362743 1.020621  1.184850 yoked\nconsistent
    ## 70       Cyb561d2 0.048833023 1.311286 -2.288928   yoked\nconflict
    ## 71         Cyp4v3 0.075654079 1.121168  4.991491 yoked\nconsistent
    ## 72          Daglb 0.076507339 1.116297 -1.311362   yoked\nconflict
    ## 73        Dclre1b 0.089634767 1.047524 -2.058020   yoked\nconflict
    ## 74          Ddx47 0.040773654 1.389620 -1.375397   yoked\nconflict
    ## 75          Ddx51 0.057058567 1.243679  1.699372 yoked\nconsistent
    ## 76          Decr2 0.095122027 1.021719  1.250513 yoked\nconsistent
    ## 77           Det1 0.058409341 1.233518  5.178426 yoked\nconsistent
    ## 78         Dhtkd1 0.091018152 1.040872 -1.918612   yoked\nconflict
    ## 79       Dnase1l2 0.092829786 1.032313  4.600381 yoked\nconsistent
    ## 80           Dnlz 0.093160115 1.030770  1.434436 yoked\nconsistent
    ## 81          Dock8 0.032893929 1.482884 -3.237576   yoked\nconflict
    ## 82        Dpy19l3 0.058303918 1.234302  1.313747 yoked\nconsistent
    ## 83          Dscr3 0.052925801 1.276333 -1.819022   yoked\nconflict
    ## 84           Dus2 0.080020091 1.096801  4.599887 yoked\nconsistent
    ## 85          Dus3l 0.047270078 1.325414 -1.227621   yoked\nconflict
    ## 86         Dusp16 0.088913633 1.051032  3.905808 yoked\nconsistent
    ## 87  E130309D02Rik 0.068279693 1.165708  1.854734 yoked\nconsistent
    ## 88           Ece1 0.041595475 1.380954  1.798412 yoked\nconsistent
    ## 89           Edc3 0.094828933 1.023059 -1.362216   yoked\nconflict
    ## 90         Efcab6 0.064139340 1.192876 -1.736664   yoked\nconflict
    ## 91           Elp2 0.042060720 1.376123 -1.272605   yoked\nconflict
    ## 92           Eml4 0.067317382 1.171873 -1.251351   yoked\nconflict
    ## 93           Eno4 0.067317382 1.171873 -4.280584   yoked\nconflict
    ## 94         Entpd1 0.070107240 1.154237  2.294689 yoked\nconsistent
    ## 95         Entpd2 0.096330906 1.016234  4.979568 yoked\nconsistent
    ## 96           Eny2 0.027120875 1.566696  1.906160 yoked\nconsistent
    ## 97          Ep300 0.091099278 1.040485  1.079507 yoked\nconsistent
    ## 98           Eps8 0.057728549 1.238609  2.223432 yoked\nconsistent
    ## 99           Epyc 0.040773654 1.389620  5.895348 yoked\nconsistent
    ## 100        Ero1lb 0.079767068 1.098176  1.831698 yoked\nconsistent
    ## 101         Ethe1 0.068279693 1.165708 -3.024445   yoked\nconflict
    ## 102            F3 0.045319614 1.343714  2.577718 yoked\nconsistent
    ## 103       Fam133b 0.086483006 1.063069 -1.265882   yoked\nconflict
    ## 104       Fam219b 0.099493142 1.002207  1.880278 yoked\nconsistent
    ## 105        Fam60a 0.037554894 1.425333 -2.887595   yoked\nconflict
    ## 106         Farp1 0.087612753 1.057433  1.470293 yoked\nconsistent
    ## 107       Fastkd2 0.082686740 1.082564  3.051095 yoked\nconsistent
    ## 108         Fbln1 0.097445033 1.011240  4.595828 yoked\nconsistent
    ## 109         Fbxl4 0.054431476 1.264150  4.678084 yoked\nconsistent
    ## 110        Fbxo25 0.086680174 1.062080 -1.179674   yoked\nconflict
    ## 111         Fcgr3 0.056237885 1.249971  5.152940 yoked\nconsistent
    ## 112          Fgd1 0.063177120 1.199440  2.495098 yoked\nconsistent
    ## 113          Fgd3 0.068279693 1.165708  4.930926 yoked\nconsistent
    ## 114        Fkbp14 0.082173698 1.085267  3.631638 yoked\nconsistent
    ## 115         Fkbpl 0.068864645 1.162004 -2.642512   yoked\nconflict
    ## 116         Flrt3 0.069847384 1.155850 -1.980703   yoked\nconflict
    ## 117         Fstl4 0.080095599 1.096391 -1.936575   yoked\nconflict
    ## 118       Galnt10 0.040773654 1.389620  5.421251 yoked\nconsistent
    ## 119          Gcc2 0.070789185 1.150033 -1.098204   yoked\nconflict
    ## 120         Gdpd2 0.049566295 1.304814 -4.671692   yoked\nconflict
    ## 121        Gemin4 0.027146428 1.566287  3.510533 yoked\nconsistent
    ## 122         Gfpt2 0.043450289 1.362007 -4.474012   yoked\nconflict
    ## 123        Glyctk 0.042060720 1.376123  5.295386 yoked\nconsistent
    ## 124       Gm10146 0.067203683 1.172607 -6.823282   yoked\nconflict
    ## 125       Gm20715 0.057718203 1.238687 -3.355164   yoked\nconflict
    ## 126       Gm38393 0.067203683 1.172607 -3.177709   yoked\nconflict
    ## 127       Gm43796 0.051627647 1.287118  5.618344 yoked\nconsistent
    ## 128        Gm6741 0.099969863 1.000131  5.100647 yoked\nconsistent
    ## 129        Gm9803 0.038796683 1.411205 -4.895748   yoked\nconflict
    ## 130         Gmpr2 0.093187367 1.030643 -1.558630   yoked\nconflict
    ## 131        Gpank1 0.031082263 1.507487 -2.663846   yoked\nconflict
    ## 132         Grid1 0.038796683 1.411205  1.283537 yoked\nconsistent
    ## 133        Grin2c 0.080320367 1.095174  1.521991 yoked\nconsistent
    ## 134          Guf1 0.095675335 1.019200 -1.411535   yoked\nconflict
    ## 135          Gys1 0.068567494 1.163882 -1.633591   yoked\nconflict
    ## 136        H2-Ke6 0.091392070 1.039091 -1.530689   yoked\nconflict
    ## 137          Heg1 0.097412881 1.011384  1.772297 yoked\nconsistent
    ## 138         Hmgcl 0.096665124 1.014730 -1.479741   yoked\nconflict
    ## 139         Hmgn2 0.089988859 1.045811 -1.067040   yoked\nconflict
    ## 140        Homer3 0.058303918 1.234302  1.742975 yoked\nconsistent
    ## 141        Hrasls 0.064904893 1.187723 -2.080573   yoked\nconflict
    ## 142           Id4 0.042060720 1.376123  2.487073 yoked\nconsistent
    ## 143          Ier3 0.067205912 1.172593  5.227487 yoked\nconsistent
    ## 144            Ik 0.040696917 1.390438 -1.221337   yoked\nconflict
    ## 145        Ikbkap 0.084901312 1.071086 -1.081993   yoked\nconflict
    ## 146      Irak1bp1 0.094828933 1.023059 -1.150500   yoked\nconflict
    ## 147       Irf2bp1 0.050463698 1.297021 -1.164074   yoked\nconflict
    ## 148         Jade1 0.064139340 1.192876  1.981833 yoked\nconsistent
    ## 149         Jmjd6 0.063070454 1.200174 -1.548976   yoked\nconflict
    ## 150         Jmjd8 0.070637117 1.150967 -1.381619   yoked\nconflict
    ## 151         Kat2b 0.069847384 1.155850 -1.815229   yoked\nconflict
    ## 152        Katnb1 0.053034963 1.275438 -1.079787   yoked\nconflict
    ## 153        Kcnj12 0.059511212 1.225401 -3.219190   yoked\nconflict
    ## 154         Kcnn1 0.068279693 1.165708  2.016361 yoked\nconsistent
    ## 155         Kcnu1 0.045237693 1.344500 -3.812913   yoked\nconflict
    ## 156        Kctd15 0.082275653 1.084729  1.434824 yoked\nconsistent
    ## 157        Kdelr2 0.029360780 1.532232  3.620413 yoked\nconsistent
    ## 158       Laptm4a 0.077099653 1.112948  1.019506 yoked\nconsistent
    ## 159          Ldah 0.048833023 1.311286  1.743735 yoked\nconsistent
    ## 160         Lnpk1 0.077368302 1.111437  1.399997 yoked\nconsistent
    ## 161          Lrp5 0.057864214 1.237590  3.801397 yoked\nconsistent
    ## 162         Lypd1 0.088913633 1.051032  1.246485 yoked\nconsistent
    ## 163          Mafk 0.031796863 1.497616  1.745233 yoked\nconsistent
    ## 164        Map3k2 0.089645582 1.047471  1.591062 yoked\nconsistent
    ## 165        March9 0.060299129 1.219689  1.663745 yoked\nconsistent
    ## 166          Mdp1 0.051653076 1.286904  2.034451 yoked\nconsistent
    ## 167          Mest 0.095882891 1.018259  2.810353 yoked\nconsistent
    ## 168        Mif4gd 0.095675335 1.019200  4.511593 yoked\nconsistent
    ## 169        Mrpl16 0.077801149 1.109014  3.324596 yoked\nconsistent
    ## 170        Mrpl53 0.061919886 1.208170  2.953771 yoked\nconsistent
    ## 171         Mtmr7 0.061863668 1.208564  1.310100 yoked\nconsistent
    ## 172         Mtmr9 0.081642233 1.088085 -1.143233   yoked\nconflict
    ## 173        Mtrf1l 0.028418255 1.546403 -2.907037   yoked\nconflict
    ## 174          Mtx1 0.085297067 1.069066 -1.367930   yoked\nconflict
    ## 175           Mut 0.094828933 1.023059 -1.340435   yoked\nconflict
    ## 176         Myo10 0.040426933 1.393329  2.072325 yoked\nconsistent
    ## 177         Naa20 0.078907879 1.102880 -1.256523   yoked\nconflict
    ## 178         Naglu 0.084891880 1.071134 -2.670780   yoked\nconflict
    ## 179         Nalcn 0.078654753 1.104275  1.134995 yoked\nconsistent
    ## 180        Nanos1 0.015855592 1.799818 -2.019823   yoked\nconflict
    ## 181          Nans 0.090623777 1.042758  3.252657 yoked\nconsistent
    ## 182         Ndst2 0.042584996 1.370743 -2.008448   yoked\nconflict
    ## 183         Nedd9 0.057718203 1.238687  3.269151 yoked\nconsistent
    ## 184         Nell1 0.068279693 1.165708  1.510942 yoked\nconsistent
    ## 185        Nfkbia 0.077715664 1.109491  1.581671 yoked\nconsistent
    ## 186          Npc1 0.093836884 1.027626  1.516286 yoked\nconsistent
    ## 187        Nploc4 0.058409341 1.233518 -1.170620   yoked\nconflict
    ## 188       Nr2c2ap 0.074606969 1.127221 -1.778607   yoked\nconflict
    ## 189          Nrf1 0.004277175 2.368843 -2.102394   yoked\nconflict
    ## 190        Nudt11 0.098174114 1.008003 -1.108972   yoked\nconflict
    ## 191         Nup85 0.027188129 1.565621  2.156705 yoked\nconsistent
    ## 192         Nxph3 0.023091480 1.636548 -6.644329   yoked\nconflict
    ## 193         Oas1b 0.086156046 1.064714 -3.667190   yoked\nconflict
    ## 194          Ogg1 0.059511212 1.225401  5.383167 yoked\nconsistent
    ## 195          Optn 0.041595475 1.380954  6.025644 yoked\nconsistent
    ## 196          Orc6 0.096856343 1.013872  3.190823 yoked\nconsistent
    ## 197          P3h4 0.053867472 1.268673 -3.974249   yoked\nconflict
    ## 198         Padi2 0.072897080 1.137290  3.412527 yoked\nconsistent
    ## 199       Pcdhb16 0.056594653 1.247225  3.429710 yoked\nconsistent
    ## 200        Pcdhb6 0.073154008 1.135762  4.920790 yoked\nconsistent
    ## 201        Pcp4l1 0.085297067 1.069066  2.556690 yoked\nconsistent
    ## 202          Pdcl 0.078035802 1.107706 -2.034106   yoked\nconflict
    ## 203         Pdia4 0.075411030 1.122565  1.750810 yoked\nconsistent
    ## 204          Pdpn 0.082173698 1.085267  5.351491 yoked\nconsistent
    ## 205          Pdpr 0.073062173 1.136307  1.850812 yoked\nconsistent
    ## 206         Pds5a 0.068279693 1.165708  1.195244 yoked\nconsistent
    ## 207          Pecr 0.076507339 1.116297 -3.322210   yoked\nconflict
    ## 208        Pla2g7 0.067103028 1.173258  1.576937 yoked\nconsistent
    ## 209         Plaur 0.069303994 1.159242  4.951917 yoked\nconsistent
    ## 210         Plin2 0.077715664 1.109491 -3.909368   yoked\nconflict
    ## 211        Polr2d 0.082686740 1.082564  4.017044 yoked\nconsistent
    ## 212       Polr3gl 0.093160115 1.030770 -1.274596   yoked\nconflict
    ## 213         Ppdpf 0.076291995 1.117521  1.514740 yoked\nconsistent
    ## 214      Ppp1r13l 0.032624969 1.486450 -3.562537   yoked\nconflict
    ## 215          Ppt2 0.048923039 1.310487  1.781435 yoked\nconsistent
    ## 216         Pqlc2 0.085779463 1.066617 -1.882595   yoked\nconflict
    ## 217          Prcp 0.081887116 1.086784  4.554132 yoked\nconsistent
    ## 218         Prkdc 0.076291995 1.117521  1.504940 yoked\nconsistent
    ## 219       Prkrip1 0.092493448 1.033889 -1.564664   yoked\nconflict
    ## 220         Prmt2 0.084029666 1.075567  1.288555 yoked\nconsistent
    ## 221         Prmt3 0.087612753 1.057433  1.904318 yoked\nconsistent
    ## 222         Prpf3 0.068279693 1.165708  2.519612 yoked\nconsistent
    ## 223         Prpf4 0.070637117 1.150967 -1.604656   yoked\nconflict
    ## 224          Pwp2 0.043128105 1.365240 -1.242212   yoked\nconflict
    ## 225         Rab13 0.077715664 1.109491 -4.543270   yoked\nconflict
    ## 226         Rab31 0.097703142 1.010091  1.441188 yoked\nconsistent
    ## 227         Rab43 0.072125351 1.141912 -1.469960   yoked\nconflict
    ## 228         Rabif 0.067205912 1.172593 -1.552053   yoked\nconflict
    ## 229         Ramp2 0.082173698 1.085267  1.949395 yoked\nconsistent
    ## 230        Rangrf 0.085297067 1.069066  4.778751 yoked\nconsistent
    ## 231         Rap2b 0.010490243 1.979214 -1.026861   yoked\nconflict
    ## 232         Rccd1 0.027246333 1.564692 -3.288941   yoked\nconflict
    ## 233          Rem2 0.088913633 1.051032  1.560420 yoked\nconsistent
    ## 234          Rfc1 0.076942499 1.113834  1.933228 yoked\nconsistent
    ## 235         Rftn2 0.030219484 1.519713  3.428065 yoked\nconsistent
    ## 236         Rgs16 0.077715664 1.109491  5.093797 yoked\nconsistent
    ## 237          Rgs6 0.077801149 1.109014  5.061268 yoked\nconsistent
    ## 238          Rin2 0.058409341 1.233518  2.013528 yoked\nconsistent
    ## 239        Rnf166 0.080039064 1.096698 -1.537154   yoked\nconflict
    ## 240        Rnf207 0.039413555 1.404354  5.934581 yoked\nconsistent
    ## 241          Rnls 0.075883071 1.119855 -6.360421   yoked\nconflict
    ## 242         Rnpc3 0.094129578 1.026274  2.271311 yoked\nconsistent
    ## 243 RP23-220F20.2 0.031215943 1.505624 -2.230303   yoked\nconflict
    ## 244         Rpap3 0.077715664 1.109491  2.642317 yoked\nconsistent
    ## 245        Sacm1l 0.048833023 1.311286  1.470789 yoked\nconsistent
    ## 246         Sars2 0.095751221 1.018856 -1.876766   yoked\nconflict
    ## 247          Sbsn 0.046564355 1.331946 -2.995744   yoked\nconflict
    ## 248         Sdad1 0.023416918 1.630470 -1.495033   yoked\nconflict
    ## 249         Senp8 0.076291995 1.117521  4.756561 yoked\nconsistent
    ## 250         Serf1 0.095675335 1.019200 -3.370162   yoked\nconflict
    ## 251         Sesn3 0.030797928 1.511479  1.304499 yoked\nconsistent
    ## 252        Sfmbt1 0.077715664 1.109491  1.352141 yoked\nconsistent
    ## 253           Sfn 0.099969863 1.000131 -4.263196   yoked\nconflict
    ## 254         Sfxn5 0.067203683 1.172607  1.239155 yoked\nconsistent
    ## 255        Shisa5 0.051653076 1.286904  1.757217 yoked\nconsistent
    ## 256       Shroom4 0.097997489 1.008785  4.876928 yoked\nconsistent
    ## 257         Shtn1 0.048195366 1.316995  2.724997 yoked\nconsistent
    ## 258        Slain1 0.077715664 1.109491 -1.419156   yoked\nconflict
    ## 259       Slc12a2 0.091926263 1.036560 -1.508166   yoked\nconflict
    ## 260      Slc25a35 0.072858495 1.137520  4.718411 yoked\nconsistent
    ## 261       Slc38a9 0.061435642 1.211580  2.353099 yoked\nconsistent
    ## 262        Slc4a2 0.071313688 1.146827 -1.702967   yoked\nconflict
    ## 263         Smad5 0.070022502 1.154762  1.861525 yoked\nconsistent
    ## 264         Smdt1 0.072125351 1.141912 -1.530848   yoked\nconflict
    ## 265          Smg6 0.079567620 1.099264 -1.023022   yoked\nconflict
    ## 266        Smim17 0.068279693 1.165708 -1.727504   yoked\nconflict
    ## 267        Snrpd2 0.030797928 1.511479 -1.731473   yoked\nconflict
    ## 268          Snx8 0.071858244 1.143523 -2.594352   yoked\nconflict
    ## 269         Soat1 0.077715664 1.109491  1.520034 yoked\nconsistent
    ## 270         Socs4 0.053034963 1.275438 -2.448752   yoked\nconflict
    ## 271         Spag7 0.041595475 1.380954 -1.761508   yoked\nconflict
    ## 272          Spi1 0.097185220 1.012400  4.963456 yoked\nconsistent
    ## 273         Ssna1 0.094828933 1.023059  2.130753 yoked\nconsistent
    ## 274          Ssr1 0.092493448 1.033889  1.188571 yoked\nconsistent
    ## 275    St6galnac3 0.086087965 1.065058  4.898459 yoked\nconsistent
    ## 276         Sugp1 0.051612729 1.287243  1.352067 yoked\nconsistent
    ## 277         Supt3 0.082934817 1.081263 -2.145788   yoked\nconflict
    ## 278        Swsap1 0.086099708 1.064998 -1.775499   yoked\nconflict
    ## 279       Syndig1 0.068279693 1.165708  1.568618 yoked\nconsistent
    ## 280         Syvn1 0.048923039 1.310487  1.462197 yoked\nconsistent
    ## 281           Tbp 0.078842871 1.103238 -1.630764   yoked\nconflict
    ## 282        Tceal3 0.042060720 1.376123 -1.019052   yoked\nconflict
    ## 283          Tcta 0.097997489 1.008785 -1.656222   yoked\nconflict
    ## 284         Tdrd3 0.080020091 1.096801 -1.500826   yoked\nconflict
    ## 285         Thap4 0.031109226 1.507111 -1.972335   yoked\nconflict
    ## 286       Timm10b 0.075406340 1.122592 -1.700145   yoked\nconflict
    ## 287        Timm44 0.041595475 1.380954 -1.547182   yoked\nconflict
    ## 288         Tmco3 0.088617852 1.052479 -1.659585   yoked\nconflict
    ## 289       Tmem119 0.088913633 1.051032  2.869104 yoked\nconsistent
    ## 290       Tmem121 0.057428254 1.240874 -1.758438   yoked\nconflict
    ## 291      Tmem132e 0.087165449 1.059656  4.789403 yoked\nconsistent
    ## 292       Tmem143 0.099660863 1.001475  2.710642 yoked\nconsistent
    ## 293       Tmem186 0.087612753 1.057433  4.063019 yoked\nconsistent
    ## 294      Tmem229b 0.061919886 1.208170  2.038292 yoked\nconsistent
    ## 295      Tmem254a 0.038839854 1.410722 -3.617706   yoked\nconflict
    ## 296      Tmem254b 0.038839854 1.410722 -3.617706   yoked\nconflict
    ## 297       Tmem266 0.085297067 1.069066  4.727446 yoked\nconsistent
    ## 298       Tmem87b 0.032624969 1.486450  2.157978 yoked\nconsistent
    ## 299        Tmem94 0.049698081 1.303660 -1.056107   yoked\nconflict
    ## 300          Tmx3 0.095764023 1.018798  1.189006 yoked\nconsistent
    ## 301       Tomm40l 0.035203969 1.453408 -1.419861   yoked\nconflict
    ## 302      Tor1aip2 0.084029666 1.075567  1.394125 yoked\nconsistent
    ## 303         Trhde 0.085297067 1.069066  2.259483 yoked\nconsistent
    ## 304        Trim45 0.066572944 1.176702  3.717858 yoked\nconsistent
    ## 305         Tshz3 0.035654647 1.447884 -3.237049   yoked\nconflict
    ## 306       Tspan14 0.085297067 1.069066  2.933253 yoked\nconsistent
    ## 307         Tubb6 0.083493569 1.078347  4.690757 yoked\nconsistent
    ## 308       Twistnb 0.080938232 1.091846  4.623397 yoked\nconsistent
    ## 309          Tyw5 0.067317382 1.171873  4.732481 yoked\nconsistent
    ## 310         Ubac2 0.023273579 1.633137  1.716629 yoked\nconsistent
    ## 311       Ubash3a 0.032624969 1.486450 -5.793790   yoked\nconflict
    ## 312          Ubn2 0.045237693 1.344500 -1.361992   yoked\nconflict
    ## 313        Unc13c 0.082686740 1.082564  4.798918 yoked\nconsistent
    ## 314          Urb1 0.040576139 1.391729 -1.924196   yoked\nconflict
    ## 315          Use1 0.058409341 1.233518  1.896321 yoked\nconsistent
    ## 316         Usp36 0.057728549 1.238609  1.400985 yoked\nconsistent
    ## 317       Vipas39 0.090623777 1.042758 -1.187557   yoked\nconflict
    ## 318         Vipr1 0.064296777 1.191811 -1.536749   yoked\nconflict
    ## 319         Vma21 0.079567620 1.099264  1.289932 yoked\nconsistent
    ## 320        Vstm2b 0.077099653 1.112948 -1.591188   yoked\nconflict
    ## 321         Vti1a 0.058303918 1.234302 -1.056028   yoked\nconflict
    ## 322         Wash1 0.095764023 1.018798 -1.223746   yoked\nconflict
    ## 323         Wbp11 0.067432641 1.171130  1.110021 yoked\nconsistent
    ## 324         Wdfy2 0.085297067 1.069066  3.111324 yoked\nconsistent
    ## 325         Wdpcp 0.071291019 1.146965 -2.514469   yoked\nconflict
    ## 326         Whsc1 0.077368302 1.111437  1.227780 yoked\nconsistent
    ## 327           Wls 0.031682073 1.499186  2.900761 yoked\nconsistent
    ## 328          Zfat 0.080217839 1.095729  2.747299 yoked\nconsistent
    ## 329        Zfp146 0.058409341 1.233518 -1.552447   yoked\nconflict
    ## 330        Zfp184 0.095740251 1.018905  5.301730 yoked\nconsistent
    ## 331       Zfp385a 0.060699106 1.216818 -1.803569   yoked\nconflict
    ## 332        Zfp420 0.042060720 1.376123 -2.714897   yoked\nconflict
    ## 333        Zfp493 0.097654248 1.010309 -3.431675   yoked\nconflict
    ## 334        Zfp512 0.088913633 1.051032  1.473838 yoked\nconsistent
    ## 335       Zfp518a 0.094828933 1.023059  4.726568 yoked\nconsistent
    ## 336        Zfp526 0.092493448 1.033889  2.798728 yoked\nconsistent
    ## 337        Zfp644 0.085297067 1.069066  1.049193 yoked\nconsistent
    ## 338        Zfp707 0.077715664 1.109491  4.092779 yoked\nconsistent
    ## 339        Zfp746 0.095796025 1.018653 -1.563321   yoked\nconflict
    ## 340         Zfp84 0.092493448 1.033889  2.360325 yoked\nconsistent
    ## 341        Zfp937 0.058303918 1.234302  3.999816 yoked\nconsistent
    ## 342          Zic3 0.019987044 1.699251 -5.901697   yoked\nconflict
    ## 343       Zkscan3 0.087523147 1.057877  1.533165 yoked\nconsistent
    ## 344         Znfx1 0.091392070 1.039091 -1.183631   yoked\nconflict

    write.csv(CA1learnGO, "../data/04_CA1learnGO.csv", row.names = F)
    write.csv(CA1stressGO, "../data/04_CA1stressGO.csv", row.names = F)
