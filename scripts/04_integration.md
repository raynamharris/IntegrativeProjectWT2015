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
    head(CA1learn)

    ##      gene
    ## 1    Fosb
    ## 2 Gm13889
    ## 3    Irs1
    ## 4   Kcnc2
    ## 5   Lemd3
    ## 6   Npas4

    head(CA1trainpvals)

    ##            gene       padj  logpadj       lfc         direction
    ## 1 1810022K09Rik 0.08822746 1.054396 -2.103747 yoked\nconsistent
    ## 2 9430015G10Rik 0.05480315 1.261195 -1.942464 yoked\nconsistent
    ## 3         Aagab 0.08737294 1.058623 -1.129414 yoked\nconsistent
    ## 4          Abl1 0.04756029 1.322755  1.540018        consistent
    ## 5         Acad8 0.09923614 1.003330  1.777297        consistent
    ## 6         Acads 0.01979073 1.703538 -4.336036 yoked\nconsistent

    CA1forGO <- left_join(CA1trainpvals,CA1learn)

    ## Joining, by = "gene"

    CA1forGO

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
    CA1sGO <- left_join(CA1trainpvals,CA1stress)

    ## Joining, by = "gene"

    ## Warning: Column `gene` joining character vector and factor, coercing into
    ## character vector

    CA1sGO

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
