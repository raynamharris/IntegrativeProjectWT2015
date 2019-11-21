    library(tidyverse)
    library(forcats)
    library(cowplot) ## for some easy to use themes
    library(DESeq2) ## for gene expression analysis
    library(pheatmap)
    library(viridis)
    library(Rtsne) # for tSNE
    library(scales)


    library(BiocParallel)
    register(MulticoreParam(6))

    ## load functions 
    source("figureoptions.R")
    source("functions_RNAseq.R")

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/02f_trainedvyoked/', cache = F)

    a.countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)

    a.colData <- read.csv("../data/02a_colData.csv", header = T)
    a.colData <- a.colData 
    a.colData$training <- factor(a.colData$training, levels = c("yoked", "trained"))

    DGdds2 <- returndds2("DG") 

    ## [1] "DG"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 54 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    CA3dds2 <- returndds2("CA3") 

    ## [1] "CA3"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 53 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    CA1dds2 <- returndds2("CA1") 

    ## [1] "CA1"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 98 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    DGvsd <- returnvsds2(DGdds2)

    ##               143A-DG-1 143B-DG-1 143D-DG-3 144A-DG-2 144C-DG-2 144D-DG-2
    ## 0610007P14Rik  6.497231  6.520310  6.927199  6.622006  6.516323  6.739843
    ## 0610009B22Rik  5.933401  5.821822  5.483341  5.697426  5.971474  5.830043
    ## 0610009L18Rik  5.582826  5.757564  5.597896  5.205323  5.677909  5.676574
    ##               145A-DG-2 145B-DG-1 146A-DG-2 146B-DG-2 146C-DG-4 146D-DG-3
    ## 0610007P14Rik  6.669012  6.257587  6.419110  6.134992  7.017216  7.479295
    ## 0610009B22Rik  6.084135  5.939702  5.838619  6.134992  6.195432  5.205323
    ## 0610009L18Rik  5.443572  5.205323  5.838619  5.205323  5.205323  5.205323
    ##               147C-DG-3 147D-DG-1 148A-DG-3 148B-DG-4
    ## 0610007P14Rik  6.418946  6.592025  6.603472  6.467815
    ## 0610009B22Rik  6.054862  6.038992  5.754020  5.664124
    ## 0610009L18Rik  5.444089  5.874552  5.675946  5.530425

    CA3vsd <- returnvsds2(CA3dds2)

    ##               143A-CA3-1 144A-CA3-2 144B-CA3-1 144C-CA3-2 144D-CA3-2
    ## 0610007P14Rik   7.053774   7.561649   7.233758   6.996898   6.880131
    ## 0610009B22Rik   6.452876   6.973831   6.531009   6.908274   6.375915
    ## 0610009L18Rik   6.030451   6.670471   6.269295   6.224806   6.269568
    ##               145A-CA3-2 146A-CA3-2 146B-CA3-2 146D-CA3-3 147C-CA3-3
    ## 0610007P14Rik   6.784496   7.233784   6.757449   7.122555   7.137734
    ## 0610009B22Rik   6.879922   6.236353   6.895220   6.468228   6.353976
    ## 0610009L18Rik   5.734766   6.236353   5.991499   6.142415   6.111529
    ##               147D-CA3-1 148A-CA3-3 148B-CA3-4
    ## 0610007P14Rik   6.840862   6.923812   7.262408
    ## 0610009B22Rik   6.540744   6.532445   6.689017
    ## 0610009L18Rik   6.116202   6.169424   6.213015

    CA1vsd <- returnvsds2(CA1dds2)

    ##               143B-CA1-1 143C-CA1-1 143D-CA1-3 144A-CA1-2 144B-CA1-1
    ## 0610007P14Rik   7.415967   7.026659   7.234953   7.242002   7.278148
    ## 0610009B22Rik   6.904844   6.767590   6.124712   6.819405   6.644173
    ## 0610009L18Rik   6.666449   6.334833   6.124712   6.507772   6.264241
    ##               144C-CA1-2 145A-CA1-2 145B-CA1-1 146A-CA1-2 146B-CA1-2
    ## 0610007P14Rik   7.075661   7.291669   7.252495   7.032181   6.904449
    ## 0610009B22Rik   6.595087   6.743767   6.747191   6.782460   6.657871
    ## 0610009L18Rik   6.297128   6.599020   6.404817   6.606986   6.124712
    ##               146C-CA1-4 146D-CA1-3 147C-CA1-3 148A-CA1-3 148B-CA1-4
    ## 0610007P14Rik   7.290895   7.357361   7.149545   7.213656   7.068189
    ## 0610009B22Rik   6.636289   6.782490   6.641399   6.854020   6.124712
    ## 0610009L18Rik   6.703932   6.124712   6.302770   6.507437   6.124712

    print("DG")

    ## [1] "DG"

    res_summary_subfield(DGdds2, c("training", "trained", "yoked"))

    ## [1] "training" "trained"  "yoked"   
    ## [1] 214
    ## 
    ## out of 17006 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 177, 1%
    ## LFC < 0 (down)     : 37, 0.22%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 6929, 41%
    ## (mean count < 13)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    print("CA3")

    ## [1] "CA3"

    res_summary_subfield(CA3dds2, c("training", "trained", "yoked"))

    ## [1] "training" "trained"  "yoked"   
    ## [1] 0
    ## 
    ## out of 16497 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 27, 0.16%
    ## low counts [2]     : 5, 0.03%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    print("CA1")

    ## [1] "CA1"

    res_summary_subfield(CA1dds2, c("training", "trained", "yoked"))

    ## [1] "training" "trained"  "yoked"   
    ## [1] 16
    ## 
    ## out of 16846 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.0059%
    ## LFC < 0 (down)     : 15, 0.089%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 2619, 16%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    contrast1 <- resvals2(DGdds2, contrastvector = c("training", "trained", "yoked"), mypval = 0.1) # 206

    ## [1] 214

    contrast2 <- resvals2(CA1dds2, contrastvector = c("training", "trained", "yoked"), mypval = 0.1) # 16

    ## [1] 16

    options(digits = 3)

    listofDEGstrainedvyoked <- function(mydds, myitssue){
      res <- results(mydds, 
                     contrast = c("training", "trained", "yoked"), 
                     independentFiltering = T)
      
      print(paste(myitssue, "trained vs yoked", sep = " "))
      
      data <- data.frame(gene = row.names(res),
                         lfc = round(res$log2FoldChange,2),
                         padj = res$padj,
                         tissue = myitssue)
      data <- data %>% 
            dplyr::mutate_each(funs(toupper)) %>% 
        dplyr::filter(padj < 0.1) %>% 
        dplyr::mutate(direction = ifelse(lfc >0, "increased", "decreased"),
                      comparison = paste("trained", "yoked", sep = "-")) %>% 
        dplyr::select(gene, lfc,padj, direction, tissue, comparison )  %>% 
        dplyr::arrange(desc(direction), gene) %>% 
        droplevels()
      print(head(data))
      return(data)
    }

    DGDEGs <- listofDEGstrainedvyoked(DGdds2, "DG")

    ## [1] "DG trained vs yoked"

    ## Warning: funs() is soft deprecated as of dplyr 0.8.0
    ## please use list() instead
    ## 
    ##   # Before:
    ##   funs(name = f(.))
    ## 
    ##   # After: 
    ##   list(name = ~ f(.))
    ## This warning is displayed once per session.

    ##            gene  lfc                 padj direction tissue    comparison
    ## 1 1190002N15RIK 1.64 0.000245166241490228 increased     DG trained-yoked
    ## 2         ABHD2 0.86   0.0153326673524323 increased     DG trained-yoked
    ## 3       ADAMTS1 1.88  0.00187777546749527 increased     DG trained-yoked
    ## 4         ADRB1 0.98   0.0311012281104977 increased     DG trained-yoked
    ## 5           AHR 1.59   0.0191420911940072 increased     DG trained-yoked
    ## 6        AMIGO2 1.36   0.0254399458003005 increased     DG trained-yoked

    CA1DEGs <- listofDEGstrainedvyoked(CA1dds2, "CA1")

    ## [1] "CA1 trained vs yoked"
    ##     gene   lfc                padj direction tissue    comparison
    ## 1 GLCCI1  1.98 0.00143873897497547 increased    CA1 trained-yoked
    ## 2  AHDC1 -1.61  0.0650203882596328 decreased    CA1 trained-yoked
    ## 3   BMT2 -2.48  0.0566256621961535 decreased    CA1 trained-yoked
    ## 4  CTCFL -2.61  0.0363254357303014 decreased    CA1 trained-yoked
    ## 5 FN3KRP -1.43  0.0480138873382015 decreased    CA1 trained-yoked
    ## 6   GNAZ -3.11  0.0112083934063577 decreased    CA1 trained-yoked

    write.csv(DGDEGs, "../data/02f_DG_DEGs.csv", row.names = F)
    write.csv(CA1DEGs, "../data/02f_CA1_DEGs.csv", row.names = F)

    upDG <- DGDEGs %>% filter(lfc > 0)
    upDG <- as.vector(upDG$gene)
    upDG

    ##   [1] "1190002N15RIK" "ABHD2"         "ADAMTS1"       "ADRB1"        
    ##   [5] "AHR"           "AMIGO2"        "ANKRD13A"      "ANKRD28"      
    ##   [9] "ANKRD33B"      "APAF1"         "ARID5B"        "ARL13B"       
    ##  [13] "ARL4A"         "ARL4D"         "ARL5B"         "ARMCX5"       
    ##  [17] "ARPP21"        "ATF3"          "B3GNT2"        "BACH1"        
    ##  [21] "BTG2"          "C2CD4B"        "CCNK"          "CIART"        
    ##  [25] "CITED2"        "CLDN12"        "CNNM1"         "CPEB4"        
    ##  [29] "CTNND1"        "CUL3"          "CWC25"         "CXADR"        
    ##  [33] "CYP51"         "DNAJA1"        "DNAJB1"        "DNAJB4"       
    ##  [37] "DUSP14"        "DUSP16"        "DUSP4"         "DUSP6"        
    ##  [41] "DUSP8"         "DYRK2"         "EGR1"          "EGR3"         
    ##  [45] "EIF5"          "EPRS"          "ERF"           "FAM107B"      
    ##  [49] "FAM118A"       "FBXW7"         "FERMT2"        "FLRT3"        
    ##  [53] "FOS"           "FOSB"          "FOXG1"         "FOXO1"        
    ##  [57] "FZD4"          "GAD1"          "GM13889"       "GMEB2"        
    ##  [61] "GPR19"         "HECA"          "HMGCR"         "HS6ST1"       
    ##  [65] "HSPA1A"        "HSPH1"         "IL16"          "ING2"         
    ##  [69] "IRF2BP2"       "IRS1"          "JDP2"          "JMJD1C"       
    ##  [73] "JUNB"          "JUND"          "KCNA4"         "KDM6B"        
    ##  [77] "KDM7A"         "KLF2"          "KLF6"          "LBH"          
    ##  [81] "LCMT2"         "LEMD3"         "LMNA"          "LONRF1"       
    ##  [85] "LRRTM2"        "MARCH11"       "MED7"          "MFAP3L"       
    ##  [89] "MN1"           "MYC"           "NAF1"          "NAP1L5"       
    ##  [93] "NEDD9"         "NEFM"          "NFIL3"         "NR4A1"        
    ##  [97] "NR4A2"         "NR4A3"         "NUAK1"         "ODC1"         
    ## [101] "OLFML2B"       "PAK6"          "PEG10"         "PELI1"        
    ## [105] "PER1"          "PER2"          "PHLDA1"        "PIGA"         
    ## [109] "PLAGL1"        "PLK3"          "POU3F3"        "PPP1R15A"     
    ## [113] "PRPF38B"       "RANBP2"        "RASL11A"       "RASL11B"      
    ## [117] "RGMB"          "RGS4"          "SCG2"          "SH2D3C"       
    ## [121] "SIAH2"         "SLC25A25"      "SLC2A3"        "SLC45A4"      
    ## [125] "SLITRK5"       "SOWAHC"        "SOX9"          "SRF"          
    ## [129] "STMN4"         "SYT4"          "THBS1"         "TIPARP"       
    ## [133] "TNIP2"         "TRA2B"         "TRIB1"         "TSC22D2"      
    ## [137] "USPL1"         "ZBTB33"        "ZFAND5"        "ZFP275"       
    ## [141] "ZFP654"        "ZFP869"

    downDG <- DGDEGs %>% filter(lfc < 0)
    downDG <- as.vector(downDG$gene)
    downDG

    ##  [1] "ANKRD27"  "BC048403" "BMT2"     "CCDC32"   "CECR6"    "COQ2"    
    ##  [7] "CPNE7"    "CTCFL"    "DPYSL2"   "EEF1E1"   "GNAZ"     "GPI1"    
    ## [13] "GYG"      "IGF2BP2"  "KLKB1"    "LRRC45"   "LYSMD4"   "MC1R"    
    ## [19] "NEUROD6"  "NXF1"     "PDE6A"    "PGAM2"    "PLCH2"    "PRUNE2"  
    ## [25] "PXN"      "RBM47"    "SCOC"     "SENP8"    "SLC5A5"   "SRGAP1"  
    ## [31] "STAC2"    "SV2B"     "TMEM170"  "TSPYL3"   "TUBB4A"   "ZFP207"  
    ## [37] "ZFP668"

    upCA1 <- CA1DEGs %>% filter(lfc > 0)
    upCA1 <- as.vector(upCA1$gene)
    upCA1

    ## [1] "GLCCI1"

    downCA1 <- CA1DEGs %>% filter(lfc < 0)
    downCA1 <- as.vector(downCA1$gene)
    downCA1

    ##  [1] "AHDC1"   "BMT2"    "CTCFL"   "FN3KRP"  "GNAZ"    "IGF2BP2" "IL4RA"  
    ##  [8] "INHBB"   "KHNYN"   "KLKB1"   "MX1"     "PDE6A"   "STAC3"   "STOX2"  
    ## [15] "TMEM170"

    ## DG 
    DEGes <- assay(vst(DGdds2))
    DEGes <- cbind(DEGes, contrast1)
    DEGes <- as.data.frame(DEGes) # convert matrix to dataframe

    DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe
    DEGes$rownames <- str_to_upper(DEGes$rownames) ## uppercase gene names
    DEGes$padjmin <- with(DEGes, pmin(padjtrainingtrainedyoked)) 
    DEGes <- DEGes %>% filter(padjmin < 0.1)
    rownames(DEGes) <- DEGes$rownames
    drop.cols <-colnames(DEGes[,grep("padj|pval|rownames", colnames(DEGes))])
    DEGes <- DEGes %>% dplyr::select(-one_of(drop.cols))
    DEGes <- as.matrix(DEGes)
    DEGes <- DEGes - rowMeans(DEGes)
    DEGes <- as.matrix(DEGes)


    write.csv(DEGes, "../data/02f_DG_DEGs_vsd.csv")

    paletteLength <- 10
    myBreaks <- c(seq(min(DEGes), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(DEGes)/paletteLength, max(DEGes), length.out=floor(paletteLength/2)))
    df <- as.data.frame(colData(DGdds2)[,c("training")]) ## matrix to df
    rownames(df) <- row.names(colData(DGdds2))
    names(df) <- c("training")
    df$subfield <- "DG"
    levels(df$training) <- c("yoked", "trained")

    pheatmap(DEGes, show_colnames=F, show_rownames = T,
             annotation_col=df, 
             annotation_colors = pheatmapcolors2,
             treeheight_row = 25, treeheight_col = 25,
             annotation_row = NA, 
             annotation_legend = TRUE,
             annotation_names_row = FALSE, annotation_names_col = FALSE,
             fontsize = 8, 
             border_color = NA ,
             color = viridis(10),
             cellwidth = 6, 
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation" ,
             clustering_distance_rows="correlation" 
             )

![](../figures/02f_trainedvyoked/pheatmap-DG-1.png)

    pheatmap(DEGes, show_colnames=F, show_rownames = T,
             annotation_col=df, 
             annotation_colors = pheatmapcolors2,
             treeheight_row = 25, treeheight_col = 12.5,
             annotation_row = NA, 
             annotation_legend = TRUE,
             annotation_names_row = FALSE, annotation_names_col = FALSE,
             fontsize = 3, 
             border_color = NA ,
             color = viridis(10),
             cellwidth = 6, 
             cluster_rows = FALSE,
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation", 
             height = 7, 
             width = 2.5,
             filename = "../figures/02f_trainedvyoked/pheatmap-DG.pdf"
             )

    ### CA1

    DEGes <- assay(vst(CA1dds2))
    DEGes <- cbind(DEGes, contrast2)
    DEGes <- as.data.frame(DEGes) # convert matrix to dataframe
    DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe
    DEGes$rownames <- str_to_upper(DEGes$rownames) ## uppercase gene names
    DEGes$padjmin <- with(DEGes, pmin(padjtrainingtrainedyoked)) 
    DEGes <- DEGes %>% filter(padjmin < 0.1)
    rownames(DEGes) <- DEGes$rownames
    drop.cols <-colnames(DEGes[,grep("padj|pval|rownames", colnames(DEGes))])
    DEGes <- DEGes %>% dplyr::select(-one_of(drop.cols))
    DEGes <- as.matrix(DEGes)
    DEGes <- DEGes - rowMeans(DEGes)
    DEGes <- as.matrix(DEGes) 
    paletteLength <- 10
    myBreaks <- c(seq(min(DEGes), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(DEGes)/paletteLength, max(DEGes), length.out=floor(paletteLength/2)))
    df <- as.data.frame(colData(CA1dds2)[,c("training")]) ## matrix to df
    rownames(df) <- row.names(colData(CA1dds2))
    names(df) <- c("training")
    df$subfield <- "CA1"
    levels(df$training) <- c("yoked", "trained")

    pheatmap(DEGes, show_colnames=F, show_rownames = T,
             annotation_col=df, 
             annotation_colors = pheatmapcolors3,
             treeheight_row = 0, treeheight_col = 25,
             annotation_row = NA, 
             annotation_legend = TRUE,
             annotation_names_row = FALSE, annotation_names_col = FALSE,
             fontsize = 8, 
             border_color = NA ,
             color = viridis(10),
             cellwidth = 6, 
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation" 
             )

![](../figures/02f_trainedvyoked/pheatmap-CA1-1.png)

    pheatmap(DEGes, show_colnames=F, show_rownames = T,
             annotation_col=df, 
             annotation_colors = pheatmapcolors3,
             treeheight_row = 25, treeheight_col = 12.5,
             annotation_row = NA, 
             annotation_legend = TRUE,
             annotation_names_row = FALSE, annotation_names_col = FALSE,
             fontsize = 3, 
             border_color = NA ,
             color = viridis(10),
             cellwidth = 6, 
             cluster_rows = FALSE,
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation", 
             height = 2.5, 
             width = 2.5,
             filename = "../figures/02f_trainedvyoked/pheatmap-CA1.pdf"
             )

tSNE
----

    a <- plot.tSNE.trained(DGdds2, 2, "DG")
    b <- plot.tSNE.trained(CA1dds2, 2, "CA1")

    plot_grid(a,b)

![](../figures/02f_trainedvyoked/tSNE-1.png)

    a

![](../figures/02f_trainedvyoked/tSNE-2.png)
