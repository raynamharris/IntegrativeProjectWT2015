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
    knitr::opts_chunk$set(fig.path = '../figures/02f_trainedvyoked/', cache = T)

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
        dplyr::mutate(gene = toupper(gene)) %>% 
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
    ##            gene  lfc     padj direction tissue    comparison
    ## 1 1190002N15RIK 1.64 2.45e-04 increased     DG trained-yoked
    ## 2 A830010M20RIK 1.53 7.89e-07 increased     DG trained-yoked
    ## 3         ABHD2 0.86 1.53e-02 increased     DG trained-yoked
    ## 4          ACAN 1.97 4.54e-09 increased     DG trained-yoked
    ## 5       ADAMTS1 1.88 1.88e-03 increased     DG trained-yoked
    ## 6         ADRB1 0.98 3.11e-02 increased     DG trained-yoked

    CA1DEGs <- listofDEGstrainedvyoked(CA1dds2, "CA1")

    ## [1] "CA1 trained vs yoked"
    ##     gene   lfc    padj direction tissue    comparison
    ## 1 GLCCI1  1.98 0.00144 increased    CA1 trained-yoked
    ## 2  AHDC1 -1.61 0.06502 decreased    CA1 trained-yoked
    ## 3   BMT2 -2.48 0.05663 decreased    CA1 trained-yoked
    ## 4  CTCFL -2.61 0.03633 decreased    CA1 trained-yoked
    ## 5 FN3KRP -1.43 0.04801 decreased    CA1 trained-yoked
    ## 6   GNAZ -3.11 0.01121 decreased    CA1 trained-yoked

    write.csv(DGDEGs, "../data/02f_DG_DEGs.csv", row.names = F)
    write.csv(CA1DEGs, "../data/02f_CA1_DEGs.csv", row.names = F)

    upDG <- DGDEGs %>% filter(lfc > 0)
    upDG <- as.vector(upDG$gene)
    upDG

    ##   [1] "1190002N15RIK" "A830010M20RIK" "ABHD2"         "ACAN"         
    ##   [5] "ADAMTS1"       "ADRB1"         "AHR"           "AMIGO2"       
    ##   [9] "ANKRD13A"      "ANKRD28"       "ANKRD33B"      "APAF1"        
    ##  [13] "ARC"           "ARID5B"        "ARL13B"        "ARL4A"        
    ##  [17] "ARL4D"         "ARL5B"         "ARMCX5"        "ARPP21"       
    ##  [21] "ATF3"          "B3GNT2"        "BACH1"         "BDNF"         
    ##  [25] "BTG2"          "C2CD4B"        "CCNK"          "CIART"        
    ##  [29] "CITED2"        "CLDN12"        "CNNM1"         "CPEB4"        
    ##  [33] "CTNND1"        "CUL3"          "CWC25"         "CXADR"        
    ##  [37] "CYP51"         "DBPHT2"        "DNAJA1"        "DNAJB1"       
    ##  [41] "DNAJB4"        "DUSP14"        "DUSP16"        "DUSP4"        
    ##  [45] "DUSP6"         "DUSP8"         "DYRK2"         "EGR1"         
    ##  [49] "EGR3"          "EGR4"          "EIF5"          "EPRS"         
    ##  [53] "ERF"           "ERRFI1"        "FAM107B"       "FAM118A"      
    ##  [57] "FBXO33"        "FBXW7"         "FERMT2"        "FLRT3"        
    ##  [61] "FOS"           "FOSB"          "FOSL2"         "FOXG1"        
    ##  [65] "FOXO1"         "FRMD6"         "FZD4"          "FZD5"         
    ##  [69] "GAD1"          "GADD45G"       "GM13889"       "GMEB2"        
    ##  [73] "GPR19"         "HECA"          "HMGCR"         "HOMER1"       
    ##  [77] "HS6ST1"        "HSPA1A"        "HSPH1"         "IL16"         
    ##  [81] "ING2"          "IRF2BP2"       "IRS1"          "IRS2"         
    ##  [85] "JDP2"          "JMJD1C"        "JUN"           "JUNB"         
    ##  [89] "JUND"          "KCNA4"         "KCNF1"         "KCNJ2"        
    ##  [93] "KDM6B"         "KDM7A"         "KITL"          "KLF2"         
    ##  [97] "KLF6"          "LBH"           "LCMT2"         "LEMD3"        
    ## [101] "LMNA"          "LONRF1"        "LRRTM2"        "MARCH11"      
    ## [105] "MED7"          "MEST"          "MFAP3L"        "MN1"          
    ## [109] "MYC"           "NAF1"          "NAP1L5"        "NEDD9"        
    ## [113] "NEFM"          "NFIL3"         "NPAS4"         "NPTX2"        
    ## [117] "NR4A1"         "NR4A2"         "NR4A3"         "NUAK1"        
    ## [121] "ODC1"          "OLFML2B"       "OTUD1"         "PAK6"         
    ## [125] "PCDH8"         "PEG10"         "PELI1"         "PER1"         
    ## [129] "PER2"          "PHLDA1"        "PIGA"          "PLAGL1"       
    ## [133] "PLK2"          "PLK3"          "POU3F3"        "PPP1R15A"     
    ## [137] "PRPF38B"       "PTGS2"         "RANBP2"        "RASD1"        
    ## [141] "RASL11A"       "RASL11B"       "RFX2"          "RGMB"         
    ## [145] "RGS2"          "RGS4"          "SCG2"          "SGK1"         
    ## [149] "SH2D3C"        "SIAH2"         "SLC16A1"       "SLC25A25"     
    ## [153] "SLC2A3"        "SLC45A4"       "SLITRK5"       "SMAD7"        
    ## [157] "SNX18"         "SOWAHC"        "SOX9"          "SPTY2D1"      
    ## [161] "SRF"           "STMN4"         "SYT4"          "THBS1"        
    ## [165] "TIPARP"        "TNIP2"         "TRA2B"         "TRIB1"        
    ## [169] "TSC22D2"       "UBC"           "USPL1"         "ZBTB33"       
    ## [173] "ZDBF2"         "ZFAND5"        "ZFP275"        "ZFP654"       
    ## [177] "ZFP869"

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

Go terms
--------

    GO_response <- read.table("../data/goterms/GO_term_summary_20191121_150656.txt", sep = "\t", row.names = NULL,  fill=TRUE)
    GO_response$GO <- "1. response to stimulus"


    GO_translation <- read.table("../data/goterms/GO_term_summary_20191121_141252.txt", sep = "\t", row.names = NULL)
    GO_translation$GO <- "2. translation"

    GO_synapse <- read.table("../data/goterms/GO_term_summary_20191121_145034.txt", sep = "\t", row.names = NULL)
    GO_synapse$GO <- "3. synapse organization"

    GO_learningormemory <- read.table("../data/goterms/GO_term_summary_20191121_142500.txt", sep = "\t", row.names = NULL)
    GO_learningormemory$GO <- "4. learning or memory"

    GOterms <- rbind(GO_learningormemory, GO_response, GO_translation, GO_synapse)

    GOterms <- GOterms %>%
      dplyr::mutate(gene = toupper(MGI.Gene.Marker.ID)) %>% 
      dplyr::select(gene, GO) %>% 
      dplyr::distinct(gene, GO) %>% 
     group_by(gene) 

    GOtermsDEGs <- inner_join(GOterms, DGDEGs) %>%
      arrange(GO,gene) %>%
      dplyr::select(gene, GO) %>% 
     group_by(GO) %>%
     summarize(genes = str_c(gene, collapse = ", "))

    ## Joining, by = "gene"

    GOtermsDEGs

    ## # A tibble: 4 x 2
    ##   GO                 genes                                                 
    ##   <chr>              <chr>                                                 
    ## 1 1. response to st… ABHD2, ADRB1, AHR, ANKRD27, APAF1, ARC, ARID5B, ARL13…
    ## 2 2. translation     CPEB4, EIF5                                           
    ## 3 3. synapse organi… AMIGO2, ARC, BDNF, FLRT3, FZD5, HOMER1, LRRTM2, NPAS4…
    ## 4 4. learning or me… ADRB1, ARC, BDNF, BTG2, EGR1, HMGCR, JUN, NPAS4, NPTX…

    write.csv(GOtermsDEGs, "../data/goterms/GOtermsDEGs2.csv", row.names = F)

    # what to konw which of our favorite genes were in those list but perhaps not captured by analysis


    sanesLichtman <- read.csv("../data/02i_sanesLichtman.csv", header = T, stringsAsFactors = F)
    head(sanesLichtman)

    ##    genes
    ## 1   ACHE
    ## 2  ADCY1
    ## 3 ADRA2A
    ## 4 ADRA2B
    ## 5 ADRA2C
    ## 6  ADRB1

    GOtermsSanesLichtman <- GOterms %>%
      filter(gene %in% sanesLichtman$genes) %>%
      dplyr::select(gene, GO) %>% 
     group_by(GO) %>%
     summarize(genes = str_c(gene, collapse = ", "))

    GOtermsSanesLichtman

    ## # A tibble: 3 x 2
    ##   GO                 genes                                                 
    ##   <chr>              <chr>                                                 
    ## 1 1. response to st… ADCY1, ADRA2A, ADRA2B, ADRA2C, ADRB1, ADRB2, ADRB3, C…
    ## 2 3. synapse organi… ACHE, BDNF, CACNA1A, CACNA1S, CAMK1, CDH1, CDH2, CHRN…
    ## 3 4. learning or me… ADCY1, ADRB1, ADRB2, BDNF, CACNA1C, CACNA1E, CALB1, C…

    write.csv(GOtermsSanesLichtman, "../data/goterms/GOtermsSanesLichtman.csv", row.names = F)
