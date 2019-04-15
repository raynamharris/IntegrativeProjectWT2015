Subfield analysis
-----------------

This script is used to identify treatement differences within each
subfield, generate volcano plots, venn diagrams, and tables for
subsequent GO analyses. The final mutlipanel figures for the manuscript
have been inserted just below the subheadings.

    library(tidyverse)
    library(cowplot) ## for some easy to use themes
    library(VennDiagram) ## venn diagrams
    library(DESeq2) ## for gene expression analysis
    library(UpSetR)
    #devtools::install_github("clauswilke/ggtextures")
    library(ggtextures)
    library(magick)

    ## load functions 
    source("figureoptions.R")
    source("functions_RNAseq.R")

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/02c_rnaseqSubfield/', cache = T)

Get varience stabilized gene expression for each tissue
-------------------------------------------------------

    a.colData <- read.csv("../data/02a_colData.csv", header = T)
    a.countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)

    returndds <- function(mytissue){
      print(mytissue)
      colData <- a.colData %>% 
        filter(Punch %in% c(mytissue))  %>% 
      droplevels()
      
      savecols <- as.character(colData$RNAseqID) 
      savecols <- as.vector(savecols) 
      countData <- a.countData %>% dplyr::select(one_of(savecols)) 

      ## create DESeq object using the factors Punch and APA
      dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ APA2)

      dds <- dds[ rowSums(counts(dds)) > 1, ]  # Pre-filtering genes with 0 counts
       # Differential expression analysis
      return(DESeq(dds))
    }

    DGdds <- returndds("DG") 

    ## [1] "DG"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    CA1dds <- returndds("CA1") 

    ## [1] "CA1"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    CA3dds <- returndds("CA3") 

    ## [1] "CA3"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

Consistent versus yoked-consistent
----------------------------------

    plot.cons.yokcons <- function(mydds, mytissue, mytitle){
      print(mytissue)
      
      res <- results(mydds, contrast =c("APA2", "consistent", "yoked_consistent"),
                     independentFiltering = T, alpha = 0.1)
      print(summary(res))

      data <- data.frame(gene = row.names(res),
                       padj = res$padj, 
                       logpadj = -log10(res$padj),
                       lfc = res$log2FoldChange)
      data <- na.omit(data)
      data <- data %>%
        dplyr::mutate(direction = ifelse(data$lfc > 1 & data$padj < 0.1, 
                            yes = "consistent", 
                            no = ifelse(data$lfc < -1 & data$padj < 0.1, 
                                        yes = "yoked\nconsistent", 
                                        no = "NS")))


      write.csv(data, file = paste0("../data/02c_", mytissue, "_consyokcons.csv", sep = ""))
      
      volcano <- ggplot(data, aes(x = lfc, y = logpadj)) + 
        geom_point(aes(color = factor(direction)), size = 0.5, alpha = 0.75, na.rm = T) + 
        theme_cowplot(font_size = 7, line_size = 0.25) +
        geom_hline(yintercept = 1,  size = 0.25, linetype = 2) + 
        scale_color_manual(values = volcano1,
                          name = "higher in",
                          breaks = c("yoked\nconsistent", "NS", "consistent"))  + 
        scale_y_continuous(limits=c(0, 12.5)) +
        scale_x_continuous(limits=c(-10, 10),
                            name="Log fold difference")+
        ylab(paste0("log10 p-value")) +  
        labs(subtitle = mytitle) +
        theme(legend.position = "bottom",
              legend.spacing.x = unit(0.1, 'cm'),
              #legend.text=element_text(size=4),
              legend.title = element_text(size=6),
              legend.key.size = unit(0.2, "cm"),
              legend.margin=margin(t=-0.1, r=0, b=0, l=-0.1, unit="cm")) 
      plot(volcano)
    }


    DGconsyokcons <-  plot.cons.yokcons(DGdds, "DG", "DG train") + theme(legend.position = "none")  

    ## [1] "DG"
    ## 
    ## out of 17011 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 119, 0.7%
    ## LFC < 0 (down)     : 6, 0.035%
    ## outliers [1]       : 20, 0.12%
    ## low counts [2]     : 4608, 27%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/02c_rnaseqSubfield/consyokcons-1.png)

    CA3consyokcons <-  plot.cons.yokcons(CA3dds, "CA3", "CA3 train")  

    ## [1] "CA3"
    ## 
    ## out of 16502 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.0061%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 11, 0.067%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/02c_rnaseqSubfield/consyokcons-2.png)

    CA1consyokcons <-  plot.cons.yokcons(CA1dds, "CA1", "CA1 train")  + theme(legend.position = c(0.6,0.9))  

    ## [1] "CA1"
    ## 
    ## out of 16852 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 522, 3.1%
    ## LFC < 0 (down)     : 360, 2.1%
    ## outliers [1]       : 32, 0.19%
    ## low counts [2]     : 4892, 29%
    ## (mean count < 5)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/02c_rnaseqSubfield/consyokcons-3.png)

    training <- plot_grid(DGconsyokcons,CA1consyokcons, nrow = 1,
                          labels = "AUTO",
                          label_size = 7,
                          rel_widths = c(0.5,0.5))
    training

![](../figures/02c_rnaseqSubfield/consyokcons-4.png)

Confict versus Consistent
-------------------------

    plot.conf.cons <- function(mydds, mytissue){
      
      print(mytissue)
      
      res <- results(mydds, contrast =c("APA2", "conflict", "consistent"),
                     independentFiltering = T, alpha = 0.1)
      print(summary(res))

      data <- data.frame(gene = row.names(res),
                       padj = res$padj, 
                       logpadj = -log10(res$padj),
                       lfc = res$log2FoldChange)
      data <- na.omit(data)
      
      data <- data %>%
        dplyr::mutate(direction = ifelse(data$lfc > 1 & data$padj < 0.1, 
                            yes = "conflict", 
                            no = ifelse(data$lfc < -1 & data$padj < 0.1, 
                                        yes = "consistent", 
                                        no = "NS")))
      
      write.csv(data, file = paste0("../data/02c_", mytissue, "_confcons.csv", sep = ""))
      
      volcano <- ggplot(data, aes(x = lfc, y = logpadj)) + 
        geom_point(aes(color = factor(direction)), size = 0.5, alpha = 0.75, na.rm = T) + 
        theme_cowplot(font_size = 7, line_size = 0.25) +
        geom_hline(yintercept = 1,  size = 0.25, linetype = 2) + 
        scale_color_manual(values = volcano2,
                           breaks = c("consistent", "NS", "conflict"),
                          name = "higher in")  + 
        scale_y_continuous(limits=c(0, 12.5)) +
        scale_x_continuous(limits=c(-10, 10),
                            name="Log fold difference")+
        ylab(paste0("log10 p-value")) +  
        labs(subtitle = mytissue) +
        theme(panel.grid.minor=element_blank(),
              legend.position = "bottom",
              legend.spacing.x = unit(-0.1, 'cm'),
              panel.grid.major=element_blank(),
              legend.margin=margin(t=-0.25, r=0, b=0, l=0, unit="cm")) 
      plot(volcano)
    }


    DGconflict <-  plot.conf.cons(DGdds, "DG")

    ## [1] "DG"
    ## 
    ## out of 17011 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 1, 0.0059%
    ## outliers [1]       : 20, 0.12%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/02c_rnaseqSubfield/confcons-1.png)

    CA3conflict <-  plot.conf.cons(CA3dds, "CA3")

    ## [1] "CA3"
    ## 
    ## out of 16502 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 11, 0.067%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/02c_rnaseqSubfield/confcons-2.png)

    CA1conflict <-  plot.conf.cons(CA1dds, "CA1")

    ## [1] "CA1"
    ## 
    ## out of 16852 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 32, 0.19%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/02c_rnaseqSubfield/confcons-3.png)

    plot_grid(DGconflict, CA3conflict, CA1conflict, nrow = 1)

![](../figures/02c_rnaseqSubfield/confcons-4.png)

Yoked confict versus yoked consistent
-------------------------------------

    plot.yokconf.yokcons <- function(mydds, mytissue, mytitle){
      
      print(mytissue)
      
      res <- results(mydds, contrast =c("APA2", "yoked_conflict", "yoked_consistent"),
                     independentFiltering = T, alpha = 0.1)
      print(summary(res))

      data <- data.frame(gene = row.names(res),
                       padj = res$padj, 
                       logpadj = -log10(res$padj),
                       lfc = res$log2FoldChange)
      data <- na.omit(data)
      data <- data %>%
        dplyr::mutate(direction = ifelse(data$lfc > 1 & data$padj < 0.1, 
                            yes = "yoked\nconsistent", 
                            no = ifelse(data$lfc < -1 & data$padj < 0.1, 
                                        yes = "yoked\nconflict", 
                                        no = "NS")))
      
      data$direction <- factor(data$direction, levels = c("yoked\nconsistent", "NS", "yoked\nconflict"))
      
      write.csv(data, file = paste0("../data/02c_", mytissue, "_yokeconfyokcons.csv", sep = ""))
      
      volcano <- ggplot(data, aes(x = lfc, y = logpadj, color = direction)) + 
        geom_point(size = 0.5, alpha = 0.75, na.rm = T) + 
        theme_cowplot(font_size = 7, line_size = 0.25) +
        geom_hline(yintercept = 1,  size = 0.25, linetype = 2) + 
        scale_color_manual(values = volcano3,
                           name = "higher in")  + 
        scale_y_continuous(limits=c(0, 12.5)) +
        scale_x_continuous(limits=c(-10, 10),
                            name="Log fold difference")+
        ylab(paste0("log10 p-value")) +  
        labs(subtitle = mytitle) +
        theme(legend.position = "bottom",
              legend.spacing.x = unit(0.1, 'cm'),
              #legend.text=element_text(size=4),
              legend.key.size = unit(0.2, "cm"),
              legend.margin=margin(t=-0.1, r=0, b=0, l=-0.1, unit="cm")) 
      plot(volcano)  
    }

    DGyoked <-  plot.yokconf.yokcons(DGdds, "DG", "DG stress")

    ## [1] "DG"
    ## 
    ## out of 17011 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 3, 0.018%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 20, 0.12%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/02c_rnaseqSubfield/yokeconfyokcons-1.png)

    CA3yoked <-  plot.yokconf.yokcons(CA3dds, "CA3", "CA3 stress")

    ## [1] "CA3"
    ## 
    ## out of 16502 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.0061%
    ## LFC < 0 (down)     : 1, 0.0061%
    ## outliers [1]       : 11, 0.067%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/02c_rnaseqSubfield/yokeconfyokcons-2.png)

    CA1yoked <-  plot.yokconf.yokcons(CA1dds, "CA1", "CA1 stress")  

    ## [1] "CA1"
    ## 
    ## out of 16852 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 545, 3.2%
    ## LFC < 0 (down)     : 372, 2.2%
    ## outliers [1]       : 32, 0.19%
    ## low counts [2]     : 4892, 29%
    ## (mean count < 5)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

![](../figures/02c_rnaseqSubfield/yokeconfyokcons-3.png)

    plot_grid(DGyoked, CA3yoked, CA1yoked, nrow = 1)

![](../figures/02c_rnaseqSubfield/yokeconfyokcons-4.png)

pkmz
====

    plotCounts(DGdds, "Prkcz", intgroup = "APA2", normalized = TRUE, main="Prkcz in DG")

![](../figures/02c_rnaseqSubfield/pkmz-1.png)

    plotCounts(CA3dds, "Prkcz", intgroup = "APA2", normalized = TRUE, main="Prkcz in CA3")

![](../figures/02c_rnaseqSubfield/pkmz-2.png)

    plotCounts(CA1dds, "Prkcz", intgroup = "APA2", normalized = TRUE, main="Prkcz in CA1")

![](../figures/02c_rnaseqSubfield/pkmz-3.png)

Upset plots
-----------

What genes overlap within cetain comparisons?

    a.colData <- read.csv("../data/02a_colData.csv", header = T)
    a.countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)

    eachsubfield <- levels(a.colData$Punch)

    listofDEGs <- function(group1, group2){
      res <- results(dds, contrast = c("APA2", group1, group2), independentFiltering = T)
      data <- data.frame(gene = row.names(res),
                         lfc = res$log2FoldChange,
                         padj = res$padj,
                         tissue = i,
                         comparison = paste(group1, group2, sep = "-"))
      data <- data %>% dplyr::filter(padj < 0.1) %>% droplevels()
      return(data)
    }

    for(i in eachsubfield){
      
      colData <- a.colData %>% 
        dplyr::filter(Punch == i)  %>%
        droplevels()
      print(i)
      
      savecols <- as.character(colData$RNAseqID) 
      savecols <- as.vector(savecols) 
      countData <- a.countData %>% dplyr::select(one_of(savecols)) 

    ## create DESeq object using the factors Punch and APA
    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ APA2)

    dds # view the DESeq object - note numnber of genes
    dds <- dds[ rowSums(counts(dds)) > 1, ]  # Pre-filtering genes with 0 counts
    dds <- DESeq(dds) # Differential expression analysis

    A <- listofDEGs("consistent","yoked_consistent")
    B <- listofDEGs("conflict","consistent")
    C <- listofDEGs("conflict","yoked_conflict")
    D <- listofDEGs("yoked_conflict","yoked_consistent")


    all <- rbind(A,B,C,D)

    write.csv(all, file = paste("../data/02c_",i,"forupset.csv", sep = ""), row.names = F)
    }

    ## [1] "CA1"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## [1] "CA3"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## [1] "DG"

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    DG <- read.csv("../data/02c_DGforupset.csv")  
    CA1<- read.csv("../data/02c_CA1forupset.csv")  
    CA3 <- read.csv("../data/02c_CA3forupset.csv") 


    # upset plot without direction
    all <- rbind(DG,CA1,CA3)

    levels(all$comparison) <- c("conf-cons",
                                "conf-yconf", 
                                "cons-ycons", 
                                "yconf-ycons")

    all$significant <- paste(all$tissue, all$comparison, sep = "-")

    myupsetdf <- all %>%
      select(gene,significant) %>%
      mutate(yesno = 1) %>%
      distinct %>%
      spread(significant, yesno, fill = 0)
    head(myupsetdf)

    ##            gene CA1-conf-yconf CA1-cons-ycons CA1-yconf-ycons
    ## 1 A830010M20Rik              0              0               0
    ## 2          Acan              0              0               0
    ## 3       Adamts1              0              1               1
    ## 4        Amigo2              0              0               0
    ## 5       Ankrd28              0              0               0
    ## 6      Ankrd33b              0              0               0
    ##   CA3-cons-ycons CA3-yconf-ycons DG-conf-cons DG-conf-yconf DG-cons-ycons
    ## 1              0               0            0             0             1
    ## 2              0               0            0             1             1
    ## 3              0               0            0             0             1
    ## 4              0               0            0             0             1
    ## 5              0               0            0             0             1
    ## 6              0               0            0             1             0
    ##   DG-yconf-ycons
    ## 1              0
    ## 2              0
    ## 3              0
    ## 4              0
    ## 5              0
    ## 6              0

    write.csv(myupsetdf, "../data/02c_upsetdf.csv")


    # upset plot with direction, only CA1learn, CA1 stress, and DGlearn
    head(all) 

    ##            gene      lfc         padj tissue comparison   significant
    ## 1 A830010M20Rik 1.991294 6.235077e-05     DG cons-ycons DG-cons-ycons
    ## 2          Acan 2.459128 1.411053e-05     DG cons-ycons DG-cons-ycons
    ## 3       Adamts1 2.370955 3.493214e-02     DG cons-ycons DG-cons-ycons
    ## 4        Amigo2 1.982597 6.941977e-02     DG cons-ycons DG-cons-ycons
    ## 5       Ankrd28 1.642931 1.247682e-02     DG cons-ycons DG-cons-ycons
    ## 6           Arc 2.891323 3.577881e-05     DG cons-ycons DG-cons-ycons

    all$direction <- ifelse(all$lfc > 0, "up", "down")
    all$sigdir <- paste(all$significant, all$direction, sep = "-")

    myupsetslim  <- all %>%
      filter(significant %in% c("DG-cons-ycons", "CA1-cons-ycons", "CA1-yconf-ycons")) %>%
      select(gene,sigdir) %>%
      mutate(yesno = 1) %>%
      distinct %>%
      spread(sigdir, yesno, fill = 0)
    head(myupsetslim)

    ##            gene CA1-cons-ycons-down CA1-cons-ycons-up CA1-yconf-ycons-down
    ## 1 A830010M20Rik                   0                 0                    0
    ## 2          Acan                   0                 0                    0
    ## 3       Adamts1                   0                 1                    0
    ## 4        Amigo2                   0                 0                    0
    ## 5       Ankrd28                   0                 0                    0
    ## 6           Arc                   0                 0                    0
    ##   CA1-yconf-ycons-up DG-cons-ycons-down DG-cons-ycons-up
    ## 1                  0                  0                1
    ## 2                  0                  0                1
    ## 3                  1                  0                1
    ## 4                  0                  0                1
    ## 5                  0                  0                1
    ## 6                  0                  0                1

    row.names(myupsetslim) <- myupsetslim$gene
    myupsetslim$gene <- NULL
    colSums(myupsetslim)

    ##  CA1-cons-ycons-down    CA1-cons-ycons-up CA1-yconf-ycons-down 
    ##                  360                  522                  372 
    ##   CA1-yconf-ycons-up   DG-cons-ycons-down     DG-cons-ycons-up 
    ##                  545                    6                  119

    myupsetslim$gene <- row.names(myupsetslim)

    upset(myupsetdf, keep.order = T)

![](../figures/02c_rnaseqSubfield/upsetplot-1.png)

    trained <- myupsetdf %>%
      select(gene, 'DG-cons-ycons', 'CA3-cons-ycons', 'CA1-cons-ycons')

    colnames(trained) <- c("gene", "DG", "CA3", "CA1")


    pdf(file="../figures/02c_rnaseqSubfield/upsettraining.pdf",  onefile=FALSE) # or other device
    upset(trained, keep.order = T,
          sets = c("CA1", "CA3", "DG"),
          sets.bar.color=c("#7570b3","#1b9e77", "#d95f02"),
          queries = list(list(query = intersects, params = list("CA1"), color = "#ca0020", active = T),
                         list(query = intersects, params = list("DG"), color = "#ca0020", active = T),
                         list(query = intersects, params = list("CA3"), color = "#ca0020", active = T),
                         list(query = intersects, params = list("CA1", "DG"), 
                              color = "#ca0020", active = T)),
          text.scale = 2,
          sets.x.label = NULL,
          #point.size = 1, line.size = 1,
          mb.ratio = c(0.6, 0.4)
    )
    dev.off()

    ## quartz_off_screen 
    ##                 2

    shared <- myupsetdf %>%
      select(gene, "CA1-yconf-ycons", "CA1-cons-ycons", "DG-cons-ycons")
    colnames(shared) <- c("gene", "CA1stress", "CA1learn", "DGlearn")

    # ca1 DG learning
    shared %>% filter(CA1learn == 1 & DGlearn == 1 & CA1stress == 0) %>%
      select(gene)  

    ##      gene
    ## 1    Fosb
    ## 2 Gm13889
    ## 3    Irs1
    ## 4   Lemd3
    ## 5   Npas4
    ## 6    Rgs2

    # CA1 and DG learning only (none in stress)
    shared %>% filter(CA1learn == 1 & DGlearn == 1 & CA1stress == 0) %>%
      select(gene) 

    ##      gene
    ## 1    Fosb
    ## 2 Gm13889
    ## 3    Irs1
    ## 4   Lemd3
    ## 5   Npas4
    ## 6    Rgs2

    # learning only CA1
    shared %>% filter(CA1learn == 1 & CA1stress == 0 & DGlearn == 0) %>%
      select(gene)

    ##              gene
    ## 1           Kcnc2
    ## 2   1810022K09Rik
    ## 3   3110035E14Rik
    ## 4   9430015G10Rik
    ## 5           Aagab
    ## 6            Abl1
    ## 7             Abr
    ## 8           Acad8
    ## 9           Acads
    ## 10          Acot9
    ## 11          Acta1
    ## 12         Actl6b
    ## 13          Adam8
    ## 14          Adat1
    ## 15         Afg3l2
    ## 16          Aftph
    ## 17          Agfg2
    ## 18          Ahdc1
    ## 19            Aip
    ## 20          Akap9
    ## 21          Aldoa
    ## 22          Alg11
    ## 23          Alms1
    ## 24         Amigo1
    ## 25         Amigo3
    ## 26         Ankfy1
    ## 27        Ankrd12
    ## 28          Ap1s1
    ## 29           Apln
    ## 30          Apol8
    ## 31           Aptx
    ## 32           Arf5
    ## 33        Arfgef3
    ## 34         Arglu1
    ## 35         Arl13b
    ## 36          Arl5a
    ## 37         Arntl2
    ## 38         Arrdc3
    ## 39          Asap1
    ## 40          Asxl3
    ## 41           Atf5
    ## 42          Atf6b
    ## 43         Atp2b2
    ## 44          Atp5e
    ## 45       Atp6v1g2
    ## 46         Atpaf2
    ## 47        Atxn7l3
    ## 48          Baz2a
    ## 49          Baz2b
    ## 50       BC030500
    ## 51        Bloc1s5
    ## 52           Bop1
    ## 53         Btbd10
    ## 54          C2cd3
    ## 55          Cabp1
    ## 56         Camkk2
    ## 57        Camsap1
    ## 58         Capn10
    ## 59           Cbx3
    ## 60           Cbx4
    ## 61        Ccdc190
    ## 62         Ccdc43
    ## 63           Ccny
    ## 64           Ccr5
    ## 65           Cct7
    ## 66          Cdh20
    ## 67           Cdh7
    ## 68          Cdk15
    ## 69          Cend1
    ## 70         Cep131
    ## 71         Cfap20
    ## 72           Chpf
    ## 73           Chrd
    ## 74           Chuk
    ## 75          Clock
    ## 76          Cntrl
    ## 77       Colgalt2
    ## 78           Coq2
    ## 79            Cpq
    ## 80           Crem
    ## 81          Crim1
    ## 82          Cstf1
    ## 83          Ctcfl
    ## 84          Cyb5b
    ## 85     D1Ertd622e
    ## 86            Dbt
    ## 87         Ddx19b
    ## 88          Ddx31
    ## 89        Dennd1b
    ## 90         Dnajb6
    ## 91           Dner
    ## 92           Dok6
    ## 93           Dpf1
    ## 94           Drp2
    ## 95         Dusp19
    ## 96         Eif2s1
    ## 97         Eif4a1
    ## 98         Eif4a2
    ## 99           Emc4
    ## 100          Eml2
    ## 101         Enkd1
    ## 102          Eri3
    ## 103         Esrrg
    ## 104         Esyt1
    ## 105        Exoc6b
    ## 106          Ext1
    ## 107       Fam110b
    ## 108       Fam131a
    ## 109       Fam19a2
    ## 110        Fam63a
    ## 111        Fam63b
    ## 112         Fanci
    ## 113          Fbn1
    ## 114         Fbxl6
    ## 115          Fgd4
    ## 116          Fgd5
    ## 117        Filip1
    ## 118          Flt1
    ## 119          Fmn1
    ## 120        Fn3krp
    ## 121        Fndc3a
    ## 122         Foxj3
    ## 123        Frmpd3
    ## 124           Fry
    ## 125          Fzd3
    ## 126          Get4
    ## 127        Glcci1
    ## 128       Gm10053
    ## 129       Gm21887
    ## 130        Gm4631
    ## 131         Gm527
    ## 132        Gm9821
    ## 133        Gpr161
    ## 134        Gpr180
    ## 135         Grem2
    ## 136         Grik3
    ## 137          Grm1
    ## 138         Gstt3
    ## 139       Gucy1a2
    ## 140        Gucy2e
    ## 141          Guk1
    ## 142         Hcfc1
    ## 143         Hdac6
    ## 144         Hecw2
    ## 145          Helz
    ## 146         Herc3
    ## 147       Herpud2
    ## 148          Hexb
    ## 149         Hmox2
    ## 150        Hspa1b
    ## 151         Hspa2
    ## 152        Hspbp1
    ## 153         Igbp1
    ## 154       Igf2bp2
    ## 155         Impa1
    ## 156         Inhbb
    ## 157         Ino80
    ## 158        Inpp4b
    ## 159          Jag2
    ## 160          Jak2
    ## 161         Josd1
    ## 162         Kat6a
    ## 163        Kbtbd3
    ## 164         Kcna1
    ## 165         Kcna4
    ## 166         Kcnc1
    ## 167         Kcnd3
    ## 168         Kctd5
    ## 169         Kdm5d
    ## 170         Khnyn
    ## 171         Kif1a
    ## 172        Kif26b
    ## 173         Klkb1
    ## 174         Kmt2e
    ## 175         Kmt5a
    ## 176          Kri1
    ## 177        Lefty1
    ## 178         Letm2
    ## 179         Lnpep
    ## 180        Lpcat4
    ## 181         Lrfn4
    ## 182        Lrrc58
    ## 183         Lsamp
    ## 184          Ltv1
    ## 185         Lzts3
    ## 186        Man1a2
    ## 187         Manba
    ## 188       Map3k12
    ## 189      Mapkapk2
    ## 190         Marf1
    ## 191         Mark3
    ## 192          Mbd5
    ## 193           Mcc
    ## 194          Med8
    ## 195       Mettl16
    ## 196      Mettl21e
    ## 197       Mfsd13a
    ## 198        Mmadhc
    ## 199      Mphosph9
    ## 200        Mpped1
    ## 201        Mrgpre
    ## 202        Mrpl28
    ## 203        Mrpl48
    ## 204        Msl3l2
    ## 205         Mtfmt
    ## 206         Myo5a
    ## 207          Nat9
    ## 208         Ncoa1
    ## 209        Ndufs7
    ## 210          Nefm
    ## 211         Neto2
    ## 212           Nf2
    ## 213          Ngef
    ## 214          Nkrf
    ## 215          Nle1
    ## 216          Nme1
    ## 217         Nolc1
    ## 218           Nov
    ## 219          Nrgn
    ## 220         Ntng1
    ## 221         Ntpcr
    ## 222        Nudt19
    ## 223         Nudt6
    ## 224          Nxt2
    ## 225         Olfm3
    ## 226         Ovca2
    ## 227         Patz1
    ## 228        Pcdh17
    ## 229         Pde6a
    ## 230         Pebp1
    ## 231          Pgk1
    ## 232         Phka1
    ## 233         Pias4
    ## 234        Plagl2
    ## 235          Plau
    ## 236        Podxl2
    ## 237        Polr2h
    ## 238        Polrmt
    ## 239       Ppp1r3f
    ## 240        Ppp3r1
    ## 241         Ppp5c
    ## 242         Prpf8
    ## 243         Psg28
    ## 244         Psmc4
    ## 245         Psmg2
    ## 246        Ptpn11
    ## 247         Ptpn4
    ## 248         Ptprm
    ## 249        Ptprn2
    ## 250         Pygo1
    ## 251         Rab3c
    ## 252          Rbak
    ## 253         Rbbp4
    ## 254         Rcan2
    ## 255          Rdh1
    ## 256          Rest
    ## 257         Rfesd
    ## 258         Rfwd3
    ## 259        Rhbdl3
    ## 260       Rnaseh1
    ## 261        Rnf165
    ## 262        Rnf180
    ## 263        Rnf216
    ## 264         Rnf25
    ## 265         Rpl10
    ## 266         Rpl36
    ## 267         Rplp2
    ## 268          Rps6
    ## 269        Rsc1a1
    ## 270        Rsph3a
    ## 271        S100a1
    ## 272        Samhd1
    ## 273        Sema4a
    ## 274         Sf3b2
    ## 275         Sgms2
    ## 276         Shoc2
    ## 277       Slc24a2
    ## 278      Slc25a38
    ## 279      Slc25a46
    ## 280       Slc26a2
    ## 281       Slc2a13
    ## 282       Slc35b4
    ## 283       Slc44a1
    ## 284        Slc4a3
    ## 285        Slc8a1
    ## 286         Slx1b
    ## 287         Smek2
    ## 288         Smim3
    ## 289        Snap29
    ## 290          Snx2
    ## 291         Snx24
    ## 292        Sorcs1
    ## 293         Sox10
    ## 294          Sox5
    ## 295         Spast
    ## 296       Specc1l
    ## 297        Sptbn1
    ## 298         Srek1
    ## 299           Srr
    ## 300          Srrd
    ## 301    St6galnac4
    ## 302       St8sia4
    ## 303       Stard13
    ## 304         Stk25
    ## 305          Stk3
    ## 306         Stox2
    ## 307        Strip1
    ## 308          Stx3
    ## 309       Tbc1d30
    ## 310        Tbc1d9
    ## 311         Tdrd7
    ## 312          Tet3
    ## 313          Tfrc
    ## 314          Thra
    ## 315       Tmem50b
    ## 316        Tmem57
    ## 317        Tmem65
    ## 318       Tmem88b
    ## 319        Tmem8b
    ## 320          Tns2
    ## 321          Tox2
    ## 322         Traf6
    ## 323           Ttn
    ## 324        Tuba4a
    ## 325        Tvp23a
    ## 326         Uckl1
    ## 327         Uqcrh
    ## 328         Usmg5
    ## 329         Usp45
    ## 330        Usp6nl
    ## 331         Xrcc3
    ## 332         Xrcc6
    ## 333         Ylpm1
    ## 334           Zak
    ## 335       Zc3h12a
    ## 336        Zc3h13
    ## 337        Zfp114
    ## 338        Zfp395
    ## 339        Zfp414
    ## 340        Zfp446
    ## 341        Zfp580
    ## 342        Zfp617
    ## 343        Zfp711
    ## 344        Zfp738
    ## 345        Zfp821
    ## 346        Zfp831
    ## 347       Zfyve16
    ## 348         Zmym4
    ## 349         Znrd1

    # stress only CA1
    shared %>% filter(CA1learn == 0 & CA1stress == 1 & DGlearn == 0) %>%
      select(gene)

    ##              gene
    ## 1   1110032F04Rik
    ## 2   1600002K03Rik
    ## 3   2010107G23Rik
    ## 4   2210013O21Rik
    ## 5   2210016L21Rik
    ## 6   2310009B15Rik
    ## 7   5730409E04Rik
    ## 8           Abcc4
    ## 9           Abhd4
    ## 10           Ache
    ## 11          Adcy5
    ## 12          Adcy6
    ## 13         Adgrb3
    ## 14          Adpgk
    ## 15         Afg3l1
    ## 16           Alg9
    ## 17            Ank
    ## 18        Ankrd40
    ## 19         Apcdd1
    ## 20        Arfgap3
    ## 21        Arhgap4
    ## 22      Arhgef10l
    ## 23         Arid3b
    ## 24          Armt1
    ## 25          Arvcf
    ## 26           Asph
    ## 27          Astn1
    ## 28          Atg4b
    ## 29         Atp1a2
    ## 30         Atp2a2
    ## 31          Atp5d
    ## 32            Axl
    ## 33        B4galt7
    ## 34          Bbof1
    ## 35           Bbs9
    ## 36           Bcl6
    ## 37         Bhlhb9
    ## 38           Bin3
    ## 39          Bnip2
    ## 40           Bod1
    ## 41            Bok
    ## 42           Braf
    ## 43  C130074G19Rik
    ## 44          C1ql3
    ## 45        C1qtnf6
    ## 46         C77370
    ## 47         Cacfd1
    ## 48        Cacna1c
    ## 49         Cacnb1
    ## 50         Camk1g
    ## 51          Carm1
    ## 52        Carnmt1
    ## 53         Ccdc53
    ## 54         Ccdc59
    ## 55          Ccng1
    ## 56           Ccr2
    ## 57           Cd63
    ## 58         Cdadc1
    ## 59         Cdc123
    ## 60       Cdc42bpg
    ## 61           Cdk8
    ## 62         Cdkal1
    ## 63         Cep135
    ## 64          Cers1
    ## 65          Ces2b
    ## 66            Cfp
    ## 67           Chd6
    ## 68          Chmp7
    ## 69          Ciao1
    ## 70         Cirh1a
    ## 71          Clcn7
    ## 72        Clec11a
    ## 73          Cnbd2
    ## 74         Cntrob
    ## 75         Commd9
    ## 76          Crocc
    ## 77          Csmd3
    ## 78           Ctsh
    ## 79         Cuedc1
    ## 80          Cxcl9
    ## 81         Cyb561
    ## 82       Cyb561d2
    ## 83         Cyp4v3
    ## 84          Daglb
    ## 85        Dclre1b
    ## 86          Ddx24
    ## 87          Ddx47
    ## 88          Ddx51
    ## 89          Decr2
    ## 90           Det1
    ## 91         Dhtkd1
    ## 92       Dnase1l2
    ## 93           Dnlz
    ## 94          Dock8
    ## 95        Dpy19l3
    ## 96          Dscr3
    ## 97           Dus2
    ## 98          Dus3l
    ## 99       Dync1li1
    ## 100 E130309D02Rik
    ## 101          Ece1
    ## 102          Edc3
    ## 103        Efcab6
    ## 104        Elmod1
    ## 105          Elp2
    ## 106          Eml4
    ## 107          Eno4
    ## 108        Entpd2
    ## 109          Eny2
    ## 110         Ep300
    ## 111          Eps8
    ## 112          Epyc
    ## 113        Ero1lb
    ## 114         Ethe1
    ## 115         Extl2
    ## 116            F3
    ## 117       Fam133b
    ## 118       Fam134a
    ## 119       Fam219b
    ## 120        Fam60a
    ## 121         Farp1
    ## 122       Fastkd2
    ## 123         Fbln1
    ## 124         Fbxl4
    ## 125        Fbxo25
    ## 126         Fcgr3
    ## 127          Fgd1
    ## 128          Fgd3
    ## 129        Fkbp14
    ## 130         Fkbpl
    ## 131         Flrt3
    ## 132         Fstl4
    ## 133       Galnt10
    ## 134         Garem
    ## 135          Gcc2
    ## 136         Gdpd2
    ## 137        Gemin4
    ## 138         Gfpt2
    ## 139        Glyctk
    ## 140       Gm10146
    ## 141       Gm20715
    ## 142       Gm38393
    ## 143       Gm43796
    ## 144        Gm6741
    ## 145        Gm9803
    ## 146         Gmpr2
    ## 147         Gosr2
    ## 148        Gpank1
    ## 149         Gpm6a
    ## 150        Gpr158
    ## 151         Grid1
    ## 152        Grin2c
    ## 153          Guf1
    ## 154          Gys1
    ## 155        H2-Ke6
    ## 156          Heg1
    ## 157         Hmgcl
    ## 158         Hmgn2
    ## 159        Homer3
    ## 160        Hrasls
    ## 161         Htr1a
    ## 162           Id4
    ## 163            Ik
    ## 164        Ikbkap
    ## 165          Ipo5
    ## 166      Irak1bp1
    ## 167       Irf2bp1
    ## 168         Jade1
    ## 169         Jmjd6
    ## 170         Jmjd8
    ## 171         Kat2b
    ## 172        Katnb1
    ## 173        Kcnj12
    ## 174         Kcnn1
    ## 175         Kcnu1
    ## 176        Kctd15
    ## 177        Kdelr2
    ## 178       Laptm4a
    ## 179          Ldah
    ## 180         Lnpk1
    ## 181          Lrp5
    ## 182         Lypd1
    ## 183          Mafk
    ## 184        Map3k2
    ## 185        March9
    ## 186          Mcl1
    ## 187          Mdp1
    ## 188         Mgrn1
    ## 189        Mif4gd
    ## 190        Mrpl16
    ## 191        Mrpl53
    ## 192          Msl2
    ## 193         Mtmr7
    ## 194         Mtmr9
    ## 195        Mtrf1l
    ## 196          Mtx1
    ## 197           Mut
    ## 198         Myo10
    ## 199         Naa20
    ## 200         Naglu
    ## 201         Nalcn
    ## 202        Nanos1
    ## 203          Nans
    ## 204         Ndst2
    ## 205         Nedd9
    ## 206         Nell1
    ## 207        Nfkbia
    ## 208          Npc1
    ## 209        Nploc4
    ## 210       Nr2c2ap
    ## 211          Nrf1
    ## 212        Nudt11
    ## 213         Nup85
    ## 214         Nxph3
    ## 215         Oas1b
    ## 216          Ogg1
    ## 217          Optn
    ## 218          Orc6
    ## 219          P3h4
    ## 220         Padi2
    ## 221       Pcdhb16
    ## 222        Pcdhb6
    ## 223        Pcp4l1
    ## 224          Pdcl
    ## 225         Pdia4
    ## 226          Pdpn
    ## 227          Pdpr
    ## 228         Pds5a
    ## 229          Pecr
    ## 230        Pgrmc1
    ## 231        Pla2g7
    ## 232         Plaur
    ## 233       Plekha1
    ## 234         Plin2
    ## 235        Pnpla8
    ## 236        Polr2d
    ## 237       Polr3gl
    ## 238         Ppdpf
    ## 239        Ppp1cb
    ## 240      Ppp1r13l
    ## 241        Ppp6r2
    ## 242          Ppt2
    ## 243         Pqlc2
    ## 244          Prcp
    ## 245         Prkdc
    ## 246       Prkrip1
    ## 247         Prmt2
    ## 248         Prmt3
    ## 249         Prpf3
    ## 250         Prpf4
    ## 251         Psmc5
    ## 252          Pwp2
    ## 253         Rab13
    ## 254         Rab31
    ## 255         Rab43
    ## 256       Rabgap1
    ## 257         Rabif
    ## 258         Ramp2
    ## 259        Rangrf
    ## 260         Rap2b
    ## 261         Rccd1
    ## 262          Rem2
    ## 263          Rfc1
    ## 264         Rftn2
    ## 265         Rgs16
    ## 266          Rgs6
    ## 267          Rin2
    ## 268        Rnf166
    ## 269        Rnf207
    ## 270          Rnls
    ## 271         Rnpc3
    ## 272 RP23-220F20.2
    ## 273         Rpap3
    ## 274          Rrp1
    ## 275          Rtn4
    ## 276         Rufy3
    ## 277        Sacm1l
    ## 278         Sars2
    ## 279          Sbsn
    ## 280         Sdad1
    ## 281       Sec23ip
    ## 282         Senp2
    ## 283         Senp8
    ## 284         Serf1
    ## 285         Sesn3
    ## 286           Sf1
    ## 287        Sfmbt1
    ## 288           Sfn
    ## 289         Sfxn5
    ## 290        Shisa5
    ## 291       Shroom4
    ## 292         Shtn1
    ## 293        Slain1
    ## 294       Slc12a2
    ## 295       Slc24a3
    ## 296      Slc25a35
    ## 297       Slc38a9
    ## 298        Slc4a2
    ## 299       Slitrk1
    ## 300         Smad5
    ## 301         Smdt1
    ## 302          Smg6
    ## 303        Smim17
    ## 304        Snrpd2
    ## 305          Snx8
    ## 306         Soat1
    ## 307         Socs4
    ## 308         Spag7
    ## 309         Spag9
    ## 310          Spi1
    ## 311         Ssna1
    ## 312          Ssr1
    ## 313    St6galnac3
    ## 314         Sugp1
    ## 315         Supt3
    ## 316        Swsap1
    ## 317       Syndig1
    ## 318         Synpo
    ## 319         Syvn1
    ## 320           Tbp
    ## 321        Tceal3
    ## 322          Tcta
    ## 323         Tdrd3
    ## 324         Thap4
    ## 325       Timm10b
    ## 326        Timm44
    ## 327         Tmco3
    ## 328       Tmem119
    ## 329       Tmem121
    ## 330      Tmem132e
    ## 331       Tmem143
    ## 332       Tmem186
    ## 333      Tmem229b
    ## 334      Tmem254a
    ## 335      Tmem254b
    ## 336       Tmem266
    ## 337       Tmem87b
    ## 338        Tmem94
    ## 339          Tmx3
    ## 340       Tomm40l
    ## 341      Tor1aip2
    ## 342         Trhde
    ## 343        Trim45
    ## 344         Tshz3
    ## 345       Tspan14
    ## 346         Tubb6
    ## 347       Twistnb
    ## 348          Tyw5
    ## 349         Ubac2
    ## 350       Ubash3a
    ## 351         Ube3c
    ## 352          Ubn2
    ## 353          Ulk2
    ## 354        Unc13c
    ## 355          Urb1
    ## 356          Use1
    ## 357         Usp36
    ## 358       Vipas39
    ## 359         Vipr1
    ## 360         Vma21
    ## 361        Vstm2b
    ## 362         Vti1a
    ## 363         Wash1
    ## 364         Wbp11
    ## 365         Wdfy2
    ## 366         Wdpcp
    ## 367         Whsc1
    ## 368           Wls
    ## 369          Zfat
    ## 370        Zfp146
    ## 371        Zfp184
    ## 372       Zfp385a
    ## 373        Zfp420
    ## 374        Zfp493
    ## 375        Zfp512
    ## 376       Zfp518a
    ## 377        Zfp526
    ## 378        Zfp644
    ## 379        Zfp707
    ## 380        Zfp746
    ## 381         Zfp84
    ## 382        Zfp937
    ## 383          Zic3
    ## 384       Zkscan3
    ## 385         Znfx1

    # stress only CA1 and DG
    shared %>% filter(CA1learn == 0 & CA1stress == 1 & DGlearn == 1) %>%
      select(gene)

    ##      gene
    ## 1  Dusp16
    ## 2  Entpd1
    ## 3    Ier3
    ## 4    Mest
    ## 5 Rasl11b

    # All degs
    updown <- read.csv("../data/02c_setsize_updown.csv") 
    updown$direction <- factor(updown$direction,  levels = c("yoked_consistent", "yoked_conflict", "consistent"))
    updown$set <- factor(updown$set,  levels = c("DGtrain", "CA1train", "CA1stress"))
    updown$status <- factor(updown$status,  levels = c("train", "stress"))

    levels(updown$set)[levels(updown$set)=="DGtrain"] <- "DG train"
    levels(updown$set)[levels(updown$set)=="CA1train"]   <- "CA1 train"
    levels(updown$set)[levels(updown$set)=="CA1stress"]   <- "CA1 stress"
    levels(updown$direction)[levels(updown$direction)=="yoked_consistent"] <- "yoked\nconsistent"
    levels(updown$direction)[levels(updown$direction)=="yoked_conflict"]   <- "yoked\nconflict"

    # unique and shared
    shared2 <- read.csv("../data/02c_intersect_updown_shared2.csv") 
    shared2$group <- paste(shared2$set, shared2$direction, sep = " ")
    shared2$set <- factor(shared2$set,  levels = c("DGtrain", "CA1train", "CA1stress"))
    shared2$direction <- factor(shared2$direction,  
                                levels = c("yoked_consistent", "yoked_conflict", "consistent"))
    shared2$group <-factor(shared2$group,  
                           levels = c("DGtrain yoked_consistent", "DGtrain consistent", 
                                      "CA1train yoked_consistent","CA1train consistent",
                                      "CA1stress yoked_consistent", "CA1stress yoked_conflict" ))
    shared2$status <- factor(shared2$status,  levels = c(  "stress","train","unique"))

    levels(shared2$direction)[levels(shared2$direction)=="yoked_consistent"] <- "yoked\nconsistent"
    levels(shared2$direction)[levels(shared2$direction)=="yoked_conflict"]   <- "yoked\nconflict"
    levels(shared2$status)[levels(shared2$status)=="stress"]   <- "shared"
    levels(shared2$status)[levels(shared2$status)=="train"]   <- "shared"

    d1 <- ggplot(updown, aes(x=direction, y=setsize, fill = direction)) +
      geom_bar(stat="identity", position=position_dodge()) +
      theme_cowplot(font_size = 6, line_size = 0.25) +
      scale_fill_manual(values = c("#404040", "#bababa", "#ca0020"),
                        name = NULL)  +
      labs(x = NULL, y = "Total DEGs") +
      scale_y_continuous(limits = c(0, 560),
                         breaks = c(0,125,250,375,500)) +
        facet_wrap(~set, scales = "free_x") +
        theme(axis.text.x = element_blank(),
              legend.position = "none",
              legend.text=element_text(size=4),
              legend.key.size = unit(0.2, "cm"),
              #legend.margin=margin(t=-0.25, r=0, b=0, l=0, unit="cm"),
              strip.text.x = element_text(size = 5),
              strip.background = element_rect(colour=NA, fill=NA),
            panel.border = element_rect(fill = NA, color = "black"),
            legend.margin=margin(t=-0.1, r=0, b=-0.1, l=0, unit="cm"))



    images = c(
      unique = "../figures/00_schematics/patterns_blank.png",
      shared = "../figures/00_schematics/patterns_crosed-lines.png")
    images

    ##                                               unique 
    ##        "../figures/00_schematics/patterns_blank.png" 
    ##                                               shared 
    ## "../figures/00_schematics/patterns_crosed-lines.png"

    head(shared2)

    ##         set setsize         direction subfield status
    ## 1  CA1train     168 yoked\nconsistent      CA1 unique
    ## 2  CA1train     181        consistent      CA1 unique
    ## 3 CA1stress     179 yoked\nconsistent      CA1 unique
    ## 4 CA1stress     206   yoked\nconflict      CA1 unique
    ## 5   DGtrain       6 yoked\nconsistent       DG unique
    ## 6   DGtrain     103        consistent       DG unique
    ##                        group
    ## 1  CA1train yoked_consistent
    ## 2        CA1train consistent
    ## 3 CA1stress yoked_consistent
    ## 4   CA1stress yoked_conflict
    ## 5   DGtrain yoked_consistent
    ## 6         DGtrain consistent

    d2 <- ggplot(shared2, aes(x=direction, y=setsize, image = status)) +
      geom_textured_bar(stat = "identity") +
      theme_cowplot(font_size = 6, line_size = 0.25) +
      scale_image_manual(values = images,
                         name = NULL) +
      labs(x = "subfield * treatment", y = "Shared DEGs") +
        scale_y_continuous(limits = c(0, 560),
                         breaks = c(0,125,250,375,500)) +
      theme(axis.text.x=element_text(angle=60, vjust = 1, hjust = 1),
            legend.position = c(0.05, 0.9),
            legend.text=element_text(size=4),
            legend.key.size = unit(0.2, "cm"),
            strip.text.x = element_text(size = 0),
            #axis.text=element_text(size=7),
            strip.background = element_rect(colour=NA, fill=NA),
            panel.border = element_rect(fill = NA, color = "black"),
            legend.margin=margin(t=-0.1, r=0, b=-0.1, l=-0.1, unit="cm")) +
      facet_wrap(~set, scales = "free_x") 
    d2

![](../figures/02c_rnaseqSubfield/newbarplot-1.png)

    newbarplots <- plot_grid(d1,d2, nrow = 2, rel_heights =c(0.4,0.6))
    newbarplots

![](../figures/02c_rnaseqSubfield/newbarplot-2.png)

    pdf(file="../figures/02c_rnaseqSubfield/barplots.pdf", width=1.9, height=2.15)
    plot(newbarplots)    
    dev.off() 

    ## quartz_off_screen 
    ##                 2

    bottomplots <- plot_grid(CA1yoked, newbarplots,
                        labels = c("C", "D"),
                        label_size = 7)

    volcanos <- plot_grid(training, bottomplots, nrow = 2,
                          rel_heights = c(0.475,0.525)) 
    volcanos

![](../figures/02c_rnaseqSubfield/combo-1.png)

    #althought this technically now contains volcanos and bar plots

    pdf(file="../figures/02c_rnaseqSubfield/volcanos.pdf", width=3.15, height=4)
    plot(volcanos)    
    dev.off()     

    ## quartz_off_screen 
    ##                 2
