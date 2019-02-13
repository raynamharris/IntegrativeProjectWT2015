Subfield analysis
-----------------

This script is used to identify treatement differences within each
subfield, generate volcano plots, venn diagrams, and tables for
subsequent GO analyses. The final mutlipanel figures for the manuscript
have been inserted just below the subheadings.

    library(ggplot2) ## for awesome plots!
    library(cowplot) ## for some easy to use themes
    library(dplyr) ## for filtering and selecting rows

    ## Warning: package 'dplyr' was built under R version 3.5.1

    library(car) ## stats
    library(VennDiagram) ## venn diagrams
    library(pheatmap) ## awesome heatmaps
    library(viridis) # for awesome color pallette
    library(reshape2) ## for melting dataframe
    library(DESeq2) ## for gene expression analysis

    ## Warning: package 'IRanges' was built under R version 3.5.1

    ## Warning: package 'GenomicRanges' was built under R version 3.5.1

    ## Warning: package 'DelayedArray' was built under R version 3.5.1

    library(edgeR)  ## for basic read counts status

    ## Warning: package 'edgeR' was built under R version 3.5.1

    ## Warning: package 'limma' was built under R version 3.5.1

    library(magrittr) ## to use the weird pipe
    library(genefilter)  ## for PCA fuction
    library(xtable) # for latex or html tables

    ## load functions 
    source("figureoptions.R")
    source("functions_RNAseq.R")

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/02c_rnaseqSubfield/')

DG
--

The most notable comparison within DG is the consistent verses
yoked-consistent.

![fig4](../figures/figures2-01.png)

    colData <- read.csv("../data/02a_colData.csv", header = T)
    countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)

    colData <- colData %>% 
      filter(Punch %in% c("DG"))  %>% 
      droplevels()

    savecols <- as.character(colData$RNAseqID) 
    savecols <- as.vector(savecols) 
    countData <- countData %>% dplyr::select(one_of(savecols)) 

    colData %>% select(APA2,Punch)  %>%  summary()

    ##                APA2   Punch  
    ##  conflict        :5   DG:16  
    ##  consistent      :3          
    ##  yoked_conflict  :4          
    ##  yoked_consistent:4

    ## create DESeq object using the factors Punch and APA
    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ APA2)

    dds # view the DESeq object - note numnber of genes

    ## class: DESeqDataSet 
    ## dim: 22485 16 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(22485): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(16): 143A-DG-1 143B-DG-1 ... 148A-DG-3 148B-DG-4
    ## colData names(8): RNAseqID Mouse ... ID APA2

    dds <- dds[ rowSums(counts(dds)) > 1, ]  # Pre-filtering genes with 0 counts
    dds # view number of genes afternormalization and the number of samples

    ## class: DESeqDataSet 
    ## dim: 17011 16 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(17011): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(16): 143A-DG-1 143B-DG-1 ... 148A-DG-3 148B-DG-4
    ## colData names(8): RNAseqID Mouse ... ID APA2

    dds <- DESeq(dds) # Differential expression analysis
    rld <- rlog(dds, blind=FALSE) ## log transformed data
    vsd <- getVarianceStabilizedData(dds)

    write.csv(colData, file = "../data/02c_DGcolData.csv", row.names = T)
    write.csv(vsd, file = "../data/02c_DGvsd.csv", row.names = T)

    ###  "consistent", "yoked_consistent"
    res <- results(dds, contrast =c("APA2", "consistent", "yoked_consistent"), independentFiltering = T, alpha = 0.1)
    summary(res)

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

    resOrdered <- res[order(res$padj),]
    head(resOrdered, 10)

    ## log2 fold change (MLE): APA2 consistent vs yoked_consistent 
    ## Wald test p-value: APA2 consistent vs yoked_consistent 
    ## DataFrame with 10 rows and 6 columns
    ##                baseMean   log2FoldChange             lfcSE
    ##               <numeric>        <numeric>         <numeric>
    ## Smad7  171.392871064045 3.53587630546276 0.418085263001235
    ## Sgk1   341.089572273562 2.52942417594796 0.361917138618265
    ## Lmna   127.261228543472 2.38190944527576 0.360880330414225
    ## Tiparp 146.843901753813 3.00078251228731 0.456274430076876
    ## Fzd5   26.8401177227407 4.05654356592169 0.655253175606779
    ## Acan   50.8597490321187 2.45912773236628 0.428543544588545
    ## Egr4   683.770839985962 3.23264668960927 0.562954537617405
    ## Errfi1  196.30327794802 2.16723775140178 0.378461151913533
    ## Rasd1  72.8100929443534 3.11825818416375 0.542898092373717
    ## Per1   512.774508684412 1.82273176771818 0.322402773641593
    ##                    stat               pvalue                 padj
    ##               <numeric>            <numeric>            <numeric>
    ## Smad7  8.45730911460593 2.73617917722616e-17 3.38821067515915e-13
    ## Sgk1      6.98895936679 2.76932604930945e-12 1.71462822342994e-08
    ## Lmna   6.60027506221178 4.10395627378058e-11 1.48901359411528e-07
    ## Tiparp 6.57670540902702 4.80986382658575e-11 1.48901359411528e-07
    ## Fzd5    6.1908033672103 5.98583433033886e-10 1.48245173025172e-06
    ## Acan    5.7383380602019 9.56101231319715e-09 1.41105309237205e-05
    ## Egr4   5.74228729604137 9.34061706254104e-09 1.41105309237205e-05
    ## Errfi1 5.72644706185573 1.02555744418545e-08 1.41105309237205e-05
    ## Rasd1  5.74372654457076 9.26153120045455e-09 1.41105309237205e-05
    ## Per1   5.65358587685248 1.57134545200805e-08 1.94579707322157e-05

    data <- data.frame(gene = row.names(res),
                       padj = res$padj, 
                       logpadj = -log10(res$padj),
                       lfc = res$log2FoldChange)
    data <- na.omit(data)
    data <- data %>%
      mutate(direction = ifelse(data$lfc > 1 & data$padj < 0.1, 
                            yes = "consistent", 
                            no = ifelse(data$lfc < -1 & data$padj < 0.1, 
                                        yes = "yoked_consistent", 
                                        no = "neither")))
    DGvolcano <- ggplot(data, aes(x = lfc, y = logpadj)) + 
      geom_point(aes(color = factor(direction)), size = 1, alpha = 0.5, na.rm = T) + # add gene points
      theme_cowplot(font_size = 8, line_size = 0.25) +
      geom_hline(yintercept = 1,  size = 0.25, linetype = 2) + 
      scale_color_manual(values = volcano1)  + 
      scale_y_continuous(limits=c(0, 8)) +
      scale_x_continuous( limits=c(-3, 3),
                          name="Log fold change")+
      ylab(paste0("log10 p-value")) +       
      theme(panel.grid.minor=element_blank(),
            legend.position = "none", # remove legend 
            panel.grid.major=element_blank())
    DGvolcano

![](../figures/02c_rnaseqSubfield/DG-1.png)

    pdf(file="../figures/02c_rnaseqSubfield/DGvolcano.pdf", width=1.5, height=2)
    plot(DGvolcano)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    # save DEGs
    DGvolcanoDEGs <- data %>%
      filter(direction != "neither") %>%
      arrange(desc(lfc))
    head(DGvolcanoDEGs,10)

    ##       gene         padj   logpadj      lfc  direction
    ## 1  Col10a1 2.958539e-04  3.528923 6.862672 consistent
    ## 2   Lrrc32 9.769907e-03  2.010110 6.049415 consistent
    ## 3    Thbs1 3.275529e-04  3.484719 5.491177 consistent
    ## 4     Fzd5 1.482452e-06  5.829019 4.056544 consistent
    ## 5    Nlrp3 9.260166e-04  3.033381 3.762818 consistent
    ## 6     Ier3 1.588050e-02  1.799136 3.758436 consistent
    ## 7    Npas4 2.205531e-04  3.656487 3.537464 consistent
    ## 8    Smad7 3.388211e-13 12.470030 3.535876 consistent
    ## 9    Nr4a3 1.114124e-02  1.953066 3.372927 consistent
    ## 10    Rfx2 2.636311e-04  3.579003 3.356056 consistent

    write.csv(DGvolcanoDEGs, "../data/DG-consistent-yokedconsistent.csv", row.names = F)

    # are any protein kinases differentially expressed?
    pkcs <- data[grep("Prkc", data$gene), ]
    pkcs # no pkcs are differentially expressed

    ##         gene      padj      logpadj         lfc direction
    ## 8139   Prkca 0.9999357 2.793108e-05 -0.06147923   neither
    ## 8140   Prkcb 0.9999357 2.793108e-05 -0.23086376   neither
    ## 8141   Prkcd 0.9999357 2.793108e-05 -1.74799258   neither
    ## 8142 Prkcdbp 0.9999357 2.793108e-05  0.89326487   neither
    ## 8143   Prkce 0.9999357 2.793108e-05 -0.08024100   neither
    ## 8144   Prkcg 0.9999357 2.793108e-05 -0.35258124   neither
    ## 8145   Prkci 0.9999357 2.793108e-05  0.14227925   neither
    ## 8146  Prkcsh 0.9999357 2.793108e-05 -0.12210803   neither
    ## 8147   Prkcz 0.9999357 2.793108e-05 -0.07894751   neither

    ## go setup
    table(res$padj<0.1)

    ## 
    ## FALSE  TRUE 
    ## 12258   125

    logs <- data.frame(cbind("gene"=row.names(res),"logP"=round(-log(res$pvalue+1e-10,10),1)))
    logs$logP=as.numeric(as.character(logs$logP))
    sign <- rep(1,nrow(logs))
    sign[res$log2FoldChange<0]=-1  ##change to correct model
    table(sign)

    ## sign
    ##   -1    1 
    ## 8853 8158

    logs$logP <- logs$logP*sign
    write.csv(logs, file = "./02e_GO_MWU/DGconsistentyoked.csv", row.names = F)

    ## yoked yoked
    res <- results(dds, contrast =c("APA2", "yoked_conflict", "yoked_consistent"), independentFiltering = T, alpha = 0.1)
    summary(res)

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

    resOrdered <- res[order(res$padj),]
    head(resOrdered, 10)

    ## log2 fold change (MLE): APA2 yoked_conflict vs yoked_consistent 
    ## Wald test p-value: APA2 yoked_conflict vs yoked_consistent 
    ## DataFrame with 10 rows and 6 columns
    ##                baseMean   log2FoldChange             lfcSE
    ##               <numeric>        <numeric>         <numeric>
    ## Nlrp3  20.9424446578681 4.09843296366574  0.77923374127775
    ## Kcnc2  22.0989714545398  4.0853572173812 0.828862390929107
    ## Gm2115 18.9055674115327 3.48302773153726 0.740349228112595
    ## Rnase4 15.5820343950564 3.55390988651979 0.842613820491949
    ## Cxcl14 58.4698286328215 1.83345102740602 0.452447198842891
    ## Sst    5.79777379074028 6.11617852820523  1.54412893086366
    ## Cnr1   121.554806190537 3.88279836743752  1.00562221165236
    ## Dner    48.134823231393 1.78996353992257 0.466270200246526
    ## Itga5  18.4196792457967   2.978384440916 0.770052016771093
    ## Myb    6.92940033536349 5.85347258329016  1.53376995782851
    ##                    stat               pvalue                padj
    ##               <numeric>            <numeric>           <numeric>
    ## Nlrp3  5.25956814568287 1.44394111343109e-07 0.00245340034583077
    ## Kcnc2    4.928872708052 8.27054357503709e-07 0.00702624029417276
    ## Gm2115 4.70457400275367 2.54396579921199e-06  0.0144081742981369
    ## Rnase4 4.21772085870238 2.46784058046291e-05   0.104827698256613
    ## Cxcl14 4.05229832805898 5.07169366376427e-05   0.172346294082037
    ## Sst    3.96092476862301 7.46600544114441e-05   0.211424830750808
    ## Cnr1   3.86109049944074 0.000112882079380255   0.228300316931863
    ## Dner   3.83889757264389 0.000123587970402978   0.228300316931863
    ## Itga5  3.86777045712401 0.000109834980276631   0.228300316931863
    ## Myb     3.8163953814674 0.000135415452559867   0.228300316931863

    ## go setup
    table(res$padj<0.1)

    ## 
    ## FALSE  TRUE 
    ## 16988     3

    logs <- data.frame(cbind("gene"=row.names(res),"logP"=round(-log(res$pvalue+1e-10,10),1)))
    logs$logP=as.numeric(as.character(logs$logP))
    sign <- rep(1,nrow(logs))
    sign[res$log2FoldChange<0]=-1  ##change to correct model
    table(sign)

    ## sign
    ##   -1    1 
    ## 7543 9468

    logs$logP <- logs$logP*sign
    write.csv(logs, file = "./02e_GO_MWU/DGyokedyoked.csv", row.names = F)


    ##  "conflict", "yoked_conflict"
    res <- results(dds, contrast =c("APA2", "conflict", "yoked_conflict"), independentFiltering = T, alpha = 0.1)
    summary(res)

    ## 
    ## out of 17011 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 8, 0.047%
    ## LFC < 0 (down)     : 2, 0.012%
    ## outliers [1]       : 20, 0.12%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    resOrdered <- res[order(res$padj),]
    head(resOrdered, 10)

    ## log2 fold change (MLE): APA2 conflict vs yoked_conflict 
    ## Wald test p-value: APA2 conflict vs yoked_conflict 
    ## DataFrame with 10 rows and 6 columns
    ##                  baseMean    log2FoldChange             lfcSE
    ##                 <numeric>         <numeric>         <numeric>
    ## Rps12    18.1519752295379 -21.4660372235101  3.98668946399401
    ## Smad7    171.392871064045  1.75395877649941 0.367789349169671
    ## Dbpht2   179.735977990203  1.43077032426827  0.32133350449148
    ## Insm1    9.54062645923209 -4.64275029977526  1.03675452742903
    ## Slc16a1  51.2786850104172  1.93718965194453 0.446154808676806
    ## Ankrd33b  209.87162559181  1.14655990310839 0.276132373178601
    ## Nptx2    287.767100019705  1.36387651791528 0.326912129485928
    ## Sgk1     341.089572273562  1.32296838420707  0.31899157276228
    ## Fzd5     26.8401177227407  2.51256156493201 0.616467582890519
    ## Acan     50.8597490321187  1.57030753918313 0.390249971345474
    ##                       stat               pvalue               padj
    ##                  <numeric>            <numeric>          <numeric>
    ## Rps12    -5.38442670726719 7.26759094316345e-08 0.0012348363771529
    ## Smad7     4.76892215736857 1.85214205533671e-06 0.0157348728311131
    ## Dbpht2    4.45260237189555 8.48357704857436e-06 0.0360361144080817
    ## Insm1    -4.47815772870406  7.5289965568597e-06 0.0360361144080817
    ## Slc16a1     4.341967438813 1.41212464833709e-05  0.047986819799791
    ## Ankrd33b  4.15221109321652 3.29278313316402e-05 0.0714369056620668
    ## Nptx2     4.17199728887387 3.01941230562506e-05 0.0714369056620668
    ## Sgk1      4.14734587735638 3.36351742273282e-05 0.0714369056620668
    ## Fzd5       4.0757399653539 4.58682534733302e-05 0.0865941660850393
    ## Acan        4.023850491953 5.72542716941817e-05 0.0972807330355842

    ## go setup
    table(res$padj<0.1)

    ## 
    ## FALSE  TRUE 
    ## 16981    10

    logs <- data.frame(cbind("gene"=row.names(res),"logP"=round(-log(res$pvalue+1e-10,10),1)))
    logs$logP=as.numeric(as.character(logs$logP))
    sign <- rep(1,nrow(logs))
    sign[res$log2FoldChange<0]=-1  ##change to correct model
    table(sign)

    ## sign
    ##   -1    1 
    ## 9674 7337

    logs$logP <- logs$logP*sign
    write.csv(logs, file = "./02e_GO_MWU/DGconflictyoked.csv", row.names = F)

    #### "conflict", "consistent"
    res <- results(dds, contrast =c("APA2", "conflict", "consistent"), independentFiltering = T, alpha = 0.1)
    summary(res)

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

    resOrdered <- res[order(res$padj),]
    head(resOrdered, 10)

    ## log2 fold change (MLE): APA2 conflict vs consistent 
    ## Wald test p-value: APA2 conflict vs consistent 
    ## DataFrame with 10 rows and 6 columns
    ##                       baseMean      log2FoldChange             lfcSE
    ##                      <numeric>           <numeric>         <numeric>
    ## Rps12         18.1519752295379   -20.9602515412509  4.31657511842665
    ## 0610007P14Rik 41.0311096797814  -0.147853585358429 0.401936989088285
    ## 0610009B22Rik 8.74131849081542  -0.757911677387021 0.749958107311731
    ## 0610009L18Rik 2.61860548764425   0.604880308193334  1.47685825843458
    ## 0610009O20Rik 48.3693444464377   0.245771257518523 0.338998927574654
    ## 0610010F05Rik 60.2560857800178  -0.366265579236547 0.337328940633006
    ## 0610010K14Rik 21.0476159827708   0.408843299351544 0.633957053036292
    ## 0610012G03Rik 54.2812123171702   0.107773629510678 0.399888012498432
    ## 0610030E20Rik 44.9946840589216 -0.0960001442131723 0.418155894506165
    ## 0610037L13Rik  88.714181322051   0.113994422485766 0.427197151641241
    ##                             stat               pvalue               padj
    ##                        <numeric>            <numeric>          <numeric>
    ## Rps12          -4.85575970907479 1.19926058988033e-06 0.0203766366826567
    ## 0610007P14Rik -0.367852646987792    0.712983110459147  0.999877921572714
    ## 0610009B22Rik  -1.01060535248269    0.312205353000397  0.999877921572714
    ## 0610009L18Rik  0.409572350453244    0.682119683062293  0.999877921572714
    ## 0610009O20Rik  0.724991253738995    0.468457396118792  0.999877921572714
    ## 0610010F05Rik  -1.08578166625502    0.277575599518543  0.999877921572714
    ## 0610010K14Rik  0.644906933984594    0.518987497960177  0.999877921572714
    ## 0610012G03Rik  0.269509528023425    0.787537611086085  0.999877921572714
    ## 0610030E20Rik -0.229579794221355    0.818418309399155  0.999877921572714
    ## 0610037L13Rik  0.266842655780387    0.789590310071293  0.999877921572714

    ## go setup
    table(res$padj<0.1)

    ## 
    ## FALSE  TRUE 
    ## 16990     1

    logs <- data.frame(cbind("gene"=row.names(res),"logP"=round(-log(res$pvalue+1e-10,10),1)))
    logs$logP=as.numeric(as.character(logs$logP))
    sign <- rep(1,nrow(logs))
    sign[res$log2FoldChange<0]=-1  ##change to correct model
    table(sign)

    ## sign
    ##   -1    1 
    ## 8670 8341

    logs$logP <- logs$logP*sign
    write.csv(logs, file = "./02e_GO_MWU/DGconflictconsistent.csv", row.names = F)

    ## plot of pkmz
    plotCounts(dds, "Prkcz", intgroup = "APA2", normalized = TRUE, main="Prkcz in DG")

![](../figures/02c_rnaseqSubfield/DG-2.png)

    # order results table by the smallest adjusted p value:
    res <- res[order(res$padj),]

    results = as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res))
    head(results)

    ##                baseMean log2FoldChange     lfcSE       stat       pvalue
    ## Rps12         18.151975    -20.9602515 4.3165751 -4.8557597 1.199261e-06
    ## 0610007P14Rik 41.031110     -0.1478536 0.4019370 -0.3678526 7.129831e-01
    ## 0610009B22Rik  8.741318     -0.7579117 0.7499581 -1.0106054 3.122054e-01
    ## 0610009L18Rik  2.618605      0.6048803 1.4768583  0.4095724 6.821197e-01
    ## 0610009O20Rik 48.369344      0.2457713 0.3389989  0.7249913 4.684574e-01
    ## 0610010F05Rik 60.256086     -0.3662656 0.3373289 -1.0857817 2.775756e-01
    ##                     padj      sig
    ## Rps12         0.02037664 FDR<0.05
    ## 0610007P14Rik 0.99987792  Not Sig
    ## 0610009B22Rik 0.99987792  Not Sig
    ## 0610009L18Rik 0.99987792  Not Sig
    ## 0610009O20Rik 0.99987792  Not Sig
    ## 0610010F05Rik 0.99987792  Not Sig

    ## [1] 125

    ## [1] 10

    ## [1] 1

    ## [1] 3

![](../figures/02c_rnaseqSubfield/DGvenndiagrams-1.png)

CA3
---

There are so few differences in the CA3 that I don’t make any figures
for the manuscript.

    colData <- read.csv("../data/02a_colData.csv", header = T)
    countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)

    colData <- colData %>% 
      filter(Punch %in% c("CA3"))  %>% 
      droplevels()
    savecols <- as.character(colData$RNAseqID) 
    savecols <- as.vector(savecols) 
    countData <- countData %>% dplyr::select(one_of(savecols)) 
    colData %>% select(APA2,Punch)  %>%  summary()

    ##                APA2   Punch   
    ##  conflict        :5   CA3:13  
    ##  consistent      :2           
    ##  yoked_conflict  :3           
    ##  yoked_consistent:3

    ## create DESeq object using the factors Punch and APA
    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ APA2)

    dds # view the DESeq object - note numnber of genes

    ## class: DESeqDataSet 
    ## dim: 22485 13 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(22485): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(13): 143A-CA3-1 144A-CA3-2 ... 148A-CA3-3 148B-CA3-4
    ## colData names(8): RNAseqID Mouse ... ID APA2

    dds <- dds[ rowSums(counts(dds)) > 1, ]  # Pre-filtering genes with 0 counts
    dds # view number of genes afternormalization and the number of samples

    ## class: DESeqDataSet 
    ## dim: 16502 13 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(16502): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(13): 143A-CA3-1 144A-CA3-2 ... 148A-CA3-3 148B-CA3-4
    ## colData names(8): RNAseqID Mouse ... ID APA2

    dds <- DESeq(dds) # Differential expression analysis
    rld <- rlog(dds, blind=FALSE) ## log transformed data

    res <- results(dds, contrast =c("APA2", "consistent", "yoked_consistent"), independentFiltering = T, alpha = 0.1)
    summary(res)

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

    res <- results(dds, contrast =c("APA2", "conflict", "yoked_conflict"), independentFiltering = T, alpha = 0.1)
    summary(res)

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

    res <- results(dds, contrast =c("APA2", "yoked_conflict", "yoked_consistent"), independentFiltering = T, alpha = 0.1)
    summary(res)

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

    res <- results(dds, contrast =c("APA2", "conflict", "consistent"), independentFiltering = T, alpha = 0.1)
    summary(res)

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

    ## go setup
    table(res$padj<0.1)

    ## 
    ## FALSE 
    ## 16491

    logs <- data.frame(cbind("gene"=row.names(res),"logP"=round(-log(res$pvalue+1e-10,10),1)))
    logs$logP=as.numeric(as.character(logs$logP))
    sign <- rep(1,nrow(logs))
    sign[res$log2FoldChange<0]=-1  ##change to correct model
    table(sign)

    ## sign
    ##   -1    1 
    ## 7619 8883

    logs$logP <- logs$logP*sign
    write.csv(logs, file = "./02e_GO_MWU/CA3conflictconsistent.csv", row.names = F)

    plotCounts(dds, "Prkcz", intgroup = "APA2", normalized = TRUE, main="Prkcz in CA3")

![](../figures/02c_rnaseqSubfield/CA3-1.png)

CA1
---

Two comparisons within CA1 are noteable

![fig5](../figures/figures2-02.png)

![fig6](../figures/figures2-03.png)

    colData <- read.csv("../data/02a_colData.csv", header = T)
    countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)
    colData <- colData %>% 
      filter(Punch %in% c("CA1"))  %>% 
      droplevels()
    savecols <- as.character(colData$RNAseqID) 
    savecols <- as.vector(savecols) 
    countData <- countData %>% dplyr::select(one_of(savecols)) 
    colData %>% select(APA2,Punch)  %>%  summary()

    ##                APA2   Punch   
    ##  conflict        :4   CA1:15  
    ##  consistent      :4           
    ##  yoked_conflict  :5           
    ##  yoked_consistent:2

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ APA2)

    dds # view the DESeq object - note numnber of genes

    ## class: DESeqDataSet 
    ## dim: 22485 15 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(22485): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(15): 143B-CA1-1 143C-CA1-1 ... 148A-CA1-3 148B-CA1-4
    ## colData names(8): RNAseqID Mouse ... ID APA2

    dds <- dds[ rowSums(counts(dds)) > 1, ]  # Pre-filtering genes with 0 counts
    dds # view number of genes afternormalization and the number of samples

    ## class: DESeqDataSet 
    ## dim: 16852 15 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(16852): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(15): 143B-CA1-1 143C-CA1-1 ... 148A-CA1-3 148B-CA1-4
    ## colData names(8): RNAseqID Mouse ... ID APA2

    dds <- DESeq(dds) # Differential expression analysis
    rld <- rlog(dds, blind=FALSE) ## log transformed data


    res <- results(dds, contrast =c("APA2", "consistent", "yoked_consistent"), independentFiltering = T, alpha = 0.1)
    summary(res)

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

    resOrdered <- res[order(res$padj),]
    head(resOrdered, 10)

    ## log2 fold change (MLE): APA2 consistent vs yoked_consistent 
    ## Wald test p-value: APA2 consistent vs yoked_consistent 
    ## DataFrame with 10 rows and 6 columns
    ##                 baseMean    log2FoldChange             lfcSE
    ##                <numeric>         <numeric>         <numeric>
    ## Agap1   141.608762023137  2.77875046323028 0.429199296190091
    ## Mga     103.422400328703  2.88512466619299 0.476961395194488
    ## Adamts1  114.05510912055  3.03089016770375 0.542440999554141
    ## Gpd1    249.757427758137 -1.22793747874702 0.227810028910739
    ## Sdhaf2  77.6262495970029 -1.88444201930518 0.348951742801298
    ## Scn4b   126.513273389751  2.75983601093444  0.53291136892103
    ## Lhfpl4  214.234684832542  1.68458658679185 0.339420793178537
    ## Ncoa4   101.428981706278  2.35081948916179 0.479986115716983
    ## Lats2   79.0118946108766  2.50181097858582 0.521225803364738
    ## Gad2    132.201113296669  2.85189408073533 0.598342010017182
    ##                      stat               pvalue                 padj
    ##                 <numeric>            <numeric>            <numeric>
    ## Agap1    6.47426612274679 9.52738560798276e-11 1.13642655532018e-06
    ## Mga      6.04896894226951 1.45775755684507e-09 8.69406606902402e-06
    ## Adamts1   5.5875019959682 2.30359111553791e-08 9.15907827537872e-05
    ## Gpd1    -5.39018183096829 7.03864293094537e-08 0.000167913865760633
    ## Sdhaf2  -5.40029404689985 6.65317533489317e-08 0.000167913865760633
    ## Scn4b    5.17878989243971 2.23329823700518e-07 0.000443979689516631
    ## Lhfpl4   4.96312135451805 6.93692079381607e-07  0.00118205130326626
    ## Ncoa4    4.89768227076785 9.69736686637254e-07  0.00144587739977615
    ## Lats2    4.79986018043532 1.58776440997872e-06  0.00210431709802513
    ## Gad2     4.76632767378883 1.87614050241445e-06  0.00223786039127996

    topGene <- rownames(res)[which.min(res$padj)]
    plotCounts(dds, gene = topGene, intgroup=c("APA2"))

![](../figures/02c_rnaseqSubfield/CA1-1.png)

    data <- data.frame(gene = row.names(res),
                       padj = res$padj, 
                       logpadj = -log10(res$padj),
                       lfc = res$log2FoldChange)
    data <- na.omit(data)
    data <- data %>%
      mutate(direction = ifelse(data$lfc > 1 & data$padj < 0.1, 
                            yes = "consistent", 
                            no = ifelse(data$lfc < -1 & data$padj < 0.1, 
                                        yes = "yoked_consistent", 
                                        no = "neither")))
    top_labelled <- top_n(data, n = 5, wt = lfc)

    CA1volcano <- ggplot(data, aes(x = lfc, y = logpadj)) + 
      geom_point(aes(color = factor(direction)), size = 1, alpha = 0.5, na.rm = T) + # add gene points
        theme_cowplot(font_size = 8, line_size = 0.25) +
      scale_color_manual(values = volcano1)  + 
      geom_hline(yintercept = 1,  size = 0.25, linetype = 2) + 
      scale_y_continuous(limits=c(0, 8)) +
      scale_x_continuous( limits=c(-3, 3),
                          name="Log fold change")+
      ylab(paste0("log10 p-value")) +       
      theme(panel.grid.minor=element_blank(),
            legend.position = "none", # remove legend 
            panel.grid.major=element_blank())
    CA1volcano

![](../figures/02c_rnaseqSubfield/CA1-2.png)

    pdf(file="../figures/02c_rnaseqSubfield/CA1volcano.pdf",  width=1.5, height=2)
    plot(CA1volcano)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    # save DEGs
    CA1volcanoDEGs <- data %>%
      filter(direction != "neither") %>%
      arrange(desc(lfc))
    head(CA1volcanoDEGs,10)

    ##       gene        padj  logpadj      lfc  direction
    ## 1    Srprb 0.009269714 2.032934 6.362542 consistent
    ## 2   Zfp536 0.016971128 1.770289 6.278678 consistent
    ## 3    Ahnak 0.019268528 1.715151 6.183278 consistent
    ## 4    Uvssa 0.023546445 1.628075 6.091268 consistent
    ## 5  Gm43951 0.019842316 1.702408 6.014951 consistent
    ## 6     Vcan 0.018269414 1.738275 5.994322 consistent
    ## 7    Pex26 0.025320503 1.596528 5.857468 consistent
    ## 8    Med26 0.030274380 1.518925 5.805185 consistent
    ## 9  Gm16485 0.035396074 1.451045 5.797461 consistent
    ## 10  Adgrf5 0.028566992 1.544135 5.791507 consistent

    write.csv(CA1volcanoDEGs, "../data/CA1-consistent-yokedconsistent.csv", row.names = F)

    ## go setup
    table(res$padj<0.1)

    ## 
    ## FALSE  TRUE 
    ## 11046   882

    logs <- data.frame(cbind("gene"=row.names(res),"logP"=round(-log(res$pvalue+1e-10,10),1)))
    logs$logP=as.numeric(as.character(logs$logP))
    sign <- rep(1,nrow(logs))
    sign[res$log2FoldChange<0]=-1  ##change to correct model
    table(sign)

    ## sign
    ##   -1    1 
    ## 7694 9158

    logs$logP <- logs$logP*sign
    write.csv(logs, file = "./02e_GO_MWU/CA1consistentyoked.csv", row.names = F)

    #conflict yoked conflict

    res <- results(dds, contrast =c("APA2", "conflict", "yoked_conflict"), independentFiltering = T, alpha = 0.1)
    summary(res)

    ## 
    ## out of 16852 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.0059%
    ## LFC < 0 (down)     : 3, 0.018%
    ## outliers [1]       : 32, 0.19%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    resOrdered <- res[order(res$padj),]
    head(resOrdered, 3)

    ## log2 fold change (MLE): APA2 conflict vs yoked_conflict 
    ## Wald test p-value: APA2 conflict vs yoked_conflict 
    ## DataFrame with 3 rows and 6 columns
    ##                 baseMean    log2FoldChange             lfcSE
    ##                <numeric>         <numeric>         <numeric>
    ## Gm20390  47.871821328339  3.00058138879586 0.573131448872674
    ## Il4ra   21.8885208760384 -5.08358235078372  1.03917819120446
    ## Gm21949 21.1841188999619 -17.0604245457745  3.93136294526155
    ##                      stat               pvalue                padj
    ##                 <numeric>            <numeric>           <numeric>
    ## Gm20390  5.23541570559054  1.6461394013611e-07 0.00276880647308937
    ## Il4ra   -4.89192555599306 9.98542153556464e-07 0.00839773951140986
    ## Gm21949 -4.33956996169416 1.42761812252115e-05  0.0800417894026856

    # volcano plots

    res <- results(dds, contrast =c("APA2", "yoked_conflict", "yoked_consistent"), independentFiltering = T, alpha = 0.1)
    summary(res)

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

    resOrdered <- res[order(res$padj),]
    head(resOrdered, 10)

    ## log2 fold change (MLE): APA2 yoked_conflict vs yoked_consistent 
    ## Wald test p-value: APA2 yoked_conflict vs yoked_consistent 
    ## DataFrame with 10 rows and 6 columns
    ##                 baseMean    log2FoldChange             lfcSE
    ##                <numeric>         <numeric>         <numeric>
    ## Agap1   141.608762023137  2.68968384470947 0.425950937249105
    ## Pcdhb12  28.059174402023 -4.44059921537497 0.793522420606138
    ## Sdc3    215.373472922264  1.85322249231914 0.337204018822692
    ## Gm20390  47.871821328339 -3.80408224976329 0.702192890495312
    ## Adamts1  114.05510912055  2.88221063805472 0.537733929661711
    ## Ncoa4   101.428981706278   2.4713781982471 0.475325965110602
    ## Lats2   79.0118946108766   2.6661801320419 0.517454551359232
    ## Scn4b   126.513273389751  2.63717217475313 0.525235880320185
    ## Tgoln1  253.956891941181  1.57640904952372 0.317332343314282
    ## Mga     103.422400328703  2.31524201605359   0.4748411371402
    ##                      stat               pvalue                 padj
    ##                 <numeric>            <numeric>            <numeric>
    ## Agap1    6.31453909241308 2.70967824020662e-10 3.23210420491846e-06
    ## Pcdhb12  -5.5960601743086 2.19277621297338e-08 0.000130777173341732
    ## Sdc3      5.4958493638049 3.88834321833904e-08  0.00015460052636116
    ## Gm20390 -5.41743202082261 6.04611256628736e-08 0.000180295076726689
    ## Adamts1  5.35991961650611 8.32589898656502e-08 0.000198622646223495
    ## Ncoa4    5.19933346723872 2.00004427353209e-07 0.000397608801578179
    ## Lats2    5.15249141212203 2.57048392342877e-07 0.000438010460552263
    ## Scn4b     5.0209292121199 5.14221090642523e-07 0.000766703646148002
    ## Tgoln1   4.96769107447225  6.7754785642546e-07 0.000897976759049209
    ## Mga      4.87582442834981 1.08354951670463e-06  0.00129245786352528

    topGene <- rownames(res)[which.min(res$padj)]
    plotCounts(dds, gene = topGene, intgroup=c("APA2"))

![](../figures/02c_rnaseqSubfield/CA1-3.png)

    data <- data.frame(gene = row.names(res),
                       padj = res$padj, 
                       logpadj = -log10(res$padj),
                       lfc = res$log2FoldChange)
    data <- na.omit(data)
    data <- data %>%
      mutate(direction = ifelse(data$lfc > 1 & data$padj < 0.1, 
                            yes = "yoked_conflict", 
                            no = ifelse(data$lfc < -1 & data$padj < 0.1, 
                                        yes = "yoked_consistent", 
                                        no = "neither")))
    top_labelled <- top_n(data, n = 5, wt = lfc)


    CA1volcano2 <- ggplot(data, aes(x = lfc, y = logpadj)) + 
      geom_point(aes(color = factor(direction)), size = 1, alpha = 0.5, na.rm = T) + # add gene points
      theme_cowplot(font_size = 8, line_size = 0.25) +
      scale_color_manual(values = volcano3)  + 
      geom_hline(yintercept = 1,  size = 0.25, linetype = 2) + 
      scale_y_continuous(limits=c(0, 8)) +
      scale_x_continuous( limits=c(-3, 3),
                          name="Log fold change")+
      ylab(paste0("log10 p-value")) +       
      theme(panel.grid.minor=element_blank(),
            legend.position = "none", # remove legend 
            panel.grid.major=element_blank())
    CA1volcano2

![](../figures/02c_rnaseqSubfield/CA1-4.png)

    pdf(file="../figures/02c_rnaseqSubfield/CA1volcano2.pdf",  width=1.5, height=2)
    plot(CA1volcano2)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    # save DEGs
    CA1volcano2DEGs <- data %>%
      filter(direction != "neither") %>%
      arrange(desc(lfc))
    head(CA1volcano2DEGs,10)

    ##       gene        padj  logpadj      lfc      direction
    ## 1    Srprb 0.002523184 2.598051 6.962723 yoked_conflict
    ## 2   Adgrf5 0.014393982 1.841819 6.323867 yoked_conflict
    ## 3    Ahnak 0.014393982 1.841819 6.288004 yoked_conflict
    ## 4   Fndc10 0.014393982 1.841819 6.187727 yoked_conflict
    ## 5  Slco1a4 0.017582699 1.754914 6.093251 yoked_conflict
    ## 6  Ppfibp1 0.031682073 1.499186 6.031033 yoked_conflict
    ## 7     Optn 0.041595475 1.380954 6.025644 yoked_conflict
    ## 8   Rnf207 0.039413555 1.404354 5.934581 yoked_conflict
    ## 9     Mdc1 0.020331763 1.691825 5.901609 yoked_conflict
    ## 10    Epyc 0.040773654 1.389620 5.895348 yoked_conflict

    write.csv(CA1volcano2DEGs, "../data/CA1-yokedconflict-yokedconsistent.csv", row.names = F)


    res <- results(dds, contrast =c("APA2", "conflict", "consistent"), independentFiltering = T, alpha = 0.1)
    summary(res)

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

    plotCounts(dds, "Prkcz", intgroup = "APA2", normalized = TRUE, main="Prkcz in CA1")

![](../figures/02c_rnaseqSubfield/CA1-5.png)
