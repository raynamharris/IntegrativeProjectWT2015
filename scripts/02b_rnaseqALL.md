The figures made from this script were compiled in Adobe.

<img src="../figures/figures-02.png" width="1370" />

    library(ggplot2) ## for awesome plots!
    library(cowplot) ## for some easy to use themes
    library(dplyr) ## for filtering and selecting rows
    library(car) ## stats
    library(VennDiagram) ## venn diagrams
    library(pheatmap) ## awesome heatmaps
    library(viridis) # for awesome color pallette
    library(reshape2) ## for melting dataframe
    library(DESeq2) ## for gene expression analysis
    library(edgeR)  ## for basic read counts status
    library(magrittr) ## to use the weird pipe
    library(genefilter)  ## for PCA fuction
    library(ggrepel) ## for labeling volcano plot

    ## load functions 
    source("figureoptions.R")
    source("functions_RNAseq.R")

    ## set output file for figures 
    knitr::opts_chunk$set(fig.path = '../figures/02b_RNAseqAll/')

Design
------

The two two catagorical variables are

-   Hippocampal subfield: DG, CA3, CA1
-   Treatment: yoked\_consistent, consistent, yoked\_conflict, conflict

<!-- -->

    colData <- read.csv("../data/02a_colData.csv", header = T)
    countData <- read.csv("../data/02a_countData.csv", header = T, check.names = F, row.names = 1)
    colData %>% select(APA2,Punch)  %>%  summary()

    ##                APA2    Punch   
    ##  conflict        :14   CA1:15  
    ##  consistent      : 9   CA3:13  
    ##  yoked_conflict  :12   DG :16  
    ##  yoked_consistent: 9

    head(colData)

    ##     RNAseqID   Mouse Punch      Group   Conflict Treatment     ID
    ## 1 143A-CA3-1 15-143A   CA3   conflict   Conflict  conflict 15143A
    ## 2  143A-DG-1 15-143A    DG   conflict   Conflict  conflict 15143A
    ## 3 143B-CA1-1 15-143B   CA1    control   Conflict   shocked 15143B
    ## 4  143B-DG-1 15-143B    DG    control   Conflict   shocked 15143B
    ## 5 143C-CA1-1 15-143C   CA1 consistent NoConflict   trained 15143C
    ## 6 143D-CA1-3 15-143D   CA1    control NoConflict     yoked 15143D
    ##               APA2
    ## 1         conflict
    ## 2         conflict
    ## 3   yoked_conflict
    ## 4   yoked_conflict
    ## 5       consistent
    ## 6 yoked_consistent

    totalCounts=colSums(countData)
    ### on average 1 million gene counts per sample 
    summary((colSums(countData)/1000000))

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    ##  0.09042  1.08461  2.17911  2.52574  3.30608 11.70070

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ Punch + APA2 + Punch*APA2)

    dds$Punch <- factor(dds$Punch, levels=c("DG","CA3", "CA1")) ## specify the factor levels
    dds$APA2 <- factor(dds$APA2, levels=c("yoked_consistent", "consistent", "yoked_conflict" , "conflict")) ## specify the factor levels

    dds # view the DESeq object - note numnber of genes

    ## class: DESeqDataSet 
    ## dim: 22485 44 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(22485): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(44): 143A-CA3-1 143A-DG-1 ... 148B-CA3-4 148B-DG-4
    ## colData names(8): RNAseqID Mouse ... ID APA2

    dds <- dds[ rowSums(counts(dds)) > 1, ]  # Pre-filtering genes with 0 counts
    dds # view number of genes afternormalization and the number of samples

    ## class: DESeqDataSet 
    ## dim: 17929 44 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(17929): 0610007P14Rik 0610009B22Rik ... Zzef1 Zzz3
    ## rowData names(0):
    ## colnames(44): 143A-CA3-1 143A-DG-1 ... 148B-CA3-4 148B-DG-4
    ## colData names(8): RNAseqID Mouse ... ID APA2

    dds <- DESeq(dds) # Differential expression analysis
    rld <- rlog(dds, blind=FALSE) ## log transformed data
    vsd <- getVarianceStabilizedData(dds)
    write.csv(vsd, file = "../data/02b_vsd.csv", row.names = T)

this is for CA1 DG
------------------

    res <- results(dds, contrast =c("Punch", "CA1", "DG"), independentFiltering = T)
    summary(res)

    ## 
    ## out of 17929 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1223, 6.8%
    ## LFC < 0 (down)     : 1542, 8.6%
    ## outliers [1]       : 8, 0.045%
    ## low counts [2]     : 5210, 29%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    head((res[order(res$padj),]), 10)

    ## log2 fold change (MLE): Punch CA1 vs DG 
    ## Wald test p-value: Punch CA1 vs DG 
    ## DataFrame with 10 rows and 6 columns
    ##                  baseMean    log2FoldChange             lfcSE
    ##                 <numeric>         <numeric>         <numeric>
    ## Pou3f1   219.411952155873  5.98811582532105 0.529427167656916
    ## Prkcg    1597.44305102598   2.9761049905205 0.285040386563396
    ## Wfs1     558.096580361287  6.41720721777274 0.633476148881828
    ## Synj2    114.149797568097  6.05491769796562 0.603560024461183
    ## St8sia5  68.4953636213458  6.85817710097697 0.695949666130441
    ## Pex5l    417.241798395133  3.75565247595177 0.387901968935755
    ## Tmem200a 51.0216626372409  7.74247836767241 0.800495709544138
    ## Tiam1    321.235486821507 -5.24975765311302 0.550613429163663
    ## Man1a    138.573627062308  5.71253818001133 0.603077877738542
    ## Zfp462   93.4387577077342  5.67696987081699 0.608847677804295
    ##                       stat               pvalue                 padj
    ##                  <numeric>            <numeric>            <numeric>
    ## Pou3f1     11.310556373264 1.16349388204921e-29 1.47891707347276e-25
    ## Prkcg     10.4409940864944 1.61114601446348e-25 1.02396384949226e-21
    ## Wfs1      10.1301481185992 4.06033915011148e-24 1.72036569790223e-20
    ## Synj2      10.032005852891 1.10253681838088e-23 3.50358637460984e-20
    ## St8sia5   9.85441539056871 6.55973177434624e-23  1.6676150116743e-19
    ## Pex5l     9.68196290999953 3.59741650689522e-22 7.19343776710422e-19
    ## Tmem200a  9.67210476628482 3.96145577607816e-22 7.19343776710422e-19
    ## Tiam1    -9.53437997523411 1.50783878004262e-21 2.39576734164022e-18
    ## Man1a     9.47230596723686 2.73732874947874e-21 3.86602063718048e-18
    ## Zfp462    9.32412174304419 1.11906519693207e-20 1.42244377182035e-17

this is for CA1 CA3
-------------------

    res <- results(dds, contrast =c("Punch", "CA1", "CA3"), independentFiltering = T)
    summary(res)

    ## 
    ## out of 17929 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 877, 4.9%
    ## LFC < 0 (down)     : 1288, 7.2%
    ## outliers [1]       : 8, 0.045%
    ## low counts [2]     : 4863, 27%
    ## (mean count < 3)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    head((res[order(res$padj),]), 10)

    ## log2 fold change (MLE): Punch CA1 vs CA3 
    ## Wald test p-value: Punch CA1 vs CA3 
    ## DataFrame with 10 rows and 6 columns
    ##                baseMean    log2FoldChange             lfcSE
    ##               <numeric>         <numeric>         <numeric>
    ## Doc2b  349.417572153451  7.18133518109858 0.485599785928926
    ## Itpka  710.071023144903  3.09537793188557 0.232315140416282
    ## Pou3f1 219.411952155873  6.52776905240328 0.548908573891511
    ## C1ql3  285.386723855134  6.43352476856748 0.609253058871802
    ## Syn2   1094.47006534942 -2.23244918862592 0.211235038450167
    ## Fibcd1 360.508601669652   7.6909143460939 0.794190465956115
    ## Nptxr  1873.89195088131 -2.47431772210074 0.270208264146182
    ## Bcr    341.475318016323  2.32744268995085 0.256737538464715
    ## Ahi1   479.309563295498  1.79803940606164 0.204797358974687
    ## Homer3 182.933674454584 -4.51664918141573 0.517315559977438
    ##                     stat               pvalue                 padj
    ##                <numeric>            <numeric>            <numeric>
    ## Doc2b   14.7885880290517 1.73556525225735e-49 2.26630110639765e-45
    ## Itpka   13.3240473536895 1.67760521698626e-40 1.09530844617033e-36
    ## Pou3f1  11.8922701573495 1.29829032027963e-32 5.65102500073714e-29
    ## C1ql3   10.5596921917485 4.58161703975939e-26 1.19653510610356e-22
    ## Syn2   -10.5685553164163 4.16866681051938e-26 1.19653510610356e-22
    ## Fibcd1  9.68396710332566  3.5275607705144e-22 7.67714809022951e-19
    ## Nptxr  -9.15707641259313 5.33223150962408e-20 9.94689700752446e-17
    ## Bcr     9.06545534349558 1.24084882330284e-19 2.02537549183606e-16
    ## Ahi1    8.77960250592819 1.64053414176282e-18 2.38023275812654e-15
    ## Homer3 -8.73093626182965 2.52572992932246e-18 3.29809814170926e-15

this is for CA3 DG
------------------

    res <- results(dds, contrast =c("Punch", "CA3", "DG"), independentFiltering = T)
    summary(res)

    ## 
    ## out of 17929 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1571, 8.8%
    ## LFC < 0 (down)     : 1377, 7.7%
    ## outliers [1]       : 8, 0.045%
    ## low counts [2]     : 3824, 21%
    ## (mean count < 2)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    head((res[order(res$padj),]), 10)

    ## log2 fold change (MLE): Punch CA3 vs DG 
    ## Wald test p-value: Punch CA3 vs DG 
    ## DataFrame with 10 rows and 6 columns
    ##                 baseMean    log2FoldChange             lfcSE
    ##                <numeric>         <numeric>         <numeric>
    ## Fam163b 628.819601522936 -5.70013539670791 0.336089093685364
    ## Doc2b   349.417572153451 -7.02195793637976 0.443508811683749
    ## C1ql3   285.386723855134 -8.02894919282852 0.547842656156047
    ## Pter    159.610534988909 -8.21109827989322 0.656588791993648
    ## Cpe     3054.96826641212  2.95874463126907 0.240400534705513
    ## Adcy1   2607.63934923601 -3.33229530549559 0.279435782261391
    ## Pitpnm2 935.979067712815 -2.94903024539698 0.256559415494396
    ## Pde1b   258.300187916992 -3.36029952761273 0.293356997935408
    ## Kcnc2   144.358588449054  6.22171240084158 0.558906647495304
    ## Me1     87.7852191159365  4.42363434729768 0.412345208515744
    ##                      stat               pvalue                 padj
    ##                 <numeric>            <numeric>            <numeric>
    ## Fam163b -16.9601915200622 1.61816164273402e-64 2.28112246776215e-60
    ## Doc2b   -15.8327360164986 1.84998786056728e-56 1.30396394352084e-52
    ## C1ql3   -14.6555751046548 1.24101383873535e-48 5.83152402821739e-45
    ## Pter     -12.505693639639  6.9490587180363e-36 2.44902201870394e-32
    ## Cpe      12.3075626054389 8.24764911040337e-35 2.32534219018712e-31
    ## Adcy1   -11.9250844631576 8.75951409278765e-33 2.05804783610046e-29
    ## Pitpnm2  -11.494531353348 1.40541696884345e-30 2.83030900139802e-27
    ## Pde1b    -11.454642470648 2.22881011954894e-30 3.92744203191018e-27
    ## Kcnc2    11.1319348744977  8.7710653025308e-29 1.37384119521974e-25
    ## Me1      10.7279877538065 7.52162176228496e-27 1.06032301982931e-23

this is for consistent yoked-consistent
---------------------------------------

    res <- results(dds, contrast =c("APA2", "consistent", "yoked_consistent"), independentFiltering = T)
    summary(res)

    ## 
    ## out of 17929 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 106, 0.59%
    ## LFC < 0 (down)     : 8, 0.045%
    ## outliers [1]       : 8, 0.045%
    ## low counts [2]     : 4863, 27%
    ## (mean count < 3)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    head((res[order(res$padj),]), 10)

    ## log2 fold change (MLE): APA2 consistent vs yoked_consistent 
    ## Wald test p-value: APA2 consistent vs yoked consistent 
    ## DataFrame with 10 rows and 6 columns
    ##                baseMean   log2FoldChange             lfcSE
    ##               <numeric>        <numeric>         <numeric>
    ## Plk2   690.892982346694 2.33695480383662 0.287228730093021
    ## Frmd6  115.436102973485 3.28456907410619 0.460199207192014
    ## Sgk1   243.126788527358 2.51900483921696 0.352679536541189
    ## Smad7  101.648353958591 3.55571966960197 0.497558063998369
    ## Arc    865.501371190336 2.88376455109463 0.419805445423536
    ## Fbxo33 220.541754552169 3.00911903213762  0.44446108800189
    ## Lmna   104.873184036613 2.37485429509161 0.350610025411292
    ## Junb   509.171327352875 2.38648428642696 0.353936652165955
    ## Fosl2  294.146709990136 2.56630291348006 0.405437600249251
    ## Dnaja1 528.535428970117 1.01336750965102 0.163895125829566
    ##                    stat               pvalue                 padj
    ##               <numeric>            <numeric>            <numeric>
    ## Plk2   8.13621535380455 4.07826108076359e-16 5.32539331926109e-12
    ## Frmd6  7.13727668969176 9.51980996986564e-13 3.10774196466264e-09
    ## Sgk1   7.14247518844285 9.16650670302539e-13 3.10774196466264e-09
    ## Smad7  7.14634115469512 8.91214084671919e-13 3.10774196466264e-09
    ## Arc    6.86928810126616 6.45230809755564e-12 1.68508478275763e-08
    ## Fbxo33 6.77026428942375 1.28547354747413e-11 2.39795908327389e-08
    ## Lmna   6.77349226481964  1.2571029992111e-11 2.39795908327389e-08
    ## Junb   6.74268763018071 1.55483254346742e-11 2.53787541907469e-08
    ## Fosl2  6.32971118589489 2.45620507490642e-10 3.56368065201422e-07
    ## Dnaja1  6.1830240803184 6.28851172544568e-10 7.47676388243095e-07

this is for consistent yoked-conflict yoked-consistent DG
---------------------------------------------------------

    res <- results(dds, contrast =c("APA2", "yoked_conflict", "yoked_consistent"), independentFiltering = T)
    summary(res)

    ## 
    ## out of 17929 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 38, 0.21%
    ## LFC < 0 (down)     : 2, 0.011%
    ## outliers [1]       : 8, 0.045%
    ## low counts [2]     : 2434, 14%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    head((res[order(res$padj),]), 10)

    ## log2 fold change (MLE): APA2 yoked_conflict vs yoked_consistent 
    ## Wald test p-value: APA2 yoked conflict vs yoked consistent 
    ## DataFrame with 10 rows and 6 columns
    ##                 baseMean   log2FoldChange             lfcSE
    ##                <numeric>        <numeric>         <numeric>
    ## Kcnc2   144.358588449054  3.8023112290593 0.556925143128729
    ## St8sia5 68.4953636213458 3.78869790825454 0.673902508231409
    ## Gm2115  249.257751832989 3.48213722650287 0.658706293164926
    ## Cnr1    350.970175675373 3.92670174540615 0.762127039058243
    ## Dner    169.806900335675 1.75376437273897 0.352584944840884
    ## Camk1g  40.6858606603921 3.07441617926896 0.651972053515784
    ## Tenm4   155.064731202057 2.09080515700774  0.45595322720424
    ## Klk8    34.8871312628441 5.00760720075993   1.1276964055125
    ## Pou3f1  219.411952155873 2.14125662429894 0.479904106088455
    ## Slc32a1 91.6157514759355 4.03971461796684 0.915087632883217
    ##                     stat               pvalue                 padj
    ##                <numeric>            <numeric>            <numeric>
    ## Kcnc2   6.82732908717039 8.65100952684687e-12 1.33978184542277e-07
    ## St8sia5 5.62202672044894 1.88730063642428e-08 0.000146143124781514
    ## Gm2115  5.28632755240888 1.24796350258549e-07  0.00064424035881805
    ## Cnr1     5.1522929172784 2.57320693829982e-07 0.000996281396336233
    ## Dner    4.97401944808056 6.55787255813261e-07  0.00203123544615599
    ## Camk1g  4.71556435999067 2.41041682634293e-06  0.00622168756492884
    ## Tenm4   4.58556937918368 4.52750896214222e-06   0.0100167901852424
    ## Klk8     4.4405632369504 8.97237334987668e-06   0.0154394606743933
    ## Pou3f1  4.46184268301358 8.12578728841866e-06   0.0154394606743933
    ## Slc32a1 4.41456585446214 1.01212902788901e-05   0.0156748422549171

this is for consistent-conflict vs yoked-conflict
-------------------------------------------------

    res <- results(dds, contrast =c("APA2", "conflict", "yoked_conflict"), independentFiltering = T)
    summary(res)

    ## 
    ## out of 17929 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 17, 0.095%
    ## LFC < 0 (down)     : 44, 0.25%
    ## outliers [1]       : 8, 0.045%
    ## low counts [2]     : 6945, 39%
    ## (mean count < 10)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    head((res[order(res$padj),]), 10)

    ## log2 fold change (MLE): APA2 conflict vs yoked_conflict 
    ## Wald test p-value: APA2 conflict vs yoked_conflict 
    ## DataFrame with 10 rows and 6 columns
    ##                 baseMean    log2FoldChange             lfcSE
    ##                <numeric>         <numeric>         <numeric>
    ## Camk1g  40.6858606603921 -3.40364428501139 0.620423076444881
    ## Insm1   21.1075753227574 -4.60154941196191 0.836635136992648
    ## Kcnc2   144.358588449054 -2.52302302621629 0.471272220797739
    ## Neurod6 331.806267643941 -3.56020520569808 0.703632664264456
    ## Sv2b    568.239350545566 -3.75707558827698 0.743875755350764
    ## Khdrbs3 327.062082451813 -1.63815194028196 0.355869901580467
    ## Stmn2    416.33755501586 -2.26179598177661 0.497962984351558
    ## Stox2    4694.9508569463 -5.24078244314957  1.17277513950018
    ## Slc6a7  90.3559406351636 -1.81763187373056 0.410369551227979
    ## Plk2    690.892982346694    1.105671692536 0.251850117127218
    ##                      stat               pvalue                 padj
    ##                 <numeric>            <numeric>            <numeric>
    ## Camk1g  -5.48600529902079 4.11124568870954e-08  0.00022562516339638
    ## Insm1   -5.50006712424552 3.79646693599839e-08  0.00022562516339638
    ## Kcnc2    -5.3536425761431 8.62010539514554e-08 0.000315380922723725
    ## Neurod6 -5.05974976221399 4.19807053469746e-07 0.000966434064873371
    ## Sv2b    -5.05067622012413  4.4024875404217e-07 0.000966434064873371
    ## Khdrbs3 -4.60323262239009 4.15983148344033e-06  0.00760971839370685
    ## Stmn2   -4.54209660728477 5.56975048478885e-06  0.00873336876014892
    ## Stox2    -4.4687018565069 7.86957276846453e-06   0.0107970538383333
    ## Slc6a7   -4.4292561869942 9.45586449223749e-06   0.0115319520740887
    ## Plk2     4.39019725362086 1.13247884371376e-05    0.011619191029025

this is for consistent conflict yoked-conflict
----------------------------------------------

    res <- results(dds, contrast =c("APA2", "conflict", "consistent"), independentFiltering = T)
    summary(res)

    ## 
    ## out of 17929 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 8, 0.045%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

    head((res[order(res$padj),]), 10)

    ## log2 fold change (MLE): APA2 conflict vs consistent 
    ## Wald test p-value: APA2 conflict vs consistent 
    ## DataFrame with 10 rows and 6 columns
    ##                       baseMean     log2FoldChange             lfcSE
    ##                      <numeric>          <numeric>         <numeric>
    ## 0610007P14Rik 43.5653432205546  -0.13431443318451 0.381085088890507
    ## 0610009B22Rik  13.157668511865 -0.748166636504363 0.660643972055877
    ## 0610009L18Rik 3.86186075059207  0.565189744483188  1.23716707518475
    ## 0610009O20Rik 48.2470462620772  0.262824867865574 0.402802206893118
    ## 0610010F05Rik 64.5504551982064 -0.361872413635425 0.332203886323633
    ## 0610010K14Rik 19.7011887603264   0.41586826561794  0.68358362003477
    ## 0610012G03Rik 61.3515272894271  0.118168136843711 0.369940445873239
    ## 0610030E20Rik 39.7023920740888 -0.110873849749036 0.507992878780045
    ## 0610037L13Rik 70.3041268921024  0.117013727428709 0.365489732352497
    ## 0610040J01Rik 3.44839145828991 -0.344274193984147  1.26364621724818
    ##                             stat            pvalue              padj
    ##                        <numeric>         <numeric>         <numeric>
    ## 0610007P14Rik -0.352452607304981 0.724498854846802 0.999821112533798
    ## 0610009B22Rik  -1.13248083408091 0.257432337627188 0.999821112533798
    ## 0610009L18Rik  0.456841889684776 0.647784694333414 0.999821112533798
    ## 0610009O20Rik  0.652491131795891 0.514084391925527 0.999821112533798
    ## 0610010F05Rik  -1.08930818853482 0.276018002979688 0.999821112533798
    ## 0610010K14Rik  0.608364878018562  0.54294550025247 0.999821112533798
    ## 0610012G03Rik  0.319424756503109 0.749404440357756 0.999821112533798
    ## 0610030E20Rik -0.218258669324857 0.827227575239532 0.999821112533798
    ## 0610037L13Rik  0.320155990909904 0.748850082871842 0.999821112533798
    ## 0610040J01Rik -0.272445079394031 0.785279814758617 0.999821112533798

    contrast1 <- resvals(contrastvector = c("Punch", "CA1", "DG"), mypval = 0.1) # 2765

    ## [1] 2765

    contrast2 <- resvals(contrastvector = c("Punch", "CA1", "CA3"), mypval = 0.1) # 2165

    ## [1] 2165

    contrast3 <- resvals(contrastvector = c("Punch", "CA3", "DG"), mypval = 0.1) # 2948

    ## [1] 2948

    contrast4 <- resvals(contrastvector = c("APA2", "consistent", "yoked_consistent"), mypval = 0.1) #  114

    ## [1] 114

    contrast5 <- resvals(contrastvector = c("APA2", "conflict", "yoked_conflict"), mypval = 0.1) # 61

    ## [1] 61

    contrast6 <- resvals(contrastvector = c("APA2", "conflict", "consistent"), mypval = 0.1) # 1 or 0

    ## [1] 0

    contrast7 <- resvals(contrastvector = c("APA2", "yoked_conflict", "yoked_consistent"), mypval = 0.1) # 40

    ## [1] 40

    DEGes <- assay(rld)
    DEGes <- cbind(DEGes, contrast1, contrast2, contrast3, contrast4, contrast5, contrast6, contrast7)
    DEGes <- as.data.frame(DEGes) # convert matrix to dataframe
    DEGes$rownames <- rownames(DEGes)  # add the rownames to the dataframe
    DEGes$padjmin <- with(DEGes, pmin(padjPunchCA1DG, padjPunchCA1CA3, padjPunchCA3DG, padjAPA2conflictconsistent, padjAPA2yoked_conflictyoked_consistent,padjAPA2conflictyoked_conflict, padjAPA2consistentyoked_consistent)) 
    DEGes <- DEGes %>% filter(padjmin < 0.000000001)
    rownames(DEGes) <- DEGes$rownames
    drop.cols <-colnames(DEGes[,grep("padj|pval|rownames", colnames(DEGes))])
    DEGes <- DEGes %>% dplyr::select(-one_of(drop.cols))
    DEGes <- as.matrix(DEGes)
    DEGes <- DEGes - rowMeans(DEGes)


    df <- as.data.frame(colData(dds)[,c("APA2", "Punch")]) ## matrix to df
    rownames(df) <- names(countData)
    ann_colors <- ann_colors4 # see color options 
    DEGes <- as.matrix(DEGes) 
    paletteLength <- 30
    myBreaks <- c(seq(min(DEGes), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(DEGes)/paletteLength, max(DEGes), length.out=floor(paletteLength/2)))

    pheatmap(DEGes, show_colnames=T, show_rownames = F,
             annotation_col=df, annotation_colors = ann_colors,
             treeheight_row = 0, treeheight_col = 25,
             annotation_row = NA, 
             annotation_legend = FALSE,
             annotation_names_row = FALSE, annotation_names_col = FALSE,
             fontsize = 8, 
             border_color = "grey60" ,
             color = viridis(30),
             cellwidth = 6, 
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation" 
             )

![](../figures/02b_RNAseqAll/heatmap-1.png)

    pheatmap(DEGes, show_colnames=F, show_rownames = F,
             annotation_col=df, annotation_colors = ann_colors, 
             annotation_row = NA, 
             annotation_legend = FALSE,
             annotation_names_row = FALSE, 
             annotation_names_col = FALSE,
             treeheight_row = 0, treeheight_col = 25,
             fontsize = 8, 
             border_color = "grey60" ,
             color = viridis(30),
             height = 3.75, 
             width = 3.5,
             clustering_method="average",
             breaks=myBreaks,
             clustering_distance_cols="correlation", 
             filename = "../figures/02b_RNAseqALL/pheatmap1.pdf"
             )

Principle component analysis
----------------------------

    ## [1] 54 22  5  3  2  1  1  1  1

    ##             Df Sum Sq Mean Sq F value Pr(>F)    
    ## Punch        2  18619    9309   287.6 <2e-16 ***
    ## Residuals   41   1327      32                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ Punch, data = pcadata)
    ## 
    ## $Punch
    ##               diff       lwr       upr     p adj
    ## CA3-DG  -42.063885 -47.22958 -36.89819 0.0000000
    ## CA1-DG  -43.343230 -48.31529 -38.37117 0.0000000
    ## CA1-CA3  -1.279345  -6.52166   3.96297 0.8244072

    ##             Df Sum Sq Mean Sq F value Pr(>F)    
    ## Punch        2   8050    4025    1049 <2e-16 ***
    ## Residuals   41    157       4                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ Punch, data = pcadata)
    ## 
    ## $Punch
    ##              diff       lwr       upr p adj
    ## CA3-DG  -18.71023 -20.48834 -16.93211     0
    ## CA1-DG   15.28271  13.57125  16.99418     0
    ## CA1-CA3  33.99294  32.18845  35.79743     0

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Punch        2   26.5   13.24   0.306  0.738
    ## Residuals   41 1774.6   43.28

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC3 ~ Punch, data = pcadata)
    ## 
    ## $Punch
    ##               diff       lwr      upr     p adj
    ## CA3-DG  -1.5452815 -7.518690 4.428126 0.8050705
    ## CA1-DG  -1.6651830 -7.414683 4.084317 0.7623489
    ## CA1-CA3 -0.1199015 -6.181910 5.942107 0.9987255

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## APA2         3  255.8   85.25   2.207  0.102
    ## Residuals   40 1545.3   38.63

    ##             Df Sum Sq Mean Sq F value  Pr(>F)    
    ## APA2         3  443.7   147.9   11.83 1.1e-05 ***
    ## Residuals   40  499.9    12.5                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## APA2         3    9.8   3.277    0.24  0.868
    ## Residuals   40  547.0  13.675

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## APA2         3   23.4   7.811   0.771  0.517
    ## Residuals   40  405.5  10.137

![](../figures/02b_RNAseqAll/pca-1.png)

    ## quartz_off_screen 
    ##                 2

Number of differentially expressed genes per two-way contrast
=============================================================

venn diagrams
-------------

    # uses contrast 1-6 from heat map prep

    rldpadjs <- assay(rld)
    rldpadjs <- cbind(rldpadjs, contrast1, contrast2, contrast3, contrast4, contrast5, contrast6)
    rldpadjs <- as.data.frame(rldpadjs)
    rldpadjs <- rldpadjs[ , grepl( "padj" , names( rldpadjs ) ) ]

    venn1 <- row.names(rldpadjs[rldpadjs[1] <0.05 & !is.na(rldpadjs[1]),]) 
    venn2 <- row.names(rldpadjs[rldpadjs[2] <0.05 & !is.na(rldpadjs[2]),]) 
    venn3 <- row.names(rldpadjs[rldpadjs[3] <0.05 & !is.na(rldpadjs[3]),]) 
    venn4 <- row.names(rldpadjs[rldpadjs[4] <0.05 & !is.na(rldpadjs[4]),]) 
    venn5 <- row.names(rldpadjs[rldpadjs[5] <0.05 & !is.na(rldpadjs[5]),]) 
    venn6 <- row.names(rldpadjs[rldpadjs[6] <0.05 & !is.na(rldpadjs[6]),]) 

    candidates1 <- list("CA1-DG" = venn1, "CA3-DG" = venn3, "CA1-CA3" = venn2) 
    candidates2 <- list("Yoked-Consistent" = venn4, "Yoked-Conflict" = venn5, "Yoked-Yoked" = venn6)

    prettyvenn <- venn.diagram(
      scaled=T,
      x = candidates1, filename=NULL, 
      col = "black",
      alpha = 0.5,
      cex = 1, fontfamily = "sans", #fontface = "bold",
      cat.default.pos = "text",
      cat.cex = 1, cat.fontfamily = "sans")
    #dev.off()
    grid.draw(prettyvenn)

![](../figures/02b_RNAseqAll/venndiagram_prep-1.png)

    prettyvenn <- venn.diagram(
      scaled=T,
      x = candidates2, filename=NULL, 
      col = "black",
      alpha = 0.5,
      cex = 1, fontfamily = "sans", #fontface = "bold",
      cat.default.pos = "text",
      cat.cex = 1, cat.fontfamily = "sans")
    #dev.off()
    grid.draw(prettyvenn)

![](../figures/02b_RNAseqAll/venndiagram-1.png)

Volcanos plots and and gene lists
---------------------------------

    makevolcanodf <- function(mycontrast, myup, mydown, filename){
      res <- results(dds, contrast = mycontrast, independentFiltering = T)

      data <- data.frame(gene = row.names(res), pvalue = (res$padj), 
                         lfc = res$log2FoldChange)
      data <- na.omit(data)

      data <- data %>%
      mutate(direction = ifelse(data$lfc > 0 & data$pvalue < 0.05, 
                            yes = myup, 
                            no = ifelse(data$lfc < 0 & data$pvalue < 0.05, 
                                        yes = mydown, 
                                        no = "neither")))
      data$logp <- -log10(data$pvalue)
      data <- dplyr::arrange(data, logp)
      write.csv(data, filename, row.names = F)
      return(data)
    }

    DGvCA3 <- makevolcanodf(c("Punch", "CA3", "DG"), "CA3", "DG", "../data/DGvCA3.csv")
    DGvCA1 <- makevolcanodf(c("Punch", "CA1", "DG"),"CA1", "DG", "../data/DGvCA1.csv")
    CA3vCA1 <- makevolcanodf(c("Punch", "CA1", "CA3"),"CA1", "CA3", "../data/CA3vCA1.csv")

    volcanoplot <- function(mydata, mycolors, xlabname, filename){
      
      myvolcano <- ggplot(mydata, aes(x = lfc, y = logp)) + 
      geom_point(aes(color = direction), size = 0.5, alpha = 0.5, na.rm = T) + 
      scale_color_manual(values = mycolors) + 
      theme_cowplot(font_size = 8, line_size = 0.25) +
      geom_hline(yintercept = 1.3,  size = 0.25, linetype = 2) + 
      scale_y_continuous(limits=c(0, 60)) +
      scale_x_continuous( limits=c(-10, 10)) +
      xlab(paste0(xlabname)) +
      ylab(paste0("log10 p-value")) +       
      theme(panel.grid.minor=element_blank(),
            legend.position = "none", # remove legend 
            panel.grid.major=element_blank())
    myvolcano
    return(myvolcano)

    pdf(file=filename, width=1.75, height=2.25)
    plot(myvolcano)
    dev.off()
    }
      
    volcanoplot(DGvCA3, volcanoDGvCA3, "CA3 / DG", "../figures/02b_RNAseqAll/AllDGCA3.pdf")  

![](../figures/02b_RNAseqAll/volcanoCA3GDG-1.png)

    volcanoplot(DGvCA1, volcanoDGvCA1, "CA1 / DG", "../figures/02b_RNAseqAll/AllCA1DG.pdf")  

![](../figures/02b_RNAseqAll/volcanoCA3GDG-2.png)

    volcanoplot(CA3vCA1, volcanoCA3vCA1, "CA1 / CA3", "../figures/02b_RNAseqAll/AllCA1CA3.pdf") 

![](../figures/02b_RNAseqAll/volcanoCA3GDG-3.png)

plot single gene counts
-----------------------

    plotCounts(dds, "Prkcz", intgroup = "Punch", normalized = TRUE)

![](../figures/02b_RNAseqAll/Prkcz-1.png)

    plotCounts(dds, "Prkcz", intgroup = "APA2", normalized = TRUE)

![](../figures/02b_RNAseqAll/Prkcz-2.png)

Observed versus expected ration of DEGs
---------------------------------------

    # chisq.test equal expression of increased versus decreased expression
    chisq.test(c(1099,  1427), p = c(0.45, 0.55))$expected

    ## [1] 1136.7 1389.3

    chisq.test(c(1099,  1427), p = c(0.45, 0.55))

    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  c(1099, 1427)
    ## X-squared = 2.2734, df = 1, p-value = 0.1316

    prop.table(c(1099,  1427))

    ## [1] 0.4350752 0.5649248

    chisq.test(c(850,   1172), p = c(0.4, 0.6))$expected

    ## [1]  808.8 1213.2

    chisq.test(c(850,   1172), p = c(0.4, 0.6))

    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  c(850, 1172)
    ## X-squared = 3.4979, df = 1, p-value = 0.06145

    prop.table(c(850,   1172))

    ## [1] 0.4203759 0.5796241

    chisq.test(c(1585,  1560), p = c(0.5, 0.5))$expected

    ## [1] 1572.5 1572.5

    chisq.test(c(1585,  1560), p = c(0.5, 0.5))

    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  c(1585, 1560)
    ## X-squared = 0.19873, df = 1, p-value = 0.6557

    prop.table(c(1585,  1560))

    ## [1] 0.5039746 0.4960254

    chisq.test(c(1, 0), p = c(0.5, 0.5))$expected

    ## Warning in chisq.test(c(1, 0), p = c(0.5, 0.5)): Chi-squared approximation
    ## may be incorrect

    ## [1] 0.5 0.5

    chisq.test(c(1, 0), p = c(0.5, 0.5))

    ## Warning in chisq.test(c(1, 0), p = c(0.5, 0.5)): Chi-squared approximation
    ## may be incorrect

    ## 
    ##  Chi-squared test for given probabilities
    ## 
    ## data:  c(1, 0)
    ## X-squared = 1, df = 1, p-value = 0.3173

    prop.table(c(1, 0))

    ## [1] 1 0
