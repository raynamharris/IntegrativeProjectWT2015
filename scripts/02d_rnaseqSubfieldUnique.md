    DG <- read.csv("../data/02c_DGforupset.csv")  
    CA1 <- read.csv("../data/02c_CA1forupset.csv")  
    CA3 <- read.csv("../data/02c_CA3forupset.csv") 

    # make df for upset plots without direction
    myupsetdf <- rbind(DG,CA1,CA3)

    myupsetdf <- myupsetdf %>%
      mutate(comparison = fct_recode(comparison,
        "ConfT.ConfY"    = "conflict.trained-conflict.yoked",
        "ConfY.StdY"      = "conflict.yoked-standard.yoked",
        "StdT.StdY" = "standard.trained-standard.yoked",
        "HC.ConfT" = "home.cage-conflict.trained" ,
        "HC.ConfY" =  "home.cage-conflict.yoked" ,
        "HC.StdT" = "home.cage-standard.trained" ,
        "HC.StdY" =  "home.cage-standard.yoked")) %>%
      mutate(tissue.contrast = paste(tissue, comparison, sep = ".")) %>%
      select(gene,tissue.contrast) %>%
      mutate(yesno = 1) %>%
      distinct %>%
      spread(tissue.contrast, yesno, fill = 0)

    write.csv(myupsetdf, "../data/02d_upsetdf.csv")

    upset1 <- upset(myupsetdf, keep.order = F,
         sets.bar.color=c(
           "#d95f02","#d95f02","#d95f02","#d95f02",
           "#7570b3","#7570b3","#7570b3","#7570b3","#7570b3",
           "#1b9e77",
           "#7570b3",
           "#1b9e77",
           "#d95f02",
           "#1b9e77"),
          nsets = 14,
          order.by = "freq",
          sets.x.label = NULL,
          point.size = 1.5, 
          mb.ratio = c(0.5, 0.5))
    upset1

![](../figures/02c_rnaseqSubfield/upsetplot-1.png)

    pdf(file="../figures/02c_rnaseqSubfield/upsetall.pdf",  onefile=FALSE, width=6.69, height=5.1) # or other device
    upset1
    dev.off()

    ## quartz_off_screen 
    ##                 2

    DGupset <- myupsetdf %>% select(gene, starts_with("DG")) 
    CA1upset <- myupsetdf %>% select(gene, starts_with("CA1")) 
    CA3upset <- myupsetdf %>% select(gene, starts_with("CA3")) 

    upset2 <- upset(DGupset, keep.order = F,
          order.by = "freq",
          sets.bar.color=c("#d95f02"),
          sets.x.label = NULL,
          point.size = 1.5, 
          mb.ratio = c(0.5, 0.5),
          nsets = 6)
    upset2

![](../figures/02c_rnaseqSubfield/upsetplot-2.png)

    upset3 <- upset(CA1upset, keep.order = F,
          order.by = "freq",
          sets.bar.color=c("#7570b3"),
          sets.x.label = NULL,
          point.size = 1.5, 
          mb.ratio = c(0.5, 0.5),
          nsets = 6)
    upset3

![](../figures/02c_rnaseqSubfield/upsetplot-3.png)

    upset4 <- upset(CA3upset, keep.order = F,
          order.by = "freq",
          sets.bar.color=c("#1b9e77"),
          sets.x.label = NULL,
          point.size = 1.5, 
          mb.ratio = c(0.5, 0.5),
          nsets = 5)
    upset4

![](../figures/02c_rnaseqSubfield/upsetplot-4.png)

Create a list of genes afected by stress and learning. Save. THen use to
filter out nonspecific gene expression responses

    # DG learning but not homecage or yoked yoked
    DGtrainingonlygenes <- DGupset %>% 
      filter(DG.StdT.StdY == 1 &  DG.ConfY.StdY == 0 &DG.HC.ConfT == 0 & DG.HC.ConfY == 0 & DG.HC.StdT == 0 & DG.HC.StdY == 0 ) %>% 
      select(gene)  %>% droplevels()

    print("DG training but not homecage or yoked")

    ## [1] "DG training but not homecage or yoked"

    as.list(DGtrainingonlygenes)

    ## $gene
    ##  [1] Dnaja1  Dnajb1  Dusp16  Fos     Gm1818  Hes5    Ier3    Lbh    
    ##  [9] Lrrc32  Rasl11a Rtl1    Siah2   Slc16a1 Thbs1   Tiparp  Zfand5 
    ## 16 Levels: Dnaja1 Dnajb1 Dusp16 Fos Gm1818 Hes5 Ier3 Lbh ... Zfand5

    # CA1 learning  but not homecage or yoked yoked
    CA1trainingonlygenes <- CA1upset %>% 
      filter(CA1.StdT.StdY == 1 &  CA1.ConfY.StdY == 0 & CA1.HC.ConfT == 0 & CA1.HC.ConfY == 0 & CA1.HC.StdT == 0 & CA1.HC.StdY == 0 ) %>% 
      select(gene) %>% droplevels()

    print("CA1 training but not homecage or yoked")

    ## [1] "CA1 training but not homecage or yoked"

    as.list(CA1trainingonlygenes)

    ## $gene
    ##  [1] Entpd7   Fanci    Grin2b   Hr       Kcna4    Prrt2    Rab1b   
    ##  [8] Rbak     Sptb     Abhd14b  Abl1     Alms1    Amotl2   Ankrd26 
    ## [15] Ap1s1    Atf6b    Atp6v1b2 Atp6v1g2 Bahcc1   Cic      Colgalt2
    ## [22] Cpne6    Dcdc2b   Diaph2   Emc4     Enkd1    Erich3   Fam131a 
    ## [29] Fam171a2 Fam19a2  Fbxo32   Fndc3a   Foxj3    Gm21685  Gm21887 
    ## [36] Gm4631   Gm9821   Grik3    Grm1     H2afv    Hcn1     Itgb1   
    ## [43] Kmt5a    Larp4b   Lhx6     Lonrf2   Mcc      Mettl3   Ncam2   
    ## [50] Ovca2    Pitpnc1  Ppil2    Pura     Rabep2   Rdh1     Rnf25   
    ## [57] Samd1    Sh3rf1   Slc35b4  Slc50a1  Sowahc   Sphkap   Srrd    
    ## [64] Stox2    Tbc1d30  Tbc1d9   Tmem175  Tmem8b   Tubgcp5  Uqcrh   
    ## [71] Usmg5    Zfp395  
    ## 72 Levels: Entpd7 Fanci Grin2b Hr Kcna4 Prrt2 Rab1b Rbak Sptb ... Zfp395

    # CA3 learning  but not homecage or yoked yoked
    CA3trainingonlygenes <- CA3upset %>% 
      filter(CA3.StdT.StdY == 1  ) %>% 
      select(gene)  %>% droplevels()

    print("CA3 training but not homecage or yoked")

    ## [1] "CA3 training but not homecage or yoked"

    as.list(CA3trainingonlygenes)

    ## $gene
    ## [1] Sco2
    ## Levels: Sco2
