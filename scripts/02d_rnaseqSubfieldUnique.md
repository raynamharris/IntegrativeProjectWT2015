    DG <- read.csv("../data/02c_DGforupset.csv", stringsAsFactors = F)  
    CA1 <- read.csv("../data/02c_CA1forupset.csv", stringsAsFactors = F)  
    CA3 <- read.csv("../data/02c_CA3forupset.csv", stringsAsFactors = F) 

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
      filter(DG.StdT.StdY == 1 &  DG.ConfY.StdY == 0 &DG.HC.ConfT == 0 & 
               DG.HC.ConfY == 0 & DG.HC.StdT == 0 & DG.HC.StdY == 0 ) %>% 
      select(gene)  %>% arrange(gene) %>% droplevels()

    print("DG training but not homecage or yoked")

    ## [1] "DG training but not homecage or yoked"

    as.list(DGtrainingonlygenes)

    ## $gene
    ##  [1] "Dnaja1"  "Dnajb1"  "Dusp16"  "Fos"     "Gm1818"  "Hes5"    "Ier3"   
    ##  [8] "Lbh"     "Lrrc32"  "Rasl11a" "Rtl1"    "Siah2"   "Slc16a1" "Thbs1"  
    ## [15] "Tiparp"  "Zfand5"

    # CA1 learning  but not homecage or yoked yoked
    CA1trainingonlygenes <- CA1upset %>% 
      filter(CA1.StdT.StdY == 1 &  CA1.ConfY.StdY == 0 & CA1.HC.ConfT == 0 & 
               CA1.HC.ConfY == 0 & CA1.HC.StdT == 0 & CA1.HC.StdY == 0 ) %>% 
      select(gene) %>% arrange(gene) %>% droplevels()

    print("CA1 training but not homecage or yoked")

    ## [1] "CA1 training but not homecage or yoked"

    as.list(CA1trainingonlygenes)

    ## $gene
    ##  [1] "Abhd14b"  "Abl1"     "Alms1"    "Amotl2"   "Ankrd26"  "Ap1s1"   
    ##  [7] "Atf6b"    "Atp6v1b2" "Atp6v1g2" "Bahcc1"   "Cic"      "Colgalt2"
    ## [13] "Cpne6"    "Dcdc2b"   "Diaph2"   "Emc4"     "Enkd1"    "Entpd7"  
    ## [19] "Erich3"   "Fam131a"  "Fam171a2" "Fam19a2"  "Fanci"    "Fbxo32"  
    ## [25] "Fndc3a"   "Foxj3"    "Gm21685"  "Gm21887"  "Gm4631"   "Gm9821"  
    ## [31] "Grik3"    "Grin2b"   "Grm1"     "H2afv"    "Hcn1"     "Hr"      
    ## [37] "Itgb1"    "Kcna4"    "Kmt5a"    "Larp4b"   "Lhx6"     "Lonrf2"  
    ## [43] "Mcc"      "Mettl3"   "Ncam2"    "Ovca2"    "Pitpnc1"  "Ppil2"   
    ## [49] "Prrt2"    "Pura"     "Rab1b"    "Rabep2"   "Rbak"     "Rdh1"    
    ## [55] "Rnf25"    "Samd1"    "Sh3rf1"   "Slc35b4"  "Slc50a1"  "Sowahc"  
    ## [61] "Sphkap"   "Sptb"     "Srrd"     "Stox2"    "Tbc1d30"  "Tbc1d9"  
    ## [67] "Tmem175"  "Tmem8b"   "Tubgcp5"  "Uqcrh"    "Usmg5"    "Zfp395"

    str(CA1trainingonlygenes)

    ## 'data.frame':    72 obs. of  1 variable:
    ##  $ gene: chr  "Abhd14b" "Abl1" "Alms1" "Amotl2" ...

    # CA3 learning  but not homecage or yoked yoked
    CA3trainingonlygenes <- CA3upset %>% 
      filter(CA3.StdT.StdY == 1  ) %>% 
      select(gene)  %>% droplevels()

    print("CA3 training but not homecage or yoked")

    ## [1] "CA3 training but not homecage or yoked"

    as.list(CA3trainingonlygenes)

    ## $gene
    ## [1] "Sco2"

    write.csv(CA1trainingonlygenes, "../data/02d_CA1_trainingonlygenes.csv", row.names = F)
    write.csv(DGtrainingonlygenes, "../data/02d_DG_trainingonlygenes.csv", row.names = F)

    # to use GORILLA, I need a list of all genes used
    genesforgo <- read_csv("../data/00_geneids.csv") %>% 
      select(gene) %>% arrange(gene)

    ## Parsed with column specification:
    ## cols(
    ##   id = col_character(),
    ##   ENSMUST = col_character(),
    ##   ENSMUSG = col_character(),
    ##   OTTMUSG = col_character(),
    ##   OTTMUST = col_character(),
    ##   transcript = col_character(),
    ##   gene = col_character(),
    ##   length = col_double(),
    ##   structure1 = col_character(),
    ##   structure2 = col_character(),
    ##   structure3 = col_character(),
    ##   transcript_lenght = col_character()
    ## )

    write.csv(genesforgo, "../data/02d_genesforgo.csv")

    DGwide <- DG %>%
      select(-padj,-tissue) %>%
      spread(key = comparison, value = lfc)  %>%
      select(gene, "home.cage-standard.yoked", "home.cage-standard.trained",  "home.cage-conflict.yoked", "home.cage-conflict.trained",
             "standard.trained-standard.yoked", "conflict.yoked-standard.yoked", "conflict.trained-conflict.yoked")
    head(DGwide)

    ##            gene home.cage-standard.yoked home.cage-standard.trained
    ## 1 0610009B22Rik                 3.372990                         NA
    ## 2 0610010F05Rik                       NA                         NA
    ## 3 1110008F13Rik                 1.920996                   1.835162
    ## 4 1110032F04Rik                       NA                         NA
    ## 5 1110051M20Rik                 1.370361                         NA
    ## 6 1500011B03Rik                       NA                   1.852466
    ##   home.cage-conflict.yoked home.cage-conflict.trained
    ## 1                 3.280862                   3.216008
    ## 2                       NA                   1.468638
    ## 3                 1.917032                   2.233990
    ## 4                -4.458811                         NA
    ## 5                       NA                   1.210130
    ## 6                 1.949101                   1.936894
    ##   standard.trained-standard.yoked conflict.yoked-standard.yoked
    ## 1                              NA                            NA
    ## 2                              NA                            NA
    ## 3                              NA                            NA
    ## 4                              NA                            NA
    ## 5                              NA                            NA
    ## 6                              NA                            NA
    ##   conflict.trained-conflict.yoked
    ## 1                              NA
    ## 2                              NA
    ## 3                              NA
    ## 4                              NA
    ## 5                              NA
    ## 6                              NA
