    library(tidyverse) 

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.2.1     ✔ purrr   0.3.2
    ## ✔ tibble  2.1.3     ✔ dplyr   0.8.1
    ## ✔ tidyr   0.8.3     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

    library(corrplot)

    ## corrplot 0.84 loaded

    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    library(corrr)

    knitr::opts_chunk$set(fig.path = '../figures/02e_correlations/', cache = F)

For this analysis, I want to explor correlations between a behavioral
measure and gene expression.

    # import behavior data, create mouse id, select relvant samples
    behav <- read.csv("../data/01a_behavior.csv") 
    behav$mouse <- sapply(strsplit(as.character(behav$ID),"15"), "[", 2)
    behav <- behav %>% filter(treatment %in% c("conflict.trained", "standard.trained"),
                                      TrainSessionCombo == "Retention") %>% 
                               select(mouse,Time1stEntr,pTimeShockZone) 
    head(behav)

    ##   mouse Time1stEntr pTimeShockZone
    ## 1  140A      102.43         0.0286
    ## 2  140C      599.97         0.0021
    ## 3  141C       30.53         0.0909
    ## 4  142A      140.20         0.0203
    ## 5  142C      482.43         0.0445
    ## 6  143A       44.20         0.1122

    pcadata <- read_csv("../data/01a_pcadf.csv") %>%
      filter(treatment %in% c("conflict.trained", "standard.trained"),
             TrainSessionComboNum == 9) %>%
      select(ID,PC1,PC2) 

    ## Parsed with column specification:
    ## cols(
    ##   ID = col_character(),
    ##   treatment = col_character(),
    ##   TrainSessionComboNum = col_double(),
    ##   PC1 = col_double(),
    ##   PC2 = col_double(),
    ##   PC3 = col_double(),
    ##   PC4 = col_double(),
    ##   PC5 = col_double(),
    ##   PC6 = col_double(),
    ##   PC7 = col_double(),
    ##   PC8 = col_double(),
    ##   PC9 = col_double(),
    ##   PC10 = col_double()
    ## )

    pcadata$mouse <- sapply(strsplit(as.character(pcadata$ID),"15"), "[", 2)
    pcadata$ID <- NULL
    head(pcadata)

    ## # A tibble: 6 x 3
    ##      PC1    PC2 mouse
    ##    <dbl>  <dbl> <chr>
    ## 1 -1.97   0.880 140A 
    ## 2 -5.44  -2.13  140C 
    ## 3 -0.199 -0.760 141C 
    ## 4 -2.81   0.425 142A 
    ## 5 -5.10  -0.575 142C 
    ## 6  0.814 -1.66  143A

    DG_DEGs <- read.csv("../data/02f_DG_DEGs_vsd.csv", row.names = 1, check.names = F)
    DG_DEGs <- as.data.frame(t(DG_DEGs))
    DG_DEGs$sample <- row.names(DG_DEGs)
    DG_DEGs$mouse <- sapply(strsplit(as.character(DG_DEGs$sample),"\\-"), "[", 1)
    DG_DEGs <- DG_DEGs %>% select(mouse,`1190002N15RIK`:ZFP869)

    DEGsPCA <- left_join(DG_DEGs, pcadata)

    ## Joining, by = "mouse"

    DEGsPCAbeahv <- left_join(DEGsPCA, behav)

    ## Joining, by = "mouse"

    DEGsPCAbeahv <- as.data.frame(DEGsPCAbeahv)
    row.names(DEGsPCAbeahv) <- DEGsPCAbeahv$mouse
    DEGsPCAbeahv$mouse <- NULL
    DEGsPCAbeahv <- as.matrix(DEGsPCAbeahv)

    DGcor <- DEGsPCAbeahv %>% correlate() %>% rearrange() %>%  shave()

    ## 
    ## Correlation method: 'pearson'
    ## Missing treated using: 'pairwise.complete.obs'

    ## Registered S3 method overwritten by 'seriation':
    ##   method         from 
    ##   reorder.hclust gclus

    DGcor

    ## # A tibble: 218 x 219
    ##    rowname  SNX18  PCDH8  NPAS4  KCNJ2 ERRFI1   SGK1  NFIL3 TIPARP GADD45G
    ##    <chr>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>
    ##  1 SNX18   NA     NA     NA     NA     NA     NA     NA     NA      NA    
    ##  2 PCDH8    0.946 NA     NA     NA     NA     NA     NA     NA      NA    
    ##  3 NPAS4    0.943  0.920 NA     NA     NA     NA     NA     NA      NA    
    ##  4 KCNJ2    0.944  0.927  0.952 NA     NA     NA     NA     NA      NA    
    ##  5 ERRFI1   0.955  0.906  0.974  0.964 NA     NA     NA     NA      NA    
    ##  6 SGK1     0.961  0.920  0.953  0.925  0.940 NA     NA     NA      NA    
    ##  7 NFIL3    0.932  0.902  0.983  0.962  0.961  0.919 NA     NA      NA    
    ##  8 TIPARP   0.935  0.898  0.967  0.924  0.956  0.936  0.959 NA      NA    
    ##  9 GADD45G  0.936  0.932  0.920  0.955  0.947  0.870  0.933  0.904  NA    
    ## 10 PLK2     0.954  0.902  0.973  0.932  0.970  0.977  0.943  0.951   0.885
    ## # … with 208 more rows, and 209 more variables: PLK2 <dbl>, SLC16A1 <dbl>,
    ## #   FAM107B <dbl>, KCNF1 <dbl>, GPR19 <dbl>, KITL <dbl>, FBXO33 <dbl>,
    ## #   MEST <dbl>, PTGS2 <dbl>, APAF1 <dbl>, JUN <dbl>, CWC25 <dbl>,
    ## #   PAK6 <dbl>, FZD5 <dbl>, ZDBF2 <dbl>, SMAD7 <dbl>, CXADR <dbl>,
    ## #   A830010M20RIK <dbl>, TSC22D2 <dbl>, ARID5B <dbl>, EGR4 <dbl>,
    ## #   PIGA <dbl>, KCNA4 <dbl>, IRF2BP2 <dbl>, SIAH2 <dbl>, HMGCR <dbl>,
    ## #   AHR <dbl>, PLK3 <dbl>, RFX2 <dbl>, C2CD4B <dbl>, ARL4A <dbl>,
    ## #   BTG2 <dbl>, ARL13B <dbl>, MARCH11 <dbl>, FZD4 <dbl>, ING2 <dbl>,
    ## #   FOSB <dbl>, RASD1 <dbl>, KLF2 <dbl>, JUNB <dbl>, LONRF1 <dbl>,
    ## #   LEMD3 <dbl>, PPP1R15A <dbl>, DUSP16 <dbl>, FOS <dbl>, NR4A3 <dbl>,
    ## #   CTNND1 <dbl>, MYC <dbl>, EGR3 <dbl>, PEG10 <dbl>, SH2D3C <dbl>,
    ## #   SCG2 <dbl>, SLC2A3 <dbl>, FOXG1 <dbl>, FRMD6 <dbl>, GM13889 <dbl>,
    ## #   FOSL2 <dbl>, ABHD2 <dbl>, ERF <dbl>, LCMT2 <dbl>, MED7 <dbl>,
    ## #   CNNM1 <dbl>, TRIB1 <dbl>, GMEB2 <dbl>, STMN4 <dbl>, DYRK2 <dbl>,
    ## #   ANKRD28 <dbl>, LMNA <dbl>, ODC1 <dbl>, B3GNT2 <dbl>, RASL11A <dbl>,
    ## #   KLF6 <dbl>, LRRTM2 <dbl>, PER1 <dbl>, RGS2 <dbl>, ACAN <dbl>,
    ## #   HECA <dbl>, KDM6B <dbl>, OTUD1 <dbl>, ADAMTS1 <dbl>, FERMT2 <dbl>,
    ## #   ARC <dbl>, EIF5 <dbl>, IRS2 <dbl>, LBH <dbl>, RASL11B <dbl>,
    ## #   IL16 <dbl>, MN1 <dbl>, DUSP14 <dbl>, SLC25A25 <dbl>, ZFP654 <dbl>,
    ## #   NAF1 <dbl>, ZFP275 <dbl>, HS6ST1 <dbl>, NR4A2 <dbl>, ZFP869 <dbl>,
    ## #   `1190002N15RIK` <dbl>, POU3F3 <dbl>, DNAJB1 <dbl>, FBXW7 <dbl>, …

    DGcorSlim <- correlate(DEGsPCAbeahv) %>%  
      focus(PC1,PC2,Time1stEntr, pTimeShockZone)  

    ## 
    ## Correlation method: 'pearson'
    ## Missing treated using: 'pairwise.complete.obs'

    DGcorarranged <- fashion(DGcorSlim) %>% arrange(desc(PC1))
    DGcorarranged

    ##           rowname  PC1  PC2 Time1stEntr pTimeShockZone
    ## 1            NAF1 -.76 -.14         .62           -.74
    ## 2           TNIP2 -.68 -.04         .36           -.73
    ## 3          NAP1L5 -.67 -.04         .43           -.70
    ## 4           ADRB1 -.62 -.51         .53           -.71
    ## 5           NPAS4 -.62 -.51         .53           -.72
    ## 6          ERRFI1 -.58 -.33         .46           -.66
    ## 7           PTGS2 -.58 -.38         .50           -.65
    ## 8          PRUNE2 -.57 -.13         .58           -.50
    ## 9          LRRTM2 -.55 -.39         .63           -.55
    ## 10           PAK6 -.55 -.42         .58           -.58
    ## 11        SPTY2D1 -.54  .03         .33           -.55
    ## 12          FLRT3 -.53 -.47         .47           -.61
    ## 13        GADD45G -.53 -.39         .43           -.63
    ## 14          GPR19 -.53 -.34         .30           -.66
    ## 15           RGS2 -.53 -.40         .38           -.65
    ## 16  1190002N15RIK -.52 -.28         .42           -.57
    ## 17          KCNA4 -.52 -.47         .54           -.58
    ## 18           CCNK -.51 -.19         .50           -.51
    ## 19            LBH -.51  .38         .24           -.45
    ## 20          ABHD2 -.50 -.30         .40           -.56
    ## 21         ARMCX5 -.50 -.55         .43           -.62
    ## 22           MEST -.50  .15         .28           -.49
    ## 23          FOXO1 -.49 -.03         .30           -.52
    ## 24        SLC16A1 -.49 -.60         .60           -.56
    ## 25          GMEB2 -.48  .05         .08           -.57
    ## 26         LONRF1 -.48 -.07         .45           -.46
    ## 27           SYT4 -.48 -.12         .24           -.55
    ## 28           KLF6 -.47 -.48         .35           -.59
    ## 29           SCG2 -.47 -.31         .25           -.59
    ## 30          USPL1 -.47 -.29         .47           -.49
    ## 31  A830010M20RIK -.46 -.07         .22           -.53
    ## 32          CNNM1 -.46 -.12         .24           -.52
    ## 33          FBXW7 -.46 -.18         .22           -.54
    ## 34          CXADR -.45  .22         .03           -.50
    ## 35           FZD5 -.45 -.30         .37           -.50
    ## 36          PCDH8 -.45 -.22         .20           -.57
    ## 37         ARPP21 -.44  .19         .06           -.49
    ## 38         CITED2 -.44 -.61         .51           -.51
    ## 39        MARCH11 -.44  .19         .19           -.44
    ## 40           PIGA -.44 -.18         .23           -.54
    ## 41            AHR -.43 -.07         .16           -.51
    ## 42           ING2 -.43 -.16         .10           -.56
    ## 43         RANBP2 -.43  .13         .19           -.45
    ## 44          RASD1 -.43 -.42         .33           -.53
    ## 45          STAC2 -.43 -.35         .62           -.38
    ## 46          DUSP8 -.42 -.06         .28           -.44
    ## 47           PLK2 -.42 -.29         .22           -.54
    ## 48         SH2D3C -.42 -.04         .16           -.50
    ## 49        TMEM170 -.42  .26         .26           -.35
    ## 50         ZFAND5 -.42 -.08         .10           -.52
    ## 51           CUL3 -.41 -.06         .21           -.46
    ## 52          NFIL3 -.41 -.55         .44           -.50
    ## 53         ZBTB33 -.41  .10         .01           -.50
    ## 54          ARL4A -.40  .05         .15           -.43
    ## 55          FOSL2 -.40 -.10         .23           -.46
    ## 56         B3GNT2 -.39 -.30         .07           -.55
    ## 57         DNAJB4 -.39  .00         .10           -.47
    ## 58          ARL4D -.38 -.34         .39           -.42
    ## 59           ATF3 -.38 -.17         .20           -.46
    ## 60          CPEB4 -.38  .14         .25           -.35
    ## 61           EPRS -.38  .11         .09           -.42
    ## 62          PELI1 -.38  .16         .15           -.40
    ## 63            SRF -.38  .01         .03           -.49
    ## 64           BDNF -.37 -.20         .17           -.48
    ## 65          KCNJ2 -.37 -.37         .40           -.43
    ## 66         TSPYL3 -.37 -.22         .10           -.47
    ## 67         AMIGO2 -.36 -.45         .55           -.37
    ## 68       ANKRD33B -.36 -.24         .35           -.39
    ## 69          APAF1 -.36  .02         .21           -.37
    ## 70          ARL5B -.36  .12         .00           -.44
    ## 71         TIPARP -.36 -.16         .19           -.43
    ## 72         CLDN12 -.35  .09         .08           -.40
    ## 73          CWC25 -.35 -.32         .18           -.47
    ## 74         DUSP16 -.35 -.13         .22           -.40
    ## 75         SLC2A3 -.35  .43         .11           -.29
    ## 76         FBXO33 -.34 -.02         .07           -.42
    ## 77         ZFP869 -.34  .01         .59           -.18
    ## 78          NR4A2 -.33 -.10         .12           -.41
    ## 79           RGMB -.33 -.30         .37           -.37
    ## 80         ZFP668 -.33 -.00         .05           -.40
    ## 81            ARC -.32 -.30         .29           -.38
    ## 82           BTG2 -.32 -.09         .14           -.39
    ## 83          HMGCR -.32  .17        -.01           -.39
    ## 84           KLF2 -.32 -.33         .16           -.44
    ## 85         MFAP3L -.32 -.63         .57           -.36
    ## 86       SLC25A25 -.32 -.44         .20           -.46
    ## 87        GM13889 -.31 -.16         .05           -.42
    ## 88         ARL13B -.30  .16         .05           -.33
    ## 89          NPTX2 -.30 -.05         .04           -.39
    ## 90           PER1 -.30  .09         .01           -.37
    ## 91           RGS4 -.30 -.25         .10           -.43
    ## 92          SMAD7 -.30 -.24         .08           -.43
    ## 93           EIF5 -.29  .34         .06           -.26
    ## 94        PRPF38B -.29  .38         .01           -.27
    ## 95          FRMD6 -.28  .10         .01           -.33
    ## 96       PPP1R15A -.28 -.38         .18           -.39
    ## 97            UBC -.28 -.15         .10           -.36
    ## 98           FOSB -.27  .03         .01           -.35
    ## 99           JUNB -.27 -.07         .01           -.37
    ## 100        LYSMD4 -.27  .57        -.06           -.19
    ## 101         SIAH2 -.27 -.08        -.05           -.40
    ## 102        ZFP275 -.27  .18         .14           -.26
    ## 103          IRS1 -.26 -.09         .23           -.27
    ## 104          SCOC -.26  .22         .28           -.13
    ## 105         TRIB1 -.26 -.30         .16           -.36
    ## 106         NR4A3 -.24  .09         .01           -.29
    ## 107        ZFP654 -.24  .64         .06           -.11
    ## 108       ADAMTS1 -.23 -.28         .32           -.25
    ## 109          EGR3 -.23  .23         .02           -.22
    ## 110         CYP51 -.22  .11        -.01           -.27
    ## 111           JUN -.22 -.09         .11           -.27
    ## 112          SGK1 -.22 -.14        -.03           -.34
    ## 113        HOMER1 -.21  .26        -.04           -.22
    ## 114       SLITRK5 -.21  .14         .01           -.25
    ## 115        DBPHT2 -.20  .30        -.10           -.22
    ## 116          EGR4 -.20  .11        -.09           -.27
    ## 117           MN1 -.20 -.45         .52           -.17
    ## 118         DUSP4 -.19 -.18         .01           -.30
    ## 119         FOXG1 -.18 -.05         .17           -.17
    ## 120         PEG10 -.18 -.15         .30           -.15
    ## 121      ANKRD13A -.17 -.18         .07           -.25
    ## 122          FZD4 -.17 -.15         .31           -.14
    ## 123          IRS2 -.17  .04        -.06           -.24
    ## 124          JUND -.17 -.09        -.10           -.29
    ## 125         KCNF1 -.17 -.35         .03           -.32
    ## 126         SENP8 -.17 -.30         .20           -.19
    ## 127        C2CD4B -.16 -.29         .21           -.22
    ## 128        DNAJA1 -.16  .24         .08           -.11
    ## 129         LCMT2 -.16  .10        -.20           -.26
    ## 130         TRA2B -.16  .15        -.23           -.27
    ## 131          BMT2 -.15  .31        -.19           -.18
    ## 132         DUSP6 -.15 -.17        -.04           -.26
    ## 133         HSPH1 -.15  .24         .08           -.11
    ## 134         NR4A1 -.15 -.10        -.02           -.24
    ## 135         ZDBF2 -.15 -.30         .35           -.13
    ## 136         DYRK2 -.14 -.68         .43           -.20
    ## 137       OLFML2B -.14  .10         .27           -.03
    ## 138          RFX2 -.14  .21        -.17           -.21
    ## 139         STMN4 -.14 -.40         .17           -.22
    ## 140        DNAJB1 -.13  .47        -.11           -.09
    ## 141       IRF2BP2 -.13 -.20         .01           -.23
    ## 142          KITL -.12  .17        -.10           -.16
    ## 143         NUAK1 -.12 -.10         .28           -.08
    ## 144         SNX18 -.12  .05        -.15           -.22
    ## 145         CIART -.11  .35        -.26           -.15
    ## 146          EGR1 -.10  .12        -.10           -.14
    ## 147         CTCFL -.09  .10        -.20           -.16
    ## 148           ERF -.09  .08        -.16           -.17
    ## 149         LEMD3 -.09  .33        -.29           -.16
    ## 150         RBM47 -.09 -.25         .29           -.05
    ## 151        FERMT2 -.08 -.03         .05           -.11
    ## 152         KLKB1 -.08  .35        -.32           -.13
    ## 153          NEFM -.08 -.39        -.01           -.22
    ## 154       TSC22D2 -.08  .02        -.13           -.15
    ## 155       RASL11B -.07 -.34        -.01           -.18
    ## 156        CTNND1 -.06  .28        -.27           -.13
    ## 157        ARID5B -.05 -.25         .13           -.09
    ## 158       FAM107B -.05  .04        -.28           -.18
    ## 159           FOS -.05  .01        -.12           -.13
    ## 160        SLC5A5 -.05  .19         .20            .06
    ## 161         CPNE7 -.03  .08         .09            .02
    ## 162        HSPA1A -.03 -.07         .38            .10
    ## 163        JMJD1C -.03  .31        -.10           -.02
    ## 164        DUSP14 -.02 -.44        -.21           -.23
    ## 165          IL16 -.02  .18        -.40           -.14
    ## 166       NEUROD6 -.02  .49        -.13            .06
    ## 167          ACAN -.01 -.19        -.04           -.09
    ## 168         PGAM2 -.00 -.22         .41            .11
    ## 169          SV2B -.00 -.37         .23           -.02
    ## 170          MC1R  .77  .30        -.86            .72
    ## 171        TUBB4A  .70 -.09        -.53            .64
    ## 172          COQ2  .55  .58        -.75            .53
    ## 173          PLK3  .52  .11        -.53            .47
    ## 174          JDP2  .50 -.41        -.30            .40
    ## 175         CECR6  .42 -.22        -.06            .46
    ## 176       SLC45A4  .42  .04        -.04            .53
    ## 177        CCDC32  .32  .05        -.38            .27
    ## 178        PHLDA1  .31 -.40        -.26            .18
    ## 179          MED7  .30  .42        -.51            .27
    ## 180          ODC1  .29  .33        -.57            .22
    ## 181           PXN  .28  .53        -.22            .40
    ## 182         PLCH2  .27 -.67         .11            .20
    ## 183           GYG  .26  .46        -.62            .20
    ## 184        HS6ST1  .26  .02        -.23            .25
    ## 185       RASL11A  .22 -.06        -.20            .20
    ## 186        DPYSL2  .19  .44        -.41            .19
    ## 187          PER2  .19  .27        -.42            .13
    ## 188        ZFP207  .18  .37         .08            .36
    ## 189       IGF2BP2  .17  .17        -.37            .11
    ## 190         BACH1  .15 -.36         .02            .10
    ## 191         OTUD1  .15 -.03         .04            .19
    ## 192         PDE6A  .15 -.01        -.33            .06
    ## 193          GAD1  .14 -.07        -.30            .03
    ## 194          GPI1  .14 -.01        -.12            .16
    ## 195          SOX9  .14 -.07         .17            .22
    ## 196       ANKRD28  .13  .28        -.31            .09
    ## 197          NXF1  .13  .52        -.22            .22
    ## 198       FAM118A  .12  .42        -.49            .06
    ## 199        POU3F3  .12  .14        -.44            .02
    ## 200        SRGAP1  .12  .01        -.05            .13
    ## 201        EEF1E1  .11 -.32         .34            .19
    ## 202        SOWAHC  .11 -.34         .22            .14
    ## 203           MYC  .08 -.30        -.25           -.08
    ## 204         NEDD9  .08 -.41         .22            .06
    ## 205      BC048403  .07 -.16         .24            .16
    ## 206          GNAZ  .07  .20        -.35            .00
    ## 207          LMNA  .06  .22        -.42           -.04
    ## 208       ANKRD27  .04  .54        -.13            .12
    ## 209        PLAGL1  .04  .38        -.30            .02
    ## 210         THBS1  .04  .03         .15            .11
    ## 211          HECA  .03  .54        -.25            .09
    ## 212         KDM6B  .03  .36        -.38           -.04
    ## 213        LRRC45  .02  .34        -.11            .07
    ## 214         KDM7A  .00  .50        -.36           -.01

    DEGsPCAbeahvDf <- as.data.frame(DEGsPCAbeahv)
    a <- ggplot(DEGsPCAbeahvDf, aes(x = PC1, y = NAF1, label = rownames(DEGsPCAbeahvDf))) +
      geom_point(colour = "darkred") + geom_smooth(method='lm', colour = "darkgrey") + geom_text(vjust = -0.1)
    b <- ggplot(DEGsPCAbeahvDf, aes(x = Time1stEntr, y = NAF1, label = rownames(DEGsPCAbeahvDf))) +
      geom_point(colour = "darkred") + geom_smooth(method='lm', colour = "darkgrey") + geom_text(vjust = -0.1)

    c <- ggplot(DEGsPCAbeahvDf, aes(x = PC1, y = ARC, label = rownames(DEGsPCAbeahvDf))) +
      geom_point(colour = "darkred") + geom_smooth(method='lm', colour = "darkgrey") + geom_text(vjust = -0.1)
    d <- ggplot(DEGsPCAbeahvDf, aes(x = Time1stEntr, y = ARC, label = rownames(DEGsPCAbeahvDf))) +
      geom_point(colour = "darkred") + geom_smooth(method='lm', colour = "darkgrey") + geom_text(vjust = -0.1)


    a

    ## Warning: Removed 8 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 8 rows containing missing values (geom_point).

    ## Warning: Removed 8 rows containing missing values (geom_text).

![](../figures/02e_correlations/unnamed-chunk-1-1.png)

    b

    ## Warning: Removed 8 rows containing non-finite values (stat_smooth).

    ## Warning: Removed 8 rows containing missing values (geom_point).

    ## Warning: Removed 8 rows containing missing values (geom_text).

![](../figures/02e_correlations/unnamed-chunk-1-2.png)
