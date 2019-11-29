    library(tidyverse)

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✔ ggplot2 3.2.1     ✔ purrr   0.3.3
    ## ✔ tibble  2.1.3     ✔ dplyr   0.8.3
    ## ✔ tidyr   1.0.0     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

    library(corrr)
    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    source("./figureoptions.R")
    source("./functions_RNAseq.R")

    # read the sample data, set levels, join iwth behvior PCA data
    colData <- read.csv("../data/02a_colData.csv", row.names = 1, stringsAsFactors = T)
    colData <- colData %>% filter(subfield == "DG")
    pca.Rn <- read_csv("../data/01a_pca.all.csv") %>% filter(trialNum == 9)

    ## Parsed with column specification:
    ## cols(
    ##   ID = col_character(),
    ##   treatment = col_character(),
    ##   trialNum = col_double(),
    ##   Day = col_double(),
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

    pca.Rn <- pca.Rn %>% select(ID:PC2)
    colData <- left_join(colData, pca.Rn)

    ## Joining, by = c("ID", "treatment")

    ## Warning: Column `ID` joining factor and character vector, coercing into
    ## character vector

    ## Warning: Column `treatment` joining factor and character vector, coercing
    ## into character vector

    head(colData)

    ##       ID subfield        treatment training trialNum Day        PC1
    ## 1 15143A       DG conflict.trained  trained        9   3 -0.2275039
    ## 2 15143B       DG   conflict.yoked    yoked        9   3 -3.1436627
    ## 3 15143D       DG   standard.yoked    yoked        9   3 -2.7532719
    ## 4 15144A       DG conflict.trained  trained        9   3  6.7041815
    ## 5 15144C       DG standard.trained  trained        9   4  7.0499369
    ## 6 15144D       DG   standard.yoked    yoked        9   3 -3.3026284
    ##           PC2
    ## 1  3.03543738
    ## 2 -0.48834291
    ## 3 -0.07584809
    ## 4 -0.07853719
    ## 5 -1.78499206
    ## 6  1.17314374

    # read all count data
    vsd <- read.csv("../data/02c_DGvsd.csv", row.names = 1, check.names = F) 
    vsd$gene <- row.names(vsd)
    vsd$gene <- toupper(vsd$gene)
    vsd <- as.data.frame(vsd)
    row.names(vsd) <- vsd$gene
    vsd$gene <- NULL
    head(vsd)

    ##               143A-DG-1 143B-DG-1 143D-DG-3 144A-DG-2 144C-DG-2 144D-DG-2
    ## 0610007P14RIK  6.381145  6.405587  6.834813  6.513173  6.401366  6.637581
    ## 0610009B22RIK  5.781124  5.661825  5.299160  5.528657  5.821796  5.670620
    ## 0610009L18RIK  5.405848  5.593055  5.422004  5.000789  5.507750  5.506319
    ## 0610009O20RIK  6.736671  6.691641  6.736719  6.530386  6.799892  6.642508
    ## 0610010F05RIK  6.603589  7.093384  7.135233  6.931921  6.907778  6.587132
    ## 0610010K14RIK  6.136536  5.949913  6.015194  5.782488  5.999853  6.084773
    ##               145A-DG-2 145B-DG-1 146A-DG-2 146B-DG-2 146C-DG-4 146D-DG-3
    ## 0610007P14RIK  6.562833  6.126760  6.298333  5.996243  6.929322  7.411893
    ## 0610009B22RIK  5.942028  5.787857  5.679794  5.996243  6.060620  5.000789
    ## 0610009L18RIK  5.256496  5.000789  5.679794  5.000789  5.000789  5.000789
    ## 0610009O20RIK  6.739906  6.623046  6.747840  6.383172  6.350262  6.793559
    ## 0610010F05RIK  6.739906  6.540236  6.765447  7.085062  7.264347  6.493162
    ## 0610010K14RIK  6.218650  6.200426  6.298333  7.252267  5.971751  5.000789
    ##               147C-DG-3 147D-DG-1 148A-DG-3 148B-DG-4
    ## 0610007P14RIK  6.298159  6.481476  6.493581  6.349976
    ## 0610009B22RIK  5.910805  5.893873  5.589261  5.492981
    ## 0610009L18RIK  5.257050  5.718223  5.505647  5.349661
    ## 0610009O20RIK  6.687072  6.631927  7.108382  6.425888
    ## 0610010F05RIK  6.880090  6.759923  6.892155  6.595557
    ## 0610010K14RIK  5.761732  6.136024  5.626697  6.081280

    # prep to join with sample colData
    vsd <- as.data.frame(t(vsd))
    vsd$sample <- row.names(vsd)
    vsd$mouse <- sapply(strsplit(as.character(vsd$sample),"\\-"), "[", 1)
    vsd$ID <- paste(15, vsd$mouse, sep = "")

    vsd <- left_join(colData, vsd)

    ## Joining, by = "ID"

    head(vsd)[17015:17020]

    ##       ZXDC   ZYG11B      ZYX    ZZEF1     ZZZ3    sample
    ## 1 6.579863 8.319205 7.203704 7.829641 7.742970 143A-DG-1
    ## 2 6.626521 8.118246 7.258360 7.711706 7.569642 143B-DG-1
    ## 3 6.695140 8.715286 6.853523 7.987434 7.135233 143D-DG-3
    ## 4 6.816776 9.110958 7.532625 8.240355 7.760152 144A-DG-2
    ## 5 7.219436 8.724875 7.128516 7.955995 7.634300 144C-DG-2
    ## 6 6.617663 8.814969 6.671645 7.908857 7.526267 144D-DG-2

    forcorall <-  vsd %>% select(PC1:ZZZ3)
    corrrmat <- correlate(forcorall, diagonal = 1) 

    ## 
    ## Correlation method: 'pearson'
    ## Missing treated using: 'pairwise.complete.obs'

    head(corrrmat)

    ## # A tibble: 6 x 17,014
    ##   rowname     PC1     PC2 `0610007P14RIK` `0610009B22RIK` `0610009L18RIK`
    ##   <chr>     <dbl>   <dbl>           <dbl>           <dbl>           <dbl>
    ## 1 PC1      1      -0.606          -0.0930          0.324          -0.137 
    ## 2 PC2     -0.606   1               0.0167         -0.200           0.0170
    ## 3 061000… -0.0930  0.0167          1              -0.595          -0.173 
    ## 4 061000…  0.324  -0.200          -0.595           1               0.0266
    ## 5 061000… -0.137   0.0170         -0.173           0.0266          1     
    ## 6 061000…  0.0393 -0.163           0.164          -0.342           0.465 
    ## # … with 17,008 more variables: `0610009O20RIK` <dbl>,
    ## #   `0610010F05RIK` <dbl>, `0610010K14RIK` <dbl>, `0610012G03RIK` <dbl>,
    ## #   `0610030E20RIK` <dbl>, `0610037L13RIK` <dbl>, `0610040J01RIK` <dbl>,
    ## #   `1110002E22RIK` <dbl>, `1110004E09RIK` <dbl>, `1110004F10RIK` <dbl>,
    ## #   `1110008F13RIK` <dbl>, `1110008L16RIK` <dbl>, `1110008P14RIK` <dbl>,
    ## #   `1110012L19RIK` <dbl>, `1110017D15RIK` <dbl>, `1110032A03RIK` <dbl>,
    ## #   `1110032F04RIK` <dbl>, `1110034G24RIK` <dbl>, `1110037F02RIK` <dbl>,
    ## #   `1110038F14RIK` <dbl>, `1110051M20RIK` <dbl>, `1110059E24RIK` <dbl>,
    ## #   `1110059G10RIK` <dbl>, `1110065P20RIK` <dbl>, `1190002N15RIK` <dbl>,
    ## #   `1190005I06RIK` <dbl>, `1190007I07RIK` <dbl>, `1300017J02RIK` <dbl>,
    ## #   `1500009C09RIK` <dbl>, `1500009L16RIK` <dbl>, `1500011B03RIK` <dbl>,
    ## #   `1500011K16RIK` <dbl>, `1500015O10RIK` <dbl>, `1520401A03RIK` <dbl>,
    ## #   `1600002H07RIK` <dbl>, `1600002K03RIK` <dbl>, `1600012H06RIK` <dbl>,
    ## #   `1600014C10RIK` <dbl>, `1700001C19RIK` <dbl>, `1700001K19RIK` <dbl>,
    ## #   `1700001L19RIK` <dbl>, `1700001O22RIK` <dbl>, `1700001P01RIK` <dbl>,
    ## #   `1700003E16RIK` <dbl>, `1700003F12RIK` <dbl>, `1700006E09RIK` <dbl>,
    ## #   `1700007G11RIK` <dbl>, `1700007K13RIK` <dbl>, `1700008O03RIK` <dbl>,
    ## #   `1700010I14RIK` <dbl>, `1700011E24RIK` <dbl>, `1700011H14RIK` <dbl>,
    ## #   `1700011M02RIK` <dbl>, `1700013D24RIK` <dbl>, `1700013F07RIK` <dbl>,
    ## #   `1700014D04RIK` <dbl>, `1700015E13RIK` <dbl>, `1700016D06RIK` <dbl>,
    ## #   `1700016K19RIK` <dbl>, `1700017B05RIK` <dbl>, `1700019A02RIK` <dbl>,
    ## #   `1700019B03RIK` <dbl>, `1700019D03RIK` <dbl>, `1700019O17RIK` <dbl>,
    ## #   `1700020A23RIK` <dbl>, `1700020D05RIK` <dbl>, `1700020L24RIK` <dbl>,
    ## #   `1700020N01RIK` <dbl>, `1700021F05RIK` <dbl>, `1700022I11RIK` <dbl>,
    ## #   `1700023F06RIK` <dbl>, `1700024P16RIK` <dbl>, `1700025G04RIK` <dbl>,
    ## #   `1700027J19RIK` <dbl>, `1700028J19RIK` <dbl>, `1700028K03RIK` <dbl>,
    ## #   `1700028P14RIK` <dbl>, `1700029F12RIK` <dbl>, `1700029H14RIK` <dbl>,
    ## #   `1700029I15RIK` <dbl>, `1700029J07RIK` <dbl>, `1700030J22RIK` <dbl>,
    ## #   `1700030K09RIK` <dbl>, `1700034I23RIK` <dbl>, `1700034J05RIK` <dbl>,
    ## #   `1700037C18RIK` <dbl>, `1700037H04RIK` <dbl>, `1700040L02RIK` <dbl>,
    ## #   `1700047I17RIK2` <dbl>, `1700048O20RIK` <dbl>, `1700057G04RIK` <dbl>,
    ## #   `1700061G19RIK` <dbl>, `1700064H15RIK` <dbl>, `1700066B19RIK` <dbl>,
    ## #   `1700066M21RIK` <dbl>, `1700067K01RIK` <dbl>, `1700088E04RIK` <dbl>,
    ## #   `1700092M07RIK` <dbl>, `1700093K21RIK` <dbl>, `1700102P08RIK` <dbl>, …

    summary(lm( PC1 ~ ARC, data = vsd))

    ## 
    ## Call:
    ## lm(formula = PC1 ~ ARC, data = vsd)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -4.6903 -1.7333  0.2229  1.2413  5.5850 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -25.9124     5.3643  -4.831 0.000267 ***
    ## ARC           2.8553     0.5619   5.081 0.000167 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.565 on 14 degrees of freedom
    ## Multiple R-squared:  0.6484, Adjusted R-squared:  0.6233 
    ## F-statistic: 25.82 on 1 and 14 DF,  p-value: 0.0001674

    summary(lm( PC2 ~ ARC, data = vsd))

    ## 
    ## Call:
    ## lm(formula = PC2 ~ ARC, data = vsd)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.0654 -0.8975 -0.2179  0.5500  3.2926 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)   9.9971     2.9862   3.348  0.00479 **
    ## ARC          -0.9639     0.3128  -3.081  0.00813 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.428 on 14 degrees of freedom
    ## Multiple R-squared:  0.4041, Adjusted R-squared:  0.3616 
    ## F-statistic: 9.495 on 1 and 14 DF,  p-value: 0.008128

    cor.test(vsd$PC1, vsd$ARC, method = c("pearson"))

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  vsd$PC1 and vsd$ARC
    ## t = 5.0812, df = 14, p-value = 0.0001674
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.5151616 0.9298016
    ## sample estimates:
    ##       cor 
    ## 0.8052353

    cor.test(vsd$PC2, vsd$ARC, method = c("pearson"))

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  vsd$PC2 and vsd$ARC
    ## t = -3.0814, df = 14, p-value = 0.008128
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.8603085 -0.2044175
    ## sample estimates:
    ##        cor 
    ## -0.6357062

    a <- ggplot(vsd, aes(x = ARC, y = PC1)) +
       geom_point(aes( color = treatment)) + 
       geom_smooth(method = "lm", color = "grey") +
       scale_color_manual(values = treatmentcolors) +
      theme_ms() +
       theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank())  +
      labs(subtitle = "R2 = 0.81, p = 0.0002")
    b <- ggplot(vsd, aes(x = ARC, y = PC2)) +
       geom_point(aes( color = treatment)) + 
       geom_smooth(method = "lm", color = "grey") +
       scale_color_manual(values = treatmentcolors) +
      theme_ms() +
       theme(legend.position = "none")  +
      labs(subtitle = "R2 = -0.63, p = 0.009 ")

    left <- plot_grid(a,b, nrow = 2, labels = c("a","b"), label_size = 8)
    left

![](../figures/02d_correlations/ARC-1.png)

    corrsTop10 <- corrrmat %>% 
      focus(PC1, PC2)  %>% 
      arrange(desc(PC1)) %>% 
      filter(PC1 > 0.78) 
    corrsTop10

    ## # A tibble: 11 x 3
    ##    rowname    PC1    PC2
    ##    <chr>    <dbl>  <dbl>
    ##  1 NAF1     0.862 -0.699
    ##  2 PTGS2    0.839 -0.645
    ##  3 RGS2     0.834 -0.611
    ##  4 HIST1H1D 0.815 -0.544
    ##  5 COL10A1  0.813 -0.579
    ##  6 ARC      0.805 -0.636
    ##  7 HSPB3    0.803 -0.286
    ##  8 NPAS4    0.799 -0.641
    ##  9 FZD5     0.798 -0.670
    ## 10 ACAN     0.785 -0.546
    ## 11 AREG     0.784 -0.438

    topcorrrs <- corrsTop10$rowname

    p1 <- corrrmat %>% 
      focus(PC1, PC2, topcorrrs,  mirror = TRUE) %>% 
       focus(PC1, PC2) %>%
      mutate(rowname = reorder(rowname, PC1)) %>%
      ggplot(aes(rowname, PC1, fill = PC1)) +
        geom_col() + coord_flip() +
      scale_fill_gradient2(low = "#67a9cf",  high = "#ef8a62", midpoint = 0) +
      theme_ms() +
      theme(legend.position = "none") +
      #ylim(-1,1) +
      labs(x = NULL, y = "Correlation to PC1") +
      labs(subtitle = " ")
    p1 

![](../figures/02d_correlations/corrr-1.png)

    p2 <- corrrmat %>% 
      focus(PC1, PC2, topcorrrs,  mirror = TRUE) %>% 
      rearrange() %>% 
      network_plot(colors = c("#67a9cf", "white", "#ef8a62"),
                   min_cor = .6, curved = T, legend = T,
                   repel = TRUE) + 
      theme(legend.position = "bottom") +
      labs(subtitle = " ")

    ## Registered S3 method overwritten by 'seriation':
    ##   method         from 
    ##   reorder.hclust gclus

    p2 

![](../figures/02d_correlations/corrr-2.png)

    right <- plot_grid(p1,p2,  nrow = 1, labels = c("c","d"), label_size = 8)
    right

![](../figures/02d_correlations/corrr-3.png)

    fig4 <- plot_grid(left, right, nrow = 1, rel_widths = c(1,2))
    fig4

![](../figures/02d_correlations/correlations-1.png)

    pdf(file="../figures/02d_correlations/correlations.pdf", width=6.69, height=3.5)
    plot(fig4)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf(file="../figures/figure_4.pdf", width=6.69, height=3.5)
    plot(fig4)
    dev.off()

    ## quartz_off_screen 
    ##                 2

What genes genes correlated with PC1? How many gene correlated with PC1 and are differentially epxressed in DG?
---------------------------------------------------------------------------------------------------------------

    # 57 genes are correlated with PC1 > 0.7
    corrsPC1sig <- corrrmat %>% 
      focus(PC1)  %>% 
      arrange(desc(PC1)) %>% 
      filter(PC1 > .7) %>%
       mutate_each(funs=toupper) %>%
      mutate(gene = rowname) %>% 
      arrange(gene) %>% 
      select(gene,PC1)   
    corrsPC1sig$gene

    ##  [1] "ABRA"     "ACAN"     "ADAMTS1"  "AIM1L"    "AMIGO2"   "ARC"     
    ##  [7] "AREG"     "ARL4D"    "ARL5B"    "ARMCX5"   "ATF3"     "BDNF"    
    ## [13] "BTG2"     "CHST14"   "COL10A1"  "CPEB4"    "DUSP8"    "EGR1"    
    ## [19] "EGR4"     "ERRFI1"   "FBXO33"   "FOSL2"    "FRMD6"    "FZD5"    
    ## [25] "GM10269"  "HIST1H1D" "HIST1H3I" "HOMER1"   "HOXC4"    "HSPB3"   
    ## [31] "IGHD"     "KCNK10"   "MEST"     "NAF1"     "NEXN"     "NFIL3"   
    ## [37] "NPAS4"    "NR4A3"    "PCDH8"    "PELI1"    "PER1"     "PLK2"    
    ## [43] "PPIC"     "PTGS2"    "RASD1"    "RGS2"     "SGK1"     "SLC16A1" 
    ## [49] "SLC25A25" "SMAD7"    "SPTY2D1"  "SYT4"     "TIPARP"   "TNFRSF23"
    ## [55] "TRIB1"    "UBC"      "ZFP804B"

    # 214 DEGs in DG
    # read DG DEGs
    DG_DEGs <- read.csv("../data/02c_DEGs.DG.training.yoked.trained.csv", row.names = 1, check.names = F) 
    DG_DEGs$gene <- toupper(DG_DEGs$gene)
    dim(DG_DEGs)

    ## [1] 214   5

    # 41 found  in both correlate and deg datasets
    PC1corrsDEGS <- inner_join(corrsPC1sig, DG_DEGs) %>% 
      arrange(gene) %>% select(gene,PC1)   

    ## Joining, by = "gene"

    length(PC1corrsDEGS$gene)

    ## [1] 41

    PC1corrsDEGS$gene

    ##  [1] "ACAN"     "ADAMTS1"  "AMIGO2"   "ARC"      "ARL4D"    "ARL5B"   
    ##  [7] "ARMCX5"   "ATF3"     "BDNF"     "BTG2"     "CPEB4"    "DUSP8"   
    ## [13] "EGR1"     "EGR4"     "ERRFI1"   "FBXO33"   "FOSL2"    "FRMD6"   
    ## [19] "FZD5"     "HOMER1"   "MEST"     "NAF1"     "NFIL3"    "NPAS4"   
    ## [25] "NR4A3"    "PCDH8"    "PELI1"    "PER1"     "PLK2"     "PTGS2"   
    ## [31] "RASD1"    "RGS2"     "SGK1"     "SLC16A1"  "SLC25A25" "SMAD7"   
    ## [37] "SPTY2D1"  "SYT4"     "TIPARP"   "TRIB1"    "UBC"

    citation("corrr")

    ## 
    ## To cite package 'corrr' in publications use:
    ## 
    ##   Edgar Ruiz, Simon Jackson and Jorge Cimentada (2019). corrr:
    ##   Correlations in R. R package version 0.4.0.
    ##   https://CRAN.R-project.org/package=corrr
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {corrr: Correlations in R},
    ##     author = {Edgar Ruiz and Simon Jackson and Jorge Cimentada},
    ##     year = {2019},
    ##     note = {R package version 0.4.0},
    ##     url = {https://CRAN.R-project.org/package=corrr},
    ##   }

    write.csv(corrsTop10, "../data/02d_corrsTop10.csv", row.names = F)
    write.csv(corrsPC1sig, "../data/02d_corrsPC1sig.csv", row.names = F)
