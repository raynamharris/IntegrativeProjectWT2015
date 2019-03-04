Genes from Cembrowski sublement file 1
<a href="https://elifesciences.org/articles/14997/figures" class="uri">https://elifesciences.org/articles/14997/figures</a>

    # import cembrowski markers
    cembrowskisupp <- read.table("../data/cembrowksi_markers/elife-14997-supp1-v1_dendrogram.txt", sep="\t", header = T)

    # select just columns with gene symbol and enriched column
    # then rename gene column and convert to uppercase
    cembrowskisupp <- cembrowskisupp %>% dplyr::select(gene_short_name, enriched) 
    colnames(cembrowskisupp)[1] <- "gene"
    cembrowskisupp$gene <- str_to_upper(cembrowskisupp$gene)
    head(cembrowskisupp)

    ##       gene  enriched
    ## 1   ABLIM3 dg_d-dg_v
    ## 2    AKAP7 dg_d-dg_v
    ## 3 ARHGAP20 dg_d-dg_v
    ## 4     BTG2 dg_d-dg_v
    ## 5    C1QL2 dg_d-dg_v
    ## 6    CALD1 dg_d-dg_v

    # subset my maker

    cembrowksimarkers <- function(subfields){
      mydf <- cembrowskisupp %>% 
      dplyr::filter(enriched %in% subfields) %>% 
      select(gene)
      names(mydf) <- NULL
      mylist <- as.list(mydf[,1])
      return(mylist)
    }

    levels(cembrowskisupp$enriched)

    ##  [1] "ca1_d"                           "ca1_d-ca1_v"                    
    ##  [3] "ca1_v"                           "ca2"                            
    ##  [5] "ca3_d"                           "ca3_d-ca3_v-ca2-ca1_d-ca1_v"    
    ##  [7] "ca3_v"                           "ca4"                            
    ##  [9] "ca4-ca3_d-ca3_v-ca2-ca1_d-ca1_v" "dg_d"                           
    ## [11] "dg_d-dg_v"                       "dg_v"

    CA1_markers <- cembrowksimarkers(c("ca1_d", "ca1_v", "ca1_d-ca1_v"))
    DG_markers <- cembrowksimarkers(c("dg_d", "dg_v", "dg_d-dg_v"))
    CA3_markers <- cembrowksimarkers(c("ca3_d", "ca3_v", "ca3_d-ca3_v-ca2-ca1_d-ca1_v"))
    CA2_markers <- cembrowksimarkers(c("ca2"))
    CA4_markers <- cembrowksimarkers(c("ca4"))



    # make data frames of genes expression results for markers 
    marker_expression <- function(mydf, subfield, markers){
        MARKERS <- str_to_upper(markers)
        df <- mydf %>%
        dplyr::filter(gene %in% c(MARKERS)) %>%
        dplyr::mutate(marker = subfield) %>%
        droplevels()
        return(df)
    }

    # import subfield specific data
    wrangledata <- function(filename, mycomparison){
      mydata <- read.csv(filename, header = T)
      mydata$gene <- str_to_upper(mydata$gene)
      mydata$comparison <- mycomparison
      return(mydata)
    }

    CA1DG <- wrangledata("../data/DGvCA1.csv", "CA1-DG")
    CA1CA3 <- wrangledata("../data/CA3vCA1.csv", "CA1-CA3")
    CA3DG <- wrangledata("../data/DGvCA3.csv", "CA3-DG")

    #look for markers in each supfield

    #CA1 markers
    CA1DG_CA1 <- marker_expression(CA1DG, "CA1", CA1_markers)
    kable(CA1DG_CA1)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
direction
</th>
<th style="text-align:right;">
logp
</th>
<th style="text-align:left;">
comparison
</th>
<th style="text-align:left;">
marker
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
TOX
</td>
<td style="text-align:right;">
0.7774696
</td>
<td style="text-align:right;">
-1.1358768
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1093166
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
NDRG1
</td>
<td style="text-align:right;">
0.6335538
</td>
<td style="text-align:right;">
-0.6418147
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1982165
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
ZFP386
</td>
<td style="text-align:right;">
0.5158355
</td>
<td style="text-align:right;">
-0.7714939
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2874888
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
ADRA1A
</td>
<td style="text-align:right;">
0.5069893
</td>
<td style="text-align:right;">
-1.7962602
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2950012
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
ENPP2
</td>
<td style="text-align:right;">
0.4513398
</td>
<td style="text-align:right;">
0.7348180
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3454963
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
GRIN3A
</td>
<td style="text-align:right;">
0.2125827
</td>
<td style="text-align:right;">
1.3625266
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.6724721
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
SERPINE2
</td>
<td style="text-align:right;">
0.1410322
</td>
<td style="text-align:right;">
1.1432638
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.8506817
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
KLHL13
</td>
<td style="text-align:right;">
0.0834265
</td>
<td style="text-align:right;">
1.7412683
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.0786962
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
SATB2
</td>
<td style="text-align:right;">
0.0087898
</td>
<td style="text-align:right;">
5.1229973
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
2.0560189
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNAB1
</td>
<td style="text-align:right;">
0.0036131
</td>
<td style="text-align:right;">
1.4613243
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
2.4421177
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
LDB2
</td>
<td style="text-align:right;">
0.0030478
</td>
<td style="text-align:right;">
5.5764104
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
2.5160157
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
NDST3
</td>
<td style="text-align:right;">
0.0001106
</td>
<td style="text-align:right;">
3.0985396
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
3.9561942
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNH7
</td>
<td style="text-align:right;">
0.0000688
</td>
<td style="text-align:right;">
2.7532453
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
4.1623781
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
SCN3B
</td>
<td style="text-align:right;">
0.0000474
</td>
<td style="text-align:right;">
1.3634918
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
4.3244785
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
CDS1
</td>
<td style="text-align:right;">
0.0000018
</td>
<td style="text-align:right;">
3.0993865
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
5.7445323
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
PLEKHG1
</td>
<td style="text-align:right;">
0.0000001
</td>
<td style="text-align:right;">
3.5829567
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
6.8560137
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
FIBCD1
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
4.7255644
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
7.7169185
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
BC030500
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
7.1974313
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
8.6233201
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
MPPED1
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
6.1974917
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
9.1570208
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
NOV
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
7.3443639
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
10.9216106
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
ZDHHC2
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
3.2394718
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
11.9127303
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
MAN1A
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
5.7125382
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
17.4127358
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
PEX5L
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
3.7556525
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
18.1430635
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
WFS1
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
6.4172072
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
19.7643792
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
POU3F1
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
5.9881158
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
24.8300562
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
</tbody>
</table>

    summary(CA1DG_CA1$direction)

    ##     CA1 neither 
    ##      17       8

    CA1CA3_CA1 <- marker_expression(CA1CA3, "CA1", CA1_markers)
    kable(CA1CA3_CA1)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
direction
</th>
<th style="text-align:right;">
logp
</th>
<th style="text-align:left;">
comparison
</th>
<th style="text-align:left;">
marker
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
ZFP386
</td>
<td style="text-align:right;">
0.9339140
</td>
<td style="text-align:right;">
-0.1362492
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0296931
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
LDB2
</td>
<td style="text-align:right;">
0.8002605
</td>
<td style="text-align:right;">
0.6585487
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0967686
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
TOX
</td>
<td style="text-align:right;">
0.7629620
</td>
<td style="text-align:right;">
-1.3036435
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1174971
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
GRIN3A
</td>
<td style="text-align:right;">
0.4880650
</td>
<td style="text-align:right;">
0.8692593
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3115224
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
ADRA1A
</td>
<td style="text-align:right;">
0.3169317
</td>
<td style="text-align:right;">
-2.6677237
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4990344
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
KLHL13
</td>
<td style="text-align:right;">
0.2707597
</td>
<td style="text-align:right;">
1.2432127
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.5674160
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
SERPINE2
</td>
<td style="text-align:right;">
0.2447905
</td>
<td style="text-align:right;">
-0.9847938
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.6112055
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
ENPP2
</td>
<td style="text-align:right;">
0.0426756
</td>
<td style="text-align:right;">
-1.7808215
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.3698208
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
NDRG1
</td>
<td style="text-align:right;">
0.0148826
</td>
<td style="text-align:right;">
-2.6728472
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.8273206
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
PLEKHG1
</td>
<td style="text-align:right;">
0.0039268
</td>
<td style="text-align:right;">
2.1941445
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
2.4059591
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
BC030500
</td>
<td style="text-align:right;">
0.0029966
</td>
<td style="text-align:right;">
3.7100773
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
2.5233757
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
SATB2
</td>
<td style="text-align:right;">
0.0011222
</td>
<td style="text-align:right;">
7.6268474
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
2.9499240
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
MAN1A
</td>
<td style="text-align:right;">
0.0006248
</td>
<td style="text-align:right;">
2.4528810
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
3.2042362
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
ZDHHC2
</td>
<td style="text-align:right;">
0.0001704
</td>
<td style="text-align:right;">
1.8920977
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
3.7684152
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
CDS1
</td>
<td style="text-align:right;">
0.0001444
</td>
<td style="text-align:right;">
2.6390858
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
3.8404995
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
MPPED1
</td>
<td style="text-align:right;">
0.0000509
</td>
<td style="text-align:right;">
4.4723893
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
4.2929340
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
NDST3
</td>
<td style="text-align:right;">
0.0000006
</td>
<td style="text-align:right;">
4.1565003
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
6.2031257
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNH7
</td>
<td style="text-align:right;">
0.0000001
</td>
<td style="text-align:right;">
3.6692709
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
6.8557475
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
SCN3B
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
1.9792790
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
8.5076320
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNAB1
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
2.8158528
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
8.5397888
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
NOV
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
6.8411388
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
8.8840650
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
PEX5L
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
2.8804017
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
9.9564733
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
WFS1
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
5.4268652
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
13.9291120
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
FIBCD1
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
7.6909143
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
18.1148001
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
POU3F1
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
6.5277691
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
28.2478728
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
</tbody>
</table>

    summary(CA1CA3_CA1$direction)

    ##     CA1     CA3 neither 
    ##      16       2       7

    CA3DG_CA1 <- marker_expression(CA3DG, "CA1", CA1_markers)
    kable(CA3DG_CA1)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
direction
</th>
<th style="text-align:right;">
logp
</th>
<th style="text-align:left;">
comparison
</th>
<th style="text-align:left;">
marker
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
TOX
</td>
<td style="text-align:right;">
0.9628007
</td>
<td style="text-align:right;">
0.1677667
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0164636
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
NOV
</td>
<td style="text-align:right;">
0.8042800
</td>
<td style="text-align:right;">
0.5032251
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0945927
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
GRIN3A
</td>
<td style="text-align:right;">
0.6895603
</td>
<td style="text-align:right;">
0.4932673
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1614277
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
KLHL13
</td>
<td style="text-align:right;">
0.6804242
</td>
<td style="text-align:right;">
0.4980556
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1672203
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
ADRA1A
</td>
<td style="text-align:right;">
0.6539837
</td>
<td style="text-align:right;">
0.8714635
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1844331
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
CDS1
</td>
<td style="text-align:right;">
0.6165617
</td>
<td style="text-align:right;">
0.4603008
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2100235
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
POU3F1
</td>
<td style="text-align:right;">
0.5652385
</td>
<td style="text-align:right;">
-0.5396532
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2477683
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
DCN
</td>
<td style="text-align:right;">
0.5454322
</td>
<td style="text-align:right;">
2.5826405
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2632592
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
ZFP386
</td>
<td style="text-align:right;">
0.5341102
</td>
<td style="text-align:right;">
-0.6352448
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2723691
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
AGRP
</td>
<td style="text-align:right;">
0.4894116
</td>
<td style="text-align:right;">
-2.2469865
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3103257
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
SATB2
</td>
<td style="text-align:right;">
0.3602730
</td>
<td style="text-align:right;">
-2.5038501
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4433683
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
NDST3
</td>
<td style="text-align:right;">
0.2909527
</td>
<td style="text-align:right;">
-1.0579606
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.5361776
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
WFS1
</td>
<td style="text-align:right;">
0.2655729
</td>
<td style="text-align:right;">
0.9903421
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.5758162
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNH7
</td>
<td style="text-align:right;">
0.2553987
</td>
<td style="text-align:right;">
-0.9160256
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.5927814
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
MPPED1
</td>
<td style="text-align:right;">
0.1333669
</td>
<td style="text-align:right;">
1.7251024
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.8749521
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
SCN3B
</td>
<td style="text-align:right;">
0.0696512
</td>
<td style="text-align:right;">
-0.6157871
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.1570714
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
PEX5L
</td>
<td style="text-align:right;">
0.0649024
</td>
<td style="text-align:right;">
0.8752508
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.1877392
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
PLEKHG1
</td>
<td style="text-align:right;">
0.0518044
</td>
<td style="text-align:right;">
1.3888122
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.2856336
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
NDRG1
</td>
<td style="text-align:right;">
0.0085225
</td>
<td style="text-align:right;">
2.0310326
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
2.0694324
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
BC030500
</td>
<td style="text-align:right;">
0.0049207
</td>
<td style="text-align:right;">
3.4873540
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
2.3079733
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
ZDHHC2
</td>
<td style="text-align:right;">
0.0028177
</td>
<td style="text-align:right;">
1.3473741
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
2.5501101
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
LDB2
</td>
<td style="text-align:right;">
0.0023833
</td>
<td style="text-align:right;">
4.9178617
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
2.6228197
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNAB1
</td>
<td style="text-align:right;">
0.0016449
</td>
<td style="text-align:right;">
-1.3545285
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.7838705
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
FIBCD1
</td>
<td style="text-align:right;">
0.0004483
</td>
<td style="text-align:right;">
-2.9653499
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.3484577
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
SERPINE2
</td>
<td style="text-align:right;">
0.0001147
</td>
<td style="text-align:right;">
2.1280576
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
3.9404271
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
ENPP2
</td>
<td style="text-align:right;">
0.0000324
</td>
<td style="text-align:right;">
2.5156396
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
4.4899446
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
<tr>
<td style="text-align:left;">
MAN1A
</td>
<td style="text-align:right;">
0.0000001
</td>
<td style="text-align:right;">
3.2596572
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
6.8823158
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA1
</td>
</tr>
</tbody>
</table>

    summary(CA3DG_CA1$direction)

    ##     CA3      DG neither 
    ##       7       2      18

    #CA3 markers
    CA1DG_CA3 <- marker_expression(CA1DG, "CA3", CA3_markers)
    kable(CA1DG_CA3)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
direction
</th>
<th style="text-align:right;">
logp
</th>
<th style="text-align:left;">
comparison
</th>
<th style="text-align:left;">
marker
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
PDLIM1
</td>
<td style="text-align:right;">
0.9732320
</td>
<td style="text-align:right;">
0.1303192
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0117836
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
PRKCD
</td>
<td style="text-align:right;">
0.9721015
</td>
<td style="text-align:right;">
0.0766420
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0122884
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC1A2
</td>
<td style="text-align:right;">
0.9449609
</td>
<td style="text-align:right;">
0.0662632
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0245862
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
CD109
</td>
<td style="text-align:right;">
0.8043630
</td>
<td style="text-align:right;">
-0.9243913
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0945479
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
TRHDE
</td>
<td style="text-align:right;">
0.4054731
</td>
<td style="text-align:right;">
-1.0898797
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3920380
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
CHGB
</td>
<td style="text-align:right;">
0.1672312
</td>
<td style="text-align:right;">
0.7591654
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.7766827
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
IYD
</td>
<td style="text-align:right;">
0.1233424
</td>
<td style="text-align:right;">
-3.9572575
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.9088878
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
PTGS2
</td>
<td style="text-align:right;">
0.0370060
</td>
<td style="text-align:right;">
-1.8753955
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
1.4317278
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
EPHB1
</td>
<td style="text-align:right;">
0.0239441
</td>
<td style="text-align:right;">
-2.3850739
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
1.6208022
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
OCIAD2
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
3.7728251
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
13.8753838
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
</tbody>
</table>

    summary(CA1DG_CA3$direction)

    ##     CA1      DG neither 
    ##       1       2       7

    CA1CA3_CA3 <- marker_expression(CA1CA3, "CA3", CA3_markers)
    kable(CA1CA3_CA3)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
direction
</th>
<th style="text-align:right;">
logp
</th>
<th style="text-align:left;">
comparison
</th>
<th style="text-align:left;">
marker
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
OCIAD2
</td>
<td style="text-align:right;">
0.6987759
</td>
<td style="text-align:right;">
0.3071165
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1556621
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
PDLIM1
</td>
<td style="text-align:right;">
0.5304873
</td>
<td style="text-align:right;">
-1.8920813
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2753250
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC1A2
</td>
<td style="text-align:right;">
0.0886139
</td>
<td style="text-align:right;">
-1.1731458
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.0524979
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
PRKCD
</td>
<td style="text-align:right;">
0.0694934
</td>
<td style="text-align:right;">
-2.5972483
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.1580563
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
CD109
</td>
<td style="text-align:right;">
0.0352274
</td>
<td style="text-align:right;">
-5.9177574
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.4531190
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
IYD
</td>
<td style="text-align:right;">
0.0150522
</td>
<td style="text-align:right;">
-6.0481411
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.8223991
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
PTGS2
</td>
<td style="text-align:right;">
0.0015335
</td>
<td style="text-align:right;">
-2.7737245
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
2.8143214
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
CHGB
</td>
<td style="text-align:right;">
0.0000609
</td>
<td style="text-align:right;">
-1.9726655
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
4.2154143
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
EPHB1
</td>
<td style="text-align:right;">
0.0000401
</td>
<td style="text-align:right;">
-4.0152213
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
4.3965974
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
TRHDE
</td>
<td style="text-align:right;">
0.0000042
</td>
<td style="text-align:right;">
-4.6453464
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
5.3794538
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
</tbody>
</table>

    summary(CA1CA3_CA3$direction)

    ##     CA3 neither 
    ##       6       4

    CA3DG_CA3 <- marker_expression(CA3DG, "CA3", CA3_markers)
    kable(CA3DG_CA3)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
direction
</th>
<th style="text-align:right;">
logp
</th>
<th style="text-align:left;">
comparison
</th>
<th style="text-align:left;">
marker
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
PDLIM1
</td>
<td style="text-align:right;">
0.3388145
</td>
<td style="text-align:right;">
2.0224005
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4700380
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
PTGS2
</td>
<td style="text-align:right;">
0.2274971
</td>
<td style="text-align:right;">
0.8983291
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.6430242
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
IYD
</td>
<td style="text-align:right;">
0.1480578
</td>
<td style="text-align:right;">
2.0908836
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.8295686
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC1A2
</td>
<td style="text-align:right;">
0.0153756
</td>
<td style="text-align:right;">
1.2394090
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.8131674
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
PRKCD
</td>
<td style="text-align:right;">
0.0063352
</td>
<td style="text-align:right;">
2.6738903
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
2.1982396
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
EPHB1
</td>
<td style="text-align:right;">
0.0023402
</td>
<td style="text-align:right;">
1.6301474
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
2.6307432
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
CD109
</td>
<td style="text-align:right;">
0.0020406
</td>
<td style="text-align:right;">
4.9933661
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
2.6902488
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
TRHDE
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
3.5554668
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
8.6831707
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
CHGB
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
2.7318309
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
13.1585646
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
<tr>
<td style="text-align:left;">
OCIAD2
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
3.4657086
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
15.4400102
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA3
</td>
</tr>
</tbody>
</table>

    summary(CA3DG_CA3$direction)

    ##     CA3 neither 
    ##       7       3

    #DG markers
    CA1DG_DG <- marker_expression(CA1DG, "DG", DG_markers)
    kable(CA1DG_DG)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
direction
</th>
<th style="text-align:right;">
logp
</th>
<th style="text-align:left;">
comparison
</th>
<th style="text-align:left;">
marker
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
IL33
</td>
<td style="text-align:right;">
0.9689680
</td>
<td style="text-align:right;">
-0.0680812
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0136905
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NRGN
</td>
<td style="text-align:right;">
0.9386389
</td>
<td style="text-align:right;">
0.0428860
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0275014
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
BTG2
</td>
<td style="text-align:right;">
0.9108047
</td>
<td style="text-align:right;">
0.1347767
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0405747
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NCOR2
</td>
<td style="text-align:right;">
0.8713718
</td>
<td style="text-align:right;">
0.0756319
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0597965
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ZCCHC12
</td>
<td style="text-align:right;">
0.8636423
</td>
<td style="text-align:right;">
-0.2876254
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0636661
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
RAVER1
</td>
<td style="text-align:right;">
0.8213156
</td>
<td style="text-align:right;">
0.1473491
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0854899
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
EXT2
</td>
<td style="text-align:right;">
0.7990373
</td>
<td style="text-align:right;">
0.1984583
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0974329
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GSE1
</td>
<td style="text-align:right;">
0.7928806
</td>
<td style="text-align:right;">
-0.2011861
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1007922
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
EGR4
</td>
<td style="text-align:right;">
0.7773827
</td>
<td style="text-align:right;">
0.3359632
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1093651
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
B4GALT3
</td>
<td style="text-align:right;">
0.7411251
</td>
<td style="text-align:right;">
-0.2589210
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1301085
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
H2AFZ
</td>
<td style="text-align:right;">
0.6679720
</td>
<td style="text-align:right;">
-0.2532493
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1752418
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CACNG7
</td>
<td style="text-align:right;">
0.6335538
</td>
<td style="text-align:right;">
-0.3220101
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1982165
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
STK32B
</td>
<td style="text-align:right;">
0.5664586
</td>
<td style="text-align:right;">
-1.7028908
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2468318
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TAOK2
</td>
<td style="text-align:right;">
0.5290231
</td>
<td style="text-align:right;">
-0.2756038
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2765253
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GTF2A2
</td>
<td style="text-align:right;">
0.5267071
</td>
<td style="text-align:right;">
0.6390556
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2784308
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GLP2R
</td>
<td style="text-align:right;">
0.5139021
</td>
<td style="text-align:right;">
-1.5200803
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2891196
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PDZRN3
</td>
<td style="text-align:right;">
0.5006098
</td>
<td style="text-align:right;">
-2.0826302
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3005006
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
FAT4
</td>
<td style="text-align:right;">
0.4890830
</td>
<td style="text-align:right;">
-0.9369753
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3106174
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SH3BP1
</td>
<td style="text-align:right;">
0.4807497
</td>
<td style="text-align:right;">
-0.2437067
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3180810
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TRHR
</td>
<td style="text-align:right;">
0.4608419
</td>
<td style="text-align:right;">
-2.2442757
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3364481
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ETS2
</td>
<td style="text-align:right;">
0.4438707
</td>
<td style="text-align:right;">
-0.4344350
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3527435
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
LCORL
</td>
<td style="text-align:right;">
0.3433642
</td>
<td style="text-align:right;">
-1.4233329
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4642450
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TOX3
</td>
<td style="text-align:right;">
0.3169527
</td>
<td style="text-align:right;">
1.8921544
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4990055
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CHUK
</td>
<td style="text-align:right;">
0.2801177
</td>
<td style="text-align:right;">
0.6051016
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.5526595
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ECM2
</td>
<td style="text-align:right;">
0.2788556
</td>
<td style="text-align:right;">
-2.1379750
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.5546207
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NAPA
</td>
<td style="text-align:right;">
0.2508142
</td>
<td style="text-align:right;">
-0.4664596
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.6006479
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
FBLN1
</td>
<td style="text-align:right;">
0.1922964
</td>
<td style="text-align:right;">
-3.5919721
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.7160288
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC4A4
</td>
<td style="text-align:right;">
0.1545687
</td>
<td style="text-align:right;">
-1.3259896
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.8108783
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP2B4
</td>
<td style="text-align:right;">
0.1462349
</td>
<td style="text-align:right;">
-1.4550796
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.8349490
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
MARCKS
</td>
<td style="text-align:right;">
0.1217731
</td>
<td style="text-align:right;">
-1.0026880
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.9144485
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NPY5R
</td>
<td style="text-align:right;">
0.1202272
</td>
<td style="text-align:right;">
-2.8046686
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.9199974
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PNKP
</td>
<td style="text-align:right;">
0.1047571
</td>
<td style="text-align:right;">
-1.4159956
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.9798166
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PRDM5
</td>
<td style="text-align:right;">
0.0952666
</td>
<td style="text-align:right;">
-3.0416695
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.0210595
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ZFP57
</td>
<td style="text-align:right;">
0.0920165
</td>
<td style="text-align:right;">
-4.6128349
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.0361344
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
INPP5J
</td>
<td style="text-align:right;">
0.0920056
</td>
<td style="text-align:right;">
-0.8559112
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.0361858
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TRPC6
</td>
<td style="text-align:right;">
0.0752220
</td>
<td style="text-align:right;">
-2.8042214
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.1236551
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
FAM53B
</td>
<td style="text-align:right;">
0.0667546
</td>
<td style="text-align:right;">
-1.4018487
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.1755189
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ZFP423
</td>
<td style="text-align:right;">
0.0548335
</td>
<td style="text-align:right;">
2.6224862
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.2609543
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CD47
</td>
<td style="text-align:right;">
0.0476658
</td>
<td style="text-align:right;">
-0.8819554
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
1.3217934
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ABLIM3
</td>
<td style="text-align:right;">
0.0205912
</td>
<td style="text-align:right;">
-1.6690306
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
1.6863192
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
STRA6
</td>
<td style="text-align:right;">
0.0184625
</td>
<td style="text-align:right;">
-4.7284254
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
1.7337097
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TSHZ1
</td>
<td style="text-align:right;">
0.0151200
</td>
<td style="text-align:right;">
-2.7817807
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
1.8204492
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
VPS37B
</td>
<td style="text-align:right;">
0.0118652
</td>
<td style="text-align:right;">
-1.5595107
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
1.9257258
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CYGB
</td>
<td style="text-align:right;">
0.0117389
</td>
<td style="text-align:right;">
-1.8147956
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
1.9303736
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GLIS3
</td>
<td style="text-align:right;">
0.0069128
</td>
<td style="text-align:right;">
-3.2013136
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.1603481
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
IGFBP5
</td>
<td style="text-align:right;">
0.0063220
</td>
<td style="text-align:right;">
-2.8407580
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.1991446
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
KITL
</td>
<td style="text-align:right;">
0.0059249
</td>
<td style="text-align:right;">
-5.9808999
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.2273179
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ARHGAP20
</td>
<td style="text-align:right;">
0.0029272
</td>
<td style="text-align:right;">
-1.6294237
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.5335411
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
DGKH
</td>
<td style="text-align:right;">
0.0022362
</td>
<td style="text-align:right;">
-2.7388342
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.6504858
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NPNT
</td>
<td style="text-align:right;">
0.0018068
</td>
<td style="text-align:right;">
-6.8199667
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.7430930
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
AKAP7
</td>
<td style="text-align:right;">
0.0018022
</td>
<td style="text-align:right;">
-3.2682468
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.7441910
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
MEF2C
</td>
<td style="text-align:right;">
0.0013489
</td>
<td style="text-align:right;">
-1.8481172
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.8700094
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SIPA1L2
</td>
<td style="text-align:right;">
0.0012933
</td>
<td style="text-align:right;">
-1.8850961
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.8882943
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PLEKHG5
</td>
<td style="text-align:right;">
0.0011572
</td>
<td style="text-align:right;">
-1.0314271
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.9366002
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PDE7B
</td>
<td style="text-align:right;">
0.0009437
</td>
<td style="text-align:right;">
-5.9474625
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.0251717
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NBEAL2
</td>
<td style="text-align:right;">
0.0007408
</td>
<td style="text-align:right;">
-3.5431199
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.1302758
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNK1
</td>
<td style="text-align:right;">
0.0005719
</td>
<td style="text-align:right;">
-1.5745229
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.2426895
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SCN3A
</td>
<td style="text-align:right;">
0.0004597
</td>
<td style="text-align:right;">
-2.2927343
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.3375512
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
JUN
</td>
<td style="text-align:right;">
0.0003963
</td>
<td style="text-align:right;">
-1.8852497
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.4019653
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
DSP
</td>
<td style="text-align:right;">
0.0003746
</td>
<td style="text-align:right;">
-7.2929743
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.4264274
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
AU040320
</td>
<td style="text-align:right;">
0.0003466
</td>
<td style="text-align:right;">
-1.5092363
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.4601222
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
KLF13
</td>
<td style="text-align:right;">
0.0002770
</td>
<td style="text-align:right;">
-1.1272778
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.5575971
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
C1QL2
</td>
<td style="text-align:right;">
0.0001618
</td>
<td style="text-align:right;">
-6.4076389
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.7910726
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
LRRC16A
</td>
<td style="text-align:right;">
0.0001488
</td>
<td style="text-align:right;">
-7.4429728
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.8274212
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
MAML2
</td>
<td style="text-align:right;">
0.0001222
</td>
<td style="text-align:right;">
-5.2002803
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.9128878
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ZFP536
</td>
<td style="text-align:right;">
0.0001170
</td>
<td style="text-align:right;">
-7.2847453
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.9319331
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CALD1
</td>
<td style="text-align:right;">
0.0000922
</td>
<td style="text-align:right;">
-4.3073732
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
4.0350734
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NEUROD2
</td>
<td style="text-align:right;">
0.0000390
</td>
<td style="text-align:right;">
-1.2232543
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
4.4090788
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SHISA9
</td>
<td style="text-align:right;">
0.0000381
</td>
<td style="text-align:right;">
-6.0269146
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
4.4188299
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SV2C
</td>
<td style="text-align:right;">
0.0000265
</td>
<td style="text-align:right;">
-5.4381019
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
4.5768355
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GPRC5B
</td>
<td style="text-align:right;">
0.0000095
</td>
<td style="text-align:right;">
-2.4021016
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
5.0238459
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PCDH8
</td>
<td style="text-align:right;">
0.0000077
</td>
<td style="text-align:right;">
-4.1391565
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
5.1145086
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
LRRTM4
</td>
<td style="text-align:right;">
0.0000026
</td>
<td style="text-align:right;">
-4.3620250
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
5.5810153
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
RFX3
</td>
<td style="text-align:right;">
0.0000016
</td>
<td style="text-align:right;">
-2.6646879
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
5.8003128
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
STXBP6
</td>
<td style="text-align:right;">
0.0000008
</td>
<td style="text-align:right;">
-4.2333137
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
6.1147881
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PDZD2
</td>
<td style="text-align:right;">
0.0000003
</td>
<td style="text-align:right;">
-3.6840699
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
6.5151689
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ROBO3
</td>
<td style="text-align:right;">
0.0000003
</td>
<td style="text-align:right;">
-3.0697955
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
6.5740844
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
DOCK10
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-4.3657314
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
7.5322095
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNC3
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-1.5215015
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
7.7373947
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PITPNC1
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-3.2711751
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
7.9985607
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
DGAT2
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-2.6278393
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
8.2735698
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PROX1
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-7.1395609
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
9.0270069
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC29A4
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-5.3274631
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
9.6469553
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PLEKHA2
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-5.1522919
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
10.7749767
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
IL1RAP
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-3.4068378
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
10.9271751
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PITPNM2
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-2.6294539
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
14.5209386
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
FAM163B
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-3.4644652
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
15.5760874
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TIAM1
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-5.2497577
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
17.6205554
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
</tbody>
</table>

    summary(CA1DG_DG$direction)

    ##      DG neither 
    ##      50      38

    CA1CA3_DG <- marker_expression(CA1CA3, "DG", DG_markers)
    kable(CA1CA3_DG)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
direction
</th>
<th style="text-align:right;">
logp
</th>
<th style="text-align:left;">
comparison
</th>
<th style="text-align:left;">
marker
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
INPP5J
</td>
<td style="text-align:right;">
0.9960477
</td>
<td style="text-align:right;">
0.0048900
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0017198
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PRDM5
</td>
<td style="text-align:right;">
0.9913720
</td>
<td style="text-align:right;">
0.0345403
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0037633
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
STXBP6
</td>
<td style="text-align:right;">
0.9810477
</td>
<td style="text-align:right;">
-0.0400911
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0083099
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TAOK2
</td>
<td style="text-align:right;">
0.9788483
</td>
<td style="text-align:right;">
0.0175271
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0092846
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SCN3A
</td>
<td style="text-align:right;">
0.9658984
</td>
<td style="text-align:right;">
-0.0548579
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0150686
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ZFP423
</td>
<td style="text-align:right;">
0.9623187
</td>
<td style="text-align:right;">
-0.1044163
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0166811
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NPNT
</td>
<td style="text-align:right;">
0.9561092
</td>
<td style="text-align:right;">
-0.2321580
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0194925
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP2B4
</td>
<td style="text-align:right;">
0.9549242
</td>
<td style="text-align:right;">
-0.0933218
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0200311
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PDZD2
</td>
<td style="text-align:right;">
0.9501894
</td>
<td style="text-align:right;">
-0.0907684
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0221898
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
H2AFZ
</td>
<td style="text-align:right;">
0.9477175
</td>
<td style="text-align:right;">
0.0518395
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0233211
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ZCCHC12
</td>
<td style="text-align:right;">
0.8704441
</td>
<td style="text-align:right;">
-0.2974691
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0602591
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
FAM53B
</td>
<td style="text-align:right;">
0.8641699
</td>
<td style="text-align:right;">
0.2007073
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0634009
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
C1QL2
</td>
<td style="text-align:right;">
0.8530851
</td>
<td style="text-align:right;">
-0.5498678
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0690076
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CACNG7
</td>
<td style="text-align:right;">
0.8482817
</td>
<td style="text-align:right;">
-0.1578502
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0714599
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GLIS3
</td>
<td style="text-align:right;">
0.8463729
</td>
<td style="text-align:right;">
0.3758195
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0724383
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SCN9A
</td>
<td style="text-align:right;">
0.8284518
</td>
<td style="text-align:right;">
0.8080024
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0817327
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
KLF13
</td>
<td style="text-align:right;">
0.8201485
</td>
<td style="text-align:right;">
-0.1176743
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0861075
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NAPA
</td>
<td style="text-align:right;">
0.8079486
</td>
<td style="text-align:right;">
0.1348031
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0926163
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
DSP
</td>
<td style="text-align:right;">
0.8054653
</td>
<td style="text-align:right;">
0.8651266
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0939532
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
LCORL
</td>
<td style="text-align:right;">
0.7618253
</td>
<td style="text-align:right;">
-0.5826119
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1181446
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GSE1
</td>
<td style="text-align:right;">
0.7353232
</td>
<td style="text-align:right;">
0.2762146
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1335217
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TRPC6
</td>
<td style="text-align:right;">
0.7319333
</td>
<td style="text-align:right;">
0.7685376
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1355285
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GTF2A2
</td>
<td style="text-align:right;">
0.7220605
</td>
<td style="text-align:right;">
0.4156935
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1414264
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TRHR
</td>
<td style="text-align:right;">
0.7175505
</td>
<td style="text-align:right;">
-1.3237536
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1441475
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SV2C
</td>
<td style="text-align:right;">
0.6967302
</td>
<td style="text-align:right;">
-0.8086928
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1569354
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
IGFBP5
</td>
<td style="text-align:right;">
0.6858173
</td>
<td style="text-align:right;">
0.6265635
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1637916
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SH3BP1
</td>
<td style="text-align:right;">
0.6797519
</td>
<td style="text-align:right;">
0.1654830
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1676496
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NPY5R
</td>
<td style="text-align:right;">
0.6557918
</td>
<td style="text-align:right;">
-1.0386595
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1832340
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNK1
</td>
<td style="text-align:right;">
0.6506523
</td>
<td style="text-align:right;">
-0.3073742
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1866510
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GLP2R
</td>
<td style="text-align:right;">
0.6248674
</td>
<td style="text-align:right;">
1.3137265
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2042122
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ROBO3
</td>
<td style="text-align:right;">
0.6244817
</td>
<td style="text-align:right;">
-0.4631014
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2044803
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
MARCKS
</td>
<td style="text-align:right;">
0.6217946
</td>
<td style="text-align:right;">
0.4134640
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2063531
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC4A4
</td>
<td style="text-align:right;">
0.5670385
</td>
<td style="text-align:right;">
-0.6682564
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2463875
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PITPNM2
</td>
<td style="text-align:right;">
0.5566320
</td>
<td style="text-align:right;">
0.3195764
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2544319
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
VPS37B
</td>
<td style="text-align:right;">
0.5341154
</td>
<td style="text-align:right;">
0.5227741
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2723649
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
FBLN1
</td>
<td style="text-align:right;">
0.5266732
</td>
<td style="text-align:right;">
-2.1034832
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2784588
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
JUN
</td>
<td style="text-align:right;">
0.4792543
</td>
<td style="text-align:right;">
-0.5324696
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3194340
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNC3
</td>
<td style="text-align:right;">
0.4358144
</td>
<td style="text-align:right;">
-0.3115423
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3606984
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PDE7B
</td>
<td style="text-align:right;">
0.4320643
</td>
<td style="text-align:right;">
-1.9062244
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3644516
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ECM2
</td>
<td style="text-align:right;">
0.3959879
</td>
<td style="text-align:right;">
-1.8322880
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4023180
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PLEKHG5
</td>
<td style="text-align:right;">
0.3796605
</td>
<td style="text-align:right;">
0.3732809
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4206046
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
MEF2C
</td>
<td style="text-align:right;">
0.3705249
</td>
<td style="text-align:right;">
0.6914092
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4311825
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PNKP
</td>
<td style="text-align:right;">
0.3687803
</td>
<td style="text-align:right;">
-0.9121780
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4332323
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
LRRTM4
</td>
<td style="text-align:right;">
0.3662961
</td>
<td style="text-align:right;">
-1.1832014
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4361677
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TSHZ1
</td>
<td style="text-align:right;">
0.3620846
</td>
<td style="text-align:right;">
-1.3088795
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4411900
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SIPA1L2
</td>
<td style="text-align:right;">
0.3489717
</td>
<td style="text-align:right;">
0.7307136
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4572098
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
B4GALT3
</td>
<td style="text-align:right;">
0.3489510
</td>
<td style="text-align:right;">
0.6751344
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4572356
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
AKAP7
</td>
<td style="text-align:right;">
0.3342639
</td>
<td style="text-align:right;">
-1.3061999
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4759105
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PDZRN3
</td>
<td style="text-align:right;">
0.3241727
</td>
<td style="text-align:right;">
-3.0188387
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4892236
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
EXT2
</td>
<td style="text-align:right;">
0.3113138
</td>
<td style="text-align:right;">
0.6913527
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.5068016
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TOX3
</td>
<td style="text-align:right;">
0.2892826
</td>
<td style="text-align:right;">
2.1353173
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.5386777
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
DGAT2
</td>
<td style="text-align:right;">
0.2712343
</td>
<td style="text-align:right;">
-0.6957288
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.5666555
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CYGB
</td>
<td style="text-align:right;">
0.2696407
</td>
<td style="text-align:right;">
-0.9707975
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.5692146
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
DOCK10
</td>
<td style="text-align:right;">
0.2381598
</td>
<td style="text-align:right;">
-1.2980138
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.6231315
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PLEKHA2
</td>
<td style="text-align:right;">
0.2322239
</td>
<td style="text-align:right;">
-1.3004322
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.6340932
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
RAVER1
</td>
<td style="text-align:right;">
0.2125335
</td>
<td style="text-align:right;">
0.6851594
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.6725725
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NEUROD2
</td>
<td style="text-align:right;">
0.1732261
</td>
<td style="text-align:right;">
-0.5254213
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.7613867
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
BTG2
</td>
<td style="text-align:right;">
0.1634914
</td>
<td style="text-align:right;">
1.3175662
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.7865050
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
DGKH
</td>
<td style="text-align:right;">
0.1592407
</td>
<td style="text-align:right;">
1.5939591
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.7979460
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PROX1
</td>
<td style="text-align:right;">
0.1520947
</td>
<td style="text-align:right;">
-2.2883668
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.8178859
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GPRC5B
</td>
<td style="text-align:right;">
0.1454802
</td>
<td style="text-align:right;">
-1.0315136
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.8371961
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
EGR4
</td>
<td style="text-align:right;">
0.1371590
</td>
<td style="text-align:right;">
1.5064941
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.8627758
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
STRA6
</td>
<td style="text-align:right;">
0.1340285
</td>
<td style="text-align:right;">
-3.4117039
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.8728027
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CHUK
</td>
<td style="text-align:right;">
0.1225572
</td>
<td style="text-align:right;">
0.8688364
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.9116613
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ARHGAP20
</td>
<td style="text-align:right;">
0.1222818
</td>
<td style="text-align:right;">
1.0314514
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.9126382
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TIAM1
</td>
<td style="text-align:right;">
0.1154888
</td>
<td style="text-align:right;">
-1.3201120
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.9374602
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC29A4
</td>
<td style="text-align:right;">
0.1040846
</td>
<td style="text-align:right;">
-1.8128861
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.9826136
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
FAT4
</td>
<td style="text-align:right;">
0.0950403
</td>
<td style="text-align:right;">
2.1778648
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.0220921
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
MAML2
</td>
<td style="text-align:right;">
0.0923696
</td>
<td style="text-align:right;">
-2.7853786
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.0344711
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ZFP57
</td>
<td style="text-align:right;">
0.0910430
</td>
<td style="text-align:right;">
-4.8859909
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.0407533
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ETS2
</td>
<td style="text-align:right;">
0.0899719
</td>
<td style="text-align:right;">
0.9200101
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.0458929
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
IL1RAP
</td>
<td style="text-align:right;">
0.0835377
</td>
<td style="text-align:right;">
-1.1568181
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.0781176
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
RFX3
</td>
<td style="text-align:right;">
0.0694287
</td>
<td style="text-align:right;">
1.3062993
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.1584608
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NCOR2
</td>
<td style="text-align:right;">
0.0634713
</td>
<td style="text-align:right;">
0.6707897
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.1974229
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PITPNC1
</td>
<td style="text-align:right;">
0.0565902
</td>
<td style="text-align:right;">
-1.3943194
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.2472585
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PCDH8
</td>
<td style="text-align:right;">
0.0520789
</td>
<td style="text-align:right;">
-2.2491685
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.2833385
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CALD1
</td>
<td style="text-align:right;">
0.0411376
</td>
<td style="text-align:right;">
-2.6885426
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.3857607
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NRGN
</td>
<td style="text-align:right;">
0.0378107
</td>
<td style="text-align:right;">
0.8217013
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
1.4223848
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
IL33
</td>
<td style="text-align:right;">
0.0365066
</td>
<td style="text-align:right;">
-2.3268491
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.4376291
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
STK32B
</td>
<td style="text-align:right;">
0.0338094
</td>
<td style="text-align:right;">
-5.3113457
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.4709626
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CD47
</td>
<td style="text-align:right;">
0.0325580
</td>
<td style="text-align:right;">
-1.0064556
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.4873421
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
KITL
</td>
<td style="text-align:right;">
0.0311453
</td>
<td style="text-align:right;">
-5.1514744
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.5066072
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NBEAL2
</td>
<td style="text-align:right;">
0.0258096
</td>
<td style="text-align:right;">
-2.6798044
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.5882182
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SHISA9
</td>
<td style="text-align:right;">
0.0207047
</td>
<td style="text-align:right;">
-3.9708970
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.6839308
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
AU040320
</td>
<td style="text-align:right;">
0.0151918
</td>
<td style="text-align:right;">
-1.1594907
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.8183902
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
LRRC16A
</td>
<td style="text-align:right;">
0.0078968
</td>
<td style="text-align:right;">
-5.8568333
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
2.1025479
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ZFP536
</td>
<td style="text-align:right;">
0.0020865
</td>
<td style="text-align:right;">
-6.3163692
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
2.6805913
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ABLIM3
</td>
<td style="text-align:right;">
0.0012594
</td>
<td style="text-align:right;">
2.4621399
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
2.8998254
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
FAM163B
</td>
<td style="text-align:right;">
0.0000097
</td>
<td style="text-align:right;">
2.2356702
</td>
<td style="text-align:left;">
CA1
</td>
<td style="text-align:right;">
5.0116236
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
DG
</td>
</tr>
</tbody>
</table>

    summary(CA1CA3_DG$direction)

    ##     CA1     CA3 neither 
    ##       3      10      76

    CA3DG_DG <- marker_expression(CA3DG, "DG", DG_markers)
    kable(CA3DG_DG)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
lfc
</th>
<th style="text-align:left;">
direction
</th>
<th style="text-align:right;">
logp
</th>
<th style="text-align:left;">
comparison
</th>
<th style="text-align:left;">
marker
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
ZCCHC12
</td>
<td style="text-align:right;">
0.9948188
</td>
<td style="text-align:right;">
0.0098437
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0022560
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TOX3
</td>
<td style="text-align:right;">
0.9227969
</td>
<td style="text-align:right;">
-0.2431629
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0348939
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ZFP57
</td>
<td style="text-align:right;">
0.9203646
</td>
<td style="text-align:right;">
0.2731560
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0360401
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ECM2
</td>
<td style="text-align:right;">
0.8733296
</td>
<td style="text-align:right;">
-0.3056870
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0588218
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GTF2A2
</td>
<td style="text-align:right;">
0.8436795
</td>
<td style="text-align:right;">
0.2233622
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0738225
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CD47
</td>
<td style="text-align:right;">
0.8341661
</td>
<td style="text-align:right;">
0.1245003
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0787475
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CACNG7
</td>
<td style="text-align:right;">
0.8226318
</td>
<td style="text-align:right;">
-0.1641599
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0847945
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TRHR
</td>
<td style="text-align:right;">
0.7298578
</td>
<td style="text-align:right;">
-0.9205221
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1367617
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PDZRN3
</td>
<td style="text-align:right;">
0.6888003
</td>
<td style="text-align:right;">
0.9362085
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1619067
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CHUK
</td>
<td style="text-align:right;">
0.6567806
</td>
<td style="text-align:right;">
-0.2637348
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1825797
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SCN9A
</td>
<td style="text-align:right;">
0.6343613
</td>
<td style="text-align:right;">
1.4256019
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1976633
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
KITL
</td>
<td style="text-align:right;">
0.5687005
</td>
<td style="text-align:right;">
-0.8294255
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2451164
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
H2AFZ
</td>
<td style="text-align:right;">
0.5514649
</td>
<td style="text-align:right;">
-0.3050888
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2584821
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PNKP
</td>
<td style="text-align:right;">
0.5453818
</td>
<td style="text-align:right;">
-0.5038176
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2632994
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC4A4
</td>
<td style="text-align:right;">
0.5024077
</td>
<td style="text-align:right;">
-0.6577333
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.2989437
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
LCORL
</td>
<td style="text-align:right;">
0.4954726
</td>
<td style="text-align:right;">
-0.8407211
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3049803
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
FBLN1
</td>
<td style="text-align:right;">
0.4662145
</td>
<td style="text-align:right;">
-1.4884889
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3314142
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TAOK2
</td>
<td style="text-align:right;">
0.4481848
</td>
<td style="text-align:right;">
-0.2931309
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3485429
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GSE1
</td>
<td style="text-align:right;">
0.4364086
</td>
<td style="text-align:right;">
-0.4774007
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.3601067
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
EXT2
</td>
<td style="text-align:right;">
0.3962864
</td>
<td style="text-align:right;">
-0.4928944
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4019908
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
AU040320
</td>
<td style="text-align:right;">
0.3786259
</td>
<td style="text-align:right;">
-0.3497456
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4217896
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NBEAL2
</td>
<td style="text-align:right;">
0.3189220
</td>
<td style="text-align:right;">
-0.8633155
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.4963155
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
STRA6
</td>
<td style="text-align:right;">
0.2960985
</td>
<td style="text-align:right;">
-1.3167215
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.5285637
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ZFP536
</td>
<td style="text-align:right;">
0.2599906
</td>
<td style="text-align:right;">
-0.9683762
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.5850423
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
RAVER1
</td>
<td style="text-align:right;">
0.2275340
</td>
<td style="text-align:right;">
-0.5378103
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.6429538
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CYGB
</td>
<td style="text-align:right;">
0.1937324
</td>
<td style="text-align:right;">
-0.8439981
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.7127978
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
EGR4
</td>
<td style="text-align:right;">
0.1648591
</td>
<td style="text-align:right;">
-1.1705309
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.7828870
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NPY5R
</td>
<td style="text-align:right;">
0.1356180
</td>
<td style="text-align:right;">
-1.7660091
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.8676827
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TSHZ1
</td>
<td style="text-align:right;">
0.1204634
</td>
<td style="text-align:right;">
-1.4729011
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.9191448
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
BTG2
</td>
<td style="text-align:right;">
0.1158687
</td>
<td style="text-align:right;">
-1.1827894
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.9360337
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SH3BP1
</td>
<td style="text-align:right;">
0.1049011
</td>
<td style="text-align:right;">
-0.4091897
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.9792198
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GLP2R
</td>
<td style="text-align:right;">
0.0949612
</td>
<td style="text-align:right;">
-2.8338068
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.0224540
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
LRRC16A
</td>
<td style="text-align:right;">
0.0913203
</td>
<td style="text-align:right;">
-1.5861395
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.0394328
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP2B4
</td>
<td style="text-align:right;">
0.0833983
</td>
<td style="text-align:right;">
-1.3617578
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.0788429
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
B4GALT3
</td>
<td style="text-align:right;">
0.0632118
</td>
<td style="text-align:right;">
-0.9340554
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.1992018
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NAPA
</td>
<td style="text-align:right;">
0.0609959
</td>
<td style="text-align:right;">
-0.6012627
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
1.2146997
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NCOR2
</td>
<td style="text-align:right;">
0.0345292
</td>
<td style="text-align:right;">
-0.5951578
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
1.4618137
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
CALD1
</td>
<td style="text-align:right;">
0.0339002
</td>
<td style="text-align:right;">
-1.6188306
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
1.4697975
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
INPP5J
</td>
<td style="text-align:right;">
0.0335293
</td>
<td style="text-align:right;">
-0.8608012
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
1.4745749
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PCDH8
</td>
<td style="text-align:right;">
0.0204157
</td>
<td style="text-align:right;">
-1.8899880
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
1.6900354
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ZFP423
</td>
<td style="text-align:right;">
0.0116787
</td>
<td style="text-align:right;">
2.7269025
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
1.9326056
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PRDM5
</td>
<td style="text-align:right;">
0.0108312
</td>
<td style="text-align:right;">
-3.0762098
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
1.9653245
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NRGN
</td>
<td style="text-align:right;">
0.0093363
</td>
<td style="text-align:right;">
-0.7788153
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.0298232
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SHISA9
</td>
<td style="text-align:right;">
0.0077122
</td>
<td style="text-align:right;">
-2.0560175
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.1128212
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NEUROD2
</td>
<td style="text-align:right;">
0.0066962
</td>
<td style="text-align:right;">
-0.6978331
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.1741739
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
STK32B
</td>
<td style="text-align:right;">
0.0043598
</td>
<td style="text-align:right;">
3.6084550
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
2.3605377
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
FAM53B
</td>
<td style="text-align:right;">
0.0043390
</td>
<td style="text-align:right;">
-1.6025560
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.3626085
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
AKAP7
</td>
<td style="text-align:right;">
0.0036809
</td>
<td style="text-align:right;">
-1.9620470
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.4340470
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
MARCKS
</td>
<td style="text-align:right;">
0.0035590
</td>
<td style="text-align:right;">
-1.4161520
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.4486762
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GPRC5B
</td>
<td style="text-align:right;">
0.0031532
</td>
<td style="text-align:right;">
-1.3705880
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.5012535
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
JUN
</td>
<td style="text-align:right;">
0.0030986
</td>
<td style="text-align:right;">
-1.3527801
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.5088307
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TRPC6
</td>
<td style="text-align:right;">
0.0022415
</td>
<td style="text-align:right;">
-3.5727589
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
2.6494673
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
IL33
</td>
<td style="text-align:right;">
0.0010572
</td>
<td style="text-align:right;">
2.2587679
</td>
<td style="text-align:left;">
CA3
</td>
<td style="text-align:right;">
2.9758296
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
FAT4
</td>
<td style="text-align:right;">
0.0005130
</td>
<td style="text-align:right;">
-3.1148401
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.2898875
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNK1
</td>
<td style="text-align:right;">
0.0002495
</td>
<td style="text-align:right;">
-1.2671487
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.6028731
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ETS2
</td>
<td style="text-align:right;">
0.0001738
</td>
<td style="text-align:right;">
-1.3544451
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
3.7600325
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
MAML2
</td>
<td style="text-align:right;">
0.0000618
</td>
<td style="text-align:right;">
-2.4149016
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
4.2090335
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
KLF13
</td>
<td style="text-align:right;">
0.0000608
</td>
<td style="text-align:right;">
-1.0096035
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
4.2159365
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
IGFBP5
</td>
<td style="text-align:right;">
0.0000148
</td>
<td style="text-align:right;">
-3.4673216
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
4.8305762
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
GLIS3
</td>
<td style="text-align:right;">
0.0000084
</td>
<td style="text-align:right;">
-3.5771331
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
5.0757374
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
C1QL2
</td>
<td style="text-align:right;">
0.0000065
</td>
<td style="text-align:right;">
-5.8577710
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
5.1884486
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SCN3A
</td>
<td style="text-align:right;">
0.0000043
</td>
<td style="text-align:right;">
-2.2378763
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
5.3711561
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PITPNC1
</td>
<td style="text-align:right;">
0.0000024
</td>
<td style="text-align:right;">
-1.8768557
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
5.6194318
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
VPS37B
</td>
<td style="text-align:right;">
0.0000016
</td>
<td style="text-align:right;">
-2.0822848
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
5.7995859
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
NPNT
</td>
<td style="text-align:right;">
0.0000015
</td>
<td style="text-align:right;">
-6.5878086
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
5.8190292
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PDE7B
</td>
<td style="text-align:right;">
0.0000010
</td>
<td style="text-align:right;">
-4.0412381
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
6.0156402
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
DOCK10
</td>
<td style="text-align:right;">
0.0000004
</td>
<td style="text-align:right;">
-3.0677176
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
6.4491689
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
LRRTM4
</td>
<td style="text-align:right;">
0.0000001
</td>
<td style="text-align:right;">
-3.1788236
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
6.8587537
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SV2C
</td>
<td style="text-align:right;">
0.0000001
</td>
<td style="text-align:right;">
-4.6294091
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
7.0872255
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
KCNC3
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-1.2099592
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
7.4719548
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ROBO3
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-2.6066941
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
7.8360097
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
DSP
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-8.1581009
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
7.9380971
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PLEKHG5
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-1.4047080
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
8.2244562
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
DGAT2
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-1.9321105
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
8.2480882
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
MEF2C
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-2.5395265
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
8.2849814
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SIPA1L2
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-2.6158096
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
8.5116791
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
STXBP6
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-4.1932226
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
9.3916071
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PROX1
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-4.8511940
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
9.5093658
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
DGKH
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-4.3327932
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
9.9684896
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ARHGAP20
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-2.6608751
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
10.8886059
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PDZD2
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-3.5933015
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
11.2024418
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
IL1RAP
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-2.2500197
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
11.3592354
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
ABLIM3
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-4.1311706
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
14.8531613
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PLEKHA2
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-3.8518597
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
16.2335462
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC29A4
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-3.5145770
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
18.0259317
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
RFX3
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-3.9709872
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
19.4776869
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
TIAM1
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-3.9296457
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
21.9371765
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
PITPNM2
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-2.9490302
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
26.5481661
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
<tr>
<td style="text-align:left;">
FAM163B
</td>
<td style="text-align:right;">
0.0000000
</td>
<td style="text-align:right;">
-5.7001354
</td>
<td style="text-align:left;">
DG
</td>
<td style="text-align:right;">
59.6418514
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
DG
</td>
</tr>
</tbody>
</table>

    summary(CA3DG_DG$direction)

    ##     CA3      DG neither 
    ##       3      50      36
