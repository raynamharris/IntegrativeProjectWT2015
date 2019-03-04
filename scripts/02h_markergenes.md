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
    # make data frames of genes expression results for markers 
    marker_expression <- function(mydf, subfield, markers){
        MARKERS <- str_to_upper(markers)
        df <- mydf %>%
        dplyr::filter(gene %in% c(MARKERS)) %>%
        dplyr::mutate(marker = subfield) %>%
        droplevels()
        return(kable(head(df)))
    }

    marker_summary <- function(mydf, subfield, markers){
        MARKERS <- str_to_upper(markers)
        df <- mydf %>%
        dplyr::filter(gene %in% c(MARKERS)) %>%
        dplyr::mutate(marker = subfield) %>%
        droplevels()
        return(summary(df$direction))
    }


    #CA1 markers
    marker_expression(CA1DG, "CA1", CA1_markers)

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
</tbody>
</table>

    marker_expression(CA1CA3, "CA1", CA1_markers)

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
</tbody>
</table>

    marker_expression(CA3DG, "CA1", CA1_markers)

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
</tbody>
</table>

    marker_summary(CA1DG, "CA1", CA1_markers)

    ##     CA1 neither 
    ##      17       8

    marker_summary(CA1CA3, "CA1", CA1_markers)

    ##     CA1     CA3 neither 
    ##      16       2       7

    marker_summary(CA3DG, "CA1", CA1_markers)

    ##     CA3      DG neither 
    ##       7       2      18

    #CA3 markers
    marker_expression(CA1DG, "CA3", CA3_markers)

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
</tbody>
</table>

    marker_expression(CA1CA3, "CA3", CA3_markers)

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
</tbody>
</table>

    marker_expression(CA3DG, "CA3", CA3_markers)

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
</tbody>
</table>

    marker_summary(CA1DG, "CA3", CA3_markers)

    ##     CA1      DG neither 
    ##       1       2       7

    marker_summary(CA1CA3, "CA3", CA3_markers)

    ##     CA3 neither 
    ##       6       4

    marker_summary(CA3DG, "CA3", CA3_markers)

    ##     CA3 neither 
    ##       7       3

    #DG markers
    marker_expression(CA1DG, "DG", DG_markers)

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
</tbody>
</table>

    marker_expression(CA1CA3, "DG", DG_markers)

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
</tbody>
</table>

    marker_expression(CA3DG, "DG", DG_markers)

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
</tbody>
</table>

    marker_summary(CA1DG, "DG", DG_markers)

    ##      DG neither 
    ##      50      38

    marker_summary(CA1CA3, "DG", DG_markers)

    ##     CA1     CA3 neither 
    ##       3      10      76

    marker_summary(CA3DG, "DG", DG_markers)

    ##     CA3      DG neither 
    ##       3      50      36

    #CA2 markers
    marker_expression(CA1DG, "CA2", CA2_markers)

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
STARD5
</td>
<td style="text-align:right;">
0.9935771
</td>
<td style="text-align:right;">
0.0162932
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0027984
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
TGFB1I1
</td>
<td style="text-align:right;">
0.9898905
</td>
<td style="text-align:right;">
0.0441772
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0044128
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
CCDC3
</td>
<td style="text-align:right;">
0.9551579
</td>
<td style="text-align:right;">
-0.1950221
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0199248
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
PLCXD3
</td>
<td style="text-align:right;">
0.9351512
</td>
<td style="text-align:right;">
-0.3196008
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0291182
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
ADCY5
</td>
<td style="text-align:right;">
0.8410592
</td>
<td style="text-align:right;">
-0.1573519
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0751734
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
LMAN2
</td>
<td style="text-align:right;">
0.8361923
</td>
<td style="text-align:right;">
-0.1663754
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0776938
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
</tbody>
</table>

    marker_expression(CA1CA3, "CA2", CA2_markers)

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
S100B
</td>
<td style="text-align:right;">
0.9977255
</td>
<td style="text-align:right;">
-0.0037770
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0009889
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
SRGAP2
</td>
<td style="text-align:right;">
0.9490577
</td>
<td style="text-align:right;">
-0.0622659
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0227074
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC4A8
</td>
<td style="text-align:right;">
0.8972082
</td>
<td style="text-align:right;">
-0.1079909
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0471068
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
GSTO1
</td>
<td style="text-align:right;">
0.8119602
</td>
<td style="text-align:right;">
-0.2493835
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0904652
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
NDRG3
</td>
<td style="text-align:right;">
0.7476589
</td>
<td style="text-align:right;">
0.1336132
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1262965
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
STARD5
</td>
<td style="text-align:right;">
0.7019315
</td>
<td style="text-align:right;">
-0.6818641
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1537053
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
</tbody>
</table>

    marker_expression(CA3DG, "CA2", CA2_markers)

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
GSTO1
</td>
<td style="text-align:right;">
0.9964599
</td>
<td style="text-align:right;">
-0.0042628
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0015402
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
STXBP5L
</td>
<td style="text-align:right;">
0.9618199
</td>
<td style="text-align:right;">
-0.0517211
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0169062
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
PRSS23
</td>
<td style="text-align:right;">
0.9355109
</td>
<td style="text-align:right;">
0.1434952
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0289511
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
ARG2
</td>
<td style="text-align:right;">
0.9199638
</td>
<td style="text-align:right;">
-0.3804561
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0362293
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
FGF5
</td>
<td style="text-align:right;">
0.9025627
</td>
<td style="text-align:right;">
-0.2625428
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0445226
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
<tr>
<td style="text-align:left;">
FGF2
</td>
<td style="text-align:right;">
0.8736317
</td>
<td style="text-align:right;">
-0.3495466
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0586716
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA2
</td>
</tr>
</tbody>
</table>

    marker_summary(CA1DG, "CA2", CA2_markers)

    ##     CA1      DG neither 
    ##       4       3      29

    marker_summary(CA1CA3, "CA2", CA2_markers)

    ##     CA1     CA3 neither 
    ##       1       8      27

    marker_summary(CA3DG, "CA2", CA2_markers)

    ##     CA3 neither 
    ##      16      23

    #CA4 markers
    marker_expression(CA1DG, "CA4", CA4_markers)

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
RGS12
</td>
<td style="text-align:right;">
0.9870056
</td>
<td style="text-align:right;">
0.0341754
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0056804
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
PMP22
</td>
<td style="text-align:right;">
0.9575135
</td>
<td style="text-align:right;">
-0.1540648
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0188551
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
HCFC1R1
</td>
<td style="text-align:right;">
0.9258916
</td>
<td style="text-align:right;">
-0.0971213
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0334398
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
CAMK2N2
</td>
<td style="text-align:right;">
0.8963837
</td>
<td style="text-align:right;">
-0.1539523
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0475061
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
FXYD6
</td>
<td style="text-align:right;">
0.8766904
</td>
<td style="text-align:right;">
-0.3016679
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0571537
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
RASL11A
</td>
<td style="text-align:right;">
0.8484333
</td>
<td style="text-align:right;">
0.3478579
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0713823
</td>
<td style="text-align:left;">
CA1-DG
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
</tbody>
</table>

    marker_expression(CA1CA3, "CA4", CA4_markers)

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
CHCHD2
</td>
<td style="text-align:right;">
0.9393319
</td>
<td style="text-align:right;">
0.0638881
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0271809
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
RPL23
</td>
<td style="text-align:right;">
0.8476304
</td>
<td style="text-align:right;">
-0.1585312
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0717935
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
EPHB6
</td>
<td style="text-align:right;">
0.7611097
</td>
<td style="text-align:right;">
0.2075323
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1185528
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
RAC3
</td>
<td style="text-align:right;">
0.7403064
</td>
<td style="text-align:right;">
0.5358988
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1305885
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
STMN3
</td>
<td style="text-align:right;">
0.7217624
</td>
<td style="text-align:right;">
-0.2209464
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1416057
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP6V1F
</td>
<td style="text-align:right;">
0.6518695
</td>
<td style="text-align:right;">
-0.4202827
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1858393
</td>
<td style="text-align:left;">
CA1-CA3
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
</tbody>
</table>

    marker_expression(CA3DG, "CA4", CA4_markers)

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
RBPMS2
</td>
<td style="text-align:right;">
0.9817130
</td>
<td style="text-align:right;">
-0.0943706
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0080155
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP6V1F
</td>
<td style="text-align:right;">
0.9138559
</td>
<td style="text-align:right;">
0.1068304
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0391223
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
H2-D1
</td>
<td style="text-align:right;">
0.8664988
</td>
<td style="text-align:right;">
0.1570580
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0622320
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
RPL23
</td>
<td style="text-align:right;">
0.8579207
</td>
<td style="text-align:right;">
-0.1347459
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0665529
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
RAC3
</td>
<td style="text-align:right;">
0.8156012
</td>
<td style="text-align:right;">
0.3618010
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.0885221
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
<tr>
<td style="text-align:left;">
GAL
</td>
<td style="text-align:right;">
0.7029866
</td>
<td style="text-align:right;">
-3.8001751
</td>
<td style="text-align:left;">
neither
</td>
<td style="text-align:right;">
0.1530530
</td>
<td style="text-align:left;">
CA3-DG
</td>
<td style="text-align:left;">
CA4
</td>
</tr>
</tbody>
</table>

    marker_summary(CA1DG, "CA4", CA4_markers)

    ##     CA1 neither 
    ##       4      23

    marker_summary(CA1CA3, "CA4", CA4_markers)

    ##     CA3 neither 
    ##       1      27

    marker_summary(CA3DG, "CA4", CA4_markers)

    ##     CA3      DG neither 
    ##      11       1      18
