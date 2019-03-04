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

    mydfs <- list(CA1DG, CA1CA3, CA3DG)

    #look for markers in each supfield
    # make data frames of genes expression results for markers 

    marker_summary <- function(mydf, markers){
        MARKERS <- str_to_upper(markers)
        df <- mydf %>%
        dplyr::filter(gene %in% c(MARKERS)) 
        #return((head(df)))
        return(summary(df$direction))
    }

    for(i in mydfs){
      j <- marker_summary(i, CA1_markers)
      print(i[1, 6])
      print(j)
    }

    ## [1] "CA1-DG"
    ##     CA1      DG neither 
    ##      17       0       8 
    ## [1] "CA1-CA3"
    ##     CA1     CA3 neither 
    ##      16       2       7 
    ## [1] "CA3-DG"
    ##     CA3      DG neither 
    ##       7       2      18

    for(i in mydfs){
      j <- marker_summary(i, CA2_markers)
      print(i[1, 6])
      print(j)
    }

    ## [1] "CA1-DG"
    ##     CA1      DG neither 
    ##       4       3      29 
    ## [1] "CA1-CA3"
    ##     CA1     CA3 neither 
    ##       1       8      27 
    ## [1] "CA3-DG"
    ##     CA3      DG neither 
    ##      16       0      23

    for(i in mydfs){
      j <- marker_summary(i, CA3_markers)
      print(i[1, 6])
      print(j)
    }

    ## [1] "CA1-DG"
    ##     CA1      DG neither 
    ##       1       2       7 
    ## [1] "CA1-CA3"
    ##     CA1     CA3 neither 
    ##       0       6       4 
    ## [1] "CA3-DG"
    ##     CA3      DG neither 
    ##       7       0       3

    for(i in mydfs){
      j <- marker_summary(i, CA4_markers)
      print(i[1, 6])
      print(j)
    }

    ## [1] "CA1-DG"
    ##     CA1      DG neither 
    ##       4       0      23 
    ## [1] "CA1-CA3"
    ##     CA1     CA3 neither 
    ##       0       1      27 
    ## [1] "CA3-DG"
    ##     CA3      DG neither 
    ##      11       1      18

    for(i in mydfs){
      j <- marker_summary(i, DG_markers)
      print(i[1, 6])
      print(j)
    }

    ## [1] "CA1-DG"
    ##     CA1      DG neither 
    ##       0      50      38 
    ## [1] "CA1-CA3"
    ##     CA1     CA3 neither 
    ##       3      10      76 
    ## [1] "CA3-DG"
    ##     CA3      DG neither 
    ##       3      50      36
