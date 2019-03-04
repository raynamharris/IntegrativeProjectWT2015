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

    marker_summary <- function(mydf, subfield, markers){
        MARKERS <- str_to_upper(markers)
        df <- mydf %>%
        dplyr::filter(gene %in% c(MARKERS)) %>%
        dplyr::mutate(marker = subfield) %>%
        droplevels()
        #return(kable(head(df)))
        return(summary(df$direction))
    }


    #CA1 markers
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
    marker_summary(CA1DG, "CA4", CA4_markers)

    ##     CA1 neither 
    ##       4      23

    marker_summary(CA1CA3, "CA4", CA4_markers)

    ##     CA3 neither 
    ##       1      27

    marker_summary(CA3DG, "CA4", CA4_markers)

    ##     CA3      DG neither 
    ##      11       1      18
