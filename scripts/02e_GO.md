Hypothesised GO Terms
---------------------

    GO_response <- read.table("../data/goterms/GO_term_summary_20191121_150656.txt", sep = "\t", row.names = NULL,  fill=TRUE)
    GO_response$GO <- "1. Response to Stimulus (GO:0050896)"

    GO_translation <- read.table("../data/goterms/GO_term_summary_20191121_141252.txt", sep = "\t", row.names = NULL)
    GO_translation$GO <- "2. Translation (GO:0006412)"

    GO_synapse <- read.table("../data/goterms/GO_term_summary_20191121_145034.txt", sep = "\t", row.names = NULL)
    GO_synapse$GO <- "3. Synapse Organization (GO:0050808)"

    GO_learningormemory <- read.table("../data/goterms/GO_term_summary_20191121_142500.txt", sep = "\t", row.names = NULL)
    GO_learningormemory$GO <- "4. Learning or Memory (GO:0007611)"

    GOterms <- rbind(GO_learningormemory, GO_response, GO_translation, GO_synapse)
    GOterms <- GOterms %>%
      dplyr::mutate(gene = toupper(MGI.Gene.Marker.ID)) %>% 
      dplyr::select(gene, GO) %>% 
      dplyr::distinct(gene, GO) %>% 
     group_by(gene) 
    GOterms

    ## # A tibble: 7,156 x 2
    ## # Groups:   gene [6,717]
    ##    gene  GO                                
    ##    <chr> <chr>                             
    ##  1 AAAS  4. Learning or Memory (GO:0007611)
    ##  2 AAL   4. Learning or Memory (GO:0007611)
    ##  3 ABCA7 4. Learning or Memory (GO:0007611)
    ##  4 ABCC8 4. Learning or Memory (GO:0007611)
    ##  5 ABI2  4. Learning or Memory (GO:0007611)
    ##  6 ABL2  4. Learning or Memory (GO:0007611)
    ##  7 ADAM2 4. Learning or Memory (GO:0007611)
    ##  8 ADCY1 4. Learning or Memory (GO:0007611)
    ##  9 ADCY3 4. Learning or Memory (GO:0007611)
    ## 10 ADCY8 4. Learning or Memory (GO:0007611)
    ## # … with 7,146 more rows

Upregulated DG DEGS
-------------------

    DGDEGs <- read_csv("../data/02c_DEGs.DG.training.yoked.trained.csv") %>% filter(direction == "trained")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_double(),
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   direction = col_character()
    ## )

    DGDEGs$gene <- toupper(DGDEGs$gene)
    DGDEGs$gene

    ##   [1] "1190002N15RIK" "A830010M20RIK" "ABHD2"         "ACAN"         
    ##   [5] "ADAMTS1"       "ADRB1"         "AHR"           "AMIGO2"       
    ##   [9] "ANKRD13A"      "ANKRD28"       "ANKRD33B"      "APAF1"        
    ##  [13] "ARC"           "ARID5B"        "ARL13B"        "ARL4A"        
    ##  [17] "ARL4D"         "ARL5B"         "ARMCX5"        "ARPP21"       
    ##  [21] "ATF3"          "B3GNT2"        "BACH1"         "BDNF"         
    ##  [25] "BTG2"          "C2CD4B"        "CCNK"          "CIART"        
    ##  [29] "CITED2"        "CLDN12"        "CNNM1"         "CPEB4"        
    ##  [33] "CTNND1"        "CUL3"          "CWC25"         "CXADR"        
    ##  [37] "CYP51"         "DBPHT2"        "DNAJA1"        "DNAJB1"       
    ##  [41] "DNAJB4"        "DUSP14"        "DUSP16"        "DUSP4"        
    ##  [45] "DUSP6"         "DUSP8"         "DYRK2"         "EGR1"         
    ##  [49] "EGR3"          "EGR4"          "EIF5"          "EPRS"         
    ##  [53] "ERF"           "ERRFI1"        "FAM107B"       "FAM118A"      
    ##  [57] "FBXO33"        "FBXW7"         "FERMT2"        "FLRT3"        
    ##  [61] "FOS"           "FOSB"          "FOSL2"         "FOXG1"        
    ##  [65] "FOXO1"         "FRMD6"         "FZD4"          "FZD5"         
    ##  [69] "GAD1"          "GADD45G"       "GM13889"       "GMEB2"        
    ##  [73] "GPR19"         "HECA"          "HMGCR"         "HOMER1"       
    ##  [77] "HS6ST1"        "HSPA1A"        "HSPH1"         "IL16"         
    ##  [81] "ING2"          "IRF2BP2"       "IRS1"          "IRS2"         
    ##  [85] "JDP2"          "JMJD1C"        "JUN"           "JUNB"         
    ##  [89] "JUND"          "KCNA4"         "KCNF1"         "KCNJ2"        
    ##  [93] "KDM6B"         "KDM7A"         "KITL"          "KLF2"         
    ##  [97] "KLF6"          "LBH"           "LCMT2"         "LEMD3"        
    ## [101] "LMNA"          "LONRF1"        "LRRTM2"        "MARCH11"      
    ## [105] "MED7"          "MEST"          "MFAP3L"        "MN1"          
    ## [109] "MYC"           "NAF1"          "NAP1L5"        "NEDD9"        
    ## [113] "NEFM"          "NFIL3"         "NPAS4"         "NPTX2"        
    ## [117] "NR4A1"         "NR4A2"         "NR4A3"         "NUAK1"        
    ## [121] "ODC1"          "OLFML2B"       "OTUD1"         "PAK6"         
    ## [125] "PCDH8"         "PEG10"         "PELI1"         "PER1"         
    ## [129] "PER2"          "PHLDA1"        "PIGA"          "PLAGL1"       
    ## [133] "PLK2"          "PLK3"          "POU3F3"        "PPP1R15A"     
    ## [137] "PRPF38B"       "PTGS2"         "RANBP2"        "RASD1"        
    ## [141] "RASL11A"       "RASL11B"       "RFX2"          "RGMB"         
    ## [145] "RGS2"          "RGS4"          "SCG2"          "SGK1"         
    ## [149] "SH2D3C"        "SIAH2"         "SLC16A1"       "SLC25A25"     
    ## [153] "SLC2A3"        "SLC45A4"       "SLITRK5"       "SMAD7"        
    ## [157] "SNX18"         "SOWAHC"        "SOX9"          "SPTY2D1"      
    ## [161] "SRF"           "STMN4"         "SYT4"          "THBS1"        
    ## [165] "TIPARP"        "TNIP2"         "TRA2B"         "TRIB1"        
    ## [169] "TSC22D2"       "UBC"           "USPL1"         "ZBTB33"       
    ## [173] "ZDBF2"         "ZFAND5"        "ZFP275"        "ZFP654"       
    ## [177] "ZFP869"

    GOtermsDEGs <- inner_join(GOterms, DGDEGs) %>%
      arrange(GO,gene) %>%
      dplyr::select(gene, GO) %>% 
      group_by(GO) %>%
      summarize(genes = str_c(gene, collapse = ", "))

    ## Joining, by = "gene"

sanes
-----

    sanesLichtman <- read_csv("../data/02e_sanesLichtman.csv")  %>% rename(gene = genes)

    ## Parsed with column specification:
    ## cols(
    ##   genes = col_character()
    ## )

    sanesLichtman

    ## # A tibble: 237 x 1
    ##    gene   
    ##    <chr>  
    ##  1 ACHE   
    ##  2 ADCY1  
    ##  3 ADRA2A 
    ##  4 ADRA2B 
    ##  5 ADRA2C 
    ##  6 ADRB1  
    ##  7 ADRB2  
    ##  8 ADRB3  
    ##  9 BDNF   
    ## 10 CACNA1A
    ## # … with 227 more rows

    GOtermssanesLichtman <- inner_join(GOterms, sanesLichtman) %>%
      arrange(GO,gene) %>%
      dplyr::select(gene, GO) %>% 
     group_by(GO) %>%
     summarize(genes = str_c(gene, collapse = ", "))

    ## Joining, by = "gene"

DEGs and Sanes
--------------

    GO_DEGsSanesLichtman <- full_join(GOtermsDEGs, GOtermssanesLichtman, by = "GO") %>%
      rename(  "Differentially expressed in DG"  = genes.x,    "Molecules associated with LTP"  = genes.y)

    kable(GO_DEGsSanesLichtman)

<table>
<thead>
<tr>
<th style="text-align:left;">
GO
</th>
<th style="text-align:left;">
Differentially expressed in DG
</th>
<th style="text-align:left;">
Molecules associated with LTP
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">

1.  Response to Stimulus
    (<a href="GO:0050896" class="uri">GO:0050896</a>)
    </td>
    <td style="text-align:left;">
    ABHD2, ADRB1, AHR, APAF1, ARC, ARID5B, ARL13B, ARPP21, ATF3, B3GNT2,
    BACH1, BTG2, C2CD4B, CCNK, CITED2, FBXW7, FERMT2, FLRT3, FOS, FOSB,
    FOSL2, FOXG1, FOXO1, FZD4, FZD5, GADD45G, GPR19, HMGCR, HOMER1,
    HSPA1A, HSPH1, IL16, ING2, IRS1, IRS2, JUN, JUNB, JUND, KDM6B, KITL,
    KLF2, KLF6, LBH, LEMD3, LMNA, MEST, MYC, NEDD9, NFIL3, NPAS4, NPTX2,
    NR4A1, NR4A2, NR4A3, PPP1R15A, SLC16A1, SLC25A25, SLITRK5, SMAD7,
    SOX9, SRF, SYT4, THBS1, TIPARP, TNIP2, TRA2B, TRIB1, TSC22D2,
    ZBTB33, ZFAND5
    </td>
    <td style="text-align:left;">
    ADCY1, ADRA2A, ADRA2B, ADRA2C, ADRB1, ADRB2, ADRB3, CACNA1A,
    CACNA1B, CACNA1C, CACNA1D, CACNA1E, CACNA1F, CACNA1S, CALB1, CALM1,
    CALM2, CALM3, CAMK1, CAMK4, CAPN1, CAPN10, CAPN2, CAPN3, CCR7, CD47,
    CDH1, CDH2, CHRM1, CHRM2, CHRM3, CHRM4, CHRM5, CHRNA1, CHRNA3,
    CHRNA7, CHRNB1, CHRNB2, CHRNB3, CNGA2, FGF2, FYN, GABBR1, GABRA1,
    GABRA2, GABRA3, GABRA5, GABRA6, GABRB1, GABRB2, GABRB3, GABRR1,
    GAP43, GFAP, GRIA1, GRIA2, GRIN1, GRIN2A, GRIN2D, GRM1, GRM4, GRM5,
    GRM7, GUCY1A2, GUCY1B2, GUCY2C, GUCY2D, GUCY2E, GUCY2G, HOMER1,
    HOMER2, HOMER3, HTR1A, HTR1B, HTR1F, HTR2A, HTR2B, HTR2C, HTR3A,
    HTR3B, HTR4, HTR5A, HTR5B, HTR6, HTR7, IL1B, INHBA, ITGA1, ITGA10,
    ITGA11, ITGA2, ITGA2B, ITGA3, ITGA4, ITGA5, ITGA6, ITGA7, ITGA8,
    ITGA9, ITGAD, ITGAE, ITGAL, ITGAM, ITGAV, ITGAX, ITGB1, ITGB1BP1,
    ITGB2, ITGB2L, ITGB3, ITGB4, ITGB5, ITGB6, ITGB7, ITGB8, ITGBL1,
    ITPKB, L1CAM, MAPK1, MAPK11, MAPK12, MAPK14, MAPK3, MAPK4, MAPK6,
    MAPK7, MAPK8, MAPK9, MAS1, NCAM1, NGF, NOS1, NOS3, NRG1, NRG2, NRG3,
    NRGN, PNOC, SPTBN1, SRC, STX1B, SYP, TH, THY1, TNC, UBE3A, VAMP2,
    VAMP3, VAMP4, VAMP8
    </td>
    </tr>
    <tr>
    <td style="text-align:left;">

    1.  Translation (<a href="GO:0006412" class="uri">GO:0006412</a>)
        </td>
        <td style="text-align:left;">
        CPEB4, EIF5
        </td>
        <td style="text-align:left;">
        NA
        </td>
        </tr>
        <tr>
        <td style="text-align:left;">

        1.  Synapse Organization
            (<a href="GO:0050808" class="uri">GO:0050808</a>)
            </td>
            <td style="text-align:left;">
            AMIGO2, ARC, BDNF, FLRT3, FZD5, HOMER1, LRRTM2, NPAS4,
            PCDH8, SLITRK5
            </td>
            <td style="text-align:left;">
            ACHE, BDNF, CACNA1A, CACNA1S, CAMK1, CDH1, CDH2, CHRNA1,
            CHRNA7, CHRNB1, CHRNB2, DLG4, DRD1, EFNA5, EPHA5, ERBB4,
            FYN, GABRA1, GABRB2, GABRB3, GAP43, GRIN1, GRIN2A, GRM5,
            HOMER1, HTR1A, ITGA3, ITGAM, ITPKA, L1CAM, MAPK14, NRG1,
            NRG2, NTRK2, PTN, RAB3A, SYN1, TNC, UBE3A
            </td>
            </tr>
            <tr>
            <td style="text-align:left;">

            1.  Learning or Memory
                (<a href="GO:0007611" class="uri">GO:0007611</a>)
                </td>
                <td style="text-align:left;">
                ADRB1, ARC, BDNF, BTG2, EGR1, HMGCR, JUN, NPAS4, NPTX2,
                PAK6, PLK2, PTGS2, SGK1, SRF, SYT4
                </td>
                <td style="text-align:left;">
                ADCY1, ADRB1, ADRB2, BDNF, CACNA1C, CACNA1E, CALB1,
                CAMK4, CHRNA7, CHRNB2, CNR1, CREB1, DRD1, EGR1, GABRA5,
                GRIA1, GRIN1, GRIN2A, GRM4, GRM5, GRM7, HTR2A, HTR6,
                HTR7, IL1B, ITGA3, ITGA5, ITGA8, ITGB1, NCAM1, NGF,
                NTRK2, OPRL1, PLA2G6, PLCB1, PRKAR1B, PRKCZ, PTN, S100B,
                SNAP25, TH
                </td>
                </tr>
                </tbody>
                </table>

PC1 correlations
----------------

    corrsPC1sig <- read_csv("../data/02d_corrsPC1sig.csv")

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   PC1 = col_double()
    ## )

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

    ## these 57 genes were used for a go analysis using shinygo http://bioinformatics.sdstate.edu/go/
    ## results are stored as 
      #data/02e_DEGPC1enrichmentBP.csv
        #data/02e_DEGPC1enrichmentCC.csv
        #data/02e_DEGPC1enrichmentMF.csv

    BP <- read_csv("../data/02e_GO_PC1enrichmentBP.csv", n_max = 5) %>% mutate(Domain = "BP") 

    ## Parsed with column specification:
    ## cols(
    ##   `Enrichment FDR` = col_double(),
    ##   `Genes in list` = col_double(),
    ##   `Total genes` = col_double(),
    ##   `Functional Category` = col_character(),
    ##   Genes = col_character()
    ## )

    CC <- read_csv("../data/02e_GO_PC1enrichmentCC.csv", n_max = 5) %>% mutate(Domain = "CC") 

    ## Parsed with column specification:
    ## cols(
    ##   `Enrichment FDR` = col_double(),
    ##   `Genes in list` = col_double(),
    ##   `Total genes` = col_double(),
    ##   `Functional Category` = col_character(),
    ##   Genes = col_character()
    ## )

    MF <- read_csv("../data/02e_GO_PC1enrichmentMF.csv", n_max = 5) %>% mutate(Domain = "MF") 

    ## Parsed with column specification:
    ## cols(
    ##   `Enrichment FDR` = col_double(),
    ##   `Genes in list` = col_double(),
    ##   `Total genes` = col_double(),
    ##   `Functional Category` = col_character(),
    ##   Genes = col_character()
    ## )

    GO_PC1 <- rbind(BP, CC, MF) %>% select(Domain, `Functional Category`, `Total genes`, `Genes in list`, Genes)

    # alphabetize
    GO_PC1$Genes <- GO_PC1$Genes %>% str_split(., ' ') %>% lapply(., 'sort') %>%  lapply(., 'paste', collapse=' ') %>% unlist(.)

    kable(GO_PC1) 

<table>
<thead>
<tr>
<th style="text-align:left;">
Domain
</th>
<th style="text-align:left;">
Functional Category
</th>
<th style="text-align:right;">
Total genes
</th>
<th style="text-align:right;">
Genes in list
</th>
<th style="text-align:left;">
Genes
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BP
</td>
<td style="text-align:left;">
Memory
</td>
<td style="text-align:right;">
145
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;">
ARC BDNF EGR1 KCNK10 NPAS4 PLK2 PTGS2 SGK1 SYT4
</td>
</tr>
<tr>
<td style="text-align:left;">
BP
</td>
<td style="text-align:left;">
Tissue development
</td>
<td style="text-align:right;">
1855
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:left;">
ACAN ARC AREG ATF3 BDNF BTG2 COL10A1 EGR1 ERRFI1 FOSL2 FRMD6 FZD5 HOMER1
HOXC4 NR4A3 PCDH8 PTGS2 RGS2 SLC25A25 SMAD7 TIPARP
</td>
</tr>
<tr>
<td style="text-align:left;">
BP
</td>
<td style="text-align:left;">
Learning or memory
</td>
<td style="text-align:right;">
286
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
ARC BDNF BTG2 EGR1 KCNK10 NPAS4 PLK2 PTGS2 SGK1 SYT4
</td>
</tr>
<tr>
<td style="text-align:left;">
BP
</td>
<td style="text-align:left;">
Cognition
</td>
<td style="text-align:right;">
317
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
ARC BDNF BTG2 EGR1 KCNK10 NPAS4 PLK2 PTGS2 SGK1 SYT4
</td>
</tr>
<tr>
<td style="text-align:left;">
BP
</td>
<td style="text-align:left;">
Behavior
</td>
<td style="text-align:right;">
693
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
ARC BDNF BTG2 EGR1 HOMER1 KCNK10 NPAS4 NR4A3 PLK2 PTGS2 SGK1 SLC16A1
SYT4
</td>
</tr>
<tr>
<td style="text-align:left;">
CC
</td>
<td style="text-align:left;">
Neuron projection
</td>
<td style="text-align:right;">
1486
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
ACAN ARC BDNF CPEB4 FZD5 HOMER1 NEXN PCDH8 PLK2 PTGS2 RGS2 SGK1 SYT4
</td>
</tr>
<tr>
<td style="text-align:left;">
CC
</td>
<td style="text-align:left;">
Cell junction
</td>
<td style="text-align:right;">
1095
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
ARC CPEB4 FRMD6 FZD5 HOMER1 NEXN PCDH8 SLC16A1 SMAD7 SYT4
</td>
</tr>
<tr>
<td style="text-align:left;">
CC
</td>
<td style="text-align:left;">
Dendrite
</td>
<td style="text-align:right;">
736
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:left;">
ARC BDNF CPEB4 FZD5 HOMER1 PCDH8 PLK2 SYT4
</td>
</tr>
<tr>
<td style="text-align:left;">
CC
</td>
<td style="text-align:left;">
Dendritic tree
</td>
<td style="text-align:right;">
738
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:left;">
ARC BDNF CPEB4 FZD5 HOMER1 PCDH8 PLK2 SYT4
</td>
</tr>
<tr>
<td style="text-align:left;">
CC
</td>
<td style="text-align:left;">
Neuron part
</td>
<td style="text-align:right;">
1933
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
ACAN ARC BDNF CPEB4 FZD5 HOMER1 NEXN PCDH8 PLK2 PTGS2 RGS2 SGK1 SYT4
</td>
</tr>
<tr>
<td style="text-align:left;">
MF
</td>
<td style="text-align:left;">
Regulatory region nucleic acid binding
</td>
<td style="text-align:right;">
949
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
ATF3 EGR1 EGR4 FOSL2 HOXC4 NFIL3 NPAS4 NR4A3 PER1 SMAD7 TIPARP
</td>
</tr>
<tr>
<td style="text-align:left;">
MF
</td>
<td style="text-align:left;">
Transcription regulatory region DNA binding
</td>
<td style="text-align:right;">
946
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
ATF3 EGR1 EGR4 FOSL2 HOXC4 NFIL3 NPAS4 NR4A3 PER1 SMAD7 TIPARP
</td>
</tr>
<tr>
<td style="text-align:left;">
MF
</td>
<td style="text-align:left;">
RNA polymerase II regulatory region sequence-specific DNA binding
</td>
<td style="text-align:right;">
782
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;">
ATF3 EGR1 EGR4 FOSL2 HOXC4 NFIL3 NPAS4 NR4A3 PER1
</td>
</tr>
<tr>
<td style="text-align:left;">
MF
</td>
<td style="text-align:left;">
DNA-binding transcription factor activity, RNA polymerase II-specific
</td>
<td style="text-align:right;">
738
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;">
ATF3 BTG2 EGR1 EGR4 FOSL2 HOXC4 NFIL3 NPAS4 NR4A3
</td>
</tr>
<tr>
<td style="text-align:left;">
MF
</td>
<td style="text-align:left;">
RNA polymerase II regulatory region DNA binding
</td>
<td style="text-align:right;">
788
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;">
ATF3 EGR1 EGR4 FOSL2 HOXC4 NFIL3 NPAS4 NR4A3 PER1
</td>
</tr>
</tbody>
</table>

    write.csv(GO_DEGsSanesLichtman, "../data/02e_GO_DEGsSanesLichtman.csv", row.names = F)
    write.csv(GO_PC1, "../data/02e_GO_PC1.csv", row.names = F)
