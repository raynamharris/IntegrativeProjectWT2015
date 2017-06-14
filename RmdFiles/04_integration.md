Setup
-----

``` r
library(plyr)
library(dplyr)
library(reshape2)
```

Import Data
-----------

``` r
behaviorpca <- read.csv("../data/01a_scoresdf.csv", header = T)
behavior <- read.csv("../data/01a_behavior.csv", header = T)
ephys <- read.csv("../data/03_ephys.csv", header = T)
pcadata <- read.csv("../data/02a_pcadata.csv", header = T)
rossetta <- read.csv("../data/00_rossettastone.csv", header = F)
```

Wrangle Data
------------

``` r
# clearnup the rosetts data and filter extraneous samples
names(rossetta)[1] <- "Mouse"
names(rossetta)[2] <- "ID"
names(rossetta)[3] <- "Region"
names(rossetta)[4] <- "RNAseqID"
names(rossetta)[5] <- "R1filename"
rossetta$R1filename <- NULL
rossetta <- rossetta %>% dplyr::filter(Mouse != "15-100", Mouse != "15-101", Mouse != "15-147")
head(rossetta) # dictionary of names
```

    ##     Mouse     ID Region   RNAseqID
    ## 1 15-143A 15143A    CA3 143A-CA3-1
    ## 2 15-143A 15143A     DG  143A-DG-1
    ## 3 15-143B 15143B    CA1 143B-CA1-1
    ## 4 15-143B 15143B     DG  143B-DG-1
    ## 5 15-143C 15143C    CA1 143C-CA1-1
    ## 6 15-143D 15143D    CA1 143D-CA1-3

``` r
## slim behavior ephy to top 5 pcs and rename the columsn
behaviorpca <- behaviorpca[(c(1:5,35:36))]
names(behaviorpca)[names(behaviorpca)=="PC1"] <- "Behavior_PC1"
names(behaviorpca)[names(behaviorpca)=="PC2"] <- "Behavior_PC2"
names(behaviorpca)[names(behaviorpca)=="PC3"] <- "Behavior_PC3"
names(behaviorpca)[names(behaviorpca)=="PC4"] <- "Behavior_PC4"
names(behaviorpca)[names(behaviorpca)=="PC5"] <- "Behavior_PC5"


behaviorpca <- behaviorpca %>% dplyr::filter(ID != "15148", ID !=  "15140A", ID !=  "15140B", ID !=  "15140C", ID !=  "15140D", ID !=  "15141C", ID !=  "15141D", ID !=  "15142C", ID !=  "15142D", ID !=  "15142A", ID !=  "15142B", ID !=  "15145C", ID !=  "15145C", ID !=  "15145D", ID !=  "15147A", ID !=  "15147B", ID !=  "15148C", ID !=  "15148D")

pcadata <- pcadata[(c(1:9,11,13))]
names(pcadata)[names(pcadata)=="name"] <- "RNAseqID"

#wident then length RNAseq data so each row is an animals
pcadatabyregion <- left_join(pcadata, rossetta)
```

    ## Joining, by = "RNAseqID"

    ## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
    ## factors with different levels, coercing to character vector

``` r
pcadatabyregion <- melt(pcadatabyregion, id = c(10:14))
pcadatabyregion$RegionPC <- as.factor(paste(pcadatabyregion$Region, pcadatabyregion$variable, sep="_"))
pcadatabyregion <- dcast(pcadatabyregion, Mouse + ID ~ RegionPC)


alldata <- left_join(behaviorpca, pcadatabyregion)
```

    ## Joining, by = "ID"

    ## Warning in left_join_impl(x, y, by$x, by$y, suffix$x, suffix$y): joining
    ## factors with different levels, coercing to character vector

``` r
library(superheat)
# plot a heatmap of the numerical iris variables
# cluster by species and plot Sepal.Length on the right
# save the superheat object to access the membership vectors
sh <- superheat(X = iris[,-c(1, 5)],
                yr = iris[,1],
                yr.axis.name = "Sepal.Length",
                membership.rows = iris$Species)
```

![](04_integration_files/figure-markdown_github/superheatmap-1.png)

``` r
head(iris)
```

    ##   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
    ## 1          5.1         3.5          1.4         0.2  setosa
    ## 2          4.9         3.0          1.4         0.2  setosa
    ## 3          4.7         3.2          1.3         0.2  setosa
    ## 4          4.6         3.1          1.5         0.2  setosa
    ## 5          5.0         3.6          1.4         0.2  setosa
    ## 6          5.4         3.9          1.7         0.4  setosa

``` r
sh2 <- superheat(X = alldata[,-c(6,7,8)],
                #yr = behaviorpca[,10],
                #yr.axis.name = "PC10",
                membership.rows = alldata$APA,
                pretty.order.cols = TRUE,
                col.dendrogram = TRUE,
                bottom.label.size = 0.4,
                bottom.label.text.angle = 90,
                scale = TRUE)
```

![](04_integration_files/figure-markdown_github/superheatmap-2.png)
