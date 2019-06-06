numDEGs <- function(group1, group2){
  res <- results(dds, contrast = c("APA2", group1, group2), independentFiltering = T)
  sumpadj <- sum(res$padj < 0.1, na.rm = TRUE)
  return(sumpadj)
}


resvals <- function(contrastvector, mypval){
  res <- results(dds, contrast = c(contrastvector[1],contrastvector[2],contrastvector[3]), independentFiltering = T)
  sumpvalue <- sum(res$pvalue < mypval, na.rm = TRUE)
  #print(sumpvalue)
  sumpadj <- sum(res$padj < mypval, na.rm = TRUE)
  print(sumpadj)
  vals <- cbind(res$pvalue, res$padj)
  pvalcolname <- as.character(paste("pval",contrastvector[1],contrastvector[2],contrastvector[3], sep=""))
  padjcolname <- as.character(paste("padj",contrastvector[1],contrastvector[2],contrastvector[3], sep=""))
  colnames(vals) <- c(pvalcolname, padjcolname)
  return(vals)
}

myhistogram <- function(contrastvector, mypval){
  res <- results(dds, contrast = c(contrastvector[1],contrastvector[2],contrastvector[3]), independentFiltering = T)
  vals <- cbind(res$pvalue)
  pvalcolname <- as.character(paste(contrastvector[2], "vs",contrastvector[3], sep=" "))
  colnames(vals) <- c(pvalcolname)
  vals <- as.data.frame(vals)
  histogram <- hist(vals)
  return(histogram)
  print(histogram)
}

pcadataframe <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], 
                  PC4 = pca$x[, 4], PC5 = pca$x[, 5], PC6 = pca$x[, 6], 
                  PC7 = pca$x[, 7], PC8 = pca$x[, 8], PC9 = pca$x[, 9],
                  group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:9]
    return(d)
  }
}


plotPCs <- function(df, xcol, ycol, aescolor, colorname, aesshape, shapename, colorvalues){
  ggplot(df, aes(df[xcol], df[ycol], color=aescolor, shape=aesshape)) +
    geom_point(size=2) +
    xlab(paste0("PC", xcol, ": ", percentVar[xcol],"% variance")) +
    ylab(paste0("PC", ycol, ": ", percentVar[ycol],"% variance")) +
    theme_classic() +
    
    #stat_ellipse(level = 0.95, (aes(color=aescolor)),size=1) + 
    scale_colour_manual(name=colorname, values=c(colorvalues))+
    scale_shape_discrete(name=shapename) +
    theme(axis.text = element_text(size=8),
          axis.title.x = element_text(size=8),
          axis.title.y = element_text(size=8),
          legend.title = element_text(size=8),
          legend.text = element_text(size=8)) +
    #theme(legend.title=element_blank()) 
    background_grid(major = "y", minor = "y") +
    theme(legend.position="none")
}

plotPCwrap <- function(df, xcol, ycol, aescolor, colorname, colorvalues){
  ggplot(df, aes(df[xcol], df[ycol], color=aescolor)) +
    geom_point(size=2, alpha= 0.8) +
    xlab(paste0("PC", xcol, ": ", percentVar[xcol],"% variance")) +
    ylab(paste0("PC", ycol, ": ", percentVar[ycol],"% variance")) +
    stat_ellipse(level = 0.95, (aes(color=aescolor)),size=0.25) + 
    scale_colour_manual(name=colorname, values=c(colorvalues))+ 
    theme_cowplot(font_size = 8, line_size = 0.25)  +
    theme(legend.position="none") +
    facet_wrap(~wrap)
}


plotPCwrapnoe <- function(df, xcol, ycol, aescolor, colorname, colorvalues){
  ggplot(df, aes(df[xcol], df[ycol], color=aescolor)) +
    geom_point(size=2, alpha= 0.8) +
    xlab(paste0("PC", xcol, ": ", percentVar[xcol],"% variance")) +
    ylab(paste0("PC", ycol, ": ", percentVar[ycol],"% variance")) +
    #stat_ellipse(level = 0.95, (aes(color=aescolor)),size=0.25) + 
    scale_colour_manual(name=colorname, values=c(colorvalues))+ 
    theme_cowplot(font_size = 8, line_size = 0.25)  +
    theme(legend.position="none") +
    facet_wrap(~wrap)
}


## plot DEGs 

subsetDESeq <- function(eachgroup){
  
  # subset to look within one tissue in one sex
  colData <- colData %>%
    dplyr::filter(subfield == eachgroup) %>%
    droplevels()
  row.names(colData) <- colData$RNAseqID
  print(colData)
  
  # which counts to save
  savecols <- as.character(colData$RNAseqID) 
  savecols <- as.vector(savecols) 
  
  # save counts that match colData
  countData <- countData %>% dplyr::select(one_of(savecols)) 
  
  # check that row and col lenghts are equal
  print(ncol(countData) == nrow(colData))
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ treatment )
  
  print(dds)
  dds <- dds[rowSums(counts(dds) > 1) >= 10]  # filter more than sample with less 0 counts
  print(dim(dds))
  
  dds <- DESeq(dds) # Differential expression analysis
  return(dds)
}

numDEGs <- function(dds, group1, group2){
  res <- results(dds, contrast = c("treatment", group1, group2), independentFiltering = T)
  sumpadj <- sum(res$padj < 0.01, na.rm = TRUE)
  return(sumpadj)
}


returntotalDEGs <- function(dds){
  
  colData <- colData  
  group1 <- levels(colData$treatment)
  
  a <- group1
  b <- group1
  
  # comapre all contrasts, save to datafrmes
  totalDEGS=data.frame()
  for (i in a){
    for (j in b){
      if (i != j) {
        k <- paste(i,j, sep = ".") #assigns usique rownames
        print(k)
        totalDEGS[k,1]<-i               
        totalDEGS[k,2]<-j
        totalDEGS[k,3]<- numDEGs(dds, i,j) #caluculates number of DEGs
      }
    }
    b <- b[-1]  # drop 1st element of second string to not recalculate DEGs
  }
  
  print(totalDEGS)  
  return(totalDEGS)
  
}

plottotalDEGs <- function(myDEGS, mysubtitle){  
  totalDEGS <- myDEGS
  totalDEGS$V2 <- factor(totalDEGS$V2, levels =  c("home.cage",
                                                   "standard.yoked" ,"standard.trained", 
                                                   "conflict.yoked", "conflict.trained"))
  
  totalDEGS$V1 <- factor(totalDEGS$V1, levels =  c("home.cage",
                                                   "standard.yoked" ,"standard.trained", 
                                                   "conflict.yoked", "conflict.trained"))
  
  totalDEGS <- totalDEGS %>% dplyr::na_if(0)
  
  print(str(totalDEGS))
  
  
  allcontrasts <- totalDEGS %>%
    ggplot( aes(V1, V2)) +
    geom_tile(aes(fill = V3)) +
    theme_minimal(base_size = 8) + 
    geom_text(aes(label = round(V3, 1)), color = "black")+
    scale_fill_viridis(na.value="#bdbdbd", 
                      limits = c(0,602),
                       option = "C") +
    xlab(NULL) + ylab(NULL) +
    labs(fill = "# of DEGs",
         title = mysubtitle, subtitle = "  ", caption = "  ") +
    theme(axis.text.x = element_text(angle = 45)) +
    coord_flip()
  print(totalDEGS)
  plot(allcontrasts)
}


## new correlation heatmap

# subset col for heatmap
subsetcolData <- function(a.colData, eachgroup){
  
  # subset to look within one tissue in one sex
  colData <- a.colData %>%
    dplyr::filter(subfield == eachgroup) %>%
    droplevels()
  row.names(colData) <- colData$RNAseqID
  return(colData)
}

plotcorrelationheatmaps <- function(mydds, mycoldata, mysubtitle){
  dds <- mydds
  vsd <- vst(dds, blind=FALSE) # variance stabilized 
  
  colnames(vsd) = mycoldata$treatment # set col names to group name
  
  vsdm <- assay(vsd) # create matrix
  
  vsdmmean <-sapply(unique(colnames(vsdm)), function(i)
    rowMeans(vsdm[,colnames(vsdm) == i]))
  
  myannotations = data.frame(
    treatment = factor(c("home.cage","standard.yoked" ,"standard.trained", "conflict.yoked", "conflict.trained"))
  )

  myBreaks <-seq(0.8, 1.0, length.out = 10)

  pheatmap(cor(vsdmmean),
           #annotation_row = myannotations,
           #annotation_col = myannotations,
           #annotation_colors = myannotationscolors,
           annotation_names_row = F,
           main= mysubtitle,
           color = inferno(10),
           show_rowname= F, show_colnames = F,
           breaks = myBreaks
  )
}

