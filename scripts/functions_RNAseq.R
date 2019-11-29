

## used in 02_rnaseqQC.Rmd

plot.tSNE.trained <- function(mydds, myperplexity, mysubfield){
  vsddft <- as.data.frame(t(assay(vst(mydds))))
  euclidist <- dist(vsddft) # euclidean distances between the rows
  tsne_model <- Rtsne(euclidist, check_duplicates=FALSE, pca=TRUE, 
                      perplexity=myperplexity, theta=0.5, dims=2)
  
  tsne_df = as.data.frame(tsne_model$Y) 
  
  colData <- subsetcolData(colData, mysubfield)
  tsne_df <- cbind(colData, tsne_df)
  tsne_df$subfield <- factor(tsne_df$subfield, levels = c("DG", "CA3", "CA1"))
  tsne_df$training <- factor(tsne_df$training, levels = c("yoked", "trained"))
  
  tnseplot  <- ggplot(tsne_df, aes(x = V1, y = V2, color = training, label = ID)) +
    geom_point(size = 2) +
    scale_color_manual(guide = FALSE, values = c("black","darkred")) +
    theme_ms()  +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.spacing.x = unit(0.01, 'cm'),
          legend.spacing.y = unit(0.01, 'cm')) +
    labs(x = "tSNE 1", y = "tSNE 2") +
    geom_text(vjust = -0.1)
  return(tnseplot)
}


## used in 03_rnaseqSubfield

returnddstreatment <- function(mytissue){
  print(mytissue)
  colData <- colData %>% 
    filter(subfield %in% c(mytissue))  %>% 
    droplevels()
  
  savecols <- as.character(colData$RNAseqID) 
  savecols <- as.vector(savecols) 
  countData <- countData %>% dplyr::select(one_of(savecols)) 
  
  ## create DESeq object using the factors subfield and APA
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ treatment)
  
  dds <- dds[ rowSums(counts(dds)) > 1, ]  # Pre-filtering genes with 0 counts
  dds <- DESeq(dds, parallel = TRUE)
  return(dds)
}

returnddstraining <- function(mytissue){
  print(mytissue)
  colData <- colData %>% 
    filter(subfield %in% c(mytissue))  %>% 
    droplevels()
  
  savecols <- as.character(colData$RNAseqID) 
  savecols <- as.vector(savecols) 
  countData <- countData %>% dplyr::select(one_of(savecols)) 
  
  ## create DESeq object using the factors subfield and APA
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ training)
  
  dds <- dds[ rowSums(counts(dds)) > 1, ]  # Pre-filtering genes with 0 counts
  dds <- DESeq(dds, parallel = TRUE)
  return(dds)
}


savevsds <- function(mydds, vsdfilename){
  dds <- mydds
  vsd <- vst(dds, blind=FALSE) ## variance stabilized
  print(head(assay(vsd),3))
  return(write.csv(assay(vsd), file = vsdfilename, row.names = T))
}

returnvsds2 <- function(mydds, vsdfilename){
  dds <- mydds
  vsd <- vst(dds, blind=FALSE) ## variance stabilized
  print(head(assay(vsd),3))
  myvsd <- assay(vsd)
  myvsd <- as.data.frame(myvsd)
  return(myvsd)
}


res_summary <- function(mycontrast){
  res <- results(dds, contrast = mycontrast, independentFiltering = T)
  print(mycontrast)
  print(sum(res$padj < 0.1, na.rm=TRUE))
  print(summary(res))
  print(head((res[order(res$padj),]), 5))
  cat("\n")
}

res_summary_subfield <- function(mydds, mycontrast){
  res <- results(mydds, contrast = mycontrast, independentFiltering = T)
  print(mycontrast)
  print(sum(res$padj < 0.1, na.rm=TRUE))
  print(summary(res))
  cat("\n")
}

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


resvals2 <- function(mydds, contrastvector, mypval){
  res <- results(mydds, contrast = c(contrastvector[1],contrastvector[2],contrastvector[3]), independentFiltering = T)
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

## these next 3 funcitons are needed for volcano plots

calculateDEGs <-  function(mydds, whichtissue, whichfactor, up, down){
  
  # calculate DEG results
  res <- results(mydds, contrast =c(whichfactor, up, down),
                 independentFiltering = T, alpha = 0.1)
  
  # create dataframe with pvalues and lfc
  data <- data.frame(gene = row.names(res),
                     padj = res$padj, 
                     logpadj = -log10(res$padj),
                     lfc = res$log2FoldChange)
  data <- na.omit(data)
  data <- data %>%
    dplyr::mutate(direction = ifelse(data$lfc > 0 & data$padj < 0.1, 
                                     yes = up, no = ifelse(data$lfc < 0 & data$padj < 0.1, 
                                                           yes = down, no = "NS")))
  data$direction <- factor(data$direction, levels = c(down, "NS", up))
  return(data)
}  

saveDEGs <-  function(data, whichtissue, up, down){
  DEGs <- data %>% filter(direction != "NS")
  myfilename <-  paste("../data/03_", whichtissue, "_DEGs_", down, up, sep = "")
  filenamewithextension <- paste(myfilename, ".csv", sep = "")
  if (dim(DEGs)[1] != 0) {return(write.csv(DEGs, filenamewithextension))}
}  

plot.volcano <- function(data, mycolors, mysubtitle){
  
  volcano <- data %>%
    ggplot(aes(x = lfc, y = logpadj)) + 
    geom_point(aes(color = direction), size = 1, alpha = 0.75, na.rm = T) + 
    theme_ms() +
    scale_color_manual(values = mycolors,
                       name = " ",
                       drop = FALSE) +
    ylim(c(0,12.5)) +  
    xlim(c(-8,8)) +
    labs(y = NULL, x = NULL, subtitle = mysubtitle)  +
    theme(legend.position = "bottom",
          legend.spacing.x = unit(-0.1, 'cm'),
          legend.margin=margin(t=-0.25, r=0, b=0, l=0, unit="cm"),
          panel.grid = element_blank()) 
  return(volcano)
  
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
  sumpadj <- sum(res$padj < 0.1, na.rm = TRUE)
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
                      limits = c(1,1700),
                       option = "C",
                      alpha = 0.5) +
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
subsetcolData <- function(colData, eachgroup){
  
  # subset to look within one tissue in one sex
  colData <- colData %>%
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

  #myBreaks <-seq(0.8, 1.0, length.out = 10)

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

plotcorrelationheatmaps2 <- function(mydds, mycoldata, mysubtitle){
  dds <- mydds
  vsd <- vst(dds, blind=FALSE) # variance stabilized 
  
  colnames(vsd) = mycoldata$subtreat # set col names to group name
  
  vsdm <- assay(vsd) # create matrix
  
  vsdmmean <-sapply(unique(colnames(vsdm)), function(i)
    rowMeans(vsdm[,colnames(vsdm) == i]))
  
  # Generate column annotations
  annotation <- colnames(vsdmmean)
  annotation <- as.data.frame(annotation)
  annotation$subfield <- sapply(strsplit(as.character(annotation$annotation),'\\.'), "[", 1)
  annotation$t1 <- sapply(strsplit(as.character(annotation$annotation),'\\.'), "[", 2)
  annotation$t2 <- sapply(strsplit(as.character(annotation$annotation),'\\.'), "[", 3)
  annotation  <- unite_(annotation, "treatment", c("t1","t2"), sep = " ")
  row.names(annotation) <- annotation$annotation
  annotation$annotation <- NULL

  pheatmap(cor(vsdmmean),
           annotation_names_row = F, annotation_names_col = F,
           main= mysubtitle,
           color = viridis(20),
           show_rowname= F, show_colnames = F,
           display_numbers = T,
           legend_breaks = c(0.50, 0.85, 0.9, 0.95, 1.0),
           annotation  = annotation,
           annotation_colors = pheatmapcolors,
           treeheight_row = 0
  )
}



## mean standard deviation plot

plotmeansd <- function(mymeandev, mybehavior, myylab, mycolors){
  
  mymeandev %>% 
    filter(behavior == mybehavior) %>% 
    ggplot(aes(x=, TrainSessionComboNum, y=m, color=APA2)) + 
    geom_errorbar(aes(ymin=m-se, ymax=m+se, color=APA2), width=.1) +
    geom_point(size = 1.5) +
    geom_line() +
    theme_minimal(base_size = 8) + 
    theme(legend.position = "none",
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank()) + 
    scale_color_manual(values = mycolors) +
    scale_x_continuous( breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                        labels = c( "P", "T1", "T2", "T3",
                                    "Rt", "T4", "T5", "T6", "Rn")) +
    labs(x = NULL, y = myylab) 
}


returndds <- function(mytissue, mytreatment){
  print(mytissue)
  colData <- colData %>% 
    filter(subfield %in% c(mytissue),
           APA2 %in% c(mytreatment))  %>% 
    droplevels()
  
  savecols <- as.character(colData$RNAseqID) 
  savecols <- as.vector(savecols) 
  countData <- countData %>% dplyr::select(one_of(savecols)) 
  
  ## create DESeq object using the factors subfield and APA
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ APA2)
  
  dds <- dds[ rowSums(counts(dds)) > 1, ]  # Pre-filtering genes with 0 counts
  dds <- DESeq(dds, parallel = TRUE)
  return(dds)
}


returnvsds <- function(mydds, vsdfilename){
  dds <- mydds
  vsd <- vst(dds, blind=FALSE) ## variance stabilized
  print(head(assay(vsd),3))
  #return(write.csv(assay(vsd), file = vsdfilename, row.names = T))
  return(vsd)
}



# from corrr 

#' @export
network_plot.cor_df <- function(rdf,
                                min_cor = .30,
                                legend = TRUE,
                                colours = c("indianred2", "white", "skyblue1"),
                                repel = TRUE,
                                curved = TRUE,
                                colors) {
  
  if (min_cor < 0 || min_cor > 1) {
    stop ("min_cor must be a value ranging from zero to one.")
  }
  
  if (!missing(colors))
    colours <- colors
  
  rdf <-  as_matrix(rdf, diagonal = 1)
  distance <- sign(rdf) * (1 - abs(rdf))
  
  # Use multidimensional Scaling to obtain x and y coordinates for points.
  points <- data.frame(stats::cmdscale(abs(distance)))
  colnames(points) <-  c("x", "y")
  points$id <- rownames(points)
  
  # Create a proximity matrix of the paths to be plotted.
  proximity <- abs(rdf)
  proximity[upper.tri(proximity)] <- NA
  diag(proximity) <- NA
  proximity[proximity < min_cor] <- NA
  
  # Produce a data frame of data needed for plotting the paths.
  n_paths <- sum(!is.na(proximity))
  paths <- data.frame(matrix(nrow = n_paths, ncol = 6)) 
  colnames(paths) <- c("x", "y", "xend", "yend", "proximity", "sign")
  path <- 1
  for(row in 1:nrow(proximity)) {
    for(col in 1:ncol(proximity)) {
      path_proximity <- proximity[row, col]
      if (!is.na(path_proximity)) {
        path_sign <- sign(distance[row, col])
        x    <- points$x[row]
        y    <- points$y[row]
        xend <- points$x[col]
        yend <- points$y[col]
        paths[path, ] <- c(x, y, xend, yend, path_proximity, path_sign)
        path <- path + 1
      }
    }
  }
  
  plot_ <- list(
    # For plotting paths
    if (curved) geom_curve(data = paths,
                           aes(x = x, y = y, xend = xend, yend = yend,
                               alpha = proximity, size = proximity,
                               colour = proximity*sign)), 
    if (!curved) geom_segment(data = paths,
                              aes(x = x, y = y, xend = xend, yend = yend,
                                  alpha = proximity, size = proximity,
                                  colour = proximity*sign)), 
    scale_alpha(limits = c(0, 1)),
    scale_size(limits = c(0, 1)),
    scale_colour_gradientn(limits = c(-1, 1), colors = colours),
    # Plot the points
    geom_point(data = points,
               aes(x, y),
               size = 2, shape = 19, colour = "black"),
    # Plot variable labels
    if (repel) ggrepel::geom_text_repel(data = points,
                                        aes(x, y, label = id),
                                         size = 3,
                                        segment.size = 0.0,
                                        segment.color = "black"),
    if (!repel) geom_text(data = points,
                          aes(x, y, label = id),
                          fontface = 'italic', size = 3),
    # expand the axes to add space for curves
    expand_limits(x = c(min(points$x) - .1,
                        max(points$x) + .1),
                  y = c(min(points$y) - .1,
                        max(points$y) + .1)
    ),
    # Theme and legends
    theme_void(),
    guides(size = "none", alpha = "none"),
    if (legend)  labs(colour = NULL),
    if (!legend) theme(legend.position = "none")
  )
  
  ggplot() + plot_
  
}

