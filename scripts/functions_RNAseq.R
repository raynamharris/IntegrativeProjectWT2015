# functions for RNAseq analysis

## used in 02_rnaseqQC.Rmd

res_summary <- function(mycontrast){
  res <- results(dds, contrast = mycontrast, independentFiltering = T)
  print(mycontrast)
  print(sum(res$padj < 0.1, na.rm=TRUE))
  print(summary(res))
  print(head((res[order(res$padj),]), 5))
  cat("\n")
}

plot_tSNE <- function(myperplecity, mysubtitle){
  
  print(myperplecity)
  tsne_model <- Rtsne(euclidist, check_duplicates=FALSE, pca=TRUE, perplexity=myperplecity)
  tsne_df = as.data.frame(tsne_model$Y) 
  tsne_df <- cbind(colData, tsne_df)
  tsne_df$subfield <- factor(tsne_df$subfield, levels = c("DG", "CA3", "CA1"))
  tsne_df$treatment <- factor(tsne_df$treatment, levels = levelstreatment)
  
  tnseplot  <- tsne_df %>%
    ggplot(aes(x = V1, y = V2, shape = treatment, color = subfield, label = ID)) +
    geom_point(size = 2) +
    scale_color_manual(values = colorvalsubfield) +
    theme_ms()  +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.spacing.x = unit(0.01, 'cm'),
          legend.spacing.y = unit(0.01, 'cm')) +
    labs(x = "tSNE 1", y = "tSNE 2", subtitle = mysubtitle) +
    scale_shape_manual(aes(colour=colorvalsubfield), values=c(1, 16, 0, 15)) 
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


res_summary_subfield <- function(mydds, mycontrast){
  res <- results(mydds, contrast = mycontrast, independentFiltering = T)
  print(mycontrast)
  print(sum(res$padj < 0.1, na.rm=TRUE))
  print(summary(res))
  cat("\n")
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

# this is modified from DESeq
pcadataframe <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) {
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




# adapted from corrr, with modifications to size and color 

network_plot.cor_df <- function(rdf,
                                min_cor = .30,
                                legend = TRUE,
                                colours = c("skyblue1", "white", "indianred2"),
                                repel = TRUE,
                                curved = TRUE,
                                colors = c("skyblue1", "white", "indianred2") ) {
  
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
                               alpha = proximity, 
                               colour = proximity*sign)), 
    if (!curved) geom_segment(data = paths,
                              aes(x = x, y = y, xend = xend, yend = yend,
                                  alpha = proximity, 
                                  colour = proximity*sign)), 
    scale_alpha(limits = c(0, 1)),
    scale_size(limits = c(0, 1)),
    scale_colour_gradientn(limits = c(-1, 1), colors = colours),
    # Plot the points
    geom_point(data = points,
               aes(x, y),
               size = 1, shape = 19, colour = "black"),
    # Plot variable labels
    if (repel) ggrepel::geom_text_repel(data = points,
                                        aes(x, y, label = id),
                                         size = 2,
                                        segment.size = 0.0,
                                        segment.color = "black"),
    if (!repel) geom_text(data = points,
                          aes(x, y, label = id),
                          fontface = 'italic', size = 2),
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

