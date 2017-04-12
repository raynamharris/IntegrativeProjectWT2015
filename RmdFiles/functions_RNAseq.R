## I've taken the pca function from DESeq2 and elaborated it so that I could extract up to 6 PCs

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
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4],PC5 = pca$x[, 5],PC6 = pca$x[, 6],group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:6]
    return(d)
  }
}


plotPCs <- function(df, xcol, ycol, aescolor, colorname, aesshape, shapename, colorvalues){
  ggplot(df, aes(df[xcol], df[ycol], color=aescolor, shape=aesshape)) +
    geom_point(size=5) +
    xlab(paste0("PC", xcol, ": ", percentVar[xcol],"% variance")) +
    ylab(paste0("PC", ycol, ": ", percentVar[ycol],"% variance")) +
    theme_classic() +
    #stat_ellipse(level = 0.95, (aes(color=aescolor)),size=1) + 
    scale_colour_manual(name=colorname, values=c(colorvalues))+
    scale_shape_discrete(name=shapename) +
    theme(axis.text = element_text(size=14),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          legend.title = element_text(size=16),
          legend.text = element_text(size=14))
}
