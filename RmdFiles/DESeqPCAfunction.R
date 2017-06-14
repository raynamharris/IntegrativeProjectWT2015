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


plotPC1PC2 <- function(aescolor, colorname, aesshape, shapename, colorvalues){
  ggplot(pcadata, aes(PC1, PC2, color=aescolor, shape=aesshape)) +
    geom_point(size=4) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    theme_classic() +
    #stat_ellipse(level = 0.95, (aes(color=aescolor)),size=1) + 
    scale_colour_manual(name=colorname, values=c(colorvalues))+
    scale_shape_discrete(name=shapename) +
    theme(axis.text = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size=12)) + 
    guides(colour = guide_legend(override.aes = list(shape = 15), order = 1),
           shape = guide_legend(order = 2))
}

plotPC2PC1 <- function(aescolor, colorname, aesshape, shapename, colorvalues){
  ggplot(pcadata, aes(PC2, PC1, color=aescolor, shape=aesshape)) +
    geom_point(size=3) +
    xlab(paste0("PC2: ",percentVar[2],"% variance")) +
    ylab(paste0("PC1: ",percentVar[1],"% variance")) +
    theme_classic() +
    #stat_ellipse(level = 0.95, (aes(color=aescolor)),size=1) + 
    scale_colour_manual(name=colorname, values=c(colorvalues))+
    scale_shape_discrete(name=shapename) +
    theme(axis.text = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size=12)) + 
    guides(colour = guide_legend(override.aes = list(shape = 15), order = 1),
           shape = guide_legend(order = 2))
}


plotPC2PC3 <- function(aescolor, colorname, aesshape, shapename, colorvalues){
  ggplot(pcadata, aes(PC2, PC3, color=aescolor, shape=aesshape)) +
    geom_point(size=4) +
    xlab(paste0("PC2: ",percentVar[2],"% variance")) +
    ylab(paste0("PC3: ",percentVar[3],"% variance")) +
    theme_classic() +
    #stat_ellipse(level = 0.95, (aes(color=aescolor)),size=1) + 
    scale_colour_manual(name=colorname, values=c(colorvalues))+
    scale_shape_discrete(name=shapename) +
    theme(axis.text = element_text(size=14),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          legend.title = element_text(size=16),
          legend.text = element_text(size=14)) +
    theme(legend.position = "none")
}

plotPC1PC3 <- function(aescolor, colorname, aesshape, shapename, colorvalues){
  ggplot(pcadata, aes(PC1, PC3, color=aescolor, shape=aesshape)) +
    geom_point(size=4) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC3: ",percentVar[3],"% variance")) +
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


plotPC2PC4 <- function(aescolor, colorname, aesshape, shapename, colorvalues){
  ggplot(pcadata, aes(PC2, PC4, color=aescolor, shape=aesshape)) +
    geom_point(size=3) +
    xlab(paste0("PC2: ",percentVar[2],"% variance")) +
    ylab(paste0("PC4: ",percentVar[4],"% variance")) +
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

plotPC3PC4 <- function(aescolor, colorname, aesshape, shapename, colorvalues){
  ggplot(pcadata, aes(PC3, PC4, color=aescolor, shape=aesshape)) +
    geom_point(size=5) +
    xlab(paste0("PC3: ",percentVar[3],"% variance")) +
    ylab(paste0("PC4: ",percentVar[4],"% variance")) +
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
