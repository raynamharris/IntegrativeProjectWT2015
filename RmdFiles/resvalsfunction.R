resvals <- function(contrastvector, mypadj){
  res <- results(dds, contrast = c(contrastvector[1],contrastvector[2],contrastvector[3]), independentFiltering = F)
  sum <- sum(res$padj < mypadj, na.rm = TRUE)
  print(sum)
  vals <- cbind(res$pvalue, res$padj)
  pvalcolname <- as.character(paste("pval",contrastvector[1],contrastvector[2],contrastvector[3], sep=""))
  padjcolname <- as.character(paste("padj",contrastvector[1],contrastvector[2],contrastvector[3], sep=""))
  colnames(vals) <- c(pvalcolname, padjcolname)
  return(vals)
}

respadjfold <- function(contrastvector){
  res <- results(dds, contrast = c(contrastvector[1],contrastvector[2],contrastvector[3]), independentFiltering = F)
  vals <- cbind(res$log2FoldChange, res$pval)
  foldcolname <- as.character(paste("foldchange",contrastvector[1],contrastvector[2],contrastvector[3], sep=""))
  pvalcolname <- as.character(paste("pval",contrastvector[1],contrastvector[2],contrastvector[3], sep=""))
  colnames(vals) <- c(foldcolname, pvalcolname)
  return(vals)
}
