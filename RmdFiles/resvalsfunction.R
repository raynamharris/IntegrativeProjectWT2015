resvals <- function(contrastvector, mypval){
  res <- results(dds, contrast = c(contrastvector[1],contrastvector[2],contrastvector[3]), independentFiltering = F)
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
  res <- results(dds, contrast = c(contrastvector[1],contrastvector[2],contrastvector[3]), independentFiltering = F)
  vals <- cbind(res$pvalue)
  pvalcolname <- as.character(paste(contrastvector[2], "vs",contrastvector[3], sep=" "))
  colnames(vals) <- c(pvalcolname)
  vals <- as.data.frame(vals)
  histogram <- hist(vals)
  return(histogram)
  print(histogram)
}

