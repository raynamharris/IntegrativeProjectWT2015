## make a scale data matrix for heatmap from all sessions 
makescaledaveragedata <- function(data){
  longdata <- melt(data, id = c(1:3));  #longdata <- melt(behavior, id = c(1:3))
  longdata <- longdata %>% drop_na();
  # then widen with group averages, add row names, scale, and transpose
  longdata$APAsession <- as.factor(paste(longdata$APA2,longdata$TrainSessionCombo, sep="_"))
  averagedata <- dcast(longdata, APAsession ~ variable, value.var= "value", fun.aggregate=mean);
  rownames(averagedata) <- averagedata$APAsession;    
  averagedata[1] <- NULL;
  scaledaveragedata <- scale(averagedata)
  scaledaveragedata <- t(scaledaveragedata)
  scaledaveragedata <- scaledaveragedata[-1,]
  return(scaledaveragedata)
}  


# make a heat map annotation with only APA2
makecolumnannotations <- function(data){  
  columnannotations <- as.data.frame(colnames(data))
  names(columnannotations)[names(columnannotations)=="colnames(data)"] <- "column"
  rownames(columnannotations) <- columnannotations$column
  columnannotations$APA2 <- sapply(strsplit(as.character(columnannotations$column),'\\_'), "[", 1)
  columnannotations$column <- NULL
  return(columnannotations)
}

# make a heat map annotation with session and APA2
makecolumnannotations2 <- function(data){  
  columnannotations <- as.data.frame(colnames(data))
  names(columnannotations)[names(columnannotations)=="colnames(data)"] <- "column"
  rownames(columnannotations) <- columnannotations$column
  columnannotations$APA2 <- sapply(strsplit(as.character(columnannotations$column),'\\_'), "[", 1)
  columnannotations$Session <- sapply(strsplit(as.character(columnannotations$column),'\\_'), "[", 2)
  columnannotations$Session <- as.factor(columnannotations$Session)
  columnannotations$Session <- revalue(columnannotations$Session, c("T4" = "T4_C1")) 
  columnannotations$Session <- revalue(columnannotations$Session, c("T5" = "T5_C2")) 
  columnannotations$Session <- revalue(columnannotations$Session, c("T6" = "T6_C3")) 
  columnannotations$column <- NULL
  return(columnannotations)
}

## correlation heatmat ----
makecorrelationheatmap <- function(data, APAgroup, clusterTF){
  dataslim <- data %>% filter(APA==APAgroup)
  datacols <- dataslim[c(19:59)];
  cormat <- cor(datacols);
  plot <- pheatmap(cormat, 
                   show_colnames=FALSE, show_rownames=TRUE, border_color ="grey60", 
                   main=APAgroup, cluster_rows = clusterTF, cluster_cols = clusterTF)
  return(plot)         
} 

## PCA ----
makelongdata <- function(data){
  #first melt the data to make long
  longdata <- melt(data, id = c(1:18));
  longdata <- longdata %>% drop_na();
  longdata$bysession <- as.factor(paste(longdata$TrainSessionCombo, longdata$variable, sep="_"));
  longdata <- dcast(longdata, ID + APA2 ~ bysession, value.var= "value", fun.aggregate = mean)
  return(longdata)
}

makepcadf <- function(data){
  #first melt the data to make long
  longdata <- melt(data, id = c(1:18));
  longdata <- longdata %>% drop_na();
  longdata$bysession <- as.factor(paste(longdata$TrainSessionCombo, longdata$variable, sep="_"));
  longdata <- dcast(longdata, ID + APA2 ~ bysession, value.var= "value", fun.aggregate = mean)
  # calculate and save PCs
  Z <- longdata[,3:322]
  Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]
  pc = prcomp(Z, scale.=TRUE)
  loadings <- pc$rotation
  scores <- pc$x
  #get ready for ggplot
  scoresdf <- as.data.frame(scores)
  scoresdf$ID <-  longdata$ID
  scoresdf$APA2 <- longdata$APA2
  return(scoresdf)
}


makepcaloadingsdf <- function(data){
  #first melt the data to make long
  longdata <- melt(data, id = c(1:18));
  longdata <- longdata %>% drop_na();
  longdata$bysession <- as.factor(paste(longdata$TrainSessionCombo, longdata$variable, sep="_"));
  longdata <- dcast(longdata, ID + APA2 ~ bysession, value.var= "value", fun.aggregate = mean)
  # calculate and save PCs
  Z <- longdata[,3:322]
  Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]
  pc = prcomp(Z, scale.=TRUE)
  loadings <- pc$rotation
  return(loadings)
}

mkrotationdf <- function(data){
  #first melt the data to make long
  longdata <- melt(data, id = c(1:18));
  longdata <- longdata %>% drop_na();
  longdata$bysession <- as.factor(paste(longdata$TrainSessionCombo, longdata$variable, sep="_"));
  longdata <- dcast(longdata, ID + APA2 ~ bysession, value.var= "value", fun.aggregate = mean)
  # calculate and save PCs
  Z <- longdata[,3:322]
  Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]
  pc = prcomp(Z, scale.=TRUE)
  rotationdf <- data.frame(pc$rotation, variable=row.names(pc$rotation))
  #str(rotationdf)
  return(rotationdf)
}


makepcaplot <- function(data,xcol,ycol,colorcode){
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode, label=data$ID)) +
    geom_point(size = 5) +
    stat_ellipse(level = 0.95)+
    theme_cowplot(font_size = 20) + 
    theme(strip.background = element_blank()) +
    scale_colour_manual(name=NULL,
                        values=colorvalAPA) +
    theme(legend.position="none")
    
  return(plot)
}

makepcaplotwithpercent <- function(data,xcol,ycol,colorcode, newylab, newxlab){
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode, label=data$ID)) +
    geom_point(size = 1) +
    stat_ellipse(level = 0.95) +
    theme_cowplot(font_size = 15, line_size = 0.5) + 
    theme(strip.background = element_blank()) +
    scale_colour_manual(name=NULL,
                        values=colorvalAPA2) +
    theme(legend.position="none") +
    ylab(newylab) + xlab(newxlab)
  return(plot)
}

