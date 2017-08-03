## box plot ----
myboxplotlegendtop <- function(data, xcol, ycol, colorcode, session, yaxislabel){
  plot <- data %>% filter(TrainSessionCombo==session) %>%  droplevels() %>%
    ggplot(aes_string(x=xcol, y=ycol, fill=colorcode)) + 
    geom_boxplot() +
    scale_fill_manual(name="Group", 
                      values=colorvalAPA) +
    scale_y_continuous(name=yaxislabel) + 
    scale_x_discrete(name=NULL) +   
    theme_cowplot(font_size = 15) +
    theme(legend.justification=c(1,1), legend.position=c(1,1), 
          legend.title = element_blank(), axis.text.x = element_blank())
    return(plot)
}

## box plot no legend ----
myboxplotlegendbottom <- function(data, xcol, ycol, colorcode, session, yaxislabel){
  plot <- data %>% filter(TrainSessionCombo==session) %>%  droplevels() %>%
    ggplot(aes_string(x=xcol, y=ycol, fill=colorcode)) + 
    geom_boxplot() + 
    scale_fill_manual(name="Group", 
                      values=colorvalAPA) +
    scale_y_continuous(name=yaxislabel) + 
    #scale_x_discrete(name=NULL) +   
    theme_cowplot(font_size = 15) +
    #theme(legend.position="bottom") +
    theme(legend.text = element_text(size = 10))  +
    theme(legend.justification=c(1,0), legend.position=c(1,0), 
          legend.title = element_blank(), axis.text.x = element_blank())
  return(plot)
}

## plot function for single behavior with wide data ----
onebehavior <- function(data, xcol, ycol, yaxislabel, colorcode){
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size=1) + geom_jitter() +
    stat_smooth(alpha=0.25, method = "loess")  +
    theme_cowplot(font_size = 14, line_size = 0.5) + 
    #background_grid(major = "xy", minor = "none") + 
    scale_colour_manual(name="APA Training", values=colorvalAPA,
                        breaks = c("Yoked", "Same", "Conflict")) +
    scale_y_continuous(name=yaxislabel) + 
    scale_x_continuous(name = NULL, 
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                       labels=c("1" = "Hab.", "2" = "T1", "3" = "T2", 
                                "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                                "7" = "T5/C2", "8" = "T6/C3", "9"= "Reten.")) + 
    theme(legend.position="top") 
  return(plot)
}

onebehaviorhabtoretest <- function(data, xcol, ycol, yaxislabel, colorcode){
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size=1) + geom_jitter() +
    stat_smooth(alpha=0.5)  +
    theme_cowplot(font_size = 15, line_size = 0.5) + 
    #background_grid(major = "xy", minor = "none") + 
    scale_colour_manual(name="APA Training", values=colorvalAPA,
                        breaks = c("Control", "Consistent", "Conflict")) +
    scale_y_continuous(name=yaxislabel) + 
    scale_x_continuous(name = "Training Session", 
                       breaks = c(1, 2, 3, 4, 5),
                       labels=c("1" = "Habituation", "2" = "T1", "3" = "T2", 
                                "4" = "T3", "5" = "Retest")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1),
          axis.text.x = element_text(angle=60, vjust=0.5)) 
  return(plot)
}


onebehaviorc4toRentention <- function(data, xcol, ycol, yaxislabel, colorcode){
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size=1) + geom_jitter() +
    stat_smooth(alpha=0.5, method = "loess")  +
    theme_cowplot(font_size = 15, line_size = 0.5) + 
    #background_grid(major = "xy", minor = "none") + 
    scale_colour_manual(name="APA Training", values=colorvalAPA,
                        breaks = c("Control", "Consistent", "Conflict")) +
    scale_y_continuous(name=yaxislabel) + 
    scale_x_continuous(name = "Training Session", 
                       breaks = c(1, 2, 3, 4),
                       labels=c("1" = "T4_C1", "2" = "T5_C2", "3" = "T6_C3", 
                                "4" = "Retention")) +
    theme(legend.justification=c(1,1), legend.position=c(1,1),
          axis.text.x = element_text(angle=60, vjust=0.5)) 
  return(plot)
}

## make a heatmap from all sessions ----
makescaledaveragedata <- function(data){
  longdata <- melt(data, id = c(1:18));  #longdata <- melt(behavior, id = c(1:18))
  longdata <- longdata %>% drop_na();
  # then widen with group averages, add row names, scale, and transpose
  longdata$APAsession <- as.factor(paste(longdata$APA,longdata$TrainSessionCombo, sep="_"))
  averagedata <- dcast(longdata, APAsession ~ variable, value.var= "value", fun.aggregate=mean);
  rownames(averagedata) <- averagedata$APAsession;    
  averagedata[1] <- NULL;
  scaledaveragedata <- scale(averagedata)
  scaledaveragedata <- t(scaledaveragedata)
  #next lines create the annotations
  return(scaledaveragedata)
}  

makecolumnannotations <- function(data){  
  columnannotations <- as.data.frame(colnames(data))
  names(columnannotations)[names(columnannotations)=="colnames(data)"] <- "column"
  rownames(columnannotations) <- columnannotations$column
  columnannotations$APA <- sapply(strsplit(as.character(columnannotations$column),'\\_'), "[", 1)
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
  longdata <- dcast(longdata, ID + APA ~ bysession, value.var= "value", fun.aggregate = mean)
  return(longdata)
}

makepcadf <- function(data){
  #first melt the data to make long
  longdata <- melt(data, id = c(1:18));
  longdata <- longdata %>% drop_na();
  longdata$bysession <- as.factor(paste(longdata$TrainSessionCombo, longdata$variable, sep="_"));
  longdata <- dcast(longdata, ID + APA ~ bysession, value.var= "value", fun.aggregate = mean)
  # calculate and save PCs
  Z <- longdata[,3:371]
  Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]
  pc = prcomp(Z, scale.=TRUE)
  loadings <- pc$rotation
  scores <- pc$x
  #get ready for ggplot
  scoresdf <- as.data.frame(scores)
  scoresdf$ID <-  longdata$ID
  scoresdf$APA <- longdata$APA
  return(scoresdf)
}


makepcaloadingsdf <- function(data){
  #first melt the data to make long
  longdata <- melt(data, id = c(1:18));
  longdata <- longdata %>% drop_na();
  longdata$bysession <- as.factor(paste(longdata$TrainSessionCombo, longdata$variable, sep="_"));
  longdata <- dcast(longdata, ID + APA ~ bysession, value.var= "value", fun.aggregate = mean)
  # calculate and save PCs
  Z <- longdata[,3:371]
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
  longdata <- dcast(longdata, ID + APA ~ bysession, value.var= "value", fun.aggregate = mean)
  # calculate and save PCs
  Z <- longdata[,3:371]
  Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]
  pc = prcomp(Z, scale.=TRUE)
  rotationdf <- data.frame(pc$rotation, variable=row.names(pc$rotation))
  str(rotationdf)
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
                        values=colorvalAPA) +
    theme(legend.position="none") +
    ylab(newylab) + xlab(newxlab)
  return(plot)
}

