## box plot ----
myboxplot <- function(data, xcol, ycol, colorcode, session, yaxislabel){
  plot <- data %>% filter(TrainSessionCombo==session) %>%  droplevels() %>%
    ggplot(aes_string(x=xcol, y=ycol, fill=colorcode)) + 
    geom_boxplot() + 
    scale_fill_manual(name="Group", 
                        values=colorvalAPA,
                        breaks = c("Yoked", "Same", "Conflict")) +
    #scale_y_continuous(name=NULL) +
    scale_y_continuous(name=yaxislabel) + 
    scale_x_discrete(name=NULL) +   
    theme_cowplot(font_size = 15) +
    #theme(legend.position="bottom") +
    theme(legend.text = element_text(size = 10)) 
  return(plot)
}

## box plot no legend ----
myboxplotnolegend <- function(data, xcol, ycol, colorcode, session, yaxislabel){
  plot <- data %>% filter(TrainSessionCombo==session) %>%  droplevels() %>%
    ggplot(aes_string(x=xcol, y=ycol, fill=colorcode)) + 
    geom_boxplot() + 
    scale_fill_manual(name="Group", 
                      values=colorvalAPA,
                      breaks = c("Yoked", "Same", "Conflict")) +
    #scale_y_continuous(name=NULL) +
    scale_y_continuous(name=yaxislabel) + 
    scale_x_discrete(name=NULL) +   
    theme_cowplot(font_size = 15) +
    #theme(legend.position="bottom") +
    theme(legend.text = element_text(size = 10))  +
    theme(legend.position="none") 
  return(plot)
}

## plot function for single behavior with wide data ----
onebehavior <- function(data, xcol, ycol, yaxislabel, colorcode){
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size=1) + geom_jitter() +
    stat_smooth(alpha=0.5)  +
    theme_cowplot(font_size = 15, line_size = 0.5) + 
    #background_grid(major = "xy", minor = "none") + 
    theme(axis.text.x = element_text(angle=60, vjust=0.5)) +
    scale_colour_manual(name="APA Training", values=colorvalAPA,
                        breaks = c("Yoked", "Same", "Conflict")) +
    scale_y_continuous(name=yaxislabel) + 
    scale_x_continuous(name = "Training Session", 
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                       labels=c("1" = "Habituation", "2" = "T1", "3" = "T2", 
                                "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                                "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
    theme(legend.position="top",
          legend.title=element_blank()) 
  return(plot)
}




## make a heatmap from group averages ----
makeagroupaveheatmap <- function(data){
  #first melt the data to make long
  longdata <- melt(data, id = c(1:18));
  longdata <- longdata %>% drop_na();
  # then widen with group averages, add row names, scale, and transpose
  averagedata <- dcast(longdata, APA ~ variable, value.var= "value", fun.aggregate=mean);
  rownames(averagedata) <- averagedata$APA;    
  averagedata[1] <- NULL;
  averagedata <- scale(averagedata)
  averagedata <- t(averagedata)
  #next lines create the annotations
  columnannotations <- as.data.frame(colnames(averagedata))
  names(columnannotations)[names(columnannotations)=="colnames(averagedata)"] <- "APA"
  rownames(columnannotations) <- columnannotations$APA
  columnannotationcolors = list(
    APA =  c(Yoked = (values=c("#404040")), Same = (values=c("#f4a582")),
           Conflict = (values=c("#ca0020"))))
  #now plot the heatmap
  plot <- pheatmap(averagedata, 
         show_colnames=FALSE, show_rownames=TRUE,
         annotation_col = columnannotations, 
         annotation_colors = columnannotationcolors,
         annotation_names_col=TRUE,
         fontsize = 15, fontsize_row = 8, fontsize_col = 10,
         #cellwidth=75, 
         #height = 3,
         border_color = "grey60",
         annotation_legend = FALSE)
  return(plot)
}

## make a heatmap from all sessions ----
makesessionheatmap <- function(data){
  #first melt the data to make long
  longdata <- melt(data, id = c(1:18));  #longdata <- melt(behavior, id = c(1:18))
  longdata <- longdata %>% drop_na();
  # then widen with group averages, add row names, scale, and transpose
  longdata$APAsession <- as.factor(paste(longdata$APA,longdata$TrainSessionCombo, sep="_"))
  averagedata <- dcast(longdata, APAsession ~ variable, value.var= "value", fun.aggregate=mean);
  rownames(averagedata) <- averagedata$APAsession;    
  averagedata[1] <- NULL;
  averagedata <- scale(averagedata)
  averagedata <- t(averagedata)
  #next lines create the annotations
  columnannotations <- as.data.frame(colnames(averagedata))
  names(columnannotations)[names(columnannotations)=="colnames(averagedata)"] <- "column"
  rownames(columnannotations) <- columnannotations$column
  columnannotations$APA <- sapply(strsplit(as.character(columnannotations$column),'\\_'), "[", 1)
  columnannotations$Session <- sapply(strsplit(as.character(columnannotations$column),'\\_'), "[", 2)
  columnannotations$Session <- as.factor(columnannotations$Session)
  columnannotations$Session <- revalue(columnannotations$Session, c("T4" = "T4_C1")) 
  columnannotations$Session <- revalue(columnannotations$Session, c("T5" = "T5_C2")) 
  columnannotations$Session <- revalue(columnannotations$Session, c("T6" = "T6_C3")) 
    columnannotations$column <- NULL
    columnannotationcolors = list(
    APA =  c(Yoked = (values=c("#404040")), 
             Same = (values=c("#f4a582")),
             Conflict = (values=c("#ca0020"))),
    Session =  c(Hab = (values=c("#eff3ff")), 
             T1 = (values=c("#bdd7e7")),
             T2 = (values=c("#6baed6")),
             T3 = (values=c("#3182bd")),
             Retest = (values=c("#08519c")),
             T4_C1 = (values=c("#cbc9e2")),
             T5_C2 = (values=c("#9e9ac8")),
             T6_C3 = (values=c("#756bb1")),
             Retention = (values=c("#54278f"))))
        #now plot the heatmap
  plot <- pheatmap(averagedata, 
                   show_colnames=FALSE, show_rownames=TRUE,
                   annotation_col = columnannotations, 
                   annotation_colors = columnannotationcolors,
                   annotation_names_col=TRUE,
                   fontsize = 15, fontsize_row = 8, fontsize_col = 10,
                   #cellwidth=75, 
                   #height = 3,
                   border_color = "grey60",
                   annotation_legend = TRUE)
  return(plot)
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

mkrotationdata <- function(data){
  #first melt the data to make long
  longdata <- melt(data, id = c(1:18));
  longdata <- longdata %>% drop_na();
  longdata$bysession <- as.factor(paste(longdata$TrainSessionCombo, longdata$variable, sep="_"));
  longdata <- dcast(longdata, ID + APA ~ bysession, value.var= "value", fun.aggregate = mean)
  # calculate and save PCs
  Z <- longdata[,3:371]
  Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]
  pc = prcomp(Z, scale.=TRUE)
  rotation_data <- data.frame(pc$rotation, variable=row.names(pc$rotation))
  str(rotation_data)
  return(rotation_data)
}


makepcaplot <- function(data,xcol,ycol,colorcode){
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode, label=data$ID)) +
    geom_point(size = 5) +
    stat_ellipse(level = 0.95)+
    theme_cowplot(font_size = 20) + 
    theme(strip.background = element_blank()) +
    scale_colour_manual(name="APA Training",
                        values=colorvalAPA,
                        breaks=c("Yoked", "Same", "Conflict")) +
    theme(legend.position="none") +
    geom_text()
  return(plot)
}

colorvalAPA <- c("#404040", "#f4a582", "#ca0020")
#ca0020 ## red - conflict 
#f4a582 ## pink - same
#404040 ## darkgrey - yoked

colorvalAPA2 <- c("##bababa", "#404040","##f4a582", "#ca0020")
#ca0020 ## red - conflict 
#f4a582 ## pink - same
#bababa ## grey - yoked to same
#404040 ## darkgrey - yoked to conflict


