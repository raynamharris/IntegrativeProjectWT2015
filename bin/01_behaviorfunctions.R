## box plot ----
myboxplot <- function(data, xcol, ycol, yaxislabel, colorcode, session, title){
  plot <- data %>% filter(TrainSessionCombo==session) %>%  droplevels() %>%
    ggplot(aes_string(x=xcol, y=ycol, fill=colorcode)) + 
    geom_boxplot() +  
    scale_fill_manual(name="APA Training", 
                        values=c("#8073ac","#e08214",  "#7f3b08"),
                        breaks = c("Yoked", "Same", "Conflict")) +
    scale_y_continuous(name=yaxislabel) +
    scale_x_discrete(name=NULL) +   
    theme_cowplot(font_size = 20) +
    theme(legend.position="bottom") +
    theme(legend.text = element_text(size = 10)) + 
    theme(legend.title=element_blank()) 
  return(plot)
}


## plot function for single behavior with wide data ----
onebehavior <- function(data, xcol, ycol, yaxislabel, colorcode){
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size=1) + geom_jitter() +
    stat_smooth(alpha=0.5)  +
    theme_cowplot(font_size = 20, line_size = 0.5) + 
    background_grid(major = "xy", minor = "none") + 
    theme(axis.text.x = element_text(angle=60, vjust=0.5)) +
    scale_colour_manual(name="APA Training", values=c("#8073ac","#e08214",  "#7f3b08"),
                        breaks = c("Yoked", "Same", "Conflict")) +
    scale_y_continuous(name=yaxislabel) + 
    scale_x_continuous(name =NULL, 
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                       labels=c("1" = "Habituation", "2" = "T1", "3" = "T2", 
                                "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                                "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention"))+
    theme(legend.position="bottom") +
    theme(legend.title=element_blank()) 
  return(plot)
}


onebehaviornolegend <- function(data, xcol, ycol, yaxislabel, colorcode){
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size=1) + geom_jitter() +
    stat_smooth(alpha=0.5)  +
    theme_cowplot(font_size = 12, line_size = 0.5) + 
    background_grid(major = "xy", minor = "none") + 
    theme(axis.text.x = element_text(angle=60, vjust=0.5)) +
    scale_colour_manual(name="APA Training", values=c("#8073ac","#e08214",  "#7f3b08"),
                        breaks = c("Yoked", "Same", "Conflict")) +
    scale_y_continuous(name=yaxislabel) + 
    scale_x_continuous(name =NULL, 
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                       labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                                "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                                "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention")) +
    theme(legend.position="none") 
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
    APA =  c(Yoked = (values=c("#8073ac")), Same = (values=c("#e08214")),
           Conflict = (values=c("#7f3b08"))))
#now plot the heatmap
  plot <- pheatmap(averagedata, 
         show_colnames=FALSE, show_rownames=TRUE,
         annotation_col = columnannotations, 
         annotation_colors = columnannotationcolors,
         annotation_names_col=TRUE,
         fontsize = 15, fontsize_row = 8, fontsize_col = 10,
         cellwidth=75, 
         height = 3,
         border_color = "grey60"
         )
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

makepcaplot <- function(data,xcol,ycol,colorcode){
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size = 5) +
    stat_ellipse(level = 0.95)+
    theme_cowplot(font_size = 20) + 
    theme(strip.background = element_blank()) +
    scale_colour_manual(name="APA Training",
                        values=c("#8073ac","#e08214",  "#7f3b08"),
                        breaks=c("Yoked", "Same", "Conflict")) #+ 
    #labs(x = "Behavior xcol", y = "Behavior ycol") 
  return(plot)
}
