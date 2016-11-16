## plot function for single behavior with wide data ----
onebehavior <- function(data, xcol, ycol, yaxislabel, colorcode){
  plot <- data %>% 
    ggplot(aes_string(x=xcol, y=ycol, color=colorcode)) +
    geom_point(size=1) + geom_jitter() +
    stat_smooth(alpha=0.5)  +
    theme_cowplot(font_size = 12, line_size = 0.5) + 
    background_grid(major = "xy", minor = "none") + 
    theme(axis.text.x = element_text(angle=60, vjust=0.5)) +
    scale_colour_manual(name="APA Training", values=c( "#b2182b","#f4a582",  "#878787"),
                        breaks = c("Yoked", "Same", "Conflict")) +
    scale_y_continuous(name=yaxislabel) + 
    scale_x_continuous(name =NULL, 
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                       labels=c("1" = "Hab", "2" = "T1", "3" = "T2", 
                                "4" = "T3", "5" = "Retest", "6" = "T4/C1",
                                "7" = "T5/C2", "8" = "T6/C3", "9"= "Retention"))
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
    scale_colour_manual(name="APA Training", values=c("#b2182b","#f4a582", "#878787"),
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
    APA =  c(Yoked = (values=c("#878787")), Same = (values=c("#f4a582")),
           Conflict = (values=c("#b2182b"))))
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
makecorrelationheatmap <- function(data){
  datacols <- data[c(19:59)];
  cormat <- cor(datacols);
  plot <- pheatmap(cormat, 
                   show_colnames=FALSE, show_rownames=TRUE, border_color ="grey60")
  return(plot)         
} 
