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

## calculate a group average and make it a matrix ----
makeagroupmatrix <- function(data){
  longdata <- melt(data, id = c(1:18));
  longdata <- longdata %>% drop_na();
  averagedata <- dcast(longdata, APA ~ variable, value.var= "value", fun.aggregate=mean);
  rownames(averagedata) <- averagedata$APA;    
  averagedata[1] <- NULL;
  averagedata <- scale(averagedata)
  averagedata <- t(averagedata)
  return(averagedata)
}  

##
