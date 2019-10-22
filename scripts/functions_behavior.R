theme_ms <- function () { 
  theme_classic(base_size = 7) +
    theme(
      panel.grid.major  = element_blank(),  # remove major gridlines
      panel.grid.minor  = element_blank(),  # remove minor gridlines
      plot.title = element_text(hjust = 0.5, face = "bold") # center & bold 
    )
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

## make plot of behavior across time

meansdplots <- function(df, myylab, ybreaks, ylims){
  myplot <- ggplot(df, 
                   aes(x=, TrainSessionComboNum, y=m, color=treatment)) + 
    geom_errorbar(aes(ymin=m-se, ymax=m+se, color=treatment), width=.1) +
    geom_point(size = 1.5) +
    geom_line() +
    labs(subtitle = " ") +
    scale_y_continuous(name= myylab,
                       breaks = ybreaks,
                       limits = ylims) +
    scale_x_continuous(name= "training session", 
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                       labels = c( "P", "T1", "T2", "T3",
                                   "Rt", "T4", "T5", "T6", "Rn")) +
    theme_ms() +
    scale_color_manual(values = colorvalAPA00,
                       name  = NULL)  +
    theme(legend.position = "bottom",
          legend.justification=c(0,0),
          legend.text=element_text(size=5))
  return(myplot)
}  