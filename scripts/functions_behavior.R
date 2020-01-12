
## Mean and standard deviation plots  of avoidance behavior

meansdplots <- function(df, myylab, ybreaks, ylims){
  myplot <- ggplot(df, 
                   aes(x=, trialNum, y=m, color=treatment)) + 
    geom_errorbar(aes(ymin=m-se, ymax=m+se, color=treatment), width=.1) +
    geom_line() +
    geom_point(size = 1.5) +
    labs(subtitle = " ") +
    scale_y_continuous(name= myylab,
                       breaks = ybreaks,
                       limits = ylims) +
    scale_x_continuous(name= "trial", 
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                       labels = c( "P", "T1", "T2", "T3",
                                   "Rt", "T4", "T5", "T6", "Rn")) +
    scale_alpha_continuous( breaks = c(1, 2, 3)) +
    theme_ms() +
    scale_color_manual(values = treatmentcolors,
                       name  = NULL,
                       labels = levelstreatmentlegend)  +
    theme(legend.position = "bottom",
          legend.justification=c(0,0),
          legend.text=element_text(size=5))
  return(myplot)
}  


## make PCA dataframes

makepcadf <- function(data){
  Z <- data %>%
    select(TotalPath.Arena.:AnnularKurtosis) # columns 9 to 47
  Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]
  pc = prcomp(Z, scale.=TRUE)
  loadings <- pc$rotation
  scores <- pc$x
  #get ready for ggplot
  scoresdf <- as.data.frame(scores)
  scoresdf$ID <-  data$ID
  scoresdf$treatment <- data$treatment
  scoresdf$trialNum <- data$trialNum
  scoresdf$Day <- data$Day
  scoresdf <- scoresdf %>% select(ID, treatment, trialNum, Day, PC1:PC10)
  return(scoresdf)
}


## modified function for PCA from facto extra


fviz_contrib_rmh <- function(X, choice = c("row", "col", "var", "ind", "quanti.var", "quali.var", "group", "partial.axes"),
                         axes=1, fill="steelblue", color = "steelblue", 
                         sort.val = c("desc", "asc", "none"), top = Inf,
                         xtickslab.rt = 45, ggtheme = theme_minimal(), ...)
{
  
  sort.val <- match.arg(sort.val)
  choice = match.arg(choice)
  
  
  dd <- facto_summarize(X, element = choice, result = "contrib", axes = axes)
  contrib <- dd$contrib
  names(contrib) <-rownames(dd)
  
  # expected Average contribution
  theo_contrib <- 100/length(contrib)
  if(length(axes) > 1) {
    # Adjust variable contributions by the Dimension eigenvalues
    eig <- get_eigenvalue(X)[axes,1]
    theo_contrib <- sum(theo_contrib*eig)/sum(eig)
  }
  df <- data.frame(name = factor(names(contrib), levels = names(contrib)), contrib = contrib)
  
  # Define color if quanti.var
  if(choice == "quanti.var") {
    df$Groups <- .get_quanti_var_groups (X)
    if(missing(fill)) fill <- "Groups"
    if(missing(color)) color <- "Groups"
  }
  
  p <- ggpubr::ggbarplot(df, x = "name", y = "contrib", fill = fill, color = color,
                         sort.val = sort.val, top = top,
                         xtickslab.rt = xtickslab.rt, ggtheme = ggtheme,
                         sort.by.groups = FALSE, ...
  )
  
  #   p <- .ggbarplot(contrib, fill =fill, color = color,
  #                   sort.value = sort.val[1], top = top,
  #                   title = title, ylab ="Contributions (%)")+
  #     geom_hline(yintercept=theo_contrib, linetype=2, color="red")
  
  p
}

## ANOVAs

# twoway anova table function

twowayANOVAfor3measures <- function(mydata, mydescription){
  apa1 <- apa.aov.table(aov(NumEntrances ~ treatment * trial, data=mydata))
  apa1df <- as.data.frame(apa1$table_body)
  totaldf <- apa1df[5, 3]
  apa1df$df <- paste(apa1df$df, ", " , totaldf, sep = "")
  apa1df$ANOVA <- "NumEntrances ~ treatment * trial"
  apa1df
  
  apa2 <- apa.aov.table(aov(pTimeShockZone ~ treatment * trial, data=mydata))
  apa2df <- as.data.frame(apa2$table_body)
  apa2df$df <- paste(apa2df$df, ", " , totaldf, sep = "")
  apa2df$ANOVA <- "pTimeShockZone ~ treatment * trial"
  
  apa3 <- apa.aov.table(aov(Time1stEntr ~ treatment * trial, data=mydata))
  apa3df <- as.data.frame(apa3$table_body)
  apa3df$df <- paste(apa3df$df, ", " , totaldf, sep = "")
  apa3df$ANOVA <- "Time1stEntr ~ treatment * trial"
  apa3df
  
  apa123 <- as.data.frame(rbind(apa1df,apa2df,apa3df))
  apa123$trials <- mydescription
  apa123 <- apa123 %>%
    select(trials, ANOVA, Predictor, df, "F", p) %>%
    filter(!Predictor %in% c("(Intercept)", "Error"))
  
  return(apa123)
}


onewayANOVAfor3measures <- function(mydata, whichtrial, mydescription){
  
  mydata <- mydata %>% filter(trial == whichtrial)
  
  apa1 <- apa.aov.table(aov(NumEntrances ~ treatment , data=mydata))
  apa1df <- as.data.frame(apa1$table_body)
  totaldf <- apa1df[3, 3]
  apa1df$df <- paste(apa1df$df, ", " , totaldf, sep = "")
  apa1df$ANOVA <- "NumEntrances ~ treatment"
  apa1df
  
  apa2 <- apa.aov.table(aov(pTimeShockZone ~ treatment , data=mydata))
  apa2df <- as.data.frame(apa2$table_body)
  apa2df$df <- paste(apa2df$df, ", " , totaldf, sep = "")
  apa2df$ANOVA <- "pTimeShockZone ~ treatment"
  
  apa3 <- apa.aov.table(aov(Time1stEntr ~ treatment , data=mydata))
  apa3df <- as.data.frame(apa3$table_body)
  apa3df$df <- paste(apa3df$df, ", " , totaldf, sep = "")
  apa3df$ANOVA <- "Time1stEntr ~ treatment"
  apa3df
  
  apa123 <- as.data.frame(rbind(apa1df,apa2df,apa3df))
  apa123$trials <- mydescription
  apa123 <- apa123 %>%
    select(trials, ANOVA, Predictor, df, "F", p) %>%
    filter(!Predictor %in% c("(Intercept)", "Error"))
  
  return(apa123)
}
