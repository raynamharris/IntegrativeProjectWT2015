## PCA ----

makepcadf <- function(data){
  Z <- data %>%
    select(`TotalPath.Arena.`:`ShockPerEntrance`) # columns 9 to 47
  Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]
  pc = prcomp(Z, scale.=TRUE)
  loadings <- pc$rotation
  scores <- pc$x
  #get ready for ggplot
  scoresdf <- as.data.frame(scores)
  scoresdf$ID <-  data$ID
  scoresdf$treatment <- data$treatment
  scoresdf$TrainSessionComboNum <- data$TrainSessionComboNum
  scoresdf <- scoresdf %>% select(ID, treatment,TrainSessionComboNum, PC1:PC10)
  return(scoresdf)
}

getpcsd <- function(data){
  Z <- rentention %>%
    select(SdevSpeedArena:Speed2)
  Z <- Z[,apply(Z, 2, var, na.rm=TRUE) != 0]
  pc = prcomp(Z, scale.=TRUE)
  pcsummary <- summary(pc)
  return(pcsummary)
}

## make plot of behavior across time

meansdplots <- function(df, myylab, ybreaks, ylims){
  myplot <- ggplot(df, 
                   aes(x=, TrainSessionComboNum, y=m, color=treatment)) + 
    geom_errorbar(aes(ymin=m-se, ymax=m+se, color=treatment), width=.1) +
    geom_line() +
    geom_point(size = 1.5, aes(alpha = TrainSessionComboNum)) +
    labs(subtitle = " ") +
    scale_y_continuous(name= myylab,
                       breaks = ybreaks,
                       limits = ylims) +
    scale_x_continuous(name= "training session", 
                       breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                       labels = c( "P", "T1", "T2", "T3",
                                   "Rt", "T4", "T5", "T6", "Rn")) +
    scale_alpha_continuous( breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
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


## pca from facto extra


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