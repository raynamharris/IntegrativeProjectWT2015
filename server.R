#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

source("./scripts/figureoptions.R")

# table 1
table1 <- read_csv("data/table-1.csv")
 
# supp table 1
sup1 <- read_csv(file = "./data/suppltable-1.csv") %>%
  mutate(treatment = factor(treatment, levels = levelstreatment)) 


# supp table 2
sup2 <- read_csv(file = "./data/suppltable-2.csv")

# supp table 3
sup3 <- read_csv("./data/suppltable-3.csv")


behaviors <- sup1 %>%
  select(TotalPath:ShockPerEntrance)  


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  output$distPlot <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2] 
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')
    
  })
  
  
  output$behaviorplot <- renderPlot({
  
    sup1 %>%
      dplyr::group_by(treatment, trial, trialNum) %>%
      dplyr::summarise(value = mean(!!input$variable), 
                       se = sd(!!input$variable)/sqrt(length(!!input$variable)))  %>%
      ggplot(aes(x= trialNum, y=value, color=treatment)) + 
      geom_errorbar(aes(ymin=value-se, ymax=value+se, color=treatment), width=.1) +
      geom_line() +
      geom_point(size = 1.5) +
      labs(y = input$variable) +
      scale_x_continuous(name= "trial", 
                         breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                         labels = c( "P", "T1", "T2", "T3",
                                     "Rt", "T4", "T5", "T6", "Rn")) +
      scale_alpha_continuous( breaks = c(1, 2, 3)) +
      theme_minimal() +
      scale_color_manual(values = allcolors,
                         name  = NULL,
                         labels = levelstreatmentlegend)  +
      theme(legend.position = "bottom") 
  })
  
 
  
})
