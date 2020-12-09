#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Genes and Behavior"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      varSelectInput("variable", "Behavioral estimate:", behaviors)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      
      img(src='fig-1.png', align = "center", width = 600),
      plotOutput("behaviorplot")
    )
  )
))
