#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#


# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("What genes regulate learning and memory in the mouse hippocampus?"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "gene1",
        label = "Choose a gene.",
        choices = c(gene_names),
        selected = c("Grin1"),
        multiple = FALSE
      ),
      selectInput(
        inputId = "gene2",
        label = "Choose another gene.",
        choices = c(gene_names),
        selected = c("Fos"),
        multiple = FALSE
      )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      h3("Correlation matrix of candidate memory genes, immediate early genes, and behavior"),
      plotOutput("cormat"),
      h3("Correlation between two genes"),
      plotOutput("scatterplot"),
      verbatimTextOutput("corene1gene2"),
      h3("Relationship of those two genes to behavioral PC1"),
     plotOutput("pc1")
       
    )
  )
))
