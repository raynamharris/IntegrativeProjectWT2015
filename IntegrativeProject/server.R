#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#



# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  
  output$scatterplot <- renderPlot({
    
    mytitle <- paste("Correlation between ", input$gene1,  " and ",
                     input$gene2, sep = "")
    
      ggplot(df, aes_string(x = input$gene1 , 
                 y = input$gene2)) +
      geom_point(aes(color = training, shape = subfield),
                 size = 3) +
      theme_minimal(base_size = 16) +
      theme(axis.title.x = element_text(face = "italic"),
            axis.title.y = element_text(face = "italic")) +
        labs(subtitle = mytitle) +
        geom_smooth(method = "lm", color = "darkgrey")
    
  })
  
  output$pc1 <- renderPlot({
    
    
    mytitle <- paste("Correlation between ", input$gene1,  " and ",
                     input$gene2, " and behavioral PC1 ", sep = "")
    
    p1 <- ggplot(df, aes_string(y = input$gene1)) +
      geom_point(aes(x = PC1, color = training, shape = subfield),
                 size = 3) +
      theme_minimal(base_size = 16) +
      theme(axis.title.y = element_text(face = "italic"),
            legend.position = "none") +
      labs(subtitle = mytitle) +
      geom_smooth(aes(x = PC1), method = "lm", color = "darkgrey")
    
    p2 <- ggplot(df, aes_string(y = input$gene2)) +
      geom_point(aes(x = PC1, color = training, shape = subfield),
                 size = 3) +
      theme_minimal(base_size = 16) +
      theme(axis.title.y = element_text(face = "italic")) +
      labs(subtitle = " ") +
      geom_smooth(aes(x = PC1), method = "lm", color = "darkgrey")
    
    plot_grid(p1,p2, rel_widths = c(0.4,0.6))
    
    
    })
      
  output$corene1gene2 <- renderPrint({
    x <- df[input$gene1] %>% pull(input$gene1)
    y <- df[input$gene2] %>% pull(input$gene2)
    
    cor.test(x, y, method="pearson")
  })
  

  output$cormat <- renderPlot({
    
    genesonly <- df[6:20] %>%
      drop_na()
    cormat <- cor(genesonly)
    
    col3 <- colorRampPalette(c("blue", "white", "red"))
    
   corrplot(cormat, type = "upper", order = "hclust", 
             tl.col = "black", tl.srt = 45, 
             col = col3(50), diag = F)

  })
  
})
