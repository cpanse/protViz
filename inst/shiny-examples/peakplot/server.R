#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(bibliospec)
library(specL)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  
  
  #dataInput <- reactive({
    
  # Bibliospec(dbfile=paste(input$DATAROOT, input$file, sep='/'))
    
  #})
  
  
  output$distPlot <- renderPlot({
    #S <- NULL
    
    if (input$load){
    BS <- Bibliospec(dbfile=paste(input$DATAROOT, input$file, sep='/'))
    S <- BS$getPsmSet()
    
    idx <- which(sapply(S, function(x){-10 * log((1E-6 + x$score)) / log(10) > input$mascotScoreCutOff}))
    }
  
    
   SS <- do.call('rbind', lapply(S[idx], function(x){
                    data.frame(rt.predicted=ssrc(as.character(x$peptideSequence)), 
                               rt=x$rt)
                    }))
   
    op <- par(mfrow=c(4, 1))
    
    hist(sapply(S, function(x){-10 * log((1E-6 + x$score)) / log(10)}), main='score')
    plot(S[[input$id]])
    plot(SS$rt , SS$rt.predicted, main=input$file, 
         sub=paste('number of psm =', length(SS$rt))
         )
    S.lm <- lm(SS$rt.predicted ~ SS$rt)
    abline(S.lm, col='cornflowerblue', lwd=4)
    plot(S)
  })
  
})
