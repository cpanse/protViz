#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(protViz)
library(parallel)

summary.PTM_MarkerFinder <- function(data, itol_ppm = 10, 
                                     ppm.error=20, 
                                     mZmarkerIons=sort(c(428.0367, 348.0704, 250.0935, 136.0618)),
                                     minNumberIons = 2, 
                                     minMarkerIntensityRatio = 10){
  S <- mclapply(data, function(x){ 
    idx <- findNN(mZmarkerIons, x$mZ); 
    
    ppm.error <- 1e-06 * itol_ppm * x$mZ[idx]
    
    b <- (abs(mZmarkerIons - (x$mZ[idx])) < ppm.error)
    
    sum.mZmarkerIons.intensity <- sum(x$intensity[idx[b]])
    
    sum.intensity <- sum(x$intensity)
    
    percent.mZmarkerIons <- round(100 * sum.mZmarkerIons.intensity / sum.intensity, 1)
    
    if (sum.mZmarkerIons.intensity > 0 
        & sum(b) >= minNumberIons 
        & percent.mZmarkerIons > minMarkerIntensityRatio){
      
      data.frame( query=x$id, 
                  percent.mZmarkerIons=percent.mZmarkerIons, 
                  sum.intensity=sum.intensity,
                  markerIonIntensity=x$intensity[idx[b]], 
                  markerIonMZ=mZmarkerIons[b], 
                  peptideSequence=x$peptideSequence,
                  scans=x$scans,
                  markerIonPpmError =ppm.error[b],
                  mZ=x$mZ[idx[b]],
                  pepmass=x$pepmass,
                  score=x$mascotScore
      )
      
    }else{
      
      NULL
    }
  }
  
  )
  
  do.call('rbind', S)  
}


shiny_PTM_MarkerFinder  <- function(DATAROOT="/Users/cp/__projects/2016/20160525-PTM-markerFinder-shiny/data/"){
  
  #load("/Users/cp/__projects/2016/20160525-PTM-markerFinder-shiny/data/F230466.RData")
  # Define  for application that draws a histogram
  
  ui <- shinyUI(fluidPage(
     
     # Application title
     titlePanel("PTM Marker Finder"),
     
     # Sidebar with a slider input for number of bins 
     sidebarLayout(
        sidebarPanel(
           selectInput('file', 'file', list.files(path=DATAROOT)),
           sliderInput("minMarkerIntensityRatio",
                       "minMarkerIntensityRatio",
                       min = 1,
                       max = 100,
                       value = 10),
           
           sliderInput("minNumberIons",
                       "minNumberIons",
                       min = 1,
                       max = 5,
                       value = 2)
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("PTM_MarkerFinder")
        )
     )
  ))
  
  # Define server logic required to draw a histogram
  server <- shinyServer(function(input, output) {
     
     output$PTM_MarkerFinder <- renderPlot({
       
        #input <- list(minMarkerIntensityRatio=40,minNumberIons=2)
       load(paste(DATAROOT, input$file,sep='/'))
            
        S <- summary.PTM_MarkerFinder(get(unlist(strsplit(input$file, split="[.]"))[1]), 
                                      minMarkerIntensityRatio=input$minMarkerIntensityRatio,
                                      minNumberIons=input$minNumberIons)
        #x    <- faithful[, 2] 
        #bins <- seq(min(x), max(x), length.out = input$bins + 1)
        op<-par(mfrow=c(1,1))
        if (!is.null(S)){
        boxplot(markerIonIntensity ~ markerIonMZ,
                data=S,
                log = "y",
                main = input$file,
                xlab = "markerIon m/z",
                ylab = "log10 based marker ion intensity",
                sub=paste(input$minNumberIons, input$minMarkerIntensityRatio, nrow(S)))
        }else{
          plot(0,0,
               type='n',
               sub=paste(input$minNumberIons, input$minMarkerIntensityRatio, nrow(S)))
        text(0,0, "no data", cex=3)
          }
        # draw the histogram with the specified number of bins
        #hist(x, breaks = bins, col = 'darkgray', border = 'white')
     })
  })

  # Run the application 
  shinyApp(ui = ui, server = server)
}
