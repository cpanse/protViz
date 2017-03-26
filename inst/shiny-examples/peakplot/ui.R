#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(jsonlite)

DATAROOT <- "/Users/cp/data/20160628/"
DATAROOT <- "/scratch/cpanse/p1352/20160701/"


#WORKUNITID <- 123456
#WORKUNITID <- WORKUNITID.g


# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("protViz"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      textInput("WORKUNITID", 
                "Workunit Id", 
                value=WORKUNITID),
      htmlOutput("WUControls"),
      uiOutput("selectID"),
      textInput("PROJECT", 
                "Project", 
                value="123"),
      textInput("STORAGEPATH",
                "dir",
                value="p1352/bfabric/Proteomics/PTM_MarkerFinder_protViz_ADP_Ribosyl/2014/2014-12/2014-12-09/workunit_130109/"),
       textInput("DATAROOT", 
                 "data directory", 
                 value=DATAROOT),                   
       selectInput('file', 'file', list.files(path=DATAROOT)),
       actionButton("load", "load"),
       sliderInput("id", 
                   "psm id", 
                   min = 1,
                   max = 150,
                   value = 1),
       
       sliderInput("mascotScoreCutOff",
                   "score cut off:",
                   min = 1,
                   max = 150,
                   value = 1)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
       plotOutput("distPlot", height = "1200")
    )
  )
))
