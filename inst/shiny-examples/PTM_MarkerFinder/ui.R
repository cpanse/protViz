#
# This is a Shiny web application. 
#
#    https://cran.r-project.org/package=protViz
#
#
# $HeadURL: svn+ssh://cp@fgcz-148.uzh.ch/home/cp/__SUBVERSION_REPOSITORY__/__projects/2016/20160704_pptm_shiny/ui.R $
# $Id: ui.R 915 2017-04-11 12:36:53Z cp $
# $Date: 2017-04-11 14:36:53 +0200 (Tue, 11 Apr 2017) $

 
library(bfabricShiny)

# source("ptm_marker_finder.R")

  
shinyUI(fluidPage(
     # Application title
     titlePanel(paste("PTM Marker Finder -- https://CRAN.R-project.org/package=protViz Version:", packageVersion('protViz'))),
     
     
     sidebarLayout(
        sidebarPanel(
     #bfabricInput("bfabric8"),
     
     hr(),
	   htmlOutput("mZmarkerIons"),
	   hr(),
	   p('publication:'),
	   p('Nanni, P., Panse, C., Gehrig, P., Mueller, S., Grossmann, J. and Schlapbach, R. (2013), PTM MarkerFinder, a software tool to detect and validate spectra from peptides carrying post-translational modifications. Proteomics, 13: 2251–2255. ', a('DOI: 10.1002/pmic.201300036', href='http://onlinelibrary.wiley.com/doi/10.1002/pmic.201300036/abstract;jsessionid=717FB314BBA9A6BD2E722BD257D3D2A9.f01t04'))
        ),
        
        # Show a plot of the generated distribution
          mainPanel(
            tabsetPanel(
              tabPanel("bfabric", {
                list(bfabricInput("bfabric8"), 
                htmlOutput("load"))
                }),
              tabPanel("boxplot", plotOutput("PTM_MarkerFinder"))) 
	  #p('please wait some seconds until the data are processed...'),
     #      plotOutput("PTM_MarkerFinder", height = "700px")
        )
     )
  ))

