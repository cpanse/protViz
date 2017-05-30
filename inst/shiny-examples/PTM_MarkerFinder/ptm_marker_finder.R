#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# $HeadURL: svn+ssh://cp@fgcz-148.uzh.ch/home/cp/__SUBVERSION_REPOSITORY__/__projects/2016/20160704_pptm_shiny/ptm_marker_finder.R $
# $Id: ptm_marker_finder.R 915 2017-04-11 12:36:53Z cp $
# $Date: 2017-04-11 14:36:53 +0200 (Tue, 11 Apr 2017) $


 processMgf <- function(input) {

       load(paste(input$DATAROOT, input$file, sep='/'))
       S <- get(strsplit(input$file, ".RData")[[1]])

       list(summary = summary.PTM_MarkerFinder(S, 
           minMarkerIntensityRatio=input$minMarkerIntensityRatio,
           minNumberIons=input$minNumberIons),
	   PsmSet=S)
 }


my.test <- function(){
	input <- list()

	input$DATAROOT <- "/scratch/cpanse/p1352/20160701/"
	
	input$file <- "F230201.blib" 
	input$minMarkerIntensityRatio <- 15
	input$minNumberIons <- 3

	processMgf(input)
}

