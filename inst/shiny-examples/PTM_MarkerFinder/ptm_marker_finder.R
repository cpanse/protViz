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



summary.PTM_MarkerFinder <- function(data, itol_ppm = 10, 
                                     mZmarkerIons=sort(c(428.0367, 348.0704, 250.0935, 136.0618, 524.057827, 542.068392, 560.078957, 559.094941, 584.090190)),
                                     minNumberIons = 2, 
                                     minMarkerIntensityRatio = 10){
  S <- mclapply(data, function(x){ 
    idx <- findNN(mZmarkerIons, x$mZ); 
    
    ppm.error <- 1e-06 * itol_ppm * x$mZ[idx]
    
    b <- (abs(mZmarkerIons - (x$mZ[idx])) < ppm.error)
    
    sum.mZmarkerIons.intensity <- sum(x$intensity[idx[b]])
    
    sum.intensity <- sum(x$intensity)
    
    (percent.mZmarkerIons <- round(100 * sum.mZmarkerIons.intensity / sum.intensity, 1))
    
    if (sum.mZmarkerIons.intensity > 0 
        & sum(b) >= minNumberIons 
        & percent.mZmarkerIons > minMarkerIntensityRatio){
      
      data.frame( query=x$id, 
                  percent.mZmarkerIons=percent.mZmarkerIons, 
                  sum.intensity=sum.intensity,
                  markerIonIntensity=x$intensity[idx[b]], 
                  markerIonMZ=mZmarkerIons[b], 
                  peptideSequence=x$peptideSequence,
                  #scans=x$scans,
                  markerIonPpmError = ppm.error[b],
                  mZ = x$mZ[idx[b]],
                  pepmass=x$pepmass,
                  modification = x$modification,
                  score=x$mascotScore
      )
      
    }else{
      NULL
    }
  }
  )
  do.call('rbind', S)  
}

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

