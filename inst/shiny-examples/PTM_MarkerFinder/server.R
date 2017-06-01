#
# This is a Shiny web application. 
#
#
#    https://cran.r-project.org/package=protViz
#
# $HeadURL: svn+ssh://cp@fgcz-148.uzh.ch/home/cp/__SUBVERSION_REPOSITORY__/__projects/2016/20160704_pptm_shiny/server.R $
# $Id: server.R 915 2017-04-11 12:36:53Z cp $
# $Date: 2017-04-11 14:36:53 +0200 (Tue, 11 Apr 2017) $
library(protViz)
library(parallel)
library(bfabricShiny)

source("./ptm_marker_finder.R")

.ssh_load_RData <- function(host = 'fgcz-r-021.uzh.ch', user = 'cpanse', file = NULL){
  e <- new.env()
  
  cmd <- paste('cat ',  file)
  
  ssh_cmd <- paste("ssh ", user, "@", host, " '", cmd, "'", sep="")
  message(ssh_cmd)
  
  S <- load(pipe(ssh_cmd))
  
  for (x in S){
    assign(x, get(x), e)
  }
  e
}

.load_RData <- function(file = NULL){
  e <- new.env()
  
  S <- load(file)
  
  for (x in S){
    assign(x, get(x), e)
  }
  e
}

shinyServer(function(input, output, session) {
    bf <- callModule(bfabric, "bfabric8",  applicationid = c(155))
    
    
      output$mZmarkerIons <- renderUI({
        
     		markerIons <- sort(c(428.0367, 348.0704, 250.0935, 136.0618, 524.057827,
     		                     542.068392, 560.078957, 559.094941, 584.090190))

     		e <- getRDataEnv()
     		
     		if (length(ls(e)) > 0){
     		tagList(
     		  
     		 hr(),
     		    selectInput('mascot_object', 'mascot_object:',  ls(e), multiple = FALSE),
     		 hr(),
     		  sliderInput("minMarkerIntensityRatio",
     		              "minMarkerIntensityRatio",
     		              min = 1,
     		              max = 100,
     		              value = 10),
     		  
     		  sliderInput("minNumberIons",
     		              "minNumberIons",
     		              min = 1,
     		              max = 6,
     		              value = 2),
     		  
     		  sliderInput("score_cutoff",
     		              "mascot score cut-off",
     		              min = 0,
     		              max = 100,
     		              value = 25),
		      selectInput('mZmarkerIons', 'marker ions:',  markerIons, multiple = TRUE,
		              selected = markerIons[1:5]),
		      checkboxGroupInput("plotLines", label = h3("plot lines"), choices = list("yes" = 1), selected = 1),
		      sliderInput("alpha", "alpha blending %", min=1, max=100, value=40),
     		 
		      downloadButton('downloadData', 'Download'),
		      downloadButton('downloadDataWide', 'Download (wide)'),
		      downloadButton('downloadMgf', 'Generate MGF(under construction!)')
     		)}
      })


    
       
     output$load <- renderUI({
            if(length(input$relativepath) > 0){
                actionButton("load", "load selected RData", icon("upload"))
            }
     })

  # return an env
	getRDataEnv <- eventReactive(input$load, {
	  message("eventReactive(input$load")
	  message(input$relativepath)
	  
	  filename <- file.path('/srv/www/htdocs/', input$relativepath)
	
	  if (file.exists(filename)){
	    .load_RData(file=filename)
	  }else{
	   .ssh_load_RData(file = filename, host = 'fgcz-r-021.uzh.ch')
	  }
		})
	
	getData <- eventReactive(input$mascot_object, {
	  
	  message(input$mascot_object)
	  as.psmSet(get(input$mascot_object, getRDataEnv()))
	})

 processedData <- reactive({
       S <- getData()
       
       mZmarkerIons <- sapply(input$mZmarkerIons, as.numeric)
       
       return(findMz(S, 
       				itol_ppm = 10,
       				mZmarkerIons=mZmarkerIons,
                                minMarkerIntensityRatio=input$minMarkerIntensityRatio,
                                minNumberIons=input$minNumberIons
				))
  })


     output$PTM_MarkerFinder <- renderPlot({
        S <- processedData()
        op <- par(mfrow=c(1,1))
        if (!is.null(S)){
        b <- boxplot(markerIonIntensity ~ markerIonMZ,
                data=S,
                log = "y",
                main = input$file,
		xlab = "markerIon m/z",
                ylab = "log10 based marker ion intensity")
	text(1:length(b$n), b$stats[3,], b$n, col="darkblue", pos=1)
	legend('topright', legend=sapply(input$mZmarkerIons, as.numeric))


	if (1 %in% input$plotLines){
		lines(as.factor(S$markerIonMZ), S$markerIonIntensity, 
			col=rgb(0.1,0.1,0.1, alpha=input$alpha/100),
			lwd=6)
	}

	legend("topleft", c(paste("minNumberIons =", input$minNumberIons),
		paste("minMarkerIntensityRatio =", input$minMarkerIntensityRatio),
		paste("Number of queries =",  length(unique(S$query))),
		paste("Number of psm (score >=", input$score_cutoff, ") =", 
			sum(S$score >= input$score_cutoff, na.rm=TRUE)),
		paste("Number of psm (score <", input$score_cutoff, ") =", 
			sum(S$score < input$score_cutoff, na.rm=TRUE))
		))
        }else{
          plot(0,0,
               type='n',
               sub=paste(input$minNumberIons, input$minMarkerIntensityRatio, nrow(S)))
        text(0,0, "no data", cex=3)
          }
        # draw the histogram with the specified number of bins
        #hist(x, breaks = bins, col = 'darkgray', border = 'white')
     })
       output$downloadData <- downloadHandler(
           filename = function() { paste(unlist(strsplit(input$file, 
                                                         split="[.]"))[1], "csv", sep=".")  },
	   content = function(file) {
	   	write.csv(processedData(), file, row.names = FALSE)
	   }
	   )

       output$downloadDataWide <- downloadHandler(
           filename = function() { paste(unlist(strsplit(input$file, 
                                                         split="[.]"))[1], "wide", "csv", sep=".")  },
	   content = function(file) {
	   	S <- processedData()
		SS <- reshape(S[,c(1,4,5,6)], direction='wide', timevar='markerIonMZ', idvar=c('query', 'peptideSequence'))
	   	write.csv(SS, file, row.names = FALSE)
	   }

	   )
       output$downloadMgf <- downloadHandler(
           filename = function() { paste(unlist(strsplit(input$file, 
                                                         split="[.]"))[1], "mgf", sep=".")  },
	   content = function(file) {
		res <- processMgf(input)
		protViz:::.PTM_MarkerFinder_writeMGF_Header(file)
       		dump <- lapply(unique(res$summary$query), 
       		function(idx){
			pattern <- res$summary[res$summary$query==idx, c(5, 4)]
			psm <- res$PsmSet[[idx]]
			psm$title <- psm$fileName
			psm$rtinseconds <- psm$rt
			psm$scans <- psm$id
       			protViz:::.PTM_MarkerFinder_writeMGF(psm, file, pattern = pattern)
          	})
	   	#write.csv(processData(input), file, row.names = FALSE)
	   }
	)
  })
