#R

# Christian Panse <cp@fgcz.ethz.ch> 2021-01-28
# http://fgcz-ms-shiny.uzh.ch:8080/queue_generator10/
# see also
# https://doi.org/10.1021/pr8010099
# https://doi.org/10.1021/acs.jproteome.0c00536


.isBlockRandomFeasibible <- function(S, x){
	if (nrow(S) < 1) return (FALSE)
	min(table(S[, names(S) == x])) == max(table(S[, names(S) == x]))
}


#' Derive a randomization of a table.
#'
#' @param S a \code{data.frame}
#' @param x the column name of the group used for the block randomization.
#' @param check check if the size of each group is equal.
#'
#' @return returns a randomization order.
#' 
#' @export
#'
#' @examples
#' 
#' set.seed(1)
#' iris[c(1:2, 51:52, 101:103), ] |>
#'   blockRandom(x = "Species", check=FALSE)
blockRandom <- function(S, x = NA, check = TRUE){
  
  # if no group is provided just randomize
  if (is.na(x)){
    warning("no group provided switch to random.")
    S$dummy <- rep("group1", nrow(S))
    x <- 'dummy'
  }else
    if (check){
      # sanity check - all groups have same number of treatment
      stopifnot(.isBlockRandomFeasibible(S, x))
    }

	  
  if(isFALSE(.isBlockRandomFeasibible(S, x))){
    warning("Unequal cardinality of blocks.")
  }

	n <- max(table(S[, names(S) == x]))

	# split into groups
	rv <- lapply(unique(S[, x]),
		function(l, Y){treatment <- S[S[, x] == l,]; treatment[sample(nrow(treatment)),]}, Y=S)

	# compose block random data.frame
	m <- length(rv)
	base::Reduce(rbind, lapply(1:n,
	    function(i){base::Reduce(rbind, lapply(1:m, 
	        function(j){
	            if (i <= nrow(rv[[j]]))
	                rv[[j]][i, ]
	            else{
	                z <- rv[[j]][1, ]
	                z[1, ] <- NA
	                z[1, ]
	            }
	            }))[sample(m),]}))
}


#' Assign an instrument queue configuration to a plate
#' 
#' @description The function implements a space-filling curve mapping 1D to 2D.
#' This function aims to assign a sequence of samples to an instrument plate,
#' e.g., 48 well plate 85.4x127.5mm.
#'
#' @param S input data frame
#' @param x a vector of possible x-coordinates of the plate
#' @param y a vector of possible y-coordinates of the plate
#' @param plate a vector of plates
#' @param volume injection volume
#' @param reserve block plate positions
#'
#' @return a \code{data.frame}
#' @export
#'
#' @examples
#' iris[c(1:15, 51:65, 101:115), ] |>
#'   blockRandom(x = "Species", check=FALSE) |>
#'   assignPlatePosition()
assignPlatePosition <- function(S,
          x = as.character(1:8),
          y = c('A', 'B', 'C', 'D', 'E', 'F'),
          plate=1:4,
          volume=1,
          reserve = 46:48){
  
    stopifnot(is.data.frame(S))
    n <-  nrow(S)
    
    platePosition <- length(x) * length(y) 
    
    if (platePosition * length(plate) < nrow(S)){
        stop("More samples than available plate positions!")
    }
    
    S$run <- 1:n
    S$x <- x[((S$run -1) %% length(x)) + 1]
    S$y <- y[(floor((S$run - 1) / length(x)) %% length(y))  + 1]
    S$plate <- plate[floor((S$run - 1 ) / platePosition) + 1]
    
    if ("volume" %in% names(S))
      S$volume[is.na(S$volume)] <- volume
    else
      S$volume <- volume
    S
}



.insertStandardsLoop <- function(input, howoften = 1, howmany = 1, begin = FALSE, end = FALSE, between=TRUE)
{ 
  output <- data.frame()
  
  if (isTRUE(between)){
    # yes - for readability of the code we have a foor loop!
    for (i in 1:nrow(input)){
      output <- rbind(output, input[i, ])
      if (howoften > 0 && i %% howoften == 0 && howmany > 0){
        for (j in seq(1, howmany)){
          tmp <- rep(NA, ncol(input))
          output <- rbind(output, tmp)
        }
      }
    }
  }else{
    # no inserts inbeween runs
    output <- input
  }
  
  if (begin){
    tmp <- rep(NA, ncol(input))
    output <- rbind(tmp, output)
  }
  
  if (end){
    tmp <- rep(NA, ncol(input))
    output <- rbind(output, tmp)
  }
  
  output
}


#' Insert sample on a given position
#'
#' @param S input \code{data.frame}
#' @param stdName name of the sample
#' @param stdPosX x location on the plate
#' @param stdPosY y location on the plate
#' @param plate number of the plate
#' @param volume injection volume
#' @param method a path to the method file (optional)
#' @param ... addition parameter, e.g.,  \code{howoften = 1} \code{howmany = 1}.
#'
#' @return  \code{data.frame}
#' @export
#'
#' @examples
#' iris[c(1:15,51:65,101:115), ] |>
#'   assignPlatePosition() |>
#'   insertSamples(howoften=4, begin=FALSE, end=FALSE,
#'     stdPosX='6', stdPosY='F', plate=1, stdName = "clean")
#'     
insertSamples <- function(S, stdName = "autoQC01", stdPosX='8', stdPosY='F', plate=1, volume=2, method = "", ...){
  
  input <- S
  
  if (isFALSE('type' %in% names(input)))
    input$type <- "sample"
  if (isFALSE('method' %in% names(input)))
    input$method <- ""
  
  output <- .insertStandardsLoop(input, ...)
  
  output$type[is.na(output$type )] <- stdName
  output$x[is.na(output$x )] <- stdPosX
  output$y[is.na(output$y )] <- stdPosY
  output$volume[output$type == stdName] <- volume
  output$plate[output$type == stdName] <- plate
  output$name[output$type == stdName] <- stdName
  output$method[output$type == stdName] <- method
  
  output
}


formatXCalibur <- function(x, path=""){
  stopifnot(all(c('id', 'name', 'plate', 'x', 'y', 'volume') %in% names(x)))
  n <- nrow(x)
  fileName <- sprintf("S%s_%s", x$id, x$name)
  fileName <- gsub("^SNA_", "", fileName)
  
  regularSample <- grepl("^[0-9]+$", as.character(x$id))
  
  
  df <- data.frame("File Name"= sprintf("%s_%03d_%s",
                                        format(Sys.time(), "%Y%m%d"),
                                        1:n, fileName),
                   Path = rep(path, n), 
                   Position = sprintf("\"%d:%s,%s\"", x$plate, x$y, x$x),
                   "Inj Vol" = x$volume, 
                   "L3 Laboratory" = rep("FGCZ", n),
                   "Sample ID" = rep("", n),
                   "Sample Name" = x$name,
                   "L1 Study" = rep("", n),
                   "Instrument Method" = x$method)
  
  df[regularSample, 'Sample.ID'] <- x$id[regularSample]
  
  # remove dots in column names
  names(df) <- gsub("\\.", " ", names(df))
  
  df
}

writeXCalibur <- function(x, file, path){
  x <- formatXCalibur(x, path)
  fn <- file
  message(sprintf("writing XCalibur configuration to file %s ...", file))
  cat("Bracket Type=4\r\n", file = fn, append = FALSE)
  write.table(x, file = fn,
              sep = ',', row.names = FALSE,
              append = TRUE, quote = FALSE, eol = '\r\n')
}

