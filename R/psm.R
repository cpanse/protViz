#R
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/psm.R $
# $Id: psm.R 6222 2014-03-13 14:22:34Z cpanse $
# $Date: 2014-03-13 15:22:34 +0100 (Thu, 13 Mar 2014) $



# TODO 
# compute score by sum of error div. by number of hits

psm <- function(sequence, spec, FUN = defaultIon,
    plot = TRUE, 
    fi = fragmentIon(sequence, FUN=FUN)[[1]],
    fragmentIonError = 0.6) { 

    n <- nchar(sequence)

    pim <- fi$y[nrow(fi)]

    # consider only b and y ions
    # by.mZ<-c(fi$b, fi$y)
    # by.label<-c(paste("b",1:n,sep=''), paste("y",1:n,sep=''))

    by.mZ <- numeric()
    by.label <- character()
    fi.names <- names(fi)

    for (i in 1:ncol(fi)){
        by.mZ <- c(by.mZ, fi[,i])
        by.label <- c(by.label, paste(fi.names[i],1:n,sep=''))
    }

    NN <- findNN_(q=by.mZ, vec=spec$mZ)


    mZ.error<-spec$mZ[NN] - by.mZ

    if (plot == TRUE){
	plot(mZ.error ~ spec$mZ[NN],
            pch=22,
	    ylim=c(-5 * fragmentIonError,  5 * fragmentIonError))
        abline(h=fragmentIonError,col='grey')
        abline(h=-fragmentIonError,col='grey')
        abline(h=0,col='grey',lwd=2)

        plot(mZ.error[mZ.error.idx<-order(mZ.error)],
            main=paste("Error of", sequence, "(parent ion mass =", round(pim,2) ,"Da)"),
            ylim=c(-5*fragmentIonError, 5*fragmentIonError),
            pch=22,
            sub=paste('The error cut-off is', 
                fragmentIonError, 'Da (grey line).')
            )

        abline(h=fragmentIonError,col='grey')
        abline(h=-fragmentIonError,col='grey')
        abline(h=0,col='grey',lwd=2)

        text(1:length(by.label), 
            mZ.error[mZ.error.idx],  
            by.label[mZ.error.idx],
            cex=0.75,pos=3) 

        hits=(abs(mZ.error) < fragmentIonError)
        nHits<-sum(hits)

        sumMZerror=round(sum(abs(mZ.error[hits])),2)

        avgMZerror=round(sumMZerror / nHits, 2)
        cover=round(nHits/(nrow(fi)*ncol(fi)),2)

        legend("topleft", paste(c('nHits','sumMZerror','avgMZerror','cover'),
            as.character(c(nHits, sumMZerror, avgMZerror, cover)),sep='=')) 

    }


    return (list(mZ.Da.error=mZ.error, 
        mZ.ppm.error=1E+6*mZ.error/by.mZ,
        idx=NN,
        label=by.label, 
        score=-1, 
        sequence=sequence,
        fragmentIon=fi))
}




summary.psmSet <- function (object, ...){
  cat("Summary of a \"psmSet\" object.")
  
  cat("\nNumber of precursor:\n\t")
  cat(length(object))
  
  cat("\nNumber of precursors in Filename(s)\n")
  t <- (table(unlist(lapply(object, function(x){x$fileName}))))
  
  n <- names(t)
  w <- getOption("width") / 2
  for (i in 1:length(t)){
    cat('\t')
    cat(substr(n[i], nchar(n[i])-w, nchar(n[i])))
    cat('\t')
    cat(t[i])
    cat('\n')
  }
  
  cat("Number of annotated precursor:\n\t")
  cat(sum(unlist(lapply(object, function(x){x$proteinInformation != ''}))))
  cat ("\n")
}


plot.psmSet <- function (x, iRTpeptides = iRTpeptides, ...){

  if (is.psmSet(x)){
    lcmsmap(x, ...)
  }
  
  #if (!is.null(iRTpeptides)){
  #  rt <- unlist(lapply(data, function(x){x$rt}))
  #  pepmass <- unlist(lapply(data, function(xx){xx$pepmass}))
  #  
  #  peptide <- unlist(lapply(data, function(xx){xx$peptideSequence}))
  #  idx.iRT <- which(peptide %in% iRTpeptides$peptide)
  #  
  #  points(rt[idx.iRT], pepmass[idx.iRT],
  #         pch = 16, cex = 2, 
  #         col = rgb(0.9,0.1,0.1, alpha = 0.4))
  #}
}


is.psm <- function(object){

  psm.names <- c("MonoisotopicAAmass",
                 "charge",
                 "id",
                 "intensity",
                 "mZ",
                 "mascotScore",
                 "modification",
                 "pepmass",
                 "peptideSequence",
                 "proteinInformation",
                 "rtinseconds",
                 "scans",
                 "searchEngine",
                 "title")
  
  object.names <- names(object)
  
  idx.missing <- which(!psm.names %in% names(object))

  if (length(idx.missing) > 0){
    msg <- paste("while checking psm '", psm.names[idx.missing], "' is missing.", sep=' ')
    message(msg)
    return (FALSE)
  }
  
  return(TRUE)
}

is.psmSet <- function(object){
  sum(sapply(object, is.psm)) == length(object)
}

plot.psm <- function (x, ...){
  
  if (is.na(x$peptideSequence)){
	  plot(x$mZ, x$intensity, sub='no assigned peptide sequence', type='h')
  }else{

  AAmass <- aa2mass(x$peptideSequence)[[1]]#, protViz::AA$Monoisotopic, protViz::AA$letter1)
  AAmodifiedMass <- AAmass + x$varModification
  fi <- fragmentIon(AAmodifiedMass)[[1]]
  spec <- list(mZ = x$mZ, intensity = x$intensity)
  return(peakplot(peptideSequence = x$peptideSequence, 
                  spec=spec, fi = fi, ...))
  }
}

  
# TODO(cp): rename to findMarkerIon
findMz <- function(data, mZmarkerIons, itol_ppm = 10, minNumberIons = 2, minMarkerIntensityRatio = 10){
  UseMethod("findMz")
}

findMz.mascot <- function(data, mZmarkerIons, itol_ppm = 10, minNumberIons = 2, minMarkerIntensityRatio = 10){
  findMz.psmSet(as.psmSet(data), mZmarkerIons, itol_ppm , minNumberIons , minMarkerIntensityRatio)
}

findMz.psmSet <- function(data, 
                          mZmarkerIons,
                          itol_ppm = 10, 
                          minNumberIons = 2, 
                          minMarkerIntensityRatio = 10){
  if(is.psmSet(data)){
    S <- lapply(data, function(x){ 
      
      idx <- findNN(mZmarkerIons, x$mZ) 
      
      ppm.error <- 1e-06 * itol_ppm * x$mZ[idx]
      
      b <- (abs(mZmarkerIons - (x$mZ[idx])) < ppm.error)
      
      sum.mZmarkerIons.intensity <- sum(x$intensity[idx[b]])
      
      sum.intensity <- sum(x$intensity)
      
      (percent.mZmarkerIons <- round(100 * sum.mZmarkerIons.intensity / sum.intensity, 1))
      
      if (sum.mZmarkerIons.intensity > 0 
          & sum(b) >= minNumberIons 
          & percent.mZmarkerIons > minMarkerIntensityRatio){
        
        data.frame(query = x$id, 
                   percent.mZmarkerIons = percent.mZmarkerIons, 
                   sum.intensity = sum.intensity,
                   markerIonIntensity = x$intensity[idx[b]], 
                   markerIonMZ = mZmarkerIons[b], 
                   peptideSequence = x$peptideSequence,
                   #scans=x$scans,
                   markerIonPpmError = ppm.error[b],
                   mZ = x$mZ[idx[b]],
                   pepmass = x$pepmass,
		   charge = x$charge,
		   rtinseconds = x$rtinseconds,
                   modification = as.character(paste(x$varModification, collapse = '')),
                   score = x$mascotScore
        )
      }else{
        return(NULL)
      }
    }
    )
    rv <- do.call('rbind', S)  
    row.names(rv) <- 1:nrow(rv)
    rv
  }else{NULL}
}


