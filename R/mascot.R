#R
# Christian Panse <cp@fgcz.ethz.ch> 20170525

.mascot.get.pep_seq <- function(obj){
  rv <- sapply(obj$queries, function(x) { 
    if ("q_peptide" %in% names(x)){
      if ("pep_seq" %in% names(x$q_peptide)){
        return (x$q_peptide$pep_seq)
      }
    }
    NA
  })
  as.vector(rv)
}

.mascot.get.rt <- function(obj){
  as.numeric(unlist(lapply(obj$queries, function(x){x$RTINSECONDS})))
}

.mascot.get.query <- function(query){
  L <-  .mascot.get.ms2(query)
  
  rv <- list(MonoisotopicAAmass = NA,
             charge = as.numeric(gsub("[+]", "", query$query_charge, perl=TRUE)),
             id = 0,
             intensity = L$intensity,
             mZ = L$mZ, 
             mascotScore = NA,
             modification = NA,
             pepmass = as.numeric(query$query_moverz),
             peptideSequence = NA,
             proteinInformation = NA,
             rtinseconds = as.numeric(query$RTINSECONDS),
             scans = query$SCANS, 
             searchEngine = "mascot",
             title = '',
             varModification = NA)

  # TODO: add proteinInformation
  if ('q_peptide' %in% names(query)){
    rv$title <- query$q_peptide$pep_scan_title
    rv$mascotScore <- as.numeric(query$q_peptide$pep_score)
    
    if ( is.null(query$q_peptide$pep_var_mod)){
      rv$modification <- NA
    }else{
      rv$modification <- query$q_peptide$pep_var_mod
    }
    
    rv$peptideSequence <- as.character(query$q_peptide$pep_seq)
    
    # TODO: $q_peptide$pep_var_mod_pos
    rv$varModification <- rep(0.0, nchar(query$q_peptide$pep_seq))
  }

  class(rv) <- c('psm', 'list')
  rv
}


# just define a generic S3 method
as.psmSet <- function(object, ...){
  UseMethod("as.psmSet")
}

#' transformas a mascot object into a psmSet
#'
#' @param mascot obj 
#'
#' @return a psmSet object
as.psmSet.mascot <- function(object, ...){
  
  if (is.mascot(object)){
    # todo(cp): class('psmSet')
    rv <- NULL
    #if(require(parallel)){
    #  rv <- mclapply(obj$queries, .mascot.get.query, ...)
    #} else{
    rv <- lapply(object$queries, .mascot.get.query)
    #}
    
    # assign the ``query number''
    for (idx in 1:(length(rv)-1))
      rv[[idx]]$id <- idx
    
    class(rv) <- c("psmSet", "list")
    return(rv)
  }
  
  NULL
}

.mascot.get.ms2 <- function(query){
  
  
  if (is.character(query$StringIons1)){
    S <- lapply(strsplit(query$StringIons1, ","), function(x){strsplit(x, ':')})[[1]]
  
    mZ <- as.numeric(sapply(S, function(x){x[1]}))
    intensity <- as.numeric(sapply(S, function(x){x[2]}))
    idx <- order(mZ)
    
    return(list(mZ=mZ[idx], intensity=intensity[idx]))
  }
  
  return(list(mZ=NULL, intensity=NULL))
}

.mascot.get.pep_score <- function(obj){
  as.vector(unlist(lapply(obj$queries, function(x){rv <- x$q_peptide$pep_score; if(is.null(rv)){NA}else{as.numeric(rv)}})))
}

.mascot.get.pep_expect <- function(obj){
  as.vector(unlist(lapply(obj$queries, function(x){rv <- x$q_peptide$pep_expect; if(is.null(rv)){NA}else{as.numeric(rv)}})))
}

summary.mascot <- function(object, ...){
  if (is.mascot(object)){
    cat("number of queries:\n")
    cat(paste("\t", length(object$queries), "\n"))
    cat("quantile pep_score:\n")
    quantile(.mascot.get.pep_score(object), na.rm=TRUE)
  }
  #NextMethod('summary')
}

is.mascot <- function(obj){
  if ('mascot' %in% class(obj)){
    if (is.list(obj$queries) & sum(sapply(obj$queries, is.mascot_query)) == length(obj$queries)){
      return(TRUE)
    }
  }
  
  return(FALSE)
}

is.mascot_query <- function(obj){
  mascor_query.names <- c("query_charge", "query_moverz", "SCANS", "StringIons1",
                          "RTINSECONDS", "StringTitle")
  
  if ('mascot_query' %in% class(obj) 
      & length(mascor_query.names) == sum(mascor_query.names %in% names(obj))){
    return(TRUE)
  }
  
  return(FALSE)
}

plot.mascot <- function(x, ...){
  if (is.mascot(x)){
    pep_score <- .mascot.get.pep_score(x)
    pep_expect <- .mascot.get.pep_expect(x)
    
    # peptide scores versus e-value
    plot(pep_score, 1 / log(pep_expect,10),log='x', pch=16, col=rgb(0.1, 0.1, 0.1, alpha = 0.2))
    
    # SSRC
    rt.ssrc.predicted <- as.vector(sapply(.mascot.get.pep_seq(x), function(p){if(is.na(p)){NA}else{ssrc(p)}}))
    rtinseconds <- .mascot.get.rt(x)
    plot(rt.ssrc.predicted ~ rtinseconds, pch = 16, col=rgb(0.1, 0.1, 0.1, alpha = 0.2))
    
    # LC-MS map
    lcmsmap(x, ...)
  }
}


plot.mascot_query <- function(x, ...){
 
  if (is.mascot_query(x)){
    spec <- .mascot.get.ms2(x)
    pep_seq <- ""
    
    if ("q_peptide" %in% names(x)){
      if ("pep_seq" %in% names(x$q_peptide)){
        pep_seq <- x$q_peptide$pep_seq
      }
    }
    peakplot(peptideSequence = pep_seq, 
             spec = spec, 
             sub=x$StringTitle)
  }

}
