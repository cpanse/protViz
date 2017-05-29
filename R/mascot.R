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
  
  rv <- list(mZ = L$mZ, 
             intensity = L$intensity,
             modification=NA,
             mascotScore = NA,
             peptideSequence = NA,
             proteinInformation = NA,
             id = NA,
             varModification = NA,
             pepmass = as.numeric(query$query_moverz),
             charge = as.numeric(gsub("[+]", "", query$query_charge, perl=TRUE)),
             scans = query$SCANS, 
             rtinseconds = as.numeric(query$RTINSECONDS))
             

  if ('q_peptide' %in% names(query)){
    rv$title <- query$q_peptide$pep_scan_title
    rv$mascotScore <- as.numeric(query$q_peptide$pep_score)
    rv$modification <- query$q_peptide$pep_var_mod_pos
    rv$peptideSequence <- as.character(query$q_peptide$pep_seq)
    rv$varModification <- rep(0.0, nchar(query$q_peptide$pep_seq))
  }

  class(rv) <- c('psm', 'list')
  rv
}

#' transformas a mascot object into a psmSet
#'
#' @param mascot obj 
#'
#' @return a psmSet object

.mascot.get <- function(obj){
  if (is.mascot(obj)){
    # todo(cp): class('psmSet')
    rv <- lapply(obj$queries, .mascot.get.query)
    
    # assign the ``query number''
    for (id in 1:length(rv)){
        rv[[id]]$id <- id
    }
    class(rv) <- c("psmSet", "list")
    return(rv)
  }
  
  NULL
}

.mascot.get.ms2 <- function(query){
  S <- lapply(strsplit(query$StringIons1, ","), function(x){strsplit(x, ':')})[[1]]
  
  mZ <- as.numeric(sapply(S, function(x){x[1]}))
  intensity <- as.numeric(sapply(S, function(x){x[2]}))
  idx <- order(mZ)
  
  list(mZ=mZ[idx], intensity=intensity[idx])
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
  if ('mascot' %in% class(obj) 
      & is.list(obj$queries) 
      & sum(sapply(obj$queries, is.mascot_query)) == length(obj$queries)){
    return(TRUE)
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
    plot(rt.ssrc.predicted ~ rtinseconds, pch=16, col=rgb(0.1, 0.1, 0.1, alpha = 0.2))
    
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
