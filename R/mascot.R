#R
# Christian Panse <cp@fgcz.ethz.ch> 20170525

.mascot.get.pep_score <- function(data){
  as.vector(unlist(lapply(data$queries, function(x){rv <- x$q_peptide$pep_score; if(is.null(rv)){NA}else{as.numeric(rv)}})))
}

.mascot.get.pep_expect <- function(data){
  as.vector(unlist(lapply(data$queries, function(x){rv <- x$q_peptide$pep_expect; if(is.null(rv)){NA}else{as.numeric(rv)}})))
}

plot.mascot <- function(data){
  
  pep_score <- .mascot.get.pep_score(data)
  pep_expect <- .mascot.get.pep_expect(data)
  
  plot(pep_score, 1/log(pep_expect,10),log='x', pch=16, col=rgb(0.1, 0.1, 0.1, alpha = 0.2))
  
  lcmsmap(data)
}


plot.mascot_query <- function(query){
 
  S <- lapply(strsplit(query$StringIons1, ","), function(x){strsplit(x, ':')})[[1]]
  
  mZ <- as.numeric(sapply(S, function(x){x[1]}))
  intensity <- as.numeric(sapply(S, function(x){x[2]}))
  idx <- order(mZ)
  
  spec <- list(mZ=mZ[idx], intensity=intensity[idx])
  
  peakplot(peptideSequence = query$q_peptide$pep_seq, 
           spec = spec, sub='')
  # query$q_peptide$pep_seq
}
