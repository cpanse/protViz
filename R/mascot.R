#R


plot.mascot <- function(data){
  
  pep_score <- as.numeric(unlist(lapply(data$queries, function(x){x$q_peptide$pep_score})))
  pep_expect <- as.numeric(unlist(lapply(data$queries, function(x){x$q_peptide$pep_expect})))
  plot(pep_score, 1/log(pep_expect,10),log='x', pch=16, col=rgb(0.1, 0.1, 0.1, alpha = 0.2))
  
  #lcmsmap(data)
}


plot.mascot_query <- function(query){
 
  S <- lapply(strsplit(query$StringIons1, ","), function(x){strsplit(x, ':')})[[1]]
  
  mZ <- as.numeric(sapply(S, function(x){x[1]}))
  intensity <- as.numeric(sapply(S, function(x){x[2]}))
  
  spec <- list(mZ=mZ, intensity=intensity)
  
  peakplot(peptideSequence = query$q_peptide$pep_seq, 
           spec = spec, sub='')
  # query$q_peptide$pep_seq
}
