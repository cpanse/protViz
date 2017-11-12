#R

is.FASTA <- function(x){
  if(sum(!grepl("^[WFLIMVYATPEDCSQGNRHK]+$|^>.+$", x)) == 0){
    return (TRUE)
  }
  FALSE
}

as.FASTA <- function(x){
  x <- as.character(x)
  x[grepl("^[WFLIMVYATPEDCSQGNRHK]+$|^>.+$", x)]
}

#' read FASTA
#'
#' @param filename 
#'
#' @return fasta
#' @export read.FASTA
#'
#' @examples
#' F <- read.FASTA("~/p1875_db10_20170817.fasta")
read.FASTA <- function(filename){
  FASTA <-scan(filename, what='character', sep = "\n")
  if(!is.FASTA(FASTA)){
    warning(paste(filename, "is not a FASTA. casting ..."))
    FASTA <- as.FASTA(FASTA)
  }
  attr(FASTA, "filename") <- filename
  class(FASTA) <- "FASTA"
  FASTA
}


summary.FASTA <- function(x, revpattern = "^>REV.*", conpattern = "^>.*FGCZCont.*"){
  n <- sum(grepl("^>", x))
  nrev <- sum(grepl(revpattern, x))
  nAA <- sum(nchar(x[grepl("^[WFLIMVYATPEDCSQGNRHK]+$", x)]))
  
  cat("filename:", attr(x, "filename"), sep="\t", "\n")
  cat("number of IDs:", n, sep="\t", "\n")
  cat("number of REVs:", nrev, sep="\t", "\n")
  cat("number of CONs:", nrev, sep="\t", "\n")
  cat("number of AAs:", nAA, sep="\t", "\n")
  
  
}
  

fcat.FASTA <- function(x){
 # .Call("fcatCpp", x, PACKAGE = "protViz")
}

tryptic_digest.FASTA <- function(x){
  
}
