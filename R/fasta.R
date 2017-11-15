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


summary.FASTA <- function(object, revpattern = "^>REV.*", conpattern = "^>.*FGCZCont.*", ...){
  
  n <- sum(grepl("^>", object))
  nrev <- sum(grepl(revpattern, object))
  ncon <- sum(grepl(conpattern, object))
  nAA <- sum(nchar(object[grepl("^[WFLIMVYATPEDCSQGNRHK]+$", object)]))
  
  
  list("filename" = attr(object, "filename"),
       "number_of_IDs" = n,
       "number_of_REVs" = nrev,
       "number_of_Contaminats" = ncon,
       "number_of_AAs" = nAA,
       "number_of_tryptics_peptides" = number_of_tryptic_peptides_FASTA(object),
       "object_size" = object.size(object)
       )
}

# just define a generic S3 method
tryptic_digest <- function(object, ...){
  UseMethod("tryptic_digest")
}

fcat <- function(object, ...){
  UseMethod("fcat")
}

fcat.FASTA <- function(fasta) {
  fcat_FASTA(fasta)
}

tryptic_digest.FASTA <- function(fasta) {
  tryptic_digest_FASTA(fcat_FASTA(fasta))
}
