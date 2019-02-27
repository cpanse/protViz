#R

.is.cometdecoy <- function(object){
    if(sum(c('num', 'sp_score', 'protein', 'plain_peptide') %in% names(object))==4) return(TRUE)

    FALSE
}

summary.cometdecoy <- function(object, psmFdrCutoff=0.05, decoyPattern="^REV_", ...){
  stopifnot(.is.cometdecoy(object))

  # consider only first matches
  object <- object[object$num ==1, ]
  
  idobject <-order(object$sp_score)
  n <- nrow(object)
  object <- object[rev(idobject), ]
  
  object$REV <- grepl("^REV_", object$protein)
  if(sum(object$REV) == 0) warning("no decoy hit found. consider different decoy pattern or enable decoy search in comet.")
  
  object$nREVhits <- cumsum(object$REV)
  object$FDR <- object$nREVhits / (seq(1, n) - object$nREVhits)
  
  nConfidentPSM <- which(object$FDR > psmFdrCutoff)[1]
  nDecoyPSM <- object$nREVhits[which(object$FDR > psmFdrCutoff)[1]]
  confidentPeptide <- as.character(object$plain_peptide[seq(1, nConfidentPSM)])
  decoyPeptide <- grepl(decoyPattern, as.character(object$protein[seq(1, nConfidentPSM)]))
  
  nDecoyPeptide <- length(unique(confidentPeptide[decoyPeptide]))
  nConfidentPeptide <- length(unique(confidentPeptide[!decoyPeptide]))
  
  confidentProteins <- unique(sapply(strsplit(as.character(object$protein[seq(1, nConfidentPSM)]),','), function(y)y[1]))
  nDecoyProteins <- length(grep(decoyPattern, confidentProteins))
  nConfidentProteins <- length(confidentProteins)-nDecoyProteins
  
  df <- data.frame(nPSM=n,
                   psmFdrCutoff=psmFdrCutoff,
                   nDecoyPSM=nDecoyPSM,
                   nConfidentPSM=nConfidentPSM,
                   nDecoyPeptide=nDecoyPeptide,
                   nConfidentPeptide=nConfidentPeptide,
                   nDecoyProteins=nDecoyProteins,
                   nConfidentProteins=nConfidentProteins,
                   fdrPSM=nDecoyPSM/nConfidentPSM,
                   fdrPeptide=nDecoyPeptide/nConfidentPeptide,
                   fdrProtein=nDecoyProteins/nConfidentProteins
  )
  df     
}
