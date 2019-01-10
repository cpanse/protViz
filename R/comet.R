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
  
  nConvidentPSM <- which(object$FDR > psmFdrCutoff)[1]
  nDecoyPSM <- object$nREVhits[which(object$FDR > psmFdrCutoff)[1]]
  convidentPeptide <- as.character(object$plain_peptide[seq(1, nConvidentPSM)])
  decoyPeptide <- grepl(decoyPattern, as.character(object$protein[seq(1, nConvidentPSM)]))
  
  nDecoyPeptide <- length(unique(convidentPeptide[decoyPeptide]))
  nConvidentPeptide <- length(unique(convidentPeptide[!decoyPeptide]))
  
  convidentProteins <- unique(sapply(strsplit(as.character(object$protein[seq(1, nConvidentPSM)]),','), function(y)y[1]))
  nDecoyProteins <- length(grep(decoyPattern, convidentProteins))
  nConvidentProteins <- length(convidentProteins)-nDecoyProteins
  
  df <- data.frame(nPSM=n,
                   psmFdrCutoff=psmFdrCutoff,
                   nDecoyPSM=nDecoyPSM,
                   nConvidentPSM=nConvidentPSM,
                   nDecoyPeptide=nDecoyPeptide,
                   nConvidentPeptide=nConvidentPeptide,
                   nDecoyProteins=nDecoyProteins,
                   nConvidentProteins=nConvidentProteins,
                   fdrPSM=nDecoyPSM/nConvidentPSM,
                   fdrPeptide=nDecoyPeptide/nConvidentPeptide,
                   fdrProtein=nDecoyProteins/nConvidentProteins
  )
  df     
}
