#R

.is.comet <- function(x){
    if(sum(c('sp_score', 'protein', 'plain_peptide') %in% names(x))==3) return(TRUE)

    FALSE
}

#' comet summary
#'
#' @param x a comet txt file
#' @param psmFdrCutoff peptide spectrum match (PSM) cut-off  value.
#'
#' @description computes peptide spectrum match (PSM), peptide and protein False
#' Discovery Rate (FDR) for a given PSM cut-off  value.
#' 
#' @return a \code{data.frame}
#' @author Jonas Grossmann, Christian Panse, 2019
#' @examples
#' 
#'  \dontrun{
#'  S <- read.table("~/Downloads/20190103_007_autoQC4L_autoQC.txt",
#'     header = TRUE, fill = TRUE, skip=1)
#'  
#'  do.call('rbind', lapply(c(0.001,0.005,0.01,0.02,0.05,0.1),
#'    function(cutoff)protViz:::summary.comet(S, cutoff)))
#'  }
.summary.comet <- function(x, psmFdrCutoff=0.05){
    
    stopifnot(.is.comet(x))
    
    idx <-order(x$sp_score)
    n <- nrow(x)
    x <- x[rev(idx), ]
    
    x$REV <- grepl("^REV_", x$protein)
    x$nREVhits <- cumsum(x$REV)
    x$FDR <- S$nREVhits / (seq(1,n) - x$nREVhits)
    
    nConvidentPSM <- which(x$FDR > psmFdrCutoff)[1]
    nDecoyPSM <- x$nREVhits[which(x$FDR > psmFdrCutoff)[1]]
    convidentPeptide <- as.character(x$plain_peptide[seq(1, nConvidentPSM)])
    decoyPeptide <- grepl("^REV_", as.character(x$protein[seq(1, nConvidentPSM)]))
    
    nDecoyPeptide <- length(unique(convidentPeptide[decoyPeptide]))
    nConvidentPeptide <- length(unique(convidentPeptide[!decoyPeptide]))
    
    convidentProteins <- unique(sapply(strsplit(as.character(x$protein[seq(1, nConvidentPSM)]),','), function(y)y[1]))
    nDecoyProteins <- length(grep("^REV_", convidentProteins))
    
    
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








