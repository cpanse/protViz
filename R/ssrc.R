#R


# $HeadURL: https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/specL/R/ssrc.R $
# $Id: ssrc.R 104592 2015-06-05 13:23:12Z c.panse $
# $Date: 2015-06-05 15:23:12 +0200 (Fri, 05 Jun 2015) $


# Sequence Specific Retention Calculator
# http://www.ncbi.nlm.nih.gov/pubmed/15238601
# http://hs2.proteome.ca/SSRCalc/SSRCalcX.html

.ssrc <- function(x, H){
  if (is.na(x) || nchar(x) < 3){
    return (NA)
  }
  
  if (!is.character(x)){
    x <- as.character(x)
  }
  
  if (length(H) != 20){
    H <- list()
    H[['W']] <- 11.0; 
    H[['F']] <- 10.5; 
    H[['L']] <- 9.6; 
    H[['I']] <- 8.4; 
    H[['M']] <- 5.8; 
    H[['V']] <- 5.0; 
    H[['Y']] <- 4.0; 
    H[['A']] <- 0.8; 
    H[['T']] <- 0.4; 
    H[['P']] <- 0.2; 
    H[['E']] <- 0.0; 
    H[['D']] <- -0.5; 
    H[['C']] <- -0.8; 
    H[['S']] <- -0.8; 
    H[['Q']] <- -0.9; 
    H[['G']] <- -0.9; 
    H[['N']] <- -1.2; 
    H[['R']] <- -1.3; 
    H[['H']] <- -1.3; 
    H[['K']] <- -1.9; 
}
    sumRc <- sum(unlist(H)) / 20
    seq <- strsplit(x, "")[[1]]
    
    n <- length(seq)
    Kl <- 1
    if (n < 10){
      Kl <- 1 - (0.027 * (10 - n))
    } else if (n > 20){
      Kl <- 1 - (0.014 * (n - 20))
    }
    
    Hyd <- Kl * sum(unlist(H[seq])) + (0.42 * (sumRc - unlist(H[seq[1]]))) + (0.22 * (sumRc - unlist(H[seq[2]]))) + (0.05 * (sumRc - unlist(H[seq[3]])))

    if (Hyd >= 38){
      Hyd <- Hyd - 0.3 * (Hyd - 38) 
    }
    names(Hyd) <- x
    return(Hyd)
}

# lapply(c("SCHTAVGR", "SCHTGLGR", "EDLIAYLK"), ssrc)

ssrc <- function(x, H=list()){
 if (is.vector(x)){
   rv <- vapply(x,
       FUN = function(xx){.ssrc(xx, H)}, 
       FUN.VALUE = .ssrc("ELVISR", H))
 } else{
  rv <- .ssrc(x, H)
 }
 names(rv) <- x
 rv
}
