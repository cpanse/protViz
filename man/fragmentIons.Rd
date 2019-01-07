\name{fragmentIon}

\alias{fragmentIon}
\alias{as.data.frame.fragmentIonSet}
\alias{bymatrix}
\alias{ionseries}
\alias{defaultIon}

\title{Compute the b and y Fragment Ions of a Peptide Sequence}

\description{

    The function computes the fragment ions of a peptide sequence,
    according to the rules of fragmentation in a collision cell.
    It can be either CID fragmentation (collision-induced dissociation)
    or ETD (electron transfer dissociation) which are common in proteomics
    experiments. All are in a positive mode.

    If multiple peptide sequences are given it returns a list of fragment ion table.
    
}

\usage{

    fragmentIon(sequence, FUN, modified, modification, N_term, C_term)

}

\arguments{

  \item{sequence}{peptide sequence encoded as character sequence using the 20 amino acid letters or a double vector of amino acid weights.}

  \item{FUN}{the function to be applied to compute further ions. 
  If no function is assigned \code{fragmentIon} will use \code{defaultIon}.}

  \item{modified}{a vector of interger containing the varmod id of the peptide position.}

  \item{modification}{a double vector of defining the varmod mass in Dalton.}

  \item{N_term}{N-Term mass in Dalton. Default is set to 1.007825.}

  \item{C_term}{C-Term mass in Dalton. Default is set to 17.002740.}

}

\details{

    The fragment ions of a peptide can be computed following the rules proposed in PMID:6525415.
    Beside the \code{b} and \code{y} ions the \code{FUN} argument of \code{fragmentIon} defines
    which ions are computed.
    the default ions beeing computed are defined in the function \code{defaultIon}.
    The are no limits for defining other forms of fragment ions for ETD (c and z ions) CID (b and y ions).

    NOTE:
    for simplicity, we have set a Carbamidomethyl (C) fixed modification with 
    160.030649 (Cysteine mono mass is 103.00919).
    The fixed modifications setting are not enabled in the package yet.

    If only a vector of amino acids weights is given \code{fragmentIon} computes the ions series according to weights.
    If this case applies, the function returns only a list having one data frame.
    
}

\references{

      \url{http://www.ncbi.nlm.nih.gov/pubmed/6525415}

      Protein Sequencing and Identification Using Tandem Mass Spectrometry, 
      Michael Kinter and Nicholas E. Sherman, 
      Wiley-Interscience; 1st edition (September 18, 2000)


}

\author{

   Christian Panse, Bertran Gerrits 2006. 
}

\seealso{

   first used in \code{\link{peakplot}}.

}

\examples{

# Example 1
peptide.AA<-"KINHSFLR";

peptide.AA.weights <- c(128.09496,113.08406,114.04293,
    137.05891,87.03203,147.06841,113.08406,156.10111);

fragmentIon(peptide.AA);

fragmentIon(peptide.AA.weights);

HCD_Ion <- function(b, y){
    return(cbind(b = b, y = y))
}


ETD_Ion <- function(b, y){
    Hydrogen <- 1.007825
    Oxygen <- 15.994915
    Nitrogen <- 14.003074

    y_0 <- y - Oxygen - Hydrogen - Hydrogen
    c <- b + (Nitrogen + (3 * Hydrogen))
    z <- y - (Nitrogen + (3 * Hydrogen))

    return(cbind(y_0, c, z))
}

fragmentIon(peptide.AA, FUN = ETD_Ion)
        
peptides<-c('HTLNQIDSVK', 'ALGGEDVR', 'TPIVGQPSIPGGPVR')

pim <- parentIonMass(peptides)
fi <- fragmentIon(peptides)
(df <- as.data.frame(fi))

op <- par(mfrow=c(3,1)); 
for (i in 1:length(peptides)){
    plot(0, 0,
    xlab='m/Z',
    ylab='',
    xlim=range(c(fi[[i]]$b,fi[[i]]$y)),
    ylim=c(0,1),
    type='n',
    axes=FALSE,
    sub=paste(peptides[i], "/", pim[i], "Da"));
    box()
    axis(1, fi[[i]]$b, round(fi[[i]]$b,1), las=2)
    axis(1, fi[[i]]$y, round(fi[[i]]$y,1), las=2)

    pepSeq<-strsplit(peptides[i], "")
    axis(3,fi[[i]]$b, paste("b", row.names(fi[[i]]),sep=''),las=2)
    axis(3,fi[[i]]$y, paste("y", row.names(fi[[i]]),sep=''),las=2)

    text(fi[[i]]$b, rep(0.3, nchar(peptides[i])), 
    pepSeq[[1]],pos=3,cex=4, lwd=4, col="#aaaaaaaa")

    abline(v=fi[[i]]$b, col='red') 
    abline(v=fi[[i]]$y, col='blue',lwd=2)
}
par(op)


fi <- fragmentIon(c("ATSFYK","XGXFNAGVGK"))[[2]]
fi$b[1] + fi$y[9]
fi$b[2] + fi$y[8]


ION2C <- function(b, y){
    Hydrogen <- 1.007825
    Oxygen <- 15.994915
    Nitrogen <- 14.003074

    # yo <- fi$y - Oxygen - Hydrogen - Hydrogen
    c <- b + (Nitrogen + (3 * Hydrogen))
    z <- y - (Nitrogen + (3 * Hydrogen))
    
    # compute doubly charged fragment ions
    b2 <- (b + Hydrogen) / 2
    y2 <- (y + Hydrogen) / 2

    return(cbind(b, y, b2 ,y2))
}


}
