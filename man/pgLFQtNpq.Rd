\name{pgLFQtNpq}

\alias{pgLFQtNpq}
\alias{tNpq}

\title{Label Free Quantification using the top N peptide approach}

\description{

This Function implements the recently emerged TopN strategy which 
uses only the top N intense features for calculating the 
proteinVolume. This approach should reveal a quantitative protein 
value, which should make the protein itself comparable within one 
condition. This allows to estimate protein stochiometries and 
simplifies modelling and calculations with copy numbers per cell.

}


\usage{

    pgLFQtNpq(QuantitativeValue, 
    peptide, protein, N=3, plot=TRUE, FUN=asinh)

}

\arguments{

  \item{QuantitativeValue}{a data set like 
\code{pgLFQfeature$"Normalized abundance"}.}

  \item{peptide}{a vector of peptide sequences.}

  \item{protein}{a vector of protein information.}

  \item{N}{top \code{N} peptide flyers.}

  \item{plot}{logical.  If 'TRUE' (non default), a boxplot is drawn.}

  \item{FUN}{function for doing the data transformation for the 
correlation matrix for the image plot, default transformation is 
\code{asinh}.}

}

\details{

The approach has first been described by Silva et al. in 2005 for 
Waters Q-tof instruments running in the MSe mode. Grossmann et al, 
showed in 2010 that this approach also works for more widely spread 
instruments such as Orbitrap-Velos or FTICR instruments from Thermo.

todo: 
additional columns (or additonal object) for 'protein names' 
and the total number of features assigned to protein in the master map.
The length should be the same as for how many Ns chosen in the assembly method.
Double check, if 'empty' protein names.. - basically - not assigned features.. 
are also in the list! - get rid of it.

}

\author{

   Christian Panse, Jonas Grossmann 2012

}


\references{

    \url{http://www.ncbi.nlm.nih.gov/pubmed/16219938}
    \url{http://www.ncbi.nlm.nih.gov/pubmed/20576481}

}



\examples{

    data(pgLFQfeature)
    par(mfrow=c(2,4), mar=c(4,4,4,1))
    for (i in c(1, 2, 3, 4)){
        tNpq<-pgLFQtNpq(QuantitativeValue=pgLFQfeature$"Normalized abundance",
            peptide=pgLFQfeature$peptideInfo$Sequence,
            protein=pgLFQfeature$peptideInfo$Protein, 
            N=i)
    }

    for (i in c(1, 2, 3, 4)){
        tNpq<-pgLFQtNpq(QuantitativeValue=pgLFQfeature$"Normalized abundance",
            peptide=pgLFQfeature$peptideInfo$Sequence,
            protein=pgLFQfeature$peptideInfo$Protein, 
            plot=FALSE,
            N=i)

        boxplot(t(tNpq), xlab='proteins', ylab='protein value')
    }
}


