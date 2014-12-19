\docType{data}
\name{iRTpeptides}

\alias{iRTpeptides}

\alias{irt}
\alias{iRT}

\title{iRT peptides - independent retention time peptides}

\format{contains a table}


\description{

    \code{iRTpeptides} data are used for \code{genSwathIonLib} rt normalization 
    assuming.


}

\author{

       Jonas Grossmann and Christian Panse 2013

}

\examples{

    plot(sort(iRTpeptides$rt))

    plot(pim<-parentIonMass(as.character(iRTpeptides$peptide)) ~ iRTpeptides$rt)

}

\references{

    Using iRT, a normalized retention time for more targeted measurement of peptides.
    Escher C, Reiter L, MacLean B, Ossola R, Herzog F, Chilton J, MacCoss MJ, Rinner O.
    Source
    Proteomics. 2012 Apr;12(8):1111-21. doi: 10.1002/pmic.201100463.


}

