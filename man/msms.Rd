\docType{data}

\name{msms}

\alias{msms}

\title{A data set containing tandem mass spectra of an LCMS experiment.}


\description{

    The purpose of \code{msms} is for the demonstration of the peptide 
    identification methods, e.g., \code{peakplot} or \code{psm}.

    \code{msms} contains two tandem mass spectra meassured on a linear ion 
    trap coupled to an Orbitrap mass analyzer 
    (Thermo Scientific, Bremen, Germany).

    The protein sample originates from 10% wild type
    mouse brain homogenate, prepared in sterile PBS containing protease
    inhibitors (Complete; Riche, Switzerland) by repeated extrusion through
    syringe needles of successive smaller size. To remove gross cellular
    debris the homogenate was centrifuged for 10 min at 2400 rpm
    at 4 degrees centigrade.

    The data were acquired withing the Sixth Framework Programm (FP6-37457). 
    The project was called SysProt.
}

\author{

   Dorothea Rutishauser, Jonas Grossmann, Christian Panse 2008

}

\examples{
    
    data(msms)
    peakplot("TAFDEAIAELDTLNEESYK", msms[[1]])

}

\keyword{datasets}
