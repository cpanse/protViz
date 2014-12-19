\docType{data}

\name{HexNAc}

\alias{HexNAc}

\title{HexNAc - Analysis of N-HexNAc glycopeptides by LC-MS/MS, using HCD and ETD
fragmentation techniques}

\format{A data set consists of eleven tandem mass spectra.}

\description{

\code{HexNAc} is intented for demonstration of the \code{PTM_MarkerFinder} methode
descriped in [PROTEOMICS pmic.201300036].
It contains the mass pattern (m/z) of 
\code{c(126.05495, 138.05495, 144.06552, 168.06552, 186.07608, 204.08665)}.

A sample containing enriched glycopeptides from yeast was
analysed by LC-MS/MS on a LTQ-Orbitrap Velos instrument using sequential
HCD and ETD fragmentation techniques.
The dataset contains only 11 tandem mass spectra extracted from this
experiment. Three pairs of HCD/ETD spectra correspond to peptides carrying
N-HexNAc modification identified with high confidence (6 spectra in total).
One pair of HCD/ETD spectra corresponds to peptides carrying N-HexNAc
modification, identified with very low confidence (2 spectra in total). The
remaining 3 spectra are from unmodified peptides.

}

\author{

  Paolo Nanni, Christian Panse, 2013

}

\examples{

    data(HexNAc)
    HexNAc[[1]]

    plot(HexNAc[[1]]$mZ, HexNAc[[1]]$intensity, type='h')

}

\keyword{datasets}
