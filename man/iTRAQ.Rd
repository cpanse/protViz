\docType{data}
\name{iTRAQ}

\alias{iTRAQ}

\title{iTRAQ - A small 8-plex iTRAQ  data set with confident identified peptides from 5 proteins.}

\format{A data set with approx. 260 rows and 10 variables, 40KBytes file size, derived from biological data from an iTRAQ project..}

\description{
This dataset contains quantification data from iTRAQ (isobaric tag for relative and absolute quantification) LC-MSMS identified peptides using the software ProteinPilot from ABSciex. 
The data is generated on a Thermo LTQ-Oribtrap. The mascot generic files are produced with MascotDistiller and searched with ProteinPilot V.4. The peptide summary export is used and filtered towards confident peptides (Confidence above 90) and only columns with Protein accession, peptide sequence, the 8 reporter ion channels (area) are kept as well as the protein description. The rest of the peptide summary is deleted. This data is only a small illustrative part out of the whole data produced.
The package design is done in a way, that it is not restricted to ProteinPilot exports but can be used with all other kind of iTRAQ tools, as long as each reporter channel with areas or intensities can be exported. The main point of this function is, to get a hand on what peptides and proteins are used to quantify. Also to be able to not use ratio based quantification but to clearly specify the input and the experimental design.}

\author{

  Jonas Grossmann, Christian Panse, 2012

}

\examples{

    data(iTRAQ)
    iTRAQ[1:10,]

}

\references{

      \url{http://www.ncbi.nlm.nih.gov/pubmed/23000593}

}

\keyword{datasets}
