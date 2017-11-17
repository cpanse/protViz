\name{Fasta}
\alias{FASTA}
\alias{Fasta}
\alias{fasta}

\title{
FASTA format reader
}
\description{
Implements as Rcpp module performing tryptic digests and summaries.
The support is limited for protein sequences.
}
\usage{
Fasta
}

\value{
  returns a \code{list} or \code{vector}  objects.
}
\references{

  \url{https://en.wikipedia.org/wiki/FASTA_format}
}
\author{
  Christian Panse <cp@fgcz.ethz.ch> 2006-2017
}

\examples{
# >sp|P12763|FETUA_BOVIN Alpha-2-HS-glycoprotein OS=Bos taurus
  fname <- system.file("extdata", name='P12763.fasta', package = "protViz")
  F <- Fasta$new(fname)
  F
  F$summary()
  F$getTrypticPeptides()
  
  plot(parentIonMass(F$getTrypticPeptides()), ssrc(F$getTrypticPeptides()),
    log='xy')
  text(parentIonMass(F$getTrypticPeptides()), ssrc(F$getTrypticPeptides()), 
    F$getTrypticPeptides(), pos=3, cex=0.5)

}

\keyword{FASTA}% use one of  RShowDoc("KEYWORDS")
\keyword{AA}% __ONLY ONE__ keyword per line
