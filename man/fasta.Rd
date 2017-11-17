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
  fname <- system.file("extdata", name='sample.fasta', package = "protViz")
  F <- Fasta$new(fname)
  F
  F$summary()
  F$getTrypticPeptides()
}

\keyword{FASTA}% use one of  RShowDoc("KEYWORDS")
\keyword{AA}% __ONLY ONE__ keyword per line
