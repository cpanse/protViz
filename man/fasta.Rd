\name{FASTA}
\alias{FASTA}
\alias{is.FASTA}
\alias{as.FASTA}
\alias{summary.FASTA}
\alias{fcat.FASTA}
\alias{tryptic_digest.FASTA}
\alias{read.FASTA}

\title{
FASTA 
}
\description{
support limited for protein sequences.
}
\usage{
FASTA(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a vector of amino acids or FASTA formated text.
}
}
\details{
Implements FASTA object S3 methods.
}
\value{
returns a data.frame or a S3 FASTA object.
}
\references{
\url{https://en.wikipedia.org/wiki/FASTA}
}
\author{
Christian Panse <cp@fgcz.ethz.ch> 2006-2017
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
\dontrun{
 F <- read.FASTA("~/p1875_db10_20170817.fasta")
 summary(F)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{FASTA}% use one of  RShowDoc("KEYWORDS")
\keyword{AA}% __ONLY ONE__ keyword per line
