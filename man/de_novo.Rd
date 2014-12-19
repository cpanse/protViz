\name{de_novo}

\alias{de_novo}

\title{de-novo on tandem ms}

\description{
        This function computes the mass differences of MS2 peaks and matches the        'delta masses' to the amino acid weights.
}

\usage{

    de_novo(data)

}

\arguments{
  \item{data}{ms data set.}
}

\details{

    This function should only be used for demonstration!

}

\author{

   Christian Panse 2013

}
\examples{
    data(msms)
    de_novo(msms[[1]])
}

