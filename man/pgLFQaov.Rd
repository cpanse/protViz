\name{pgLFQaov}

\alias{pgLFQaov}

\title{iTRAQ two group analysis}

\description{

    The function performes an ANOVA for data sets like \code{pgLFQ}.

}

\usage{

    pgLFQaov(data, groups, names, idx, plot, FUN)

}

\arguments{

  \item{data}{a data set like \code{pgLFq}.}
  \item{groups}{a factor.}
  \item{names}{a vector of strings for the main attribute of each plot.}
  \item{idx}{a integer vector of indice to be processed.}
  \item{plot}{logical.  If 'TRUE' (non default), a boxplot is drawn.}
  \item{FUN}{function for doing the data transformation, e.g. \code{log}, \code{asinh}.}

}

\details{

    The methode performs an oav analysis using the R \code{oav} function.
    It returns a vector of the ANOVA \code{"Pr(>F)"} values in the same order.

}

\author{

   Christian Panse, Jonas Grossmann 2012

}
\examples{

    data(pgLFQprot)

    par(mfrow=c(4,3))
    ANOVA<-pgLFQaov(pgLFQprot$"Normalized abundance",
        groups=as.factor(pgLFQprot$grouping),
        names=pgLFQprot$output$Accession,
        plot=TRUE)
    ANOVA


}


