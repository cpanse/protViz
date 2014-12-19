\name{iTRAQ2GroupAnalysis}

\alias{iTRAQ2GroupAnalysis}

\title{iTRAQ two group analysis}

\description{
    The function performes a two group analysis with a t-test for data sets like \code{iTRAQ}. Result files are generated where the 2 groups are compared for each protein, while a p-value is calculated. If plot equals TRUE, the boxplots are drawn.
}

\usage{

    iTRAQ2GroupAnalysis(data, group1, group2, 
        INDEX, FUN=function(x){return(x)}, plot)

}

\arguments{
  \item{data}{a data set like iTRAQ.}
  \item{group1}{a vector of column ids.}
  \item{group2}{a vector of column ids.}
  \item{INDEX}{list of factors, each of same length as \code{X}.  The elements are coerced to factors by \code{as.factor}.}
  \item{FUN}{function for doing the data transformation, e.g. log, asinh.}
  \item{plot}{boolean. if TRUE boxplots are drawn.}
}

\details{

    t.b.d.

}

\author{

   Christian Panse, Jonas Grossmann 2011

}

\examples{
       data(iTRAQ)
       par(mfrow=c(2,3))
       qProt<-iTRAQ2GroupAnalysis(data=iTRAQ, 
        group1=c(3,4,5,6), 
        group2=7:10, INDEX=iTRAQ$prot, plot=TRUE)
       qProt
       qPeptide<-iTRAQ2GroupAnalysis(data=iTRAQ, 
        group1=c(3,4,5,6), 
        group2=7:10, 
        INDEX=paste(iTRAQ$prot,iTRAQ$peptide), plot=FALSE)
}


