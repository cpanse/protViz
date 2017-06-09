\name{lcmsmap}

\alias{pep3d}
\alias{lcmsmap}
\alias{lcmsoverview}

\alias{plot.psmSet}


\title{LC-MS Map}

\description{

The function graphs a LC-MS map having for all charge states and each single charge state.
The rt is on the x axis and the y axis represents the peptide mass. 

}

\usage{

    lcmsmap(data, charges=2:3, score.cutoff = 30, \dots)

}

\arguments{

  \item{data}{an R data  object as it can be obtained by \code{mascotDat2RData.pl} or 
  by using mascots export function.}

  \item{charges}{a vector of charge states to be displayed.}

  \item{score.cutoff}{a numeric value taken as score cut-off.}
  
  \item{\dots}{pass arguments to \code{plot} function.}
}

\author{

   Christian Panse and Jonas Grossmann 2013; 
   
}

\references{

    idea taken from from the ISB TPP Pep3D view
    PMID: 15228367

}
