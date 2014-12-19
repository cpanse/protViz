\name{lcmsmap}

\alias{pep3d}
\alias{lcmsmap}
\alias{lcmsoverview}

\title{LC-MS Map}

\description{

The function graphs a LC-MS map having for all charge states and each single charge state.
The rt is on the x axis and the y axis represents the peptide mass. 

}

\usage{

    lcmsmap(data, sub)

}

\arguments{

  \item{data}{an R data  object as it can be obtained by \code{mascotDat2RData.pl}}

  \item{sub}{a sub title for the plot}

}

\author{

   Christian Panse and Jonas Grossmann 2013; 
   
}

\references{

    idea taken from from the ISB TPP Pep3D view
    PMID: 15228367

}
