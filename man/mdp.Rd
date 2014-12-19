\name{massDeviationPlot}

\alias{massDeviationPlot}
\alias{mdp}

\title{Mass Deviation Plot}

\description{

\code{massDeviationPlot} computes the mass deviation of a data set having 
the list attributes \code{pepmass} and \code{peptideSequence}.

}

\usage{

    mdp(data, sub)

}

\arguments{

  \item{data}{an R data  object as it can be obtained by \code{mascotDat2RData.pl}}
  \item{sub}{a sub title for the plot}

}

\details{
    
    In this version it igrnoes modified peptides.

}

\author{

   Christian Panse 2012-2013; 
   
}

\references{

    Zubarev R, Mann M.
    On the proper use of mass accuracy in proteomics.
    Mol Cell Proteomics. 2007 Mar;6(3):377-81. Epub 2006 Dec 12.
    PMID: 17164402 


}



\examples{
    
    # load("F178767.Rdata")
    # pdf("massError.pdf",12,9)
    # massDeviationPlot (F178767)
    # dev.off()

}

