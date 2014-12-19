\name{PTM_MarkerFinder_util}

\alias{PTM_MarkerFinder_util}

\title{PTM MarkerFinder util plot}

\description{

\code{PTM_MarkerFinder_util} is a utility function for \code{PTM_MarkerFinder}.

}

\usage{

PTM_MarkerFinder_util(dataFileName, 
    mZmarkerIons,
    minMarkerIntensityRatio,
    minNumberIons,
    itol_ppm,
    write_csv)

}

\arguments{

  \item{dataFileName}{RData file name without 'RData' ending.}

  \item{mZmarkerIons}{A double vector of the m/z patterns which should be searched.}

  \item{minMarkerIntensityRatio}{The marker ions 
intencity percentage compared to the sum of all peak 
intensities.}

  \item{minNumberIons}{Minimal number of marker ions to 
be found for further analysis.}

  \item{itol_ppm}{The ion tolerance of the marker ions 
in ppm. default is set to 10ppm.}


  \item{write_csv}{boolean.}

}

\details{
    The function plots summeries in form of boxplots pie charts and scatter plots of the found marker ions.
}

\author{

   Paolo Nanni, Christian Panse 2012-2013;

}
