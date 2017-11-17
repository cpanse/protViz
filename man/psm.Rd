\name{psm}

\alias{psm}
\alias{is.psm}
\alias{summary.psmSet}
\alias{as.data.frame.psmSet}

\title{Compute a matching between a peptide sequence and a MS2 spectrum}

\description{
    The function computes a matching between a given peptide sequence and a given tandem mass spectrum (MS2).
    \code{\link{psm}} determines for each fragment ion mass peak the smallest mass error to a peak in the theretical spectrum from the peptide sequence.
    If the mass error is below the given \code{fragmentIonError} the match is considered a hit. 
    \code{psm} returns a list of computed fragment ions and a vector of the mass error (Da and ppm).
    \code{\link{psm}} uses a generic \code{ANSI-C} function to determines the nearest mass peak of array of double 
    values.
    If the \code{plot} is set to \code{TRUE} an error plot is drawn.

    The function \code{psm} reqires the arguments \code{sequence} and \code{spec}. All other arguments are optional.
}

\usage{

    psm(sequence, 
        spec, 
        FUN=defaultIon,
        plot=TRUE, 
        fi=fragmentIon(sequence, FUN=FUN)[[1]],
        fragmentIonError)

 
}

\arguments{

  \item{sequence}{peptide sequence encoded as character sequence using the 20 amino acid letters.}

  \item{spec}{MS2 which is a R list having a sorted mZ vector and an intensity vector of the same size.}

  \item{FUN}{this function is passed to the \code{fragmentIon} function. the function to be applied to compute further ions. If no function is assigned \code{\link{fragmentIon}} will
    use \code{defaultIon}.}
  \item{plot}{boolean if the error plot function is to be called.}

  \item{fi}{fragment ion table, if not specified \code{fragmentIon(sequence, FUN=FUN)[[1]]} is called. }
  \item{fragmentIonError}{fragment ion error cut-off. default is 0.6 Da. }

}

\details{

    This function can be very useful to make assignements to spectra. Moreover it is used for validation or clairifying ambiguities between different sequences assigned to the same spectrum. Additinally it can be used to generate spectral libraries.

}

\value{ returns a psm match.}


\author{

   Christian Panse 2007, 2008, 2009, 2010, 2012, 2017

}

\seealso{

	\code{\link{peakplot}} and \code{\link{fragmentIon}}

}
\examples{

spec <- list(scans=1138,
    title="178: (rt=22.3807) [20080816_23_fetuin_160.RAW]",
    rtinseconds=1342.8402,
    charge=2,
    mZ=c(195.139940, 221.211970, 239.251780, 290.221750, 
    316.300770, 333.300050, 352.258420, 448.384360, 466.348830, 
    496.207570, 509.565910, 538.458310, 547.253380, 556.173940, 
    560.358050, 569.122080, 594.435500, 689.536940, 707.624790, 
    803.509240, 804.528220, 822.528020, 891.631250, 909.544400, 
    916.631600, 973.702160, 990.594520, 999.430580, 1008.583600, 
    1017.692500, 1027.605900),
    intensity=c(931.8, 322.5, 5045, 733.9, 588.8, 9186, 604.6,
    1593, 531.8, 520.4, 976.4, 410.5, 2756, 2279, 5819, 2.679e+05,
    1267, 1542, 979.2, 9577, 3283, 9441, 1520, 1310, 1.8e+04,
    587.5, 2685, 671.7, 3734, 8266, 3309)
    )

m <- psm('HTLNQIDSVK', spec,plot=TRUE)
hist(m$mZ.Da.error)
hist(m$mZ.ppm.error)

}

