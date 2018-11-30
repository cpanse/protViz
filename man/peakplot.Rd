\name{peakplot}

\alias{peakplot.label}
\alias{peakplot.putlabel}
\alias{peakplot.pie}
\alias{peakplot}
\alias{protViz}
\alias{plot.psm}


\title{Labelling of Peptide Fragment Mass Spectra}

\description{

     \code{peakplot} labels a mass spectrum from a peptide sequence  
     assignment by using the \code{psm} function with the appropriate 
     fragment ion labels.
     Using \code{peakplot} it is possible to compute the fragment ion 
     coverage and customize the fragmentation ions used for labelling.
}

\usage{

    peakplot(peptideSequence, spec, 
        FUN=defaultIon, 
        fi=fragmentIon(peptideSequence, FUN=FUN)[[1]],
        main=NULL,
        sub=paste(peptideSequence, spec$title, sep=" / "),
        xlim=range(spec$mZ, na.rm=TRUE),
        ylim=range(spec$intensity, na.rm=TRUE),
        itol=0.6,
        pattern.abc="[abc].*",
        pattern.xyz="[xyz].*",
        ion.axes=TRUE, ...) 

}

\arguments{

  \item{peptideSequence}{peptide sequence encoded as character sequence using 
the 20 amino acid letters.}

  \item{spec}{ a tandem mass spec data structure which is a list having two 
eq sized vectors calles mZ and intensity. mZ values are sorted.}

  \item{FUN}{the function to be applied to compute further ions. If no 
function is assigned \code{fragmentIon} will
    use \code{defaultIon}.}

  \item{fi}{fragment ion table, if not specified 
  \code{fragmentIon(sequence, FUN=FUN)[[1]]} is called.
  This argument is useful if you have PTM or specific ion series.}

  \item{sub}{a sub title for the plot}
  \item{main}{a main title for the plot}
  \item{xlim}{xlim of the plot; default is \code{range(spec$mZ, na.rm=TRUE)}}
  \item{ylim}{ylim of the plot}

  \item{itol}{fragment ion tolerance default is 0.6 Dalton. values below are considered as hit.}

  \item{pattern.abc}{\code{regexpr pattern} for a, b, c like ions. if the pattern does not match ion are not plotted.}

  \item{pattern.xyz}{\code{regexpr pattern} for x, y, z like ions.}

  \item{ion.axes}{boolean default is \code{TRUE}. determines whether the fragment ion labels should be plotted instead of m/z values.}

  \item{...}{passed to plot function.}

}

\details{

\code{peakplot} computes the in-silico fragment ion, using the 
\code{\link{fragmentIon}} function, and matches them again 
the MS2. If the mass error between the in-silico peak and the measured one 
is below a 0.6 Da (default setting) the peptide spectrum match 
(\code{\link{psm}}) is considered has hit.

The major objective of the labelling is avoidance of overlapping labels for 
which purpose \code{peakplot} applies filtering. A label is only drawn if the 
corresponding ion count of the m/z peak is higher than a given threshold. 
Experience from several hundred annotations shows that the 0.9 percentile is 
a good cut off value. To most efficiently use the limited screen and printout 
space and to ensure that labels representing important local peaks are also 
considered for drawing, \code{peakplot} divides the display space into a 
number of bins depending on the peptide sequence length along the m/z axis. 
From these bins, the top n labels are ordered according to abundance. For the 
visual clustering effect, the abc and xyz ions are drawn on different y-axis 
levels using different colours. Ion types considered for labelling is 
dependent on the instrument setting applied during the initial search.

}

\value{

	returns a list object containing all the peptide spectrum match relevant
	information.

}

\author{

   Bertran Gerrits, Christian Panse 
   
   2006-2017; 

}
   
\references{

  \itemize{
  \item PEAKPLOT: Visualizing Fragmented Peptide Mass Spectra in Proteomics (2009)
  Panse Christian, Gerrits Bertran, Schlapbach Ralph, useR!2009,
  
  abstract: \url{https://www.r-project.org/conferences/useR-2009/abstracts/pdf/Panse+Gerrits+Schlapbach.pdf},
  
  slides: \url{https://www.r-project.org/conferences/useR-2009/slides/Panse+Gerrits+Schlapbach.pdf}.
  }
}

\seealso{
\itemize{

   \item \link{mascot}
   
	 \item{\code{\link{fragmentIon}}}
	 
	 \item{\code{\link{psm}}}
	 
	 \item{\code{\link{PTM_MarkerFinder}}}
		
		\item \code{\link{peakplot}} has been used for generating 
   figures and supplemental material in:    \itemize{

    \item{\url{http://www.ncbi.nlm.nih.gov/pubmed/17114649}}

    \item{\url{http://www.ncbi.nlm.nih.gov/pubmed/17940529}}

    \item{\url{http://www.ncbi.nlm.nih.gov/pubmed/22693562}}

    \item{\url{http://www.ncbi.nlm.nih.gov/pubmed/23713006}}
    
    \item{\url{http://www.ncbi.nlm.nih.gov/pubmed/26944343}}
    }
  }
}


\examples{

    data(msms)
    op <- par(mfrow=c(2,1))
    peakplot("TAFDEAIAELDTLNEESYK", msms[[1]])
    peakplot("TAFDEAIAELDTLSEESYK", msms[[2]])
    par(op)


    # filter cand. fragment ions
    fi <- fragmentIon("TAFDEAIAELDTLNEESYK")
    fi.cyz <- as.data.frame(cbind(c=fi[[1]]$c, y=fi[[1]]$y, z=fi[[1]]$z))

    p <- peakplot("TAFDEAIAELDTLNEESYK", spec=msms[[1]],
        fi=fi.cyz,
        itol=0.6,
        ion.axes=FALSE)
        
          
    # customizing the ion series
    ions <- function(b, y){
      Hydrogen <- 1.007825
      Oxygen <- 15.994915
      Nitrogen <- 14.003074
  
      y0 <- fi$y - Oxygen - Hydrogen - Hydrogen
      c <- b + (Nitrogen + (3 * Hydrogen))
      z <- y - (Nitrogen + (3 * Hydrogen))
      
      return(cbind(b, y, c ,z, y0))
  }

}

