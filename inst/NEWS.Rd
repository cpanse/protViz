\name{NEWS}
\title{News for Package 'protViz'}
\newcommand{\ghpr}{\href{https://github.com/protViz/protViz/pull/#1}{##1}}
\newcommand{\ghit}{\href{https://github.com/protViz/protViz/issues/#1}{##1}}

\section{Changes in protViz version 0.6.8 (2020-04-23)}{
  \itemize{
    \item fix \ghit{8}.
    \item added method for centroiding a spectrum acquired in profile mode.
    \item fix \ghit{11}.
    \item new by labeling using legend with fragment ion ordered by m/z values \ghit{10}.
    \item fix \ghit{9}.
  }
}

\section{Changes in protViz version 0.5.0 (2019-11-18)}{
  \itemize{
    \item added comet.R helper functions.
    \item recfactored peakplot functions.
  }
}

\section{Changes in protViz version 0.4.0 (2019-01-03)}{
  \itemize{
    \item resolve \ghit{6}.
    \item new S3 method for fragmentIon \code{as.data.frame}.
  }
}

\section{Changes in protViz version 0.3.0 (2017-06-22)}{
  \itemize{
    \item refactored \code{?findNN} and \code{?AA} man pages.
    \item \code{ssrc} accepts also a \code{vector} as argument.
    \item \code{as.data.frame.pepxml} for \url{http://tools.proteomecenter.org/wiki/index.php?title=Formats:pepXML}.
  }
}

\section{Changes in protViz version 0.2.45 (2017-11-17)}{
  \itemize{
    \item implement Rcpp module \code{Fasta} \ghit{3}.
    \item added the Swiss-Prot P12763 protein for testing.
    \item added \code{C++} STL \code{lower_bound} method.
	  \item added orcid in DESCRIPTION file 
  }
}

\section{Changes in protViz version 0.2.37 (2017-10-11)}{
  \itemize{
    \item \code{initial as.data.frame} function.
    \item added unit test for \code{ssrc} function.
    \item  deleted  superfluous sigmamix data; please use \url{http://bioconductor.org/packages/msqc1/} described in \url{http://doi.org/10.1002/pmic.201500502}.
  }
}

\section{Changes in protViz version 0.2.31 (2017-05-26)}{
  \itemize{
    \item fixed ``register native routine issue'' by using \code{tools::package_native_routine_registration_skeleton(".")} and renaming .C exports to \code{__findNN_} and \code{__findNN} \ghit{1}.
    
    \item added support for Matrixce Sciense Mascot xml export; new S3classes: \code{mascot} and \code{mascot_query}.
    
    \item cosmetics on \code{protViz.Rnw} vignette file.
    
    \item added \code{ssrc} method for ``Sequence Specific Retention Calculator''.
  }
  
}

\section{Changes in protViz version 0.2.16 (2016-06-01)}{
  \itemize{
    \item PTM Marker Finder:
    \itemize{
      \item added PeptideSequence and modification to csv output.
    }
  }
}

\section{Changes in protViz version 0.2.13 (2016-04-29)}{
  \itemize{
    \item General:
    \itemize{
      \item added NEWS.Rd, CITATION.
      \item changed URL.
      \item fix \code{R CMD check} issues with R version 3.3.
    }
  }
}

\section{Changes in protViz version 0.1.97}{
  \itemize{
    \item BUG FIXES:
    \itemize{
      \item free memory \code{computeFragmentIons.c}
    }
    \item Testing:
    \itemize{
      \item added more tests for \code{deisotoper.R} and \code{fragmentIon.R}.
    }
  }
}
