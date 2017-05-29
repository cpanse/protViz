\name{NEWS}
\title{News for Package 'protViz'}
\newcommand{\ghpr}{\href{https://github.com/protViz/protViz/pull/#1}{##1}}
\newcommand{\ghit}{\href{https://github.com/protViz/protViz/issues/#1}{##1}}

\section{Changes in protViz version 0.2.32 (2017-05-29)}{
  \itemize{
    \item added unit test for \code{ssrc} function.
    \item  deleted  superfluous sigmamix data; please use \url{http://bioconductor.org/packages/msqc1/} described in \url{http://doi.org/10.1002/pmic.201500502}.
  }
}

\section{Changes in protViz version 0.2.31 (2017-05-26)}{
  \itemize{
    \item fixed ``register native routine issue'' by using \code{tools::package_native_routine_registration_skeleton(".")} and renaming .C exports to \code{__findNN_} and \code{__findNN} \ghit{1}.
    
    \item added suport for Matrixce Sciense Mascot xml export; new S3classes: \code{mascot} and \code{mascot_query}.
    
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
