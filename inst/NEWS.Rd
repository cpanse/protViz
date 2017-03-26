\name{NEWS}
\title{News for Package 'protViz'}
\newcommand{\ghpr}{\href{https://github.com/protViz/protViz/pull/#1}{##1}}
\newcommand{\ghit}{\href{https://github.com/protViz/protViz/issues/#1}{##1}}

\section{Changes in protViz version 0.2.30 (2017-03-26)}{
  \itemize{
    \item fixed register native routine issue by using \code{tools::package_native_routine_registration_skeleton(".")}
  }
}

\section{Changes in protViz version 0.2.16 (2016-06-01)}{
  \itemize{
    \item PTM Marker Finder:
    \itemize{
      \item added PeptideSequence and modification to csv output
    }
  }
}

\section{Changes in protViz version 0.2.13 (2016-04-29)}{
  \itemize{
    \item General:
    \itemize{
      \item added NEWS.Rd, CITATION
      \item changed URL
      \item fix \code{R CMD check} issues with R version 3.3
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
      \item added more tests for \code{deisotoper.R} and \code{fragmentIon.R}
    }
  }
}
