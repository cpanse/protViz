\name{ssrc}

\alias{ssrc}
\alias{hydrophobicity}


\title{Sequence Specific Retention Calculator}

\description{

    This function returns as output one hydrophobicity value 
    for a given sequence of amino acids (tryptic peptide)
    which can be used to predict the retention times.
    The calculation is based on the method described in PMID:15238601.
    
}

\usage{

    ssrc(x, H=list())
}

\arguments{

  \item{x}{sequence of amino acids, e.g., \code{x="ELIVSK"}}

  \item{H}{A list of retention coefficients. The default is set to the values of PMID:15238601 table II column 2(Rc values).}

}

\references{
Krokhin, O. V. et al. An improved model for prediction of retention times of tryptic peptides in ion
    pair reversed-phase HPLC: its application to protein peptide mapping by off-line HPLC-MALDI MS.
    Mol. Cell Proteomics 3, 908-919 (2004).
    \url{http://www.ncbi.nlm.nih.gov/pubmed/15238601}

}
\seealso{ 

  \itemize{
    \item \url{http://hs2.proteome.ca/SSRCalc/SSRCalcX.html}
    
    \item \code{\link{iRTpeptides}}
  }
}

\author{

   Christian Panse, Christian Trachsel 2015

}

\examples{

  # example of table iv [PMID:15238601]
  lapply(c("SCHTAVGR", "SCHTGLGR", "EDLIAYLK"), ssrc)
  
  plot(sapply(as.character(iRTpeptides$peptide), ssrc) ~  iRTpeptides$rt)
}
