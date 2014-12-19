\docType{data}
\name{AA}

\alias{AA}

\title{AA - amino acid table}

\format{contains a table}


\description{

    Among other attributes it contains '1-letter code', 'monoisotopic mass' and 
    'average mass' of each amino acids.

}

\author{

       Christian Panse 2013

}

\examples{

    data(AA)
    AA
    AA.lm<-lm(AA$Monoisotopic ~ AA$Average)

    plot(AA$Monoisotopic, AA$Average); 
    abline(AA.lm, col='grey')
    text(AA$Monoisotopic, AA$Average, AA$letter1, pos=3)

    plot(AA$Average-AA$Monoisotopic)  
    axis(3,1:20,AA$letter1); 
    abline(v=1:20,col='grey')

}

\references{

          \url{http://www.matrixscience.com/help/aa_help.html}

          \url{http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html}

          \url{http://www.unimod.org}

}
