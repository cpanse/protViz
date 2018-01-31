\docType{data}
\name{AA}

\alias{AA}

\title{AA - amino acid table}

\format{contains a table}

\value{returns a \code{data.frame}}

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

	
	
    # computes monoisotopic mass out of formula using the CDK package
    \dontrun{
    	if (require(rcdk)){
	plot(AA$Monoisotopic, 
	     sapply(AA$formula, function(x){
		     get.formula(as.character(x), charge = 1)@mass
	     }))
	}
    }
    \dontrun{
	    if (require(XML)){
	unimodurl <- url("http://www.unimod.org/xml/unimod_tables.xml")
	unimod.list <- XML::xmlToList(
	  XML::xmlParse(
	    scan(unimodurl, what = character())))
	unimod.AA <- data.frame(
	  do.call('rbind', unimod.list$amino_acids))
	rownames(unimod.AA) <- unimod.AA$one_letter
	    }
    }
}

\seealso{

  \itemize{
        \item{\url{http://www.matrixscience.com/help/aa_help.html}}

        \item{\url{http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html}}

        \item{\url{http://www.unimod.org}}
        }

}
