\docType{data}

\name{sigmamix}

\alias{sigmamix}

\title{MS dilution series of humen peptides.}

\description{

    The data set contains a dilution series of the MS Qual/Quant 
    QC Mix standard sample (SigmaAldrich, Buchs, Switzerland) 
    measured on five different mass spec platforms.

}

\details{

    The mass spectrometers are: 

    \itemize{

    \item \code{Qtrap} - {SRM on a Qtrap 5500}

    \item \code{TSQvantage} - {SRM on a TSQ Vantage}

    \item \code{Qexactive} - {PRM/t-MS2 on a Q-Exactive}

    \item \code{SWATH} - {SWATH on the TripleTOF 5600}

    \item \code{TripleTOF_MS1} - {MS1 quantification on the
    TripleTOF 5600 to evaluate similarities as well as differences in 
    quantitative output.}

    }
    
    For each platform there are 
    three technical replicates of a mixture of 14 SILAC 
    labled tryptic peptide pairs in known amount and a complex yeast background..
    The SILAC peptides span a concentration range over three orders
    of magnitude.


    \code{sigmamix$qc} contains the absolute measured peak areas for
    each transition for all the SILAC peptide pairs, all machines and
    all dilutions (in triplicates).
    The peak areas were computed using the Skyline software.

    \code{sigmamix$ratio} contains the ratios 
    derived from \code{log(light/heavy)/log(2)}, 
    delution concentration, peptide, and
    the mass spec devise as columns.

}


\author{

    Asa Wahlander, Christian Panse, Jonas Grossmann 2014

}

\seealso{

        \url{http://www.sigmaaldrich.com/life-science/proteomics/mass-spectrometry/ms-qual-quant-qc-mix.html}, 2014

        Skyline: \url{http://www.ncbi.nlm.nih.gov/pubmed/20147306}

}

\examples{

    data(sigmamix)

    sigmamix$composition

    sigmamix$ratio[1,]
    
    if ((require(lattice))){

        bwplot(as.character(conc) ~ ratio | peptide * msDevice,
            data = sigmamix$ratio)

    }
}

\keyword{datasets}
