\name{PTM_MarkerFinder}

\alias{PTM_MarkerFinder}
\alias{pPTM}
\alias{findMz}
\alias{findMz.mascot}
\alias{findMz.psmSet}

\title{PTM MarkerFinder}

\description{

\code{PTM_MarkerFinder} is a fucntion to identify and validate spectra from
peptides carrying post-translational modifications.

}

\usage{

PTM_MarkerFinder(data,
    modification,
    modificationName,
    mZmarkerIons,
    minNumberIons=2,
    itol_ppm=10,
    minMarkerIntensityRatio=5,
    mgfFilename=-1,
    PEAKPLOT=TRUE
    )

findMz(data, mZmarkerIons, itol_ppm = 10, minNumberIons = 2, minMarkerIntensityRatio = 10)

\method{findMz}{psmSet}(data, mZmarkerIons, itol_ppm, minNumberIons, minMarkerIntensityRatio)

\method{findMz}{mascot}(data, mZmarkerIons, itol_ppm = 10, minNumberIons = 2, minMarkerIntensityRatio = 10)

}

\arguments{

  \item{data}{A list of spectra where each list 
element contains a list of \code{mZ, intensity} 
vectors. This object can be received by using the 
\code{mascotDat2RData.pl} perl script.}

  \item{modification}{A double vector containing the 
mono mass PTMs.}

  \item{modificationName}{A character vector containing the 
Name of the PTMs. This is the string which will show up in 
the peptide sequence in square brackets.}

  \item{mZmarkerIons}{The m/z patterns which should be 
searched.}

  \item{minNumberIons}{Minimal number of marker ions to 
be found for further analysis.}

  \item{itol_ppm}{The ion tolerance of the marker ions 
in ppm. default is set to 10ppm.}

  \item{minMarkerIntensityRatio}{The marker ions 
intencity percentage compared to the sum of all peak 
intensities.}

  \item{mgfFilename}{if a mgf(mascot generic file) filename is given, 
  the function writes an mgf file containing
  only the ms2 having the \code{mZmarkerIons}.}

  \item{PEAKPLOT}{If this boolean is \code{FALSE} the \code{peakplot} function is not used.}
}

\details{

The function screens MS2 spectra for the presence of fragment ions specific Post
Translational Modifications (PTMs). The function requires an R-object containing
the mass spectrometric measurement, the peptide assignments, potential
modification information, and a list of marker ions. The R-object can be
retrieved right out of the Mascot Server search result dat files using a perl
script \code{mascotDat2RData.pl} which is included in the package's \code{exec/} 
directory.

The functions iterates over each spectrum of the mass spectrometric measurement
and searches for the as input provided marker ions. If a certain number of
marker ions (default is, that two marker ions are required) are found and the
maker ion intensity ratio is higher than a given threshold, the tandem mass
spectrum is considered as HCD scan type and the corresponding ion series are
drawn using the protViz:peakplot methode. Furthermore the function is searching
for the corresponding ETD scan having the same peptide mass
by screening the succeeding scans for ETD spectra. If such a spectrum is found
the peptide spectrum assignment containing the c, z, and y ions is drawn. Note
that the PTM MarkerFinder expects the ETD scan right after the HCD scan. If the
MS protocol changes the PTM MarkerFinder methode has to be adapted.

For each HCD scan PTM MarkerFinder plots both HCD and ETD scans, a ppm error
versus maker ion m/z scatter plot, a intensity versus marker ion m/z plot, and
two pie charts where the relative and absolute maker ion intensity are shown.

As a summery report the function returns a table containing the following column
attributes: "scans", "mZ", "markerIonMZ", "markerIonIntensity",
"markerIonMzError", "markerIonPpmError", and "query" which can be used for
statistics.

Furthermore, if \code{mgfFilename} is defined, a Mascot Generic File
(MGF) is created containig the HCD scans (having the marker ions) and the corresponding
ETD scans.

}

\author{

   Paolo Nanni, Peter Gehrig, Christian Panse 2011-2013;

}

\seealso{
\code{\link{peakplot}}
}

\references{

  \itemize{
    \item Nanni P, Panse C, Gehrig P, Mueller S, Grossmann J, Schlapbach R.(2013),
    PTM MarkerFinder, a software tool to detect and validate spectra from 
    peptides carrying post-translational modifications.
    Proteomics. 2013 Aug;13(15):2251-5.
    \doi{10.1002/pmic.201300036}.
    
    \item  \code{ADP_Ribose <- c(136.0618, 250.0935, 348.0704, 428.0367)} marker ions have been used in:
    Bilan V, Leutert M, Nanni P, Panse C, Hottiger MO.
Combining Higher-Energy Collision Dissociation and Electron-Transfer/Higher-Energy Collision Dissociation Fragmentation in a Product-Dependent Manner Confidently Assigns Proteomewide ADP-Ribose Acceptor Sites,
    Anal. Chem., 2017, 89 (3), pp 1523-1530
    \doi{10.1021/acs.analchem.6b03365}.
  }

}

\examples{

    # some marker ions
    
    Glykan_MarkerIons <- c(109.02841, 127.03897, 145.04954, 163.06010, 325.11292)
    
    HexNAc_MarkerIons <- c(126.05495, 138.05495, 144.06552, 168.06552, 186.07608, 204.08665)

    # DOI: 10.1021/acs.analchem.6b03365
    # Anal Chem 2017 Feb 13;89(3):1523-1530. Epub 2017 Jan 13.
    ADP_Ribose <- c(136.0618, 250.0935, 348.0704, 428.0367)
    
    data(HexNAc)

    # prepare modification
    ptm.0 <- cbind(AA="-",
        mono=0.0, avg=0.0, desc="unmodified", unimodAccID=NA)
      
    ptm.1 <- cbind(AA='N', 
          mono=317.122300, avg=NA, desc="HexNAc",
                  unimodAccID=2)
        
    ptm.2 <- cbind(AA='M',
            mono=147.035400, avg=NA, desc="Oxidation",
                    unimodAccID=1)

    m <- as.data.frame(rbind(ptm.0, ptm.1, ptm.2), stringsAsFactors = TRUE)

    S <- PTM_MarkerFinder(data=HexNAc, modification=m$mono, 
        modificationName=m$desc,
        minMarkerIntensityRatio=3,
        itol_ppm=20,
        mZmarkerIons=HexNAc_MarkerIons) 


    boxplot(markerIonIntensity ~ markerIonMZ,
        data=S,
        log='y',
        main='Summary plot: boxplot of marker ion intensities from all pPTM spectra',
        xlab='markerIon m/z', 
        ylab='log10 based marker ion intensity')

    # export
    w <- reshape(S[,c(1,7,3,4)],
        direction='wide', 
        timevar="markerIonMZ",
	idvar=c('scans','query'))

    write.table(w,
        file=file.path(tempdir(), "HexNAc_PTM_markerFinder.csv"),
        sep=',',
	row.names=FALSE,
	col.names=TRUE,
	quote=FALSE)
}
