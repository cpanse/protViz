\name{centroid}
\alias{centroid}
\alias{profile}
\title{
Centroid a spectrum acquired in profile mode
}
\description{
returns a centroid spectrum of a recorded profile mode spectrum.
}
\usage{
centroid(mZ, intensity, tolppm=100, debug=FALSE)
}
\arguments{

  \item{mZ}{Numerical vector of profile recorded data and sorted mZ values.}

  \item{intensity}{corresponding intensity values.}
  
  \item{tolppm}{maximal distance in Da between to profile mZ values. default is set to 100ppm.}

  \item{debug}{if true all peak profiles are plotted. default is false.}

}
\details{
the method is tested on an Orbitrap Fusion Lumos FSN20242 data set.
}

\value{
returns a \code{data.frame} with a mZ and a intensity column.
}

\author{
	Christian Panse and Jonas Grossmann, April 2020
}

\note{
	thanks to Nienke Meekel and Andrea Brunner (kwrwater.nl) and Witold E. Wolski.
}

\seealso{
    \itemize{
    \item{\href{https://github.com/lgatto/MSnbase/blob/2b6d6f5162e7464c39b48d0daa4d981c6ec30bbd/R/functions-Spectrum.R#L649}{MSnbase functions-Spectrum.R}}
    \item{\url{http://proteowizard.sourceforge.net/tools.shtml}}
    \item{\url{https://planetorbitrap.com/rawfilereader}}
    }
    
}

\examples{
   # Orbitrap Fusion Lumos FSN20242 
   # p2722/.../.../stds_pos_neg_MS_highconc_UVPD_50_300.raw
   # scan 1959
   # CC(C)(C)C(O)C(OC1=CC=C(Cl)C=C1)N1C=NC=N1	1
   # exact.mass 295.1088
   # FTMS + p ESI d Full ms2 296.1162@uvpd50.00 

    p <- protViz:::.getProfileMS2()
    
    # determine eps
    plot((diff(p$mZ)) ~ p$mZ[2:length(p$mZ)], log='y');
    abline(h=0.1, col='red')
    points(p$mZ , 1E-4 * p$mZ, col='grey', type='l')
    abline(v=296.1162, col='cyan')

    op <- par(mfrow=c(2, 1))
    plot(p$mZ, p$intensity, type='h',
      main='profile', xlim=c(100,300))
    abline(v=296.1162, col='cyan')
    
    plot(centroid(p$mZ, p$intensity),type='h',
      main="centroid",xlim=c(100,300))
    par(op)
    
    
    op <- par(mfrow=c(2, 1), ask = TRUE)
        rv <- centroid(p$mZ, p$intensity, debug = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{centroid}% use one of  RShowDoc("KEYWORDS")
\keyword{profile}% __ONLY ONE__ keyword per line
