\name{centroidv2}
\alias{centroidv2}
\alias{profile}
\title{
  Centroid a spectrum acquired in profile mode
}
\description{
  returns a centroid spectrum of a recorded profile mode spectrum.
}
\usage{
  centroid(mZ, intensity, eps=0.1, debug=FALSE)
}
\arguments{

  \item{mZ}{Numerical vector of profile recorded data and sorted mZ values.}

  \item{intensity}{corresponding intensity values.}

}
\details{
  the method is tested on an Orbitrap Fusion Lumos FSN20242 data set.
}

\value{
  returns a \code{data.frame} with a mZ and a intensity column.
}

\author{
  Witold Wolski, April 2020
}

\note{
  thanks to Nienke Meekel and Andrea Brunner (kwrwater.nl).
}



\examples{
  # Orbitrap Fusion Lumos FSN20242
  # p2722/.../.../stds_pos_neg_MS_highconc_UVPD_50_300.raw
  # scan 1959
  # CC(C)(C)C(O)C(OC1=CC=C(Cl)C=C1)N1C=NC=N1	1
  # exact.mass 295.1088
  # FTMS + p ESI d Full ms2 296.1162@uvpd50.00

    p <- protViz:::.getProfileMS2()

    names(p)
    eps = 0.1
    mZ <- p$mZ
    intensity <- p$intensity
    n <- length(mZ)




    rv2 <- centroidv2(p$mZ, p$intensity)
    rv <- protViz::centroid(p$mZ, p$intensity, debug = FALSE)

    peakidx <- protViz:::.getpeakpos(mZ, intensity)

    par(mfrow = c(2,2))
    plot(mZ, intensity, type = "h", log = "y", xlim = c(147.45,147.6), main = "peak boundaries")
    abline(v = mZ[peakidx[,"peakpos"]], col = 2, lty = 2)
    abline(v = mZ[peakidx[,"start"]], col = 3, lty = 2)
    abline(v = mZ[peakidx[,"end"]], col = 4, lty = 4)


    plot(rv2$mZ, rv2$intensity, type = 'h', log = "y", col = "green", xlim = c(100,200))
    lines(rv$mZ, rv$intensity, type = 'h')
    lines(p$mZ, p$intensity, col = "grey")
    legend("topleft", legend = c("centroid","centroidv2"), col = c("black", "green","grey"), lty = c(1,1))

    plot(rv2$mZ, rv2$intensity, type = 'h', log = "y", col = "green", xlim = c(113,114))
    lines(rv$mZ, rv$intensity, type = 'h')
    lines(p$mZ, p$intensity, col = "grey")
    legend("topleft", legend = c("centroid","centroidv2","raw"), col = c("black", "green","grey"), lty = c(1,1,1))

    plot(rv2$mZ, rv2$intensity, type = 'h', log = "y", col = "green", xlim = c(147.45,147.6))
    lines(rv$mZ, rv$intensity, type = 'h')
    lines(p$mZ, p$intensity, col = "grey")
    legend("topleft", legend = c("centroid","centroidv2","raw"), col = c("black", "green","grey"), lty = c(1,1,1))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{centroid}% use one of  RShowDoc("KEYWORDS")
\keyword{profile}% __ONLY ONE__ keyword per line
