\docType{data}
\name{averagine}
\alias{averagine}
\title{averagine - a data set conatining isotope envelopes of averagine
peptides}

\format{first column mass, column 1:8 l2 normalized isotope intensities}

\description{

generated using proteowizards \code{IsotopeEnvelopeEstimator} class.

}


\author{

Witold Wolski and Christian Panse 2013
}

\examples{

data(averagine)
r<-seq(0,1,length=200); cm<-c(rgb(r,r,r), '#555599')
image(m<-as.matrix(averagine), col=rev(cm), axes=FALSE, main='protViz
averagine data - normalized isotope intensities ',
sub='the blue color indicates fields with value 0.',
xlab='isotops',
ylab='mass');
box()
axis(1, seq(0,1,length=nrow(m)), 1:nrow(m));
axis(2, seq(0,1,length=10), colnames(m)[seq(1,ncol(m), length=10)])

}

\references{

A general approach to calculating isotopic distributions for mass spectrometry James A. Yergey Volume 52, Issues 2 1 
September 1983, Pages 3379

\url{http://www.sciencedirect.com/science/article/pii/0020738183850530}
\url{http://proteowizard.sourceforge.net}

}


