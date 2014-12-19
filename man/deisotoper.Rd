\name{deisotoper}

\alias{deisotoper}

\title{find isotop pattern in a given mass spectrum}

\description{

    \code{deisotoper} returns the indices and scores of the isotop pattern.
    The score is computed by considering the isotop L2 norm of a given iostop
    model. \code{protViz} contains the \code{averagine} data.

}

\details{

    The output of the \code{deisotoper} function is a \code{list} data structure containing 
    for each input MS data set a \code{result}, \code{score}, and \code{group} lists.

    \code{result} is a list which contains for each charge state a list of isotop groups.
    Analoge to \code{result} a score is provide in the \code{score} list.

    \code{group} provides information about the charge states of each isotop cluster.

    The algorithm for finding the isotop chains as implemented in protViz is based on the idea of 
    'Features-Based Deisotoping Method for Tandem Mass Spectra'.
    Each peak represents one node of a graph \code{G=(V,E)} with \code{V = {1, ..., n}}.
    The edges are defined by \code{(v, node_array_G[v])} iff \code{node_array_G[v] > -1}.
    In other words an edge is defined between two peaks if the mZ difference is below the given 
    massError. The isotop chains are determined by a DFS run.
    The time complexity is \code{O(n log (n))} where n is the number of peaks.
    
}

\usage{

    deisotoper(data, Z=1:4, isotopPatternDF=averagine, massError=0.005, plot=FALSE)


}

\arguments{

  \item{data}{A mass spectrometric measurement containing a list.}

  \item{Z}{Charge states to be considered. The default is \code{Z=1:4}}

  \item{isotopPatternDF}{ a data frame containing isotope envelopes of peptide averagine.}

  \item{massError}{ mass error in Dalton. default is 0.005.}

  \item{plot}{boolean if the isotops should be plotted. The default is false.}

}

\author{

   Witold Eryk Wolski, Christian Trachsel, and Christian Panse 2013 

}

\references{

        isotopic-cluster graphs:

        Features-Based Deisotoping Method for Tandem Mass Spectra,
        Zheng Yuan Jinhong Shi Wenjun Lin Bolin Chen and Fang-Xiang Wu
        Advances in Bioinformatics
        Volume 2011 (2011), Article ID 210805, 12 pages
        \url{http://www.hindawi.com/journals/abi/2011/210805/}
}

\examples{
    # example1 - tandem ms
    x <- list(mZ = c(726.068, 726.337, 726.589, 726.842, 727.343, 727.846, 728.346, 
        728.846, 729.348, 730.248, 730.336, 730.581, 730.836),
        intensity = c(6.77850e+03, 2.81688e+04, 6.66884e+04, 1.22032e+07, 
            9.90405e+06, 4.61409e+06, 1.50973e+06, 3.33996e+05, 5.09421e+04, 
            1.15869e+03, 2.14788e+03, 5.37853e+03, 5.79094e+02))

    # the plain C interface function
    out1 <- .Call("deisotoper_main", x$mZ, x$intensity, Z=1:4, averagine, 
        massError=0.01, PACKAGE="protViz")

    out<-deisotoper(data=list(x), Z=2:4, isotopPatternDF=averagine)

    # example2 - a ms from heavy labeld peptide in water background    
    x <- list(mZ=c(642.572, 643.054, 643.569, 644.062, 644.557),
        intensity=c(17000, 25000, 12000, 9000,4000))

    diff(x$mZ)
    diff(x$mZ,lag=2)
    diff(x$mZ, difference=2)

    out2<-deisotoper(data=list(x), Z=1:2, isotopPatternDF=averagine, 
        massError=0.02, plot=TRUE)

}
