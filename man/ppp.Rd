\name{pressureProfilePlot}

\alias{pressureProfilePlot}
\alias{ppp}

\title{Plotting pressure profile data from Eksigent LC pumps}

\description{
    plots the pressure profiles data on a scatter plot \code{Pc} versus \code{time} grouped by time range. 
}

\usage{

    ppp(data)

}

\arguments{

  \item{data}{A data set like \code{pressureProfile}. The data must have an attribute \code{file}, \code{time}, and \code{Pc}.}

}

\details{

    This is a useful function to generate an overview over a sequence of measurments run with the same LC gradient. It shows the behaviour of the LC pumps basically the pressures over the whole gradient are plotted to spot sources of problems. 

}

\author{

    Bernd Roschitzki, Christian Trachsel, Christian Panse 2012

}

\examples{

    data(pressureProfile)
    ppp(pressureProfile[pressureProfile$filename=="F12",])

}
