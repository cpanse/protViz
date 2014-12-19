\name{pressureProfileSummary}

\alias{pressureProfileSummary}
\alias{pps}

\title{A misc function for finding NN time slots in pressure profile data.}

\description{

    The function computes a list of Pc values for the in time provided NN using findNN. 
    It returns file, time, timediff and the Pc valyes as list.

}

\usage{

    pps(data, time)

}

\arguments{

    \item{data}{the data set to be plotted. It reqires the following attributes: \code{Pc, file, time}.}
    \item{time}{specifies the timeslots to be used.}

}

\details{

        The function is useful for levelplots. It can be used to make a quality check over a sequence of experiments. At our site, it is used in conjugtion with Eksigent pumps. These devices store a text file at the instrument pc. These files are preprocessed as specified in the file structure.

}

\author{

   Christian Panse 2012

}

\examples{

library(lattice)
data(pressureProfile)

# do the pre processing
pps.data<-pps(pressureProfile, time=seq(25,40,by=5))

print(xyplot(Pc ~ as.factor(file) | paste("time =", 
    as.character(time), "minutes"),
    panel = function(x, y){
        m<-sum(y)/length(y)
        m5<-(max(y)-min(y))*0.05
        panel.abline(h=c(m-m5,m,m+m5), 
            col=rep("#ffcccc",3),lwd=c(1,2,1))
        panel.grid(h=-1, v=0)
        panel.xyplot(x, y)
    },
    ylab='Pc [psi]',
    layout=c(1,4),
    sub='The three read lines indicate avg plus min 5\%.',
    scales = list(x = list(rot = 45)),
    data=pps.data))

}
