\docType{data}
\name{pressureProfile}

\alias{pressureProfile}

\title{NanoLC pressure profile}

\description{
    
pressureProfile is a data frame with a variable number of 
cases (rows) and two variables (columns) named time(min), Pc(psi)

}

\details{

The \code{pressureProfile} data set gives the flow profile in nL/min and 
the pressure profile in psi as a function of time, respectively, for 24 
consecutive injections on an Eksigent NanoLC-Ultra 1D plus system.

The dataset consists of 24 HPLC pressure profiles of 140 min. Data was 
acquired with a Eksigent NanoLC-Ultra 1D plus system. Samples were 
separated with a 100 1E-6m ID. column of 150 mm lenght packed with C-18 
AQ 200A beads (Dr. Maisch GmbH). The column was heated to a temperature 
of 50 grad C during the entire run. Acquisition queue was as follows:

- A block of three (3) samples with a gradient length of 140 min was 
allways followed a autocalibration run with a gradient length of 25 min.

- A total of seven (7) such blocks was acquired to cover the entire 
dataset.

 
}

\author{

   Christian Trachsel 2012

}

\examples{

data(pressureProfile)
ppp(pressureProfile[pressureProfile$filename=="F01",])

par(mfrow=c(1,1))
pps.data<-pps(pressureProfile, time=seq(1,140,by=5))
boxplot(Pc~time, data=pps.data, 
    xlab='time [in minutes]', ylab='Pc(psi)')

library(lattice)
pps.data<-pps(pressureProfile, time=seq(25,40,by=5))
xyplot(Pc ~ as.factor(file) | paste("time =", as.character(time), "minutes"),
    panel = function(x, y){
        m<-sum(y)/length(y)
        m5<-(max(y)-min(y))*0.05
        panel.abline(h=c(m-m5,m,m+m5), col=rep("#ffcccc",3),lwd=c(2,4,2))
        panel.grid(h=-1, v=0)
        panel.xyplot(x, y)
    },
    ylab='Pc [psi]',
    layout=c(1,4),
    sub='The three read lines indicate avg plus min 5\%.',
    scales = list(x = list(rot = 45)),
    data=pps.data)

}

\keyword{datasets}
