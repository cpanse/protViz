\docType{data}
\name{fetuinLFQ}

\alias{fetuinLFQ}
\alias{apex}
\alias{empai}
\alias{t3pq}

\title{fetuinLFQ - A  data set for evaluation of relative and absolute label-free quantification methods.}

\format{A data set with approx. 600 rows and four variables, 4KBytes file size, derived out of 32 single LC-MS runs, 5GBytes size, of a protein mixture.}

\description{

This data set contains the read-out abundance
for three different methods (APEX, emPAI, T3PQ) for seven proteins from
a complex yeast background (abundance levels unchanged) and the
Fetuin that was spiked into the extracts in
concentrations ranging from 20fmol to 300fmol on column (in triplicates).

Samples were analyzed on a LTQ-FTICR Ultra mass spectrometer (Thermo
Fischer Scientific, Bremen, Germany) coupled to an
Eksigent-Nano-HPLC system (Eksigent Technologies, Dublin (CA), USA).


The data can be derived out of the mzXML files:

1       186240760       2008-12-11      20080816_01_fetuin.mzXML

2       179013841       2008-12-11      20080816_02_fetuin_0.mzXML

3       178692924       2008-12-11      20080816_03_fetuin_0.mzXML

4       178608302       2008-12-11      20080816_04_fetuin_0.mzXML

5       180995463       2008-12-11      20080816_05_fetuin_20.mzXML

6       175934898       2008-12-11      20080816_06_fetuin_20.mzXML

7       178364454       2008-12-11      20080816_07_fetuin_20.mzXML

8       184376874       2008-12-11      20080816_08_fetuin_40.mzXML

9       187205031       2008-12-11      20080816_09_fetuin_40.mzXML

10      185646382       2008-12-11      20080816_10_fetuin_40.mzXML

11      183498944       2008-12-11      20080816_11_fetuin_60.mzXML

12      184752098       2008-12-11      20080816_12_fetuin_60.mzXML

13      186243524       2008-12-11      20080816_13_fetuin_60.mzXML

14      187794055       2008-12-11      20080816_14_fetuin_80.mzXML

15      183463368       2008-12-11      20080816_15_fetuin_80.mzXML

16      186331090       2008-12-11      20080816_16_fetuin_80.mzXML

17      188027950       2008-12-11      20080816_17_fetuin_100.mzXML

18      186881098       2008-12-11      20080816_18_fetuin_100.mzXML

19      187219923       2008-12-11      20080816_19_fetuin_100.mzXML

20      187157009       2008-12-11      20080816_20_fetuin_120.mzXML

21      190008885       2008-12-11      20080816_21_fetuin_120.mzXML

22      184226648       2008-12-11      20080816_22_fetuin_120.mzXML

23      190681343       2008-12-11      20080816_23_fetuin_160.mzXML

24      194653066       2008-12-11      20080816_24_fetuin_160.mzXML

25      191045349       2008-12-11      20080816_25_fetuin_160.mzXML

26      184869491       2008-12-11      20080816_26_fetuin_200.mzXML

27      185490782       2008-12-11      20080816_27_fetuin_200.mzXML

28      185635558       2008-12-11      20080816_28_fetuin_200.mzXML

29      187743192       2008-12-11      20080816_29_fetuin_300.mzXML

30      190613235       2008-12-11      20080816_30_fetuin_300.mzXML

31      189570723       2008-12-11      20080816_31_fetuin_300.mzXML

32      226404551       2008-12-11      20080819_32_fetuin.mzXML


contained in a 3GBytes compressed tar ball
\url{http://fgcz-data.uzh.ch/public/fqms.tgz} md5=af804e209844d055c0eded716ef9eea8.

The t3pq data can be derived using the code from \url{https://sourceforge.net/projects/fqms/}
 
}

\author{

   Christian Panse 2009, 2010

}

\examples{

library(lattice)
data(fetuinLFQ)

cv<-1-1:7/10
t<-trellis.par.get("strip.background")
t$col<-(rgb(cv,cv,cv))
trellis.par.set("strip.background",t)

my.xlab="Fetuin concentration spiked into experiment [fmol]"
my.ylab<-"Abundance"

xyplot(abundance~conc|prot*method, 
    data=fetuinLFQ$apex, 
    groups=prot,
    aspect=1,
    panel = function(x, y, subscripts, groups) {
        if (groups[subscripts][1] == "Fetuin")  {
            panel.fill(col="#ffcccc")
        }
                panel.grid(h=-1,v=-1)
                panel.xyplot(x, y)
                panel.loess(x,y, span=1)
            },
    xlab=my.xlab,
    ylab=my.ylab
)


xyplot(abundance~conc|prot*method,
    data=fetuinLFQ$empai,
    groups=prot,
    aspect=1,
    panel = function(x, y, subscripts, groups) {
        if (groups[subscripts][1] == "Fetuin")  {
            panel.fill(col="#ffcccc")
        }
                panel.grid(h=-1,v=-1)
                panel.xyplot(x, y)
                panel.loess(x,y, span=1)
            },
    xlab=my.xlab,
    ylab=my.ylab
)


xyplot(abundance~conc|prot*method,
    data=fetuinLFQ$t3pq,
    groups=prot,
    aspect=1,
    panel = function(x, y, subscripts, groups) {
        if (groups[subscripts][1] == "Fetuin")  {
            panel.fill(col="#ffcccc")
        }
                panel.grid(h=-1,v=-1)
                panel.xyplot(x, y)
                panel.loess(x,y, span=1)
        if (groups[subscripts][1] == "Fetuin")  {
            panel.text(min(fetuinLFQ$t3pq$conc),
                max(fetuinLFQ$t3pq$abundance),
                paste("R-squared:", 
                round(summary(lm(x~y))$r.squared,2)),
                cex=0.75,
                pos=4)
        }
            },
    xlab=my.xlab,
    ylab=my.ylab
)


}

\references{

      \url{http://www.ncbi.nlm.nih.gov/pubmed/20576481}

}

\keyword{datasets}
