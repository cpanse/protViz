\docType{data}
\name{pgLFQfeature}

\alias{pgLFQfeature}

\title{pgLFQfeature - A data set with a featuremap export ProgenesisLCMS}

\format{A data object from ProgenesisImporter.R}

\description{
This data set contains a fraction of the most abundant features from 24 LC-MSMS runs measured on 
an Orbitrap-Velos. The data contains digested proteins from human HeLa cells infected with Shigella bacteria grown on a time course.
It is structured in a way, that 6 biological replicates from 4 conditions are measured (Not_infected, Infected_1hr, Infected_2hr, Infected_3hr). 
The 24 LC-MSMS runs are aligned with each other and a so called mastermap (feature map) is generated where "normalized volumes" on MS1 are extracted for all features in the respective LC-MSMS runs. Some of the features are further annotated with peptide sequences and protein accessions (using the Mascot search algorithm and an ion score cutoff of 25).

We realized, depending on your language and keyboard setting, the separator for the FeatureData as well as for ProteinMeasurements are different (semicolon and commas are used depending on your setting). We assume, that semicolons are and the individual cells are escaped by ". If this differes, we have an option that can be switched in the pgImporter function.



More information on the commercial software can be found here:
\url{http://www.nonlinear.com/products/progenesis/lc-ms/overview/}.

}

\author{

   Christian Panse, Jonas Grossmann 2012

}

\examples{
data(pgLFQfeature)
op<-par(mfrow=c(1,1),mar=c(18,18,4,1),cex=0.5)
samples<-names(pgLFQfeature$"Normalized abundance")
image(cor(asinh(pgLFQfeature$"Normalized abundance")),
    col=gray(seq(0,1,length=20)),
    main='pgLFQfeature correlation',
    axes=FALSE)

axis(1,at=seq(from=0, to=1, 
    length.out=length(samples)), 
    labels=samples, las=2)

axis(2,at=seq(from=0, to=1, 
    length.out=length(samples)), labels=samples, las=2)
par(op)

# example assemble (and quantify) all proteins from peptides 
# using intensities from the master(feature) map
# 1239 features
pgLFQfeature$peptideInfo$Sequence[1239]
pgLFQfeature$peptideInfo$Protein[1239]
# conflicts can be produced through: 
# 1. shared peptides, 
# 2. mapping to more than one peptide, 
# 3. rank two above threshold
pgLFQfeature$output$Included[1239]
pgLFQfeature$"Normalized abundance"[1239,]

# tNpq
    par(mfrow=c(4,3), mar=c(1,1,4,1))
    for (i in 1:12)
        pgLFQtNpq(QuantitativeValue=pgLFQfeature$"Normalized abundance", 
        peptide=pgLFQfeature$peptideInfo$Sequence, 
        protein=pgLFQfeature$peptideInfo$Protein, N=i)

}

\references{

      \url{http://www.ncbi.nlm.nih.gov/pubmed/20576481}

}

\keyword{datasets}
