\docType{data}
\name{pgLFQprot}

\alias{pgLFQprot}

\title{pgLFQprot - A data set with a featuremap export ProgenesisLCMS}

\format{A data object from ProgenesisImporter.R}

\description{
This data set contain the top 250 most abundant Proteins identified from 24 LC-MSMS runs measured on 
an Orbitrap-Velos. The original data contains digested proteins from HeLa cells infected with Shigella bacteria grown in a time course.
It is structured in a way, that 6 biological replicates from 4 conditions are measured (Not_infected, Infected_1hr, Infected_2hr, Infected_3hr). 
The 24 LC-MSMS runs are aligned with each other using ProgenesisLCMS. "Normalized volumes" on MS1 are extracted for all features in the respective LC-MSMS runs. The features with an annotation such as peptide sequences and protein accessions (using Mascot search algorithm) are assembled to proteins. Features volumes are stacked (for non conflicting features) to generate the quantitative protein value for each LC-MSMS run.
The csv file is an exported Protein data map.


We realized, depending on your language and keyboard setting, the separator for the FeatureData as well as for ProteinMeasurements are different (semicolon and commas are used depending on your setting). We assume, that semicolons are and the individual cells are escaped by ". If this differes, we have an option that can be switched in the pgImporter function.


More information on the commercial software can be found here:
\url{http://www.nonlinear.com/products/progenesis/lc-ms/overview/}.

}
\author{

   Christian Panse, Jonas Grossmann 2012

}

\examples{

data(pgLFQprot)
op<-par(mfrow=c(1,1),mar=c(18,18,4,1),cex=0.5)
samples<-names(pgLFQprot$"Normalized abundance")
image(cor(asinh(pgLFQprot$"Normalized abundance")),
    main='pgLFQprot correlation',
    axes=FALSE,
    col=gray(seq(0,1,length=20)))
axis(1,at=seq(from=0, to=1, 
    length.out=length(samples)), labels=samples, las=2)
axis(2,at=seq(from=0, to=1, 
    length.out=length(samples)), labels=samples, las=2)
par(op)

}

\references{

      \url{http://www.ncbi.nlm.nih.gov/pubmed/20576481}

}

\keyword{datasets}
