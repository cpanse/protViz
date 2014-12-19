\name{aa2mass}

\alias{aa2mass}

\title{determine the weight if a fiven amino acid sequence}

\description{

    The function determines the weight of each amino acid fiven sequence.

}

\usage{

    aa2mass(peptideSequence, mass=AA$Monoisotopic, letter1=AA$letter1)
}

\arguments{

  \item{peptideSequence}{a double vector which can be considered as query objects.}

  \item{mass}{A vector of size 20 containing the weight of the AA.}

  \item{letter1}{AA letter 1 code in the same size and order as the \code{mass} attribute.}

}

\details{

    For the computation no C-Term and N-Term is considered.
    \code{aa2mass} is useful if you have AA modifications.

}

\author{

   Christian Panse 2014

}

\examples{

    peptides<-c('HTLNQIDSVK', 'ALGGEDVR', 'TPIVGQPSIPGGPVR')

    C_term <- 17.002740
    N_term <- 1.007825
    H_ <- 1.008

    unlist(lapply(aa2mass(peptides), sum)) + C_term + N_term + H_ - parentIonMass(peptides)

    # determine the fragment ions
    lapply(aa2mass(peptides), function(x){fragmentIon(x)[[1]]})

    # an example with [STY] AA modification

    peptide<-'HTLNQIDSVK'
    mod <- rep(0.0, nchar(peptide)); mod[8] <- 79.998;
    
    aa2mass(peptide)[[1]] + mod

}
