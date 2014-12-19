\name{findNN_}

\alias{NN_}
\alias{findNN_}

\title{find index of nearest neighbor}

\description{

    same as \code{findNN} but is uses \code{C++} STL \code{std::lower_bound} methode.
    this function is going to replace \code{findNN}.

}

\usage{

    findNN_(q, vec, check) 

}

\arguments{

  \item{q}{a double vector which can be considered as query objects.}

  \item{vec}{a sorted double vector which can be considered as a data base.}

  \item{check}{boolean enables test if \code{vec} is sorted. default is \code{FALSE}}

}

\details{

    please check \code{?findNN}


}

\author{

    CP

}

\examples{
    # compare ANSI-C binary search with C++ std::lower_bound
    DB<-c(rep(1.0, 3), rep(2.0, 3))
    q<-c(-1, 1.0, 1.01, 1.5, 1.9)
    abs(DB[findNN(q, DB)] - q)
    abs(DB[findNN_(q, DB)] - q)


    DB<-sort(rnorm(100, mean=100, sd=10))
    # should be 0
    unique(DB[findNN_(DB,DB)] - DB)

    q<-rnorm(100, mean=100)

    idx.NN<-findNN_(q,DB)
    hist(DB[findNN_(q,DB)] - q)

    # definition of findNN_ holds
    i<-1:5
    findNN_(3.5, i)

    i<-1:6
    findNN_(3.5, i)

}
