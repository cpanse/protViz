\name{findNN}

\alias{NN}
\alias{findNN}
\alias{findNN_}
\alias{lower_bound_}

\title{find index of nearest neighbor}

\description{

    Given a vector of sorted double values \code{vec} of size
    \code{n} and a vector of \code{m} query objects \code{q}. 

    \code{findNN} determines for each element \code{q[i]} in \code{q} 
    the nearest neighbor index \code{o} so that the following remains true:

    there is no element \code{k} with \code{1} \eqn{\le} \code{k} 
    \eqn{\le} \code{n} and \code{k} is not \code{o} so that 

    \code{abs(vec[k] - q[i])} < \code{abs(vec[o] - q[i])}. 


}

\usage{

    findNN(q, vec, check) 
    findNN_(q, vec, check) 

}

\arguments{

  \item{q}{a double vector which can be considered as query objects.}

  \item{vec}{a sorted double vector which can be considered as a data base.}

  \item{check}{boolean enables test if \code{vec} is sorted. default is
  \code{FALSE}.}
}

\details{

    The internal algorithm of \code{findNN} is implemented as binary search.
    \code{findNN} has \eqn{O(m * log_2(n))} time complexity where n is defined as 
    \code{length(vec)} and m is \code{length(m)}. 

  \code{findNN} is implemented using C library function - bsearch(),
  while \code{findNN_} uses C++11 STL function lower_bound().
}

\seealso{
	\code{\link{findInterval}}
}

\references{
\itemize{
\item{\url{http://www.cplusplus.com/reference/algorithm/lower_bound/}}
}
}

\author{

   Christian Panse 2007, 2008, 2009, 2010, 2012 , 2015
   based on the C++ STL lower_bound method.

}

\examples{

    (NNidx <- findNN(q <- c(1, 1.0001, 1.24, 1.26), DB <- seq(1, 5 , by = 0.25)))
    (NNidx == c(1, 1, 2, 2))

    DB <- sort(rnorm(100, mean=100, sd = 10))

    # should be 0
    unique(DB[findNN(DB,DB)] - DB)

    q <- rnorm(100, mean=100)

    idx.NN <- findNN(q,DB)
    hist(DB[findNN(q,DB)] - q)

    # definition of findNN holds
    i <- 1:5
    findNN(3.5, i)

    i <- 1:6
    findNN(3.5, i)
    
     # compare ANSI-C binary search with C++ std::lower_bound
    DB <- c(rep(1.0, 3), rep(2.0, 3))
    q <- c(-1, 1.0, 1.01, 1.5, 1.9)
    abs(DB[findNN(q, DB)] - q)
    abs(DB[findNN_(q, DB)] - q)


    DB <- sort(rnorm(100, mean=100, sd=10))
    # should be 0
    unique(DB[findNN_(DB,DB)] - DB)

    q <- rnorm(100, mean=100)

    idx.NN <- findNN_(q, DB)
    hist(DB[findNN_(q, DB)] - q)

    # definition of findNN_ holds
    i <- 1:5
    findNN_(3.5, i)

    i <- 1:6
    findNN_(3.5, i)
}
