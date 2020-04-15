\name{cenroid}
\alias{cenroid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
centroid a profile spectrum     
}
\description{
derives a centroid spectrum of a recorded profile mode spectrum.
}
\usage{
cenroid(mZ, intensity)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mZ}{Numerical vector of profile recorded data and sorted mZ values.}
  \item{intensity}{corresponding intensity values.}
  \item{debug}{if true all peak profiles are plotted.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns a \code{data.frame} with a mZ and a intensity column.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Christian Panse and Jonas Grossmann
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \href{https://github.com/lgatto/MSnbase/blob/2b6d6f5162e7464c39b48d0daa4d981c6ec30bbd/R/functions-Spectrum.R#L649}{MSnbase functions-Spectrum.R}
}~
}
\examples{
    p <- .getProfileMS2()
    op <- par(mfrow=c(2,1))
    plot(p$mZ, p$intensity, type='h', main='profile', xlim=c(100,300))
    plot(centroid(p$mZ, p$intensity),type='h', main="centroid",xlim=c(100,300))
    
    
    op <- par(ask = TRUE)
        rv <- centroid(p$mZ, p$intensity, debug = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
