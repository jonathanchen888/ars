\name{Implementing an Adaptive-Rejection Sampler}
\alias{Documentation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Implementing an Adaptive-Rejection Sampler
%%  Implementing an Adaptive-Rejection Sampler
}
\description{
%%  Implementing the Adaptive-Rejection Sampler based on Gilks & Wild (1992)
}
\usage{
Implementing an Adaptive-Rejection Sampler
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{density}{A probability density function
%%     ~~Describe \code{x} here~~
}\item{n}{Number of samples to generate
%%     ~~Describe \code{x} here~~
}
\item{lb}{The lower bound of the domain
%%     ~~Describe \code{x} here~~
}
\item{ub}{The upper bound of the domain
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{n samples from the probability density function
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Adaptive Rejection Sampling for Gibbs Sampling
W. R. Gilks and P. Wild (1992)
%% ~put references to the literature/web site here ~
}
\author{Rachel Hu, Lin Yang, Jonathan Chen
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{ars("dnorm(x,1,1)", 500, -1, 10)
ars("dgamma(x, 2, rate = 1)", 1000, 1, 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
