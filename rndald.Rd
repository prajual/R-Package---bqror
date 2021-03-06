% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ODRII.R
\name{rndald}
\alias{rndald}
\title{Generates random numbers from an Asymmetric Laplace Distribution}
\usage{
rndald(sigma, p, n)
}
\arguments{
\item{sigma}{scale factor, a scalar.}

\item{p}{quantile or skewness parameter, p in (0,1).}

\item{n}{number of observations}
}
\value{
Returns a vector \eqn{(n x 1)} of random numbers using an AL(0, \eqn{\sigma}, p)
}
\description{
This function generates a vector of random numbers from an asymmetric
Laplace distribution with quantile p.
}
\details{
Generates a vector of random numbers from an asymmetric Laplace distribution,
as a mixture of normal–exponential distributions.
}
\examples{
set.seed(101)
sigma <- 2.503306
p <- 0.25
n <- 1
ans <- rndald(sigma, p, n)

# ans
#    1.07328

}
\references{
Kozumi, H. and Kobayashi, G. (2011). “Gibbs Sampling Methods for Bayesian Quantile Regression.”
 Journal of Statistical Computation and Simulation, 81(11): 1565–1578.

 Koenker, R. and Machado, J. (1999). “Goodness of Fit and Related
 Inference Processes for Quantile Regression.”, Journal of
  American Statistics Association, 94(3): 1296-1309.

 Keming Yu and Jin Zhang (2005). “A Three-Parameter Asymmetric
 Laplace Distribution.” Communications in Statistics - Theory and Methods: 1867-1879.
}
\seealso{
asymmetric Laplace distribution
}
