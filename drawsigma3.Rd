% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ODRII.R
\name{drawsigma3}
\alias{drawsigma3}
\title{Samples the \eqn{\sigma} for an Ordinal Model
with 3 outcomes}
\usage{
drawsigma3(z, x, beta, nu, tau2, theta, n0, d0)
}
\arguments{
\item{z}{Gibbs draw of latent response variable, a column vector.}

\item{x}{covariate matrix of dimension \eqn{(n x k)} including a column of ones.}

\item{beta}{Gibbs draw of coeffcients of dimension \eqn{(k x 1)}.}

\item{nu}{modified scale factor, row vector.}

\item{tau2}{2/(p(1-p)).}

\item{theta}{(1-2p)/(p(1-p)).}

\item{n0}{prior hyper-parameter for \eqn{\sigma}.}

\item{d0}{prior hyper-parameter for \eqn{\sigma}.}
}
\value{
Returns a column vector of the \eqn{\sigma}
from an inverse gamma distribution.
}
\description{
This function samples the \eqn{\sigma} from an inverse-gamma distribution
for an ordinal model with 3 outcomes.
}
\details{
Function samples the \eqn{\sigma} from an inverse
gamma distribution.
}
\examples{
set.seed(101)
z <- c(21.01744, 33.54702, 33.09195, -3.677646,
 21.06553, 1.490476, 0.9618205, -6.743081, 21.02186, 0.6950479)
x <- matrix(c(
     1, -0.3010490, 0.8012506,
     1,  1.2764036, 0.4658184,
     1,  0.6595495, 1.7563655,
     1, -1.5024607, -0.8251381,
     1, -0.9733585, 0.2980610,
     1, -0.2869895, -1.0130274,
     1,  0.3101613, -1.6260663,
     1, -0.7736152, -1.4987616,
     1,  0.9961420, 1.2965952,
     1, -1.1372480, 1.7537353),
     nrow = 10, ncol = 3, byrow = TRUE)
beta <- c(-0.74441, 1.364846, 0.7159231)
nu <- c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5)
tau2 <- 10.6667
theta <- 2.6667
n0 <- 5
d0 <- 8
ans <- drawsigma3(z, x, beta, nu, tau2, theta, n0, d0)

# ans
#    3.749524

}
\references{
Albert, J. and Chib, S. (1993). “Bayesian Analysis of Binary and Polychotomous
Response Data.” Journal of the American Statistical Association, 88(422): 669–679.

Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
The American Statistician, 46(3): 167-174.

Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
Gibbs Distributions, and the Bayesian Restoration of Images."
IEEE Transactions an Pattern Analysis and Machine Intelligence,
6(6): 721-741.
}
\seealso{
\link[invgamma]{rinvgamma}, Gibbs sampling
}
