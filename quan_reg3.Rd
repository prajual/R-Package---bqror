% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ODRII.R
\name{quan_reg3}
\alias{quan_reg3}
\title{Bayesian Quantile Regression for Ordinal Models
with 3 outcomes}
\usage{
quan_reg3(y, x, mc = 15000, p)
}
\arguments{
\item{y}{dependent variable i.e. ordinal outcome values.}

\item{x}{covariate matrix of dimension \eqn{(n x k)} including a column of ones.}

\item{mc}{number of MCMC iterations, post burn-in.}

\item{p}{quantile level or skewness parameter, p in (0,1).}
}
\value{
Returns a list with components
\itemize{
\item{\code{post_mean_beta}: }{a vector with mean of sampled
 \eqn{\beta} for each covariate.}
\item{\code{post_mean_sigma}: }{a vector with mean of sampled
 \eqn{\sigma}.}
\item{\code{post_std_beta}: }{a vector with standard deviation
 of sampled \eqn{\beta} for each covariate.}
 \item{\code{post_std_sigma}: }{a vector with standard deviation
 of sampled \eqn{\sigma}.}
 \item{\code{DIC_result}: }{results of the DIC criteria.}
 \item{\code{beta_draws}: }{a matrix with all sampled values for \eqn{\beta}.}
 \item{\code{sigma_draws}: }{a matrix with all sampled values for \eqn{\sigma}.}
}
}
\description{
This function estimates the Bayesian Quantile Regression for ordinal model with
3 outcomes and reports the posterior mean and posterior standard deviations
of \eqn{(\beta, \sigma)}.
}
\details{
Function implements the Bayesian quantile regression for
ordinal models with 3 outcomes using a Gibbs sampling
procedure.

Function initializes prior and then iteratively
samples \eqn{\beta}, \eqn{\delta} and latent variable z.
Burn-in is taken as \eqn{0.25*mc} and \eqn{iter = burn}-\eqn{in + mc}.
}
\examples{
set.seed(101)
data("data25j3")
x <- data25j3$x
y <- data25j3$y
p <- 0.25
ans <- quan_reg3(y, x, mc = 50, p)

# ans$post_mean_beta
#     1.7201671 1.9562172 0.8334668
# ans$post_std_beta
#     0.2400355 0.2845326 0.2036498
# ans$post_mean_sigma
#     0.9684741
# ans$post_std_sigma
#     0.1962351
# ans$Dic_Result
# dic
#    474.4673
# pd
#    5.424001
# devpostmean
#    463.6193
# ans$beta_draws
#     0.0000000 0.000000 0.0000000
#    -3.6740670 1.499495 1.3610085
#    -1.1006076 2.410271 1.3379175
#    -0.5310387 1.604194 0.7830659
#     0.4870828 1.761879 0.6921727
#     0.9481320 1.485709 1.0251322... soon
# ans$sigma_draws
#     2.0000000
#     3.6987793
#     3.2785105
#     2.9769533
#     2.9273486
#     2.5807661
#     2.2654222... soon

}
\references{
Rahman, M. A. (2016). “Bayesian
Quantile Regression for Ordinal Models.”
Bayesian Analysis, 11(1): 1-24.

Yu, K. and Moyeed, R. A. (2001). “Bayesian Quantile Regression.” Statistics and
Probability Letters, 54(4): 437–447.

Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
The American Statistician, 46(3): 167-174.

Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
Gibbs Distributions, and the Bayesian Restoration of Images."
IEEE Transactions an Pattern Analysis and Machine Intelligence,
6(6): 721-741.
}
\seealso{
tcltk, \link[stats]{rnorm}, \link[stats]{qnorm},
\link[MASS]{ginv}, Gibbs sampling
}
