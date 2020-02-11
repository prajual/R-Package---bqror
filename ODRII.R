#' Bayesian Quantile Regression for Ordinal Models
#' with 3 outcomes
#'
#' This function estimates the Bayesian Quantile Regression for ordinal model with
#' 3 outcomes and reports the posterior mean and posterior standard deviations
#' of \eqn{(\beta, \sigma)}.
#'
#' @param y         dependent variable i.e. ordinal outcome values.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param mc        number of MCMC iterations, post burn-in.
#' @param p         quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' Function implements the Bayesian quantile regression for
#' ordinal models with 3 outcomes using a Gibbs sampling
#' procedure.
#'
#' Function initializes prior and then iteratively
#' samples \eqn{\beta}, \eqn{\delta} and latent variable z.
#' Burn-in is taken as \eqn{0.25*mc} and \eqn{iter = burn}-\eqn{in + mc}.
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{post_mean_beta}: }{a vector with mean of sampled
#'  \eqn{\beta} for each covariate.}
#' \item{\code{post_mean_sigma}: }{a vector with mean of sampled
#'  \eqn{\sigma}.}
#' \item{\code{post_std_beta}: }{a vector with standard deviation
#'  of sampled \eqn{\beta} for each covariate.}
#'  \item{\code{post_std_sigma}: }{a vector with standard deviation
#'  of sampled \eqn{\sigma}.}
#'  \item{\code{DIC_result}: }{results of the DIC criteria.}
#'  \item{\code{beta_draws}: }{a matrix with all sampled values for \eqn{\beta}.}
#'  \item{\code{sigma_draws}: }{a matrix with all sampled values for \eqn{\sigma}.}
#' }
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' Yu, K. and Moyeed, R. A. (2001). “Bayesian Quantile Regression.” Statistics and
#' Probability Letters, 54(4): 437–447.
#'
#' Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
#' The American Statistician, 46(3): 167-174.
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images."
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741.
#'
#' @importFrom "MASS" "ginv"
#' @importFrom "tcltk" "tkProgressBar" "setTkProgressBar"
#' @importFrom "stats" "sd"
#' @seealso tcltk, \link[stats]{rnorm}, \link[stats]{qnorm},
#' \link[MASS]{ginv}, Gibbs sampling
#' @examples
#' set.seed(101)
#' data("data25j3")
#' x <- data25j3$x
#' y <- data25j3$y
#' p <- 0.25
#' ans <- quan_reg3(y, x, mc = 50, p)
#'
#' # ans$post_mean_beta
#' #     1.7201671 1.9562172 0.8334668
#' # ans$post_std_beta
#' #     0.2400355 0.2845326 0.2036498
#' # ans$post_mean_sigma
#' #     0.9684741
#' # ans$post_std_sigma
#' #     0.1962351
#' # ans$Dic_Result
#' # dic
#' #    474.4673
#' # pd
#' #    5.424001
#' # devpostmean
#' #    463.6193
#' # ans$beta_draws
#' #     0.0000000 0.000000 0.0000000
#' #    -3.6740670 1.499495 1.3610085
#' #    -1.1006076 2.410271 1.3379175
#' #    -0.5310387 1.604194 0.7830659
#' #     0.4870828 1.761879 0.6921727
#' #     0.9481320 1.485709 1.0251322... soon
#' # ans$sigma_draws
#' #     2.0000000
#' #     3.6987793
#' #     3.2785105
#' #     2.9769533
#' #     2.9273486
#' #     2.5807661
#' #     2.2654222... soon
#'
#' @export
quan_reg3 <- function(y, x, mc = 15000, p) {
    if ( dim(y)[2] != 1){
        stop("parameter y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be a integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( length(mc) != 1){
        stop("parameter mc must be scalar")
    }
    if ( !is.numeric(mc)){
        stop("parameter mc must be a numeric")
    }
    if ( length(p) != 1){
        stop("parameter p must be scalar")
    }
    if ( any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    J <- dim(as.array(unique(y)))[1]
    if ( J > 3 ){
        warning("This function is for 3 outcome
                variables. We are switching to quan_regg3")
        ans <- quan_regg3(y, x, mc = mc, p, tune = 0.1)
    }
    x <- as.matrix(x)
    y <- as.matrix(y)
    n <- dim(x)[1]
    k <- dim(x)[2]
    burn <- 0.2 * mc
    iter <- burn + mc

    b0 <- array(c(0, 0, 0), dim = c(1, 3))
    B0 <- diag(k)
    invb0 <- ginv(B0)
    invb0b0 <- invb0 %*% (t(b0))
    n0 <- 5
    d0 <- 8
    lengthp <- length(p)

    beta <- array(0, dim = c(iter, k))
    sigma <- array(0, dim = c(iter, 1))
    nu <- array(0, dim = c(iter, n))
    beta_draws <- array(0, dim = c(iter, k * lengthp))
    sigma_draws <- array(0, dim = c(iter, lengthp))

    beta[1, ] <- array(c(0, 0, 0), dim = c(1, 3))
    sigma[1] <- 2
    nu[1, ] <- 5 * rep(1, n)

    gammacp <- array(c(-Inf, 0, 4, Inf), dim = c(1, 4))
    lambda <- 0.5
    theta <- (1 - 2 * p) / (p * (1 - p))
    tau <- sqrt(2 / (p * (1 - p)))
    tau2 <- tau ^ 2
    total <- iter
    pb <- tkProgressBar(title = "Simulation in Progress",
                        min = 0, max = total, width = 300)
    for (i in 2:iter) {
        z <- drawlatent3(y, x, beta[(i - 1), ],
                         sigma[(i - 1)], nu[(i - 1), ],
                         theta, tau2, gammacp)
        p1 <- drawbeta3(z, x, sigma[(i - 1)],
                        nu[(i - 1), ], tau2,
                        theta, invb0, invb0b0)
        beta[i, ] <- p1
        p2 <- drawsigma3(z, x, beta[i, ],
                         nu[(i - 1), ], tau2,
                         theta, n0, d0)
        sigma[i] <- p2
        p3 <- drawnu3(z, x, beta[i, ],
                      sigma[i], tau2,
                      theta, lambda)
        nu[i, ] <- p3
        setTkProgressBar(pb, i,
                         label = paste( round( (i / iter) * 100, 0),
                                        "% done"))
    }
    close(pb)
    beta_draws <- beta
    sigma_draws <- sigma
    post_mean_beta <- colMeans(beta_draws[(burn + 1):iter, ])
    post_std_beta <- apply(beta_draws[(burn + 1):iter, ], 2, sd)
    post_mean_sigma <- mean(sigma_draws[(burn + 1):iter, ])
    post_std_sigma <- std(sigma_draws[(burn + 1):iter, ])
    dic_result <- deviance3(y, x, gammacp, p,
                            post_mean_beta, post_std_beta,
                            post_mean_sigma, post_std_sigma,
                            beta_draws, sigma_draws, burn, iter)
    result <- list("post_mean_beta" = post_mean_beta,
                   "post_std_beta" = post_std_beta,
                   "post_mean_sigma" = post_mean_sigma,
                   "post_std_sigma" = post_std_sigma,
                   "dic_result" = dic_result,
                   "beta_draws" = beta_draws,
                   "sigma_draws" = sigma_draws)
    return(result)
}
#' Samples the Latent Variable z for an Ordinal Model
#' with 3 outcomes
#'
#' This function samples the latent variable z from a truncated
#' normal distribution for an ordinal model with 3 outcomes.
#'
#' @param y         dependent variable i.e. ordinal outcome values.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param beta      column vector of coeffcients of dimension \eqn{(k x 1)}.
#' @param sigma     scale factor, a scalar value.
#' @param nu        modified scale factor, row vector.
#' @param tau2      2/(p(1-p)).
#' @param theta     (1-2p)/(p(1-p)).
#' @param gammacp   row vector of cutpoints including -Inf and Inf.
#'
#' @details
#' Function samples the latent variable z from a truncated normal
#' distribution.
#'
#' @return Returns a column vector of values for latent variable z.
#'
#' @references Albert, J. and Chib, S. (1993). “Bayesian Analysis of Binary and Polychotomous
#' Response Data.” Journal of the American Statistical Association, 88(422): 669–679.
#'
#' Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
#' The American Statistician, 46(3): 167-174.
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images."
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741.
#'
#' @seealso Gibbs sampling, truncated normal distribution,
#' \link[truncnorm]{rtruncnorm}
#' @importFrom "truncnorm" "rtruncnorm"
#' @examples
#' set.seed(101)
#' data("data25j3")
#' x <- data25j3$x
#' y <- data25j3$y
#' beta <- c(1.7201671, 1.9562172, 0.8334668)
#' sigma <- 0.9684741
#' nu <- c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5)
#' theta <- 2.6667
#' tau2 <- 10.6667
#' gammacp <- c(-Inf, 0, 4, Inf)
#' ans <- drawlatent3(y, x, beta, sigma, nu,
#' theta, tau2, gammacp)
#'
#' # ans
#' #   12.79298 20.40747 1.557821
#' #   26.07846 17.41031 12.86016
#' #   3.364703 21.61075 2.666627 .. soon
#'
#' @export
drawlatent3 <- function(y, x, beta, sigma, nu, theta, tau2, gammacp) {
    if ( dim(y)[2] != 1){
        stop("parameter y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be a integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( length(sigma) != 1){
        stop("parameter sigma must be scalar")
    }
    if ( !all(is.numeric(nu))){
        stop("each entry in nu must be numeric")
    }
    if ( length(tau2) != 1){
        stop("parameter tau2 must be scalar")
    }
    if ( !all(is.numeric(tau2))){
        stop("each entry in tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("each entry in theta must be numeric")
    }
    n <- dim(y)[1]
    z <- array(0, dim = c(1, n))
    for (i in 1:n) {
        meanp <- (x[i, ] %*% (beta)) + (theta * nu[i])
        std <- sqrt(tau2 * sigma * nu[i])
        temp <- y[i]
        a <- gammacp[temp]
        b <- gammacp[temp + 1]
        z[1, i] <- rtruncnorm(n = 1, a = a, b = b,
                              mean = meanp, sd = std)
    }
    return(z)
}

#' Samples \eqn{\beta} for an Ordinal Model
#' with 3 outcomes
#'
#' This function samples \eqn{\beta} from its conditional
#' posterior distribution for an ordinal model with 3
#' outcomes.
#'
#' @param z         dependent variable i.e. ordinal outcome values.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param sigma     scale factor, a scalar value.
#' @param nu        modified scale factor, row vector.
#' @param tau2      2/(p(1-p)).
#' @param theta     (1-2p)/(p(1-p)).
#' @param invB0     inverse of prior covariance matrix of normal distribution.
#' @param invB0b0   prior mean pre-multiplied by invB0.
#'
#' @details
#' Function samples a vector of \eqn{\beta} from a multivariate normal distribution.
#'
#' @return Returns a column vector of \eqn{\beta}
#' from a multivariate normal distribution.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988).
#' “The New S Language. Wadsworth & Brooks/Cole.”
#'
#' Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
#' The American Statistician, 46(3): 167-174.
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images."
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741.
#'
#' @importFrom "MASS" "mvrnorm" "ginv"
#' @seealso Gibbs sampling, normal distribution
#' , \link[GIGrvg]{rgig}
#' @examples
#' set.seed(101)
#' z <- c(21.01744, 33.54702, 33.09195, -3.677646,
#'  21.06553, 1.490476, 0.9618205, -6.743081, 21.02186, 0.6950479)
#' x <- matrix(c(
#'      1, -0.3010490, 0.8012506,
#'      1,  1.2764036, 0.4658184,
#'      1,  0.6595495, 1.7563655,
#'      1, -1.5024607, -0.8251381,
#'      1, -0.9733585, 0.2980610,
#'      1, -0.2869895, -1.0130274,
#'      1,  0.3101613, -1.6260663,
#'      1, -0.7736152, -1.4987616,
#'      1,  0.9961420, 1.2965952,
#'      1, -1.1372480, 1.7537353),
#'      nrow = 10, ncol = 3, byrow = TRUE)
#' sigma <- 1.809417
#' nu <- c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5)
#' tau2 <- 10.6667
#' theta <- 2.6667
#' invB0 <- matrix(c(
#'      1, 0, 0,
#'      0, 1, 0,
#'      0, 0, 1),
#'      nrow = 3, ncol = 3, byrow = TRUE)
#' invB0b0 <- c(0, 0, 0)
#'
#' ans <- drawbeta3(z, x, sigma, nu, tau2, theta, invB0, invB0b0)
#'
#' # ans
#' #    -0.74441 1.364846 0.7159231
#'
#' @export
drawbeta3 <- function(z, x, sigma, nu, tau2, theta, invB0, invB0b0) {
    if ( !all(is.numeric(z))){
        stop("each entry in z must be numeric")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( length(sigma) != 1){
        stop("parameter sigma must be scalar")
    }
    if ( !all(is.numeric(nu))){
        stop("each entry in nu must be numeric")
    }
    if ( length(tau2) != 1){
        stop("parameter tau2 must be scalar")
    }
    if ( !all(is.numeric(tau2))){
        stop("each entry in tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("each entry in theta must be numeric")
    }
    if ( !all(is.numeric(invB0))){
        stop("each entry in invB0 must be numeric")
    }
    if ( !all(is.numeric(invB0b0))){
        stop("each entry in invB0b0 must be numeric")
    }
    n <- dim(x)[1]
    k <- dim(x)[2]
    meancomp <- array(0, dim = c(n, k))
    varcomp <- array(0, dim = c(k, k, n))
    q <- array(0, dim = c(1, k))
    eye <- diag(k)
    for (j in 1:n) {
        meancomp[j, ] <- (x[j, ] *
                              (z[j] - (theta * nu[j]))) / (tau2 * sigma * nu[j])
        varcomp[, , j] <- ( (x[j, ]) %*% (t(x[j, ]))) / (tau2 * sigma * nu[j])
    }
    Btilde <- ginv(invB0 + rowSums(varcomp, dims = 2))
    btilde <- Btilde %*% (invB0b0 + colSums(meancomp))
    L <- t(chol(Btilde))
    beta <- t(btilde) + t(L %*% ( (mvrnorm(n = 1, mu = q, Sigma = eye))))
    return(beta)
}
#' Samples the \eqn{\sigma} for an Ordinal Model
#' with 3 outcomes
#'
#' This function samples the \eqn{\sigma} from an inverse-gamma distribution
#' for an ordinal model with 3 outcomes.
#'
#' @param z         Gibbs draw of latent response variable, a column vector.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param beta      Gibbs draw of coeffcients of dimension \eqn{(k x 1)}.
#' @param nu        modified scale factor, row vector.
#' @param tau2      2/(p(1-p)).
#' @param theta     (1-2p)/(p(1-p)).
#' @param n0        prior hyper-parameter for \eqn{\sigma}.
#' @param d0        prior hyper-parameter for \eqn{\sigma}.
#'
#' @details
#' Function samples the \eqn{\sigma} from an inverse
#' gamma distribution.
#'
#' @return Returns a column vector of the \eqn{\sigma}
#' from an inverse gamma distribution.
#'
#' @importFrom "invgamma" "rinvgamma"
#'
#' @references Albert, J. and Chib, S. (1993). “Bayesian Analysis of Binary and Polychotomous
#' Response Data.” Journal of the American Statistical Association, 88(422): 669–679.
#'
#' Casella, G., George E. I. (1992). “Explaining the Gibbs Sampler."
#' The American Statistician, 46(3): 167-174.
#'
#' Geman, S., and Geman, D. (1984). “Stochastic Relaxation,
#' Gibbs Distributions, and the Bayesian Restoration of Images."
#' IEEE Transactions an Pattern Analysis and Machine Intelligence,
#' 6(6): 721-741.
#'
#' @seealso \link[invgamma]{rinvgamma}, Gibbs sampling
#' @examples
#' set.seed(101)
#' z <- c(21.01744, 33.54702, 33.09195, -3.677646,
#'  21.06553, 1.490476, 0.9618205, -6.743081, 21.02186, 0.6950479)
#' x <- matrix(c(
#'      1, -0.3010490, 0.8012506,
#'      1,  1.2764036, 0.4658184,
#'      1,  0.6595495, 1.7563655,
#'      1, -1.5024607, -0.8251381,
#'      1, -0.9733585, 0.2980610,
#'      1, -0.2869895, -1.0130274,
#'      1,  0.3101613, -1.6260663,
#'      1, -0.7736152, -1.4987616,
#'      1,  0.9961420, 1.2965952,
#'      1, -1.1372480, 1.7537353),
#'      nrow = 10, ncol = 3, byrow = TRUE)
#' beta <- c(-0.74441, 1.364846, 0.7159231)
#' nu <- c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5)
#' tau2 <- 10.6667
#' theta <- 2.6667
#' n0 <- 5
#' d0 <- 8
#' ans <- drawsigma3(z, x, beta, nu, tau2, theta, n0, d0)
#'
#' # ans
#' #    3.749524
#'
#' @export
drawsigma3 <- function(z, x, beta, nu, tau2, theta, n0, d0) {
    if ( !all(is.numeric(z))){
        stop("each entry in z must be numeric")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( !all(is.numeric(nu))){
        stop("each entry in nu must be numeric")
    }
    if ( length(tau2) != 1){
        stop("parameter tau2 must be scalar")
    }
    if ( !all(is.numeric(tau2))){
        stop("each entry in tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("each entry in theta must be numeric")
    }
    if ( length(n0) != 1){
        stop("parameter n0 must be scalar")
    }
    if ( !all(is.numeric(n0))){
        stop("each entry in n0 must be numeric")
    }
    if ( length(d0) != 1){
        stop("parameter d0 must be scalar")
    }
    if ( !all(is.numeric(d0))){
        stop("each entry in d0 must be numeric")
    }
    n <- dim(x)[1]
    ntilde <- n0 + (3 * n)
    temp <- array(0, dim = c(n, 1))
    for (j in 1:n) {
        temp[j, 1] <- ( (z[j] - x[j, ] %*% ( (beta)) -
                             (theta * nu[j])) ^ 2) / (tau2 * nu[j])
    }
    dtilde <- sum(temp) + d0 + (2 * sum(nu))
    sigma <- rinvgamma(n = 1, shape = (ntilde / 2), scale = (2 / dtilde))
    return(sigma)
}

#' Samples the scale factor \eqn{\nu} for an Ordinal Model
#' with 3 outcomes
#'
#' This function samples the \eqn{\nu} from a generalized inverse Gaussian (GIG)
#' distribution for an ordinal model with 3 outcomes.
#'
#' @param z         Gibbs draw of latent response variable, a column vector.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param beta      Gibbs draw of coefficients of dimension \eqn{(k x 1)}.
#' @param sigma     scale factor, a scalar.
#' @param tau2      2/(p(1-p)).
#' @param theta     (1-2p)/(p(1-p)).
#' @param lambda    index parameter of GIG distribution which is equal to 0.5
#'
#' @details
#' Function samples the \eqn{\nu} from a GIG
#' distribution.
#'
#' @return Returns a row vector of the \eqn{\nu}
#' from GIG distribution.
#'
#' @references  Rahman, M. A. (2016), “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1), 1-24.
#'
#' Dagpunar, J. S. (1989). "An Easily Implemented Generalised Inverse Gaussian Generator."
#' Communication Statistics Simulation, 18: 703-710.
#'
#' @importFrom "GIGrvg" "rgig"
#' @seealso GIGrvg, Gibbs sampling, \link[GIGrvg]{rgig}
#' @examples
#' set.seed(101)
#' z <- c(21.01744, 33.54702, 33.09195, -3.677646,
#'  21.06553, 1.490476, 0.9618205, -6.743081, 21.02186, 0.6950479)
#' x <- matrix(c(
#'      1, -0.3010490, 0.8012506,
#'      1,  1.2764036, 0.4658184,
#'      1,  0.6595495, 1.7563655,
#'      1, -1.5024607, -0.8251381,
#'      1, -0.9733585, 0.2980610,
#'      1, -0.2869895, -1.0130274,
#'      1,  0.3101613, -1.6260663,
#'      1, -0.7736152, -1.4987616,
#'      1, 0.9961420, 1.2965952,
#'      1, -1.1372480, 1.7537353),
#'      nrow = 10, ncol = 3, byrow = TRUE)
#' beta <- c(-0.74441, 1.364846, 0.7159231)
#' sigma <- 3.749524
#' tau2 <- 10.6667
#' theta <- 2.6667
#' lambda <- 0.5
#' ans <- drawnu3(z, x, beta, sigma, tau2, theta, lambda)
#'
#' # ans
#' #    5.177456 4.042261 8.950365
#' #    1.578122 6.968687 1.031987
#' #    4.13306 0.4681557 5.109653
#' #    0.1725333
#'
#' @export
drawnu3 <- function(z, x, beta, sigma, tau2, theta, lambda) {
    if ( !all(is.numeric(z))){
        stop("each entry in z must be numeric")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( length(sigma) != 1){
        stop("parameter sigma must be scalar")
    }
    if ( length(tau2) != 1){
        stop("parameter tau2 must be scalar")
    }
    if ( !all(is.numeric(tau2))){
        stop("each entry in tau2 must be numeric")
    }
    if ( length(theta) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(theta))){
        stop("each entry in theta must be numeric")
    }
    if ( length(lambda) != 1){
        stop("parameter theta must be scalar")
    }
    if ( !all(is.numeric(lambda))){
        stop("each entry in theta must be numeric")
    }
    n <- dim(x)[1]
    tildegamma2 <- ( (theta ^ 2) / (tau2 * sigma)) + (2 / sigma)
    tildedelta2 <- array(0, dim = c(n, 1))
    nu <- array(0, dim = c(1, n))
    for (j in 1:n) {
        tildedelta2[j, 1] <- ( (z[j] -
                                    (x[j, ] %*% ( (beta)))) ^ 2) / (tau2 * sigma)
        nu[1, j] <- rgig(lambda = lambda,
                         psi = tildegamma2,
                         chi = tildedelta2[j, 1],
                         n = 1)
    }
    return(nu)
}

#' Deviance Information Criteria for Ordinal Models
#' with 3 outcomes
#'
#' Function for computing the Deviance Information Criteria for ordinal
#' models with 3 outcomes.
#'
#' @param y                dependent variable i.e. ordinal outcome values.
#' @param x                covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param gammacp          row vector of cutpoints including -Inf and Inf.
#' @param p                quantile level or skewness parameter, p in (0,1).
#' @param post_mean_beta   mean value of \eqn{\beta} obtained from MCMC draws.
#' @param post_std_beta    standard deviation of \eqn{\beta} obtained from MCMC draws.
#' @param post_mean_sigma  mean value of \eqn{\sigma} obtained from MCMC draws.
#' @param post_std_sigma   standard deviation of \eqn{\sigma} obtained from MCMC draws.
#' @param beta_draws       MCMC draw of coeffcients, dimension is \eqn{(k x iter)}.
#' @param sigma_draws      MCMC draw of scale factor, dimension is \eqn{(iter x 1)}.
#' @param burn             number of discarded MCMC iterations.
#' @param iter             total number of MCMC iterations including the burn-in.
#'
#' @details
#' The Deviance is -2*(log likelihood) and has an important role in
#' statistical model comparision because of its relation with Kullback-Leibler
#' information criteria.
#'
#' @return Returns a list with components
#' \deqn{DIC = 2*avgdeviance - devpostmean}
#' \deqn{pd = avgdeviance - devpostmean}
#' \deqn{devpostmean = -2*(logLikelihood)}.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' Spiegelhalter, D. J., Best, N. G., Carlin B. P. and Linde A. (2002).
#' “Bayesian Measures of Model Complexity and Fit.” Journal of the
#' Royal Statistical Society B, Part 4: 583-639.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., and Rubin, D. B.
#' “Bayesian Data Analysis.” 2nd Edition, Chapman and Hall.
#'
#' @seealso  decision criteria
#' @examples
#' set.seed(101)
#' data("data25j3")
#' x <- data25j3$x
#' y <- data25j3$y
#' p <- 0.25
#' ans <- quan_reg3(y, x, mc = 50, p)
#' gammacp <- c(-Inf, 0, 4, Inf)
#' p <- 0.25
#' post_mean_beta <- ans$post_mean_beta
#' post_std_beta <- ans$post_std_beta
#' post_mean_sigma <- ans$post_mean_sigma
#' post_std_sigma <- ans$post_std_sigma
#' beta_draws <- ans$beta_draws
#' sigma_draws <- ans$sigma_draws
#' mc = 50
#' burn <- 10
#' iter <- burn + mc
#' deviance <- deviance3(y, x, gammacp, p, post_mean_beta, post_std_beta,
#' post_mean_sigma, post_std_sigma, beta_draws, sigma_draws, burn, iter)
#'
#' # deviance$dic
#' #     474.4673
#' # deviance$pd
#' #     5.424001
#' # deviance$devpostmean
#' #     463.6193
#'
#' @export
deviance3 <- function(y, x, gammacp, p, post_mean_beta, post_std_beta,
                      post_mean_sigma, post_std_sigma,
                      beta_draws, sigma_draws, burn, iter) {
    if (dim(y)[2] != 1){
        stop("parameter y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be a integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( length(p) != 1){
        stop("parameter p must be scalar")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    if ( !all(is.numeric(post_mean_beta))){
        stop("each entry in post_mean_beta must be numeric")
    }
    if ( !all(is.numeric(post_std_beta))){
        stop("each entry in post_std_beta must be numeric")
    }
    if ( !all(is.numeric(post_mean_sigma))){
        stop("each entry in post_mean_sigma must be numeric")
    }
    if ( !all(is.numeric(post_std_sigma))){
        stop("each entry in post_std_sigma must be numeric")
    }
    if ( !all(is.numeric(beta_draws))){
        stop("each entry in beta_draws must be numeric")
    }
    if ( !all(is.numeric(sigma_draws))){
        stop("each entry in sigma_draws must be numeric")
    }
    if ( length(burn) != 1){
        stop("parameter burn must be scalar")
    }
    if ( length(iter) != 1){
        stop("parameter iter must be scalar")
    }
    lengthp <- length(p)
    k <- dim(x)[2]
    devpostmean <- array(0, dim = c(lengthp))
    dic <- array(0, dim = c(lengthp))
    pd <- array(0, dim = c(lengthp))
    devpostmean <- 2 * negLoglikelihood(y, x, gammacp,
                                        post_mean_beta,
                                        post_mean_sigma, p)
    nsim <- dim(beta_draws[(burn + 1):iter, ])[1]
    deviance_val <- array(0, dim = c(nsim, 1))
    for (j in 1:nsim) {
        deviance_val[j, 1] <- 2 * negLoglikelihood(y, x, gammacp,
                                                   beta_draws[(burn + j), ],
                                                   sigma_draws[(burn + j), ],
                                                   p)
    }
    avg_deviance_val <- mean(deviance_val)
    dic <- (2 * avg_deviance_val) - devpostmean
    pd <- avg_deviance_val - devpostmean
    result <- list("dic" = dic,
                   "pd" = pd,
                   "devpostmean" = devpostmean)
    return(result)
}
#' NegLoglikelihood function for Ordinal Models with 3 outcomes
#'
#' This function computes the negative of the log-likelihood for quantile
#' ordinal model with 3 outcomes where the error is assumed to follow
#' an Asymmetric Laplace distribution.
#'
#' @param y         dependent variable i.e. ordinal outcome values.
#' @param x         covariate matrix of dimension \eqn{(n x k)} including a column of ones.
#' @param gammacp   row vector of cutpoints including -Inf and Inf.
#' @param beta      column vector of coeffcients of dimension \eqn{(k x 1)}.
#' @param sigma     scale factor, a scalar.
#' @param p         quantile level or skewness parameter, p in (0,1).
#'
#' @details
#' Computes the negative of the log-likelihood for quantile
#' ordinal model with 3 outcomes where the error is assumed to follow
#' an asymmetric Laplace distribution.
#'
#' @return Returns the negative log-likelihood value.
#'
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' @seealso likelihood maximization
#' @examples
#' set.seed(101)
#' data("data25j3")
#' x <- data25j3$x
#' y <- data25j3$y
#' p <- 0.25
#' gammacp <- c(-Inf, 0, 4, Inf)
#' beta <- c(1.7201671, 1.9562172, 0.8334668)
#' sigma <- 0.9684741
#' ans <- negLoglikelihood(y, x, gammacp, beta, sigma, p)
#'
#' # ans
#' #    231.8096
#'
#' @export
negLoglikelihood <- function(y, x, gammacp, beta, sigma, p) {
    if (dim(y)[2] != 1){
        stop("parameter y should be a column vector")
    }
    if ( any(!all(y == floor(y)))){
        stop("each entry of y must be a integer")
    }
    if ( !all(is.numeric(x))){
        stop("each entry in x must be numeric")
    }
    if ( !all(is.numeric(beta))){
        stop("each entry in beta must be numeric")
    }
    if ( length(sigma) != 1){
        stop("parameter sigma must be scalar")
    }
    if (any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    J <- dim(unique(y))[1]
    n <- dim(y)[1]
    lnpdf <- array(0, dim = c(n, 1))
    mu <- x %*% beta
    for (i in 1:n) {
        meanf <- mu[i]
        if (y[i] == 1) {
            lnpdf[i] <- log(alcdf(0, meanf, sigma, p))
        }
        else if (y[i] == J) {
            lnpdf[i] <- log(1 - alcdf(gammacp[J], meanf, sigma, p))
        }
        else {
            w <- (alcdf(gammacp[J], meanf, sigma, p) -
                      alcdf(gammacp[(J - 1)], meanf, sigma, p))
            lnpdf[i] <- log(w)
        }
    }
    negsuminpdf <- -sum(lnpdf)
    return(negsuminpdf)
}
#' Generates random numbers from an Asymmetric Laplace Distribution
#'
#' This function generates a vector of random numbers from an asymmetric
#' Laplace distribution with quantile p.
#'
#' @param sigma  scale factor, a scalar.
#' @param p      quantile or skewness parameter, p in (0,1).
#' @param n      number of observations
#'
#' @details
#' Generates a vector of random numbers from an asymmetric Laplace distribution,
#' as a mixture of normal–exponential distributions.
#'
#' @return Returns a vector \eqn{(n x 1)} of random numbers using an AL(0, \eqn{\sigma}, p)
#'
#' @references
#'  Kozumi, H. and Kobayashi, G. (2011). “Gibbs Sampling Methods for Bayesian Quantile Regression.”
#'  Journal of Statistical Computation and Simulation, 81(11): 1565–1578.
#'
#'  Koenker, R. and Machado, J. (1999). “Goodness of Fit and Related
#'  Inference Processes for Quantile Regression.”, Journal of
#'   American Statistics Association, 94(3): 1296-1309.
#'
#'  Keming Yu and Jin Zhang (2005). “A Three-Parameter Asymmetric
#'  Laplace Distribution.” Communications in Statistics - Theory and Methods: 1867-1879.
#'
#' @importFrom "stats" "rnorm" "rexp"
#' @seealso asymmetric Laplace distribution
#' @examples
#' set.seed(101)
#' sigma <- 2.503306
#' p <- 0.25
#' n <- 1
#' ans <- rndald(sigma, p, n)
#'
#' # ans
#' #    1.07328
#'
#' @export
rndald <- function(sigma, p, n){
    if ( any(p < 0 | p > 1)){
        stop("parameter p must be between 0 to 1")
    }
    if ( n != floor(n)){
        stop("parameter n must be a integer")
    }
    if ( length(sigma) != 1){
        stop("parameter sigma must be scalar")
    }
    u <- rnorm(n = n, mean = 0, sd = 1)
    w <- rexp(n = n, rate = 1)
    theta <- (1 - 2 * p) / (p * (1 - p))
    tau <- sqrt(2 / (p * (1 - p)))
    eps <- sigma * ( (theta * w) + (tau * sqrt(w) * u))
    return(eps)
}

#' Trace Plots for Ordinal Models
#' with 3 outcomes
#'
#' This function generates trace plots of
#' MCMC samples for \eqn{(\beta ,\sigma)} in the quantile
#' regression model with 3 outcomes.
#'
#' @param beta_draws      Gibbs draw of \eqn{\beta} vector of dimension \eqn{(k x iter)}.
#' @param sigma_draws     Gibbs draw of scale parameter, \eqn{\sigma}.
#'
#' @details
#' Trace plot is a visual depiction of the values generated from the Markov chain
#' versus the iteration number.
#'
#' @return Returns trace plots for each element of \eqn{\beta} and \eqn{\sigma}.
#'
#' @importFrom "graphics" "plot"
#' @references Rahman, M. A. (2016). “Bayesian
#' Quantile Regression for Ordinal Models.”
#' Bayesian Analysis, 11(1): 1-24.
#'
#' @seealso traces in MCMC simulations
#' @examples
#' set.seed(101)
#' data("data25j3")
#' x <- data25j3$x
#' y <- data25j3$y
#' p <- 0.25
#' ans <- quan_reg3(y, x, mc = 50, p)
#' beta_draws <- ans$beta_draws
#' sigma_draws <- ans$sigma_draws
#' trace_plot3(beta_draws, sigma_draws)
#'
#' @export
trace_plot3 <- function(beta_draws, sigma_draws) {
    if ( !all(is.numeric(beta_draws))){
        stop("each entry in beta_draws must be numeric")
    }
    if ( !all(is.numeric(sigma_draws))){
        stop("each entry in sigma_draws must be numeric")
    }
    iter <- dim(beta_draws)[1]
    k <- dim(beta_draws)[2]
    leng <- seq(0, iter - 1)
    for (t in 1:k) {
        plot(leng, beta_draws[, t], type = "l", col = "blue")
    }
    iter2 <- dim(sigma_draws)[1]
    k2 <- dim(sigma_draws)[2]
    leng2 <- seq(0, iter2 - 1)
    for (t in 1:k2) {
        plot(leng2, sigma_draws[, t], type = "l", col = "blue")
    }
}

#' Inefficiency Factor for Ordinal Models
#' with 3 outcomes
#'
#' This function calculates the inefficiency factor from the MCMC draws
#' of \eqn{(\beta, \sigma)} for an ordinal model with 3 outcomes. The
#' inefficiency factor is calculated using the batch-means method.
#'
#' @param beta_draws      Gibbs draw of coeffcients of dimension \eqn{(k x iter)}.
#' @param nlags           scalar variable with default = 2.
#' @param sigma_draws     Gibbs draw of scale factor.
#'
#' @details
#' Calculates the inefficiency factor of \eqn{(\beta, \sigma)} using the batch-means
#' method.
#'
#' @return Returns a list with components
#' \itemize{
#' \item{\code{inefficiency_beta}: }{a vector with inefficiency facor for each \eqn{\beta}.}
#' \item{\code{inefficiency_sigma}: }{a vector with inefficiency factor for each \eqn{\sigma}.}
#' }
#'
#' @importFrom "pracma" "Reshape" "std"
#'
#' @references Greenberg, E. (2012). “Introduction to Bayesian Econometrics.”
#'  Cambridge University Press, Cambridge.
#'
#' @seealso pracma
#' @examples
#' set.seed(101)
#' data("data25j3")
#' x <- data25j3$x
#' y <- data25j3$y
#' p <- 0.25
#' ans <- quan_reg3(y, x, mc = 50, p)
#' beta_draws <- ans$beta_draws
#' sigma_draws <- ans$sigma_draws
#'
#' inefficiency <- inefficiency_factor3(beta_draws, 2, sigma_draws)
#'
#' # inefficiency$inefficiency_beta
#' #     1.322590
#' #     1.287309
#' #     1.139322
#' # inefficiency$inefficiency_sigma
#' #     1.392045
#'
#' @export
inefficiency_factor3 <- function(beta_draws, nlags = 2, sigma_draws) {
    if ( !all(is.numeric(beta_draws))){
        stop("each entry in beta must be numeric")
    }
    if ( !all(is.numeric(sigma_draws))){
        stop("each entry in sigma must be numeric")
    }
    n <- dim(beta_draws)[1]
    k <- dim(beta_draws)[2]
    inefficiency_beta <- array(0, dim = c(k, 1))
    nbatch <- floor(n / nlags)
    nuse <- nbatch * nlags
    for (i in 1:k) {
        b <- beta_draws[1:nuse, i]
        xbatch <- Reshape(b, nlags, nbatch)
        mxbatch <- colMeans(xbatch)
        varxbatch <- sum( (t(mxbatch) - mean(b)) *
                              (t(mxbatch) - mean(b))) / (nbatch - 1)
        nse <- sqrt(varxbatch / (nbatch))
        rne <- (std(b, 1) / sqrt( nuse )) / nse
        inefficiency_beta[i, 1] <- 1 / rne
    }
    n2 <- dim(sigma_draws)[1]
    k2 <- dim(sigma_draws)[2]
    inefficiency_sigma <- array(0, dim = c(k2, 1))
    nbatch2 <- floor(n2 / nlags)
    nuse2 <- nbatch * nlags
    for (i in 1:k2) {
        b2 <- sigma_draws[1:nuse2, i]
        xbatch2 <- Reshape(b2, nlags, nbatch2)
        mxbatch2 <- colMeans(xbatch2)
        varxbatch2 <- sum( (t(mxbatch2) - mean(b2)) *
                               (t(mxbatch2) - mean(b2))) / (nbatch2 - 1)
        nse2 <- sqrt(varxbatch2 / (nbatch2))
        rne2 <- (std(b2, 1) / sqrt( nuse2 )) / nse2
        inefficiency_sigma[i, 1] <- 1 / rne2
    }
    result <- list("inefficiency_beta" = inefficiency_beta,
                   "inefficiency_sigma" = inefficiency_sigma)

    return(result)
}
