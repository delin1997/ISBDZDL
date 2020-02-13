#' Generate mixed-distribution random variables
#'
#' \code{r_gaussian_mixed}(or \code{r_t_mixed}) returns a matrix of
#' mixed-Gaussian(or mixed-t) random variables/vectors
#' with or without class labels.
#'
#' The default value for \code{n} is 100.
#'
#' If \code{tag} is \code{TRUE}, there will be an extra column added at the
#' end of the output matrix representing the component to which each
#' random variable or vector belongs.
#'
#' @param n number of random variables/vectors to be generated.
#' @param mu matrix with each column the mean vector of each Gaussian component.
#' @param sigma array consisting of the sigma matrices of Gaussian components.
#' @param weight vector of the weights for Gaussian components.
#' @param tag logical; if FALSE(default), the matrix returned will not contain
#'   an extra column of class label.
#'
#' @return A matrix of mixed-Gaussian(or mixed-t) random variables/vectors
#' with or without class labels.
#'
#' @examples
#' mu <- matrix(rep(c(-1, 0, 1), 4), byrow = T, nrow = 4)
#' sigma <- matrix(rep(0, 4 * 4), nrow = 4)
#' for (i in 1:(4 - 1)) {
#'   for (j in (i + 1):4) {
#'     sigma[i, j] <- 0.5^ (j - i)
#'   }
#' }
#' sigma <- sigma + t(sigma) + diag(1, nrow = 4)
#' sigma <- array(rep(sigma, 3), dim = c(4, 4, 3))
#' weight <- rep(1, 3) / 3
#' r_g <- r_gaussian_mixed(n = 100, mu, sigma, weight, tag = FALSE)
#' r_t <- r_t_mixed(n = 100, sigma, df = 1, mu, weight, tag = FALSE)
#' @name r_mixed
NULL

#' @rdname r_mixed
r_gaussian_mixed <- function(n = 100, mu, sigma, weight, tag = FALSE) {
  r <- runif(n, 0, 1)
  c.weight <- cumsum(weight)
  if (tag) {
    r_g_m <- matrix(0, ncol = nrow(as.matrix(mu)) + 1, nrow = n)
    i <- 1
    for (x in r) {
      order <- sum(c.weight < x) + 1
      r_g_m[i, ] <- c(MASS::mvrnorm(1, mu[, order], sigma[, , order]), order)
      i <- i + 1
    }
  } else {
    r_g_m <- matrix(0, ncol = nrow(as.matrix(mu)), nrow = n)
    i <- 1
    for (x in r) {
      order <- sum(c.weight < x) + 1
      r_g_m[i, ] <- MASS::mvrnorm(1, mu[, order], sigma[, , order])
      i <- i + 1
    }
  }
  return(r_g_m)
}

#' @rdname r_mixed
#' @param df degree of freedom for t-distribution

r_t_mixed <- function(n = 100, sigma, df, mu, weight, tag = FALSE) {
  r <- runif(n, 0, 1)
  c.weight <- cumsum(weight)
  if (tag) {
    r_t_m <- matrix(0, ncol = nrow(as.matrix(mu)) + 1, nrow = n)
    i <- 1
    for (x in r) {
      order <- sum(c.weight < x) + 1
      r_t_m[i, ] <- c(mvtnorm::rmvt(1, sigma = sigma[, , order],
        df = df, delta = mu[, order]), order)
      i <- i + 1
    }
  }else{
    r_t_m <- matrix(0, ncol = nrow(as.matrix(mu)), nrow = n)
    i <- 1
    for (x in r) {
      order <- sum(c.weight < x) + 1
      r_t_m[i, ] <- mvtnorm::rmvt(1, sigma = sigma[, , order],
        df = df, delta = mu[, order])
      i <- i + 1
    }
  }
  return(r_t_m)
}
