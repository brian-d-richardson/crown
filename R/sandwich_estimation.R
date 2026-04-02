#' Sandwich variance estimator
#'
#' @param param a numeric vector, estimated parameters
#' @param get.psi estimating function
#' @param n a positive integer, the sample size
#'
#' @return estimated covariance matrix of param
#'
#' @export
get.sand.est <- function(param, get.psi, n) {

  # D: empirical mean of derivative of Psi
  D <- numDeriv::jacobian(
    f = function(x) colSums(get.psi(x)),
    x = param,
    method = "simple") / -n
  Dinv <- solve(D)

  # B: empirical mean of outer product of Psi
  Psi <- get.psi(param)
  Omega <- matrix(rowMeans(apply(Psi, 1, function(psi) psi %*% t(psi))),
                  nrow = length(param))

  # sandwich estimator
  return(Dinv %*% Omega %*% t(Dinv) / n)
}
