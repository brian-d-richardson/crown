#' Weighted trial and auxiliary example
#'
#' A synthetic cluster-randomized trial and a weighted auxiliary sample.
#'
#' @format A data frame with 100,000 rows and 12 variables:
#' \describe{
#'   \item{id}{Individual identifier.}
#'   \item{cluster}{Cluster identifier.}
#'   \item{S}{Trial indicator.}
#'   \item{A}{Cluster-level treatment.}
#'   \item{R}{Trial response indicator.}
#'   \item{C}{Trial censoring indicator.}
#'   \item{Y}{Binary outcome.}
#'   \item{X1, X2}{Cluster-level covariates.}
#'   \item{W1, W2}{Individual-level covariates.}
#'   \item{sampling_weight}{Auxiliary sampling weight; one for trial rows.}
#' }
#' @source Synthetic data.
"crown_example"
