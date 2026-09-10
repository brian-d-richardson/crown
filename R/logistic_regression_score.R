# Logistic regression estimating function ---------------------------------

#' Compute logistic regression scores
#'
#' Computes the observation-level score contributions from a weighted
#' logistic regression model.
#'
#' @param data A data frame containing:
#' \itemize{
#'   \item{the response and predictors in \code{formula};}
#'   \item{\code{wt}: observation weights.}
#' }
#' @param beta A numeric vector of regression coefficients.
#' @param formula A logistic regression formula.
#'
#' @return A numeric matrix with one score row for each observation.
#'
#' @keywords internal
#' @noRd
.logistic_score <- function(data, beta, formula) {
  frame <- model.frame(formula, data)
  response <- model.response(frame)
  design <- model.matrix(formula, data)
  probability <- plogis(as.vector(design %*% beta))
  design * as.vector(response - probability) * data$wt
}
