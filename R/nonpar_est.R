#' Nonparametric estimation wrapper
#'
#' estimate the conditional distribution of a binary outcome using either
#' xgboost or SuperLearner
#'
#' @param x a numeric matrix of covariates
#' @param y a binary vector of outcomes
#'
#' @param method either "xgboost" or "SuperLearner", method of nonparametric
#' estimation
#'
#' @param arguments an optional list with arguments to pass to `nonpar_est`. If
#' `method` is "SuperLearner", `arguments` must contain `SL.library` and
#' `cvControl`
#'
#' @param wts an optional numeric vector of observation weights
#'
#' @return a fitted model
#'
#' @export
nonpar_est <- function(x, y, method, arguments = NULL, wts = NULL) {

  mod <- NULL

  if (method == "SuperLearner") {

    stopifnot(
      "arguments must contain SL.library" = !is.null(arguments$SL.library),
      "arguments must contain cvControl" = !is.null(arguments$cvControl))

    mod <- SuperLearner::SuperLearner(
      Y = y,
      X = x,
      family = binomial(),
      obsWeights = wts,
      SL.library = arguments$SL.library,
      cvControl = arguments$cvControl)

  } else if (method == "xgboost") {

    mod <- xgboost::xgboost(
      x = x,
      y = as.factor(y),
      objective = "binary:logistic",
      weights = wts)

  } else {
    print("unrecognized method")
  }
  return(mod)
}


#' Nonparametric prediction wrapper
#'
#' predict probabilites using a fitted model from xgboost or SuperLearner
#'
#' @param method either "xgboost" or "SuperLearner", method of nonparametric
#' estimation
#'
#' @param mod a fitted model object, output of `nonpar_est`
#'
#' @param newdata a numeric matrix, predictors from new data
#'
#' @return a numeric vector, predicted probabilities
#'
#' @export
nonpar_pred <- function(mod, newdata, method) {

  pred <- NULL

  if (method == "SuperLearner") {

    pred <- predict(
      mod,
      newdata = newdata,
      onlySL = T)$pred

  } else if (method == "xgboost") {

    pred <- predict(
      mod,
      newdata = newdata)

  } else {
    stop("unrecognized method: must be 'xgboost' or 'SuperLearner'")
  }

  return(pred)

}
