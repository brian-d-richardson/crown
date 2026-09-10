#' Fit a binary nuisance regression
#'
#' @param x Predictor data frame.
#' @param y Binary outcome vector.
#' @param arguments Optional learner settings; see [dml_fit()].
#' @param wts Optional observation weights.
#' @return A fitted XGBoost model.
#' @keywords internal
#' @noRd
nonpar_est <- function(x, y, arguments = NULL, wts = NULL) {
  settings <- list(
    x = x, y = factor(y, levels = 0:1),
    objective = "binary:logistic", weights = wts
  )
  do.call(xgboost, c(settings, arguments))
}

#' Predict nuisance probabilities
#'
#' @param mod Fitted nuisance model.
#' @param newdata Predictor data frame for held-out observations.
#' @return A numeric vector of predicted probabilities.
#' @keywords internal
#' @noRd
nonpar_pred <- function(mod, newdata) {
  as.numeric(predict(mod, newdata = newdata))
}
