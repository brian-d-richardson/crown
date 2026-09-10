#' Fit a binary nuisance regression
#'
#' @param x Predictor data frame.
#' @param y Binary outcome vector.
#' @param method Either `"xgboost"` or `"SuperLearner"`.
#' @param arguments Optional learner settings; see [dml_fit()].
#' @param wts Optional observation weights.
#' @return A fitted XGBoost or SuperLearner model.
#' @keywords internal
#' @noRd
nonpar_est <- function(x, y, method, arguments = NULL, wts = NULL) {
  if (method == "xgboost") {
    settings <- list(
      x = x, y = factor(y, levels = 0:1),
      objective = "binary:logistic", weights = wts
    )
    return(do.call(xgboost, c(settings, arguments)))
  }

  SuperLearner(
    Y = y, X = x, family = binomial(), obsWeights = wts,
    SL.library = arguments$SL.library, cvControl = arguments$cvControl,
    env = asNamespace("SuperLearner")
  )
}

#' Predict nuisance probabilities
#'
#' @param mod Fitted nuisance model.
#' @param newdata Predictor data frame for held-out observations.
#' @param method Either `"xgboost"` or `"SuperLearner"`.
#' @return A numeric vector of predicted probabilities.
#' @keywords internal
#' @noRd
nonpar_pred <- function(mod, newdata, method) {
  if (method == "xgboost") {
    return(as.numeric(predict(mod, newdata = newdata)))
  }
  as.numeric(predict(mod, newdata = newdata, onlySL = TRUE)$pred)
}
