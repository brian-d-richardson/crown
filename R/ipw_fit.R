# Parametric IPW ------------------------------------------------------------

#' Fit a parametric IPW estimator
#'
#' Fits either the naive or proposed parametric inverse-probability weighted
#' estimator for the two mean potential outcomes.
#'
#' @param data Combined data containing `S`, `RW`, `RY`, `A`, `Y`, `wt`, and the
#'   variables in `weight_formula`.
#' @param weight_formula A logistic formula for outcome observation `RY` (naive)
#' or selection `Q = S * RW * RY` (proposed). The proposed model compares
#'  trial responders with observed outcome in each treatment arm with the
#'  auxiliary sample.
#' @param version `"naive"` or `"proposed"`.
#' @return Estimated risks and their covariance entries.
#' @export
fit_ipw <- function(
    data, weight_formula, version = c("naive", "proposed")) {
  version <- match.arg(version)
  if (version == "naive") {
    data <- data[data$S == 1L & data$RW == 1L, ]
    return(.fit_ipw_naive(data, weight_formula))
  }
  .fit_ipw_proposed(data, weight_formula)
}
#' Fit the naive IPW estimator
#'
#' Fits the outcome observation model and treatment-specific Hajek estimates.
#' @param data Trial-responder data.
#' @param RY_formula A formula for the outcome observation model.
#' @return Estimated risks and their covariance entries.
#' @keywords internal
#' @noRd
.fit_ipw_naive <- function(data, RY_formula) {
  # observed outcome probabilities.
  RY_fit <- glm(
    RY_formula, family = binomial(), data = data
  )
  design <- model.matrix(RY_formula, data)
  pi_RY <- predict(RY_fit, data, type = "response")
  observed0 <- data$RY == 1L & data$A == 0L
  observed1 <- data$RY == 1L & data$A == 1L

  # Hajek estimates.
  denominator0 <- sum(1 / pi_RY[observed0])
  denominator1 <- sum(1 / pi_RY[observed1])
  eta0 <- sum(data$Y[observed0] / pi_RY[observed0]) / denominator0
  eta1 <- sum(data$Y[observed1] / pi_RY[observed1]) / denominator1

  # Sandwich covariance.
  covariance <- .sandwich_covariance(
    c(eta0, eta1, denominator0, denominator1, coef(RY_fit)),
    function(value) {
      beta <- value[-seq_len(4)]
      pi_RY_value <- plogis(as.vector(design %*% beta))
      cbind(
        .logistic_score(data, beta, RY_formula),
        ifelse(observed0, 1 / pi_RY_value, 0) - value[3] / nrow(data),
        ifelse(observed1, 1 / pi_RY_value, 0) - value[4] / nrow(data),
        ifelse(observed0, nrow(data) * data$Y / pi_RY_value / value[3], 0) - value[1],
        ifelse(observed1, nrow(data) * data$Y / pi_RY_value / value[4], 0) - value[2]
      )
    },
    nrow(data)
  )

  .eta_result(eta0, eta1, covariance)
}

#' Fit the proposed IPW estimator
#'
#' Fits treatment-specific selection models and Hajek estimates.
#' @param data Combined trial and auxiliary data.
#' @param propensity_formula formula for propensity score model
#' @return Estimated risks and their covariance entries.
#' @keywords internal
#' @noRd
.fit_ipw_proposed <- function(data, propensity_formula) {
  # Observed outcomes and auxiliary rows.
  data$Q <- data$S * data$RW * data$RY
  restricted <- data$Q == 1L | data$S == 0L
  fit0 <- restricted & (data$A == 0L | data$Q == 0L)
  fit1 <- restricted & (data$A == 1L | data$Q == 0L)
  data0 <- data[fit0, ]
  data1 <- data[fit1, ]

  # treatment-specific selection models.
  propensity0 <- do.call(glm, list(
    formula = propensity_formula, family = binomial(),
    data = data0, weights = data0$wt
  ))
  propensity1 <- do.call(glm, list(
    formula = propensity_formula, family = binomial(),
    data = data1, weights = data1$wt
  ))
  design <- model.matrix(propensity_formula, data)
  observed0 <- data$Q == 1L & data$A == 0L
  observed1 <- data$Q == 1L & data$A == 1L

  # Selection odds.
  probability <- rep(0, nrow(data))
  probability[observed0] <- predict(propensity0, data[observed0, ], type = "response")
  probability[observed1] <- predict(propensity1, data[observed1, ], type = "response")
  odds <- probability / (1 - probability)

  # Hajek estimates.
  denominator0 <- sum(1 / odds[observed0])
  denominator1 <- sum(1 / odds[observed1])
  eta0 <- sum(data$Y[observed0] / odds[observed0]) / denominator0
  eta1 <- sum(data$Y[observed1] / odds[observed1]) / denominator1
  length0 <- length(coef(propensity0))

  # Sandwich covariance.
  covariance <- .sandwich_covariance(
    c(
      eta0, eta1, denominator0, denominator1,
      coef(propensity0), coef(propensity1)
    ),
    function(value) {
      beta0 <- value[4 + seq_len(length0)]
      beta1 <- tail(value, length(coef(propensity1)))
      probability_value <- rep(0, nrow(data))
      probability_value[observed0] <- plogis(
        as.vector(design[observed0, , drop = FALSE] %*% beta0)
      )
      probability_value[observed1] <- plogis(
        as.vector(design[observed1, , drop = FALSE] %*% beta1)
      )
      odds_value <- probability_value / (1 - probability_value)
      cbind(
        .logistic_score(data, beta0, propensity_formula) * fit0,
        .logistic_score(data, beta1, propensity_formula) * fit1,
        ifelse(observed0, 1 / odds_value, 0) - value[3] / nrow(data),
        ifelse(observed1, 1 / odds_value, 0) - value[4] / nrow(data),
        ifelse(observed0, nrow(data) * data$Y / odds_value / value[3], 0) - value[1],
        ifelse(observed1, nrow(data) * data$Y / odds_value / value[4], 0) - value[2]
      )
    },
    nrow(data)
  )

  .eta_result(eta0, eta1, covariance)
}
