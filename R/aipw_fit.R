# Parametric AIPW ----------------------------------------------------------

#' Fit a parametric AIPW estimator
#'
#' Fits either the naive or proposed parametric augmented inverse-probability
#' weighted estimator for the two mean potential outcomes.
#'
#' @param data Combined data containing `S`, `R`, `C`, `A`, `Y`, `wt`, and the
#'   variables in the two formulas.
#' @param outcome_formula A logistic outcome-regression formula.
#' @param weight_formula A logistic formula for censoring `C` (naive) or
#'   selection `Q = S * R * (1 - C)` (proposed). The proposed model compares
#'   uncensored trial responders in each treatment arm with the auxiliary sample.
#' @param version `"naive"` or `"proposed"`.
#' @return Estimated risks and their covariance entries.
#' @export
fit_aipw <- function(
    data, outcome_formula, weight_formula,
    version = c("naive", "proposed")) {
  version <- match.arg(version)
  if (version == "naive") {
    data <- data[data$S == 1L & data$R == 1L, ]
    return(.fit_aipw_naive(data, outcome_formula, weight_formula))
  }
  .fit_aipw_proposed(data, outcome_formula, weight_formula)
}
#' Fit the naive parametric AIPW estimator
#'
#' Fits the outcome and censoring models among trial responders.
#' @param data Trial-responder data.
#' @param outcome_formula A logistic outcome-regression formula.
#' @param censoring_formula A logistic censoring-model formula.
#' @return Estimated risks and their covariance entries.
#' @keywords internal
#' @noRd
.fit_aipw_naive <- function(data, outcome_formula, censoring_formula) {
  n <- nrow(data)
  observed0 <- data$C == 0L & data$A == 0L
  observed1 <- data$C == 0L & data$A == 1L

  # Fit the outcome and censoring models.
  outcome_fit <- glm(
    outcome_formula, family = binomial(), data = data[data$C == 0L, ]
  )
  censoring_fit <- glm(censoring_formula, family = binomial(), data = data)

  data0 <- data1 <- data
  data0$A <- 0L
  data1$A <- 1L
  design0 <- model.matrix(outcome_formula, data0)
  design1 <- model.matrix(outcome_formula, data1)
  design_c <- model.matrix(censoring_formula, data)

  # Build each observation's AIPW contribution.
  aipw_contributions <- function(beta_mu, beta_c) {
    mu0 <- plogis(as.vector(design0 %*% beta_mu))
    mu1 <- plogis(as.vector(design1 %*% beta_mu))
    pi <- 1 - plogis(as.vector(design_c %*% beta_c))
    h0 <- sum(1 / pi[observed0])
    h1 <- sum(1 / pi[observed1])
    phi0 <- mu0
    phi1 <- mu1
    phi0[observed0] <- phi0[observed0] +
      n * (data$Y[observed0] - mu0[observed0]) / pi[observed0] / h0
    phi1[observed1] <- phi1[observed1] +
      n * (data$Y[observed1] - mu1[observed1]) / pi[observed1] / h1
    list(phi0 = phi0, phi1 = phi1)
  }

  beta_mu <- coef(outcome_fit)
  beta_c <- coef(censoring_fit)
  contribution <- aipw_contributions(beta_mu, beta_c)
  eta0 <- mean(contribution$phi0)
  eta1 <- mean(contribution$phi1)
  p_mu <- length(beta_mu)

  covariance <- .sandwich_covariance(
    c(eta0, eta1, beta_mu, beta_c),
    function(value) {
      beta_mu_value <- value[2 + seq_len(p_mu)]
      beta_c_value <- tail(value, length(beta_c))
      phi <- aipw_contributions(beta_mu_value, beta_c_value)
      cbind(
        .logistic_score(data, beta_mu_value, outcome_formula) * (1L - data$C),
        .logistic_score(data, beta_c_value, censoring_formula),
        phi$phi0 - value[1],
        phi$phi1 - value[2]
      )
    },
    n
  )

  .eta_result(eta0, eta1, covariance)
}

#' Fit the proposed parametric AIPW estimator
#'
#' Fits the outcome and arm-specific selection models.
#' @param data Combined trial and auxiliary data.
#' @param outcome_formula A logistic outcome-regression formula.
#' @param propensity_formula A logistic sample-membership formula.
#' @return Estimated risks and their covariance entries.
#' @keywords internal
#' @noRd
.fit_aipw_proposed <- function(data, outcome_formula, propensity_formula) {
  data$Q <- data$S * data$R * (1L - data$C)
  observed <- data$Q == 1L
  observed0 <- observed & data$A == 0L
  observed1 <- observed & data$A == 1L
  auxiliary <- data$S == 0L
  fit0 <- observed0 | auxiliary
  fit1 <- observed1 | auxiliary
  n <- nrow(data)

  # Fit the outcome and arm-specific selection models.
  outcome_fit <- glm(
    outcome_formula, family = binomial(), data = data[observed, ]
  )
  propensity0 <- do.call(glm, list(
    formula = propensity_formula, family = binomial(),
    data = data[fit0, ], weights = data$wt[fit0]
  ))
  propensity1 <- do.call(glm, list(
    formula = propensity_formula, family = binomial(),
    data = data[fit1, ], weights = data$wt[fit1]
  ))

  data0 <- data1 <- data
  data0$A <- 0L
  data1$A <- 1L
  design0 <- model.matrix(outcome_formula, data0)
  design1 <- model.matrix(outcome_formula, data1)
  design_q <- model.matrix(propensity_formula, data)
  auxiliary_weight <- sum(data$wt[auxiliary])

  # Build each observation's AIPW contribution.
  aipw_contributions <- function(beta_mu, beta_q0, beta_q1) {
    mu0 <- plogis(as.vector(design0 %*% beta_mu))
    mu1 <- plogis(as.vector(design1 %*% beta_mu))
    probability0 <- plogis(as.vector(
      design_q[observed0, , drop = FALSE] %*% beta_q0
    ))
    probability1 <- plogis(as.vector(
      design_q[observed1, , drop = FALSE] %*% beta_q1
    ))
    odds0 <- probability0 / (1 - probability0)
    odds1 <- probability1 / (1 - probability1)
    h0 <- sum(1 / odds0)
    h1 <- sum(1 / odds1)

    phi0 <- phi1 <- numeric(n)
    phi0[auxiliary] <- n * data$wt[auxiliary] * mu0[auxiliary] / auxiliary_weight
    phi1[auxiliary] <- n * data$wt[auxiliary] * mu1[auxiliary] / auxiliary_weight
    phi0[observed0] <- n *
      (data$Y[observed0] - mu0[observed0]) / odds0 / h0
    phi1[observed1] <- n *
      (data$Y[observed1] - mu1[observed1]) / odds1 / h1
    list(phi0 = phi0, phi1 = phi1)
  }

  beta_mu <- coef(outcome_fit)
  beta_q0 <- coef(propensity0)
  beta_q1 <- coef(propensity1)
  contribution <- aipw_contributions(beta_mu, beta_q0, beta_q1)
  eta0 <- mean(contribution$phi0)
  eta1 <- mean(contribution$phi1)
  p_mu <- length(beta_mu)
  p_q0 <- length(beta_q0)

  covariance <- .sandwich_covariance(
    c(eta0, eta1, beta_mu, beta_q0, beta_q1),
    function(value) {
      beta_mu_value <- value[2 + seq_len(p_mu)]
      beta_q0_value <- value[2 + p_mu + seq_len(p_q0)]
      beta_q1_value <- tail(value, length(beta_q1))
      phi <- aipw_contributions(beta_mu_value, beta_q0_value, beta_q1_value)
      cbind(
        .logistic_score(data, beta_mu_value, outcome_formula) * observed,
        .logistic_score(data, beta_q0_value, propensity_formula) * fit0,
        .logistic_score(data, beta_q1_value, propensity_formula) * fit1,
        phi$phi0 - value[1],
        phi$phi1 - value[2]
      )
    },
    n
  )

  .eta_result(eta0, eta1, covariance)
}
