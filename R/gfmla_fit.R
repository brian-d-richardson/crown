# Parametric G-formula ------------------------------------------------------

#' Fit a parametric G-formula estimator
#'
#' Fits either the naive or proposed parametric G-formula estimator for the two
#' mean potential outcomes.
#'
#' @param data Combined data containing `S`, `R`, `C`, `A`, `Y`, `wt`, and the
#'   variables in `outcome_formula`.
#' @param outcome_formula A logistic outcome-regression formula.
#' @param version `"naive"` or `"proposed"`.
#' @return Estimated risks and their covariance entries.
#' @export
fit_gformula <- function(data, outcome_formula, version = c("naive", "proposed")) {
  version <- match.arg(version)
  proposed <- version == "proposed"

  if (!proposed) {
    data <- data[data$S == 1L & data$R == 1L, ]
  }

  fit_rows <- if (proposed) {
    data$S == 1L & data$R == 1L & data$C == 0L
  } else {
    data$C == 0L
  }
  target_rows <- if (proposed) data$S == 0L else rep(TRUE, nrow(data))

  # Fit the outcome model.
  outcome_fit <- glm(
    outcome_formula, family = binomial(), data = data[fit_rows, ]
  )

  # Predict both potential outcomes in the target sample.
  data0 <- data1 <- data
  data0$A <- 0L
  data1$A <- 1L
  design0 <- model.matrix(outcome_formula, data0)
  design1 <- model.matrix(outcome_formula, data1)
  mu0 <- predict(outcome_fit, data0, type = "response")
  mu1 <- predict(outcome_fit, data1, type = "response")
  eta0 <- weighted.mean(mu0[target_rows], data$wt[target_rows])
  eta1 <- weighted.mean(mu1[target_rows], data$wt[target_rows])

  # Estimate covariance with the stacked scores.
  covariance <- .sandwich_covariance(
    c(eta0, eta1, coef(outcome_fit)),
    function(value) {
      beta <- value[-c(1, 2)]
      cbind(
        .logistic_score(data, beta, outcome_formula) * as.numeric(fit_rows),
        as.numeric(target_rows) * data$wt *
          (plogis(as.vector(design0 %*% beta)) - value[1]),
        as.numeric(target_rows) * data$wt *
          (plogis(as.vector(design1 %*% beta)) - value[2])
      )
    },
    nrow(data)
  )

  .eta_result(eta0, eta1, covariance)
}
