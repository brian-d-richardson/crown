#' Cross-fitted naive point estimates
#'
#' Uses the existing naive AIPW formula on trial responders. Fits separate
#' outcome models by arm and one pooled censoring model using covariates only,
#' matching the default logistic censoring specification. Auxiliary rows and
#' nonresponders do not enter the fits or averages.
#'
#' This differs from the Proposed cross-fitted estimator. The covariance is
#' left unavailable because the naive score is not generally orthogonal.
#' @inheritParams crown
#' @param data Combined data with standardized columns.
#' @return Risks, unavailable covariance, fold estimates, and held-out predictions.
#' @keywords internal
#' @noRd
.naive_xgboost <- function(data, covariates, K, arguments, random_seed,
                           estimator = c("aipw", "gformula", "ipw")) {
  estimator <- match.arg(estimator)
  data <- data[data$S == 1L & data$R == 1L, ]
  if (length(K) != 1L || is.na(K) || K < 2L || K != as.integer(K)) {
    stop("K must be an integer of at least 2.")
  }
  if (any(tabulate(data$A[data$C == 0L] + 1L, nbins = 2L) < K)) {
    stop("Naive cross-fitting needs at least K uncensored responders in each arm.")
  }

  # Split responders by treatment and censoring status.
  set.seed(random_seed)
  data$fold <- 0L
  for (a in 0:1) {
    for (c in 0:1) {
      rows <- which(data$A == a & data$C == c)
      data$fold[rows] <- sample(rep(seq_len(K), length.out = length(rows)))
    }
  }
  data$muhat_0 <- data$muhat_1 <- data$pihat <- NA_real_
  fold_estimates <- data.frame(etahat_0 = numeric(K), etahat_1 = numeric(K))

  for (k in seq_len(K)) {
    training <- data[data$fold != k, ]
    test <- data[data$fold == k, ]
    observed0 <- training$A == 0L & training$C == 0L
    observed1 <- training$A == 1L & training$C == 0L

    # Fit only the nuisance models required by the selected estimator.
    if (estimator != "ipw") {
      mu0 <- nonpar_est(training[observed0, covariates, drop = FALSE],
                       training$Y[observed0], "xgboost", arguments)
      mu1 <- nonpar_est(training[observed1, covariates, drop = FALSE],
                       training$Y[observed1], "xgboost", arguments)
    }
    if (estimator != "gformula") {
      censoring <- nonpar_est(training[covariates], training$C,
                             "xgboost", arguments)
      test$pihat <- 1 - nonpar_pred(censoring, test[covariates], "xgboost")
    }
    if (estimator != "ipw") {
      test$muhat_0 <- nonpar_pred(mu0, test[covariates], "xgboost")
      test$muhat_1 <- nonpar_pred(mu1, test[covariates], "xgboost")
    }

    # Average predictions over responders, then add Hajek residual corrections.
    control <- test$A == 0L & test$C == 0L
    treated <- test$A == 1L & test$C == 0L
    if (estimator == "gformula") {
      fold_estimates[k, ] <- c(mean(test$muhat_0), mean(test$muhat_1))
    } else if (estimator == "ipw") {
      fold_estimates[k, ] <- c(
        weighted.mean(test$Y[control], 1 / test$pihat[control]),
        weighted.mean(test$Y[treated], 1 / test$pihat[treated])
      )
    } else {
      fold_estimates$etahat_0[k] <- mean(test$muhat_0) +
        weighted.mean(test$Y[control] - test$muhat_0[control], 1 / test$pihat[control])
      fold_estimates$etahat_1[k] <- mean(test$muhat_1) +
        weighted.mean(test$Y[treated] - test$muhat_1[treated], 1 / test$pihat[treated])
    }
    data[data$fold == k, ] <- test
  }

  warning("Naive XGBoost returns point estimates only; standard errors and confidence intervals are not implemented.",
          call. = FALSE)
  list(eta_hat = colMeans(fold_estimates),
       eta_hat_cov = matrix(NA_real_, 2L, 2L),
       fold_estimates = fold_estimates, dat = data)
}
