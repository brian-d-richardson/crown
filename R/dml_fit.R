#' Proposed cross-fitted DML estimator
#'
#' Fits the Proposed DML estimator with cross-fitted XGBoost nuisance models.
#'
#' @param data Combined data frame containing:
#' \itemize{
#'   \item `S`: trial (1) or auxiliary (0);
#'   \item `A`: cluster treatment (0 or 1), in both samples;
#'   \item `RW`: response indicator;
#'   \item `RY`: outcome observation indicator (1 = observed);
#'   \item `Y`: binary outcome;
#'   \item `wt`: weight (1 for trial rows, positive sampling weights for auxiliary rows);
#'   \item the baseline covariates used in the nuisance models.
#' }
#'   Use zero placeholders for unobserved `RW`, `RY`, and `Y`. Nonresponders'
#'   covariates may be missing because those rows are not used in nuisance fits.
#' @param mu_covariates,pi_covariates Character vectors of outcome and selection
#'   predictor names, respectively; not formulas. Include relevant cluster
#'   covariates; do not include `A`, `S`, `RW`, `RY`, `Y`, or `Q`.
#' @param K Number of outer cross-fitting folds, at least two.
#' @param arguments Optional XGBoost settings, e.g.,
#'   `list(nrounds = 100L, max_depth = 6L)`.
#' @param random_seed Seed set before splitting `S`-by-`A` stratified
#'   individual-level folds and fitting the nuisance models.
#' @details Uses fold-specific Hajek normalization and the equal average of
#'   the fold estimates. Predicted selection probabilities are converted to
#'   odds as `probability / (1 - probability)`. Covariance
#'   uses the mean of fold-level individual contribution covariances divided
#'   by the full observation count.
#'   For sampling weights, contribution scaling uses observation counts, not
#'   fold weight totals, so the estimate equals the weighted Hajek formula.
#'   No probability clipping or alternative learner is silently applied.
#'
#' @return A list containing:
#' \itemize{
#'   \item `eta_hat`: risks under control and treatment;
#'   \item `eta_hat_cov`: their two-by-two covariance matrix;
#'   \item `fold_estimates`: risks and covariance entries for each fold;
#'   \item `data`: input data plus folds, `Q`, out-of-fold outcome predictions,
#'     selection probabilities, and selection odds.
#' }
#' @export
fit_dml <- function(data, mu_covariates, pi_covariates, K = 5L,
                    arguments = NULL, random_seed = 1L) {

  # Split individuals within each sample and treatment arm.
  set.seed(random_seed)
  data <- as.data.frame(data)
  data$fold <- 0L
  for (s in 0:1) {
    for (a in 0:1) {
      rows <- which(data$S == s & data$A == a)
      data$fold[rows] <- sample(rep(seq_len(K), length.out = length(rows)))
    }
  }

  data$Q <- data$S * data$RW * data$RY
  data$muhat_0 <- data$muhat_1 <- data$Q_prob <- data$pihat <- NA_real_
  fold_estimates <- data.frame(
    etahat_0 = numeric(K), etahat_1 = numeric(K),
    cov_00 = numeric(K), cov_01 = numeric(K), cov_11 = numeric(K)
  )

  for (k in seq_len(K)) {

    # Training data and held-out data.
    training <- data[data$fold != k, ]
    test <- data[data$fold == k, ]
    observed0 <- training[training$Q == 1 & training$A == 0, ]
    observed1 <- training[training$Q == 1 & training$A == 1, ]
    selection0 <- training[(training$Q == 1 & training$A == 0) | training$S == 0, ]
    selection1 <- training[(training$Q == 1 & training$A == 1) | training$S == 0, ]

    mu0 <- nonpar_est(
      observed0[mu_covariates], observed0$Y, arguments
    )
    mu1 <- nonpar_est(
      observed1[mu_covariates], observed1$Y, arguments
    )
    Q0 <- nonpar_est(
      selection0[pi_covariates], selection0$Q, arguments, selection0$wt
    )
    Q1 <- nonpar_est(
      selection1[pi_covariates], selection1$Q, arguments, selection1$wt
    )

    # Predict only on the held-out fold.
    predict_mu <- test$S == 0 | test$RW == 1
    control <- test$Q == 1 & test$A == 0
    treated <- test$Q == 1 & test$A == 1
    auxiliary <- test$S == 0
    test$muhat_0[predict_mu] <- nonpar_pred(
      mu0, test[predict_mu, mu_covariates, drop = FALSE]
    )
    test$muhat_1[predict_mu] <- nonpar_pred(
      mu1, test[predict_mu, mu_covariates, drop = FALSE]
    )
    test$Q_prob[control] <- nonpar_pred(
      Q0, test[control, pi_covariates, drop = FALSE]
    )
    test$Q_prob[treated] <- nonpar_pred(
      Q1, test[treated, pi_covariates, drop = FALSE]
    )
    test$pihat <- test$Q_prob / (1 - test$Q_prob)

    # Fold-specific Hajek denominators and auxiliary weight total.
    h0 <- sum(1 / test$pihat[control])
    h1 <- sum(1 / test$pihat[treated])
    n_aux <- sum(test$wt[auxiliary])

    # Outcome averages plus inverse-probability corrections.
    phi0 <- phi1 <- numeric(nrow(test))
    phi0[control] <- (test$Y[control] - test$muhat_0[control]) / test$pihat[control] / h0
    phi1[treated] <- (test$Y[treated] - test$muhat_1[treated]) / test$pihat[treated] / h1
    phi0[auxiliary] <- test$wt[auxiliary] * test$muhat_0[auxiliary] / n_aux
    phi1[auxiliary] <- test$wt[auxiliary] * test$muhat_1[auxiliary] / n_aux

    covariance <- cov(cbind(phi0, phi1) * nrow(test)) / nrow(data)
    fold_estimates[k, ] <- c(sum(phi0), sum(phi1),
                             covariance[1, 1], covariance[1, 2], covariance[2, 2])
    data[data$fold == k, ] <- test
  }

  # Average the fold-specific estimates and covariance entries.
  average <- colMeans(fold_estimates)
  covariance <- matrix(
    c(average["cov_00"], average["cov_01"],
      average["cov_01"], average["cov_11"]), nrow = 2
  )
  list(
    eta_hat = average[c("etahat_0", "etahat_1")],
    eta_hat_cov = covariance,
    fold_estimates = fold_estimates,
    data = data
  )
}
