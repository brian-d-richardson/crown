#' Proposed cross-fitted estimators
#'
#' The default AIPW estimator uses cross-fitted nuisance models.
#' G-formula uses outcome regressions only; IPW uses selection regressions only.
#' These two alternatives return point estimates without standard errors.
#'
#' @param dat Combined data frame containing:
#' \itemize{
#'   \item `S`: trial (1) or auxiliary (0);
#'   \item `A`: cluster treatment (0 or 1), in both samples;
#'   \item `R`: response indicator;
#'   \item `C`: censoring indicator (1 = censored);
#'   \item `Y`: binary outcome;
#'   \item `wt`: weight (1 for trial rows, positive sampling weights for auxiliary rows);
#'   \item the baseline covariates used in the nuisance models.
#' }
#'   Use zero placeholders for unobserved `R`, `C`, and `Y`. Nonresponders'
#'   covariates may be missing because those rows are not used in nuisance fits.
#' @param mu_covariates,pi_covariates Character vectors of outcome and selection
#'   predictor names, respectively; not formulas. Include relevant cluster
#'   covariates; do not include `A`, `S`, `R`, `C`, `Y`, or `Q`.
#' @param K Number of outer cross-fitting folds, at least two.
#' @param method `"xgboost"` (default) or `"SuperLearner"`.
#' @param arguments Optional named learner arguments. XGBoost uses its standard
#'   defaults (100 boosting rounds); e.g.,
#'   `list(nrounds = 100L, max_depth = 6L)`. SuperLearner requires
#'   `SL.library` and `cvControl`, e.g.,
#'   `list(SL.library = c("SL.glm", "SL.xgboost"), cvControl = list(V = 5L))`.
#'   Inner SuperLearner CV is distinct from outer cross-fitting.
#' @param random_seed Seed set before splitting `S`-by-`A` stratified
#'   individual-level folds and fitting the nuisance models.
#' @param estimator `"aipw"` (default), `"gformula"`, or `"ipw"`.
#'
#' @details Uses fold-specific Hajek normalization and the equal average of
#'   the fold estimates. Selection odds are `zeta / (1 - zeta)`. Covariance
#'   uses the mean of fold-level individual contribution covariances divided
#'   by the full observation count.
#'   For sampling weights, contribution scaling uses observation counts, not
#'   fold weight totals, so the estimate equals the weighted Hajek formula.
#'   No probability clipping or alternative learner is silently applied.
#'
#' @return A list containing:
#' \itemize{
#'   \item `eta_hat`: risks under control and treatment;
#'   \item `eta_hat_cov`: their two-by-two covariance matrix (`NA` except AIPW);
#'   \item `fold_estimates`: risks and covariance entries for each fold;
#'   \item `dat`: input data plus folds, `Q`, out-of-fold outcome predictions,
#'     selection probabilities, and selection odds.
#' }
#' @export
dml_fit <- function(dat, mu_covariates, pi_covariates, K = 5L,
                    method = "xgboost", arguments = NULL, random_seed = 1L,
                    estimator = c("aipw", "gformula", "ipw")) {

  # Split individuals within each sample and treatment arm.
  method <- match.arg(method, c("xgboost", "SuperLearner"))
  estimator <- match.arg(estimator)
  set.seed(random_seed)
  dat <- as.data.frame(dat)
  dat$fold <- 0L
  for (s in 0:1) {
    for (a in 0:1) {
      rows <- which(dat$S == s & dat$A == a)
      dat$fold[rows] <- sample(rep(seq_len(K), length.out = length(rows)))
    }
  }

  dat$Q <- dat$S * dat$R * (1 - dat$C)
  dat$muhat_0 <- dat$muhat_1 <- dat$Q_prob <- dat$pihat <- NA_real_
  fold_estimates <- data.frame(
    etahat_0 = numeric(K), etahat_1 = numeric(K),
    cov_00 = numeric(K), cov_01 = numeric(K), cov_11 = numeric(K)
  )

  for (k in seq_len(K)) {

    # Training data and held-out data.
    training <- dat[dat$fold != k, ]
    test <- dat[dat$fold == k, ]
    observed0 <- training[training$Q == 1 & training$A == 0, ]
    observed1 <- training[training$Q == 1 & training$A == 1, ]
    selection0 <- training[(training$Q == 1 & training$A == 0) | training$S == 0, ]
    selection1 <- training[(training$Q == 1 & training$A == 1) | training$S == 0, ]

    # Fit only the nuisance models required by the selected estimator.
    if (estimator != "ipw") {
      mu0 <- nonpar_est(
        observed0[mu_covariates], observed0$Y, method, arguments
      )
      mu1 <- nonpar_est(
        observed1[mu_covariates], observed1$Y, method, arguments
      )
    }
    if (estimator != "gformula") {
      Q0 <- nonpar_est(
        selection0[pi_covariates], selection0$Q, method, arguments, selection0$wt
      )
      Q1 <- nonpar_est(
        selection1[pi_covariates], selection1$Q, method, arguments, selection1$wt
      )
    }

    # Predict only on the held-out fold.
    predict_mu <- test$S == 0 | test$R == 1
    control <- test$Q == 1 & test$A == 0
    treated <- test$Q == 1 & test$A == 1
    auxiliary <- test$S == 0
    if (estimator != "ipw") {
      test$muhat_0[predict_mu] <- nonpar_pred(
        mu0, test[predict_mu, mu_covariates, drop = FALSE], method
      )
      test$muhat_1[predict_mu] <- nonpar_pred(
        mu1, test[predict_mu, mu_covariates, drop = FALSE], method
      )
    }
    if (estimator != "gformula") {
      test$Q_prob[control] <- nonpar_pred(
        Q0, test[control, pi_covariates, drop = FALSE], method
      )
      test$Q_prob[treated] <- nonpar_pred(
        Q1, test[treated, pi_covariates, drop = FALSE], method
      )
      test$pihat <- test$Q_prob / (1 - test$Q_prob)
    }

    # G-formula and IPW have their own point estimates, not AIPW variances.
    if (estimator != "aipw") {
      if (estimator == "gformula") {
        eta0 <- weighted.mean(test$muhat_0[auxiliary], test$wt[auxiliary])
        eta1 <- weighted.mean(test$muhat_1[auxiliary], test$wt[auxiliary])
      } else {
        eta0 <- weighted.mean(test$Y[control], 1 / test$pihat[control])
        eta1 <- weighted.mean(test$Y[treated], 1 / test$pihat[treated])
      }
      fold_estimates[k, ] <- c(eta0, eta1, NA_real_, NA_real_, NA_real_)
      dat[dat$fold == k, ] <- test
      next
    }

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

    covariance <- cov(cbind(phi0, phi1) * nrow(test)) / nrow(dat)
    fold_estimates[k, ] <- c(sum(phi0), sum(phi1),
                             covariance[1, 1], covariance[1, 2], covariance[2, 2])
    dat[dat$fold == k, ] <- test
  }

  # Average the fold-specific estimates and covariance entries.
  average <- colMeans(fold_estimates)
  covariance <- matrix(
    c(average["cov_00"], average["cov_01"],
      average["cov_01"], average["cov_11"]), nrow = 2
  )
  if (estimator != "aipw") {
    warning("Cross-fitted G-formula/IPW returns point estimates only; standard errors and confidence intervals are not implemented.",
            call. = FALSE)
  }
  list(
    eta_hat = average[c("etahat_0", "etahat_1")],
    eta_hat_cov = covariance,
    fold_estimates = fold_estimates,
    dat = dat
  )
}
