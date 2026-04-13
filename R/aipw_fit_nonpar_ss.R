## proposed AIPW estimator with random forest and sample splitting
aipw_fit_nonpar_ss <- function(dat, mu_covariates, pi_covariates, K = 5) {

  dat <- dat %>%
    group_by(S) %>%
    mutate(fold = sample(rep(1:K, length.out = n()))) %>%
    ungroup() %>%
    mutate(Q = S * R * (1 - C))

  # initialize predicted probabilities
  dat$muhat_0 <- NA_real_
  dat$muhat_1 <- NA_real_
  dat$Q_prob  <- NA_real_


  for (k in 1:K) {

    train_idx <- dat$fold != k
    test_idx  <- dat$fold == k

    train <- dat[train_idx, ]
    test  <- dat[test_idx, ]

    ## ---------------------------
    ## Outcome regression (mu)
    ## ---------------------------

    train_outcome <- train$S == 1 & train$R == 1 & train$C == 0

    outcome_reg <- ranger(
      x = train[train_outcome, c("A", mu_covariates)],
      y = train$Y[train_outcome],
      probability = TRUE
    )

    ## test design matrices
    idx_mu <- test$S == 0 | test$R == 1

    test0 <- test[idx_mu, ]
    test1 <- test0

    test0$A <- 0
    test1$A <- 1

    ## predictions
    mu0 <- predict(outcome_reg,
                   data = test0[, c("A", mu_covariates)])$predictions[,2]

    mu1 <- predict(outcome_reg,
                   data = test1[, c("A", mu_covariates)])$predictions[,2]

    dat$muhat_0[test_idx][idx_mu] <- mu0
    dat$muhat_1[test_idx][idx_mu] <- mu1


    ## ---------------------------
    ## Q models (pi)
    ## ---------------------------

    train_Q0 <- (train$Q == 1 & train$A == 0) | train$S == 0
    train_Q1 <- (train$Q == 1 & train$A == 1) | train$S == 0

    Q_reg_0 <- ranger(
      x = train[train_Q0, pi_covariates],
      y = train$Q[train_Q0],
      case.weights = train$wt[train_Q0],
      probability = TRUE
    )

    Q_reg_1 <- ranger(
      x = train[train_Q1, pi_covariates],
      y = train$Q[train_Q1],
      case.weights = train$wt[train_Q1],
      probability = TRUE
    )

    ## test indices
    Q_ind_test <- test$Q == 1

    idx0 <- Q_ind_test & test$A == 0
    idx1 <- Q_ind_test & test$A == 1

    ## predictions
    if (any(idx0)) {
      dat$Q_prob[test_idx][idx0] <-
        predict(Q_reg_0,
                data = test[idx0, pi_covariates])$predictions[,2]
    }

    if (any(idx1)) {
      dat$Q_prob[test_idx][idx1] <-
        predict(Q_reg_1,
                data = test[idx1, pi_covariates])$predictions[,2]
    }
  }


  ## estimated propensity scores in trial data
  Q_ind <- dat$Q == 1
  dat$pihat <- NA_real_
  dat$pihat[Q_ind] <- dat$Q_prob[Q_ind] /
    (1 - dat$Q_prob[Q_ind])


  # AIPW estimator ----------------------------------------------------------

  ## Hajek denominator
  n_trial_hat <- dat %>%
    filter(S == 1, R == 1, C == 0) %>%
    mutate(
      term_0 = (1 - A) / pihat,
      term_1 = A / pihat) %>%
    summarise(
      nhat_0 = sum(term_0),
      nhat_1 = sum(term_1)) %>%
    unlist()

  ## auxiliary sample size
  n_aux <- sum(dat$wt[dat$S == 0])

  ## AIPW estimator
  etahat <- dat %>%
    mutate(

      ## IPW terms
      ipw_0 = ifelse(S == 1 & R == 1 & C == 0 & A == 0,
                     (Y - muhat_0) / pihat,
                     0),

      ipw_1 = ifelse(S == 1 & R == 1 & C == 0 & A == 1,
                     (Y - muhat_1) / pihat,
                     0),

      ## outcome regression terms
      or_0 = ifelse(S == 0, muhat_0 * wt, 0),
      or_1 = ifelse(S == 0, muhat_1 * wt, 0),

      ## influence function
      if_0 = (S * R * (1 - C) * ipw_0 / n_trial_hat[1] +
        (1 - S) * or_0 / n_aux) * nrow(dat),
      if_1 = (S * R * (1 - C) * ipw_1 / n_trial_hat[2] +
        (1 - S) * or_1 / n_aux) * nrow(dat)) %>%

    summarise(

      ## AIPW estimator
      etahat_0 = mean(if_0),
      etahat_1 = mean(if_1),

      ## covariance estimator
      cov_00 = var(if_0) / nrow(dat),
      cov_01 = cov(if_0, if_1) / nrow(dat),
      cov_11 = var(if_1) / nrow(dat))


  # return list of results --------------------------------------------------

  res <- list(

    # causal parameter estimates and covariance
    eta_results =
      data.frame(
        etahat_0 = etahat$etahat_0,
        etahat_1 = etahat$etahat_1,
        cov_00 = etahat$cov_00,
        cov_01 = etahat$cov_01,
        cov_11 = etahat$cov_11),

    # data set
    dat = dat)

  return(res)

}
