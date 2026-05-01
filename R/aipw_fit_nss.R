## proposed AIPW estimator with nonparametric sample splitting
aipw_fit_nss_v1 <- function(dat, mu_covariates, pi_covariates, K,
                            method, arguments) {

  dat <- dat %>%
    group_by(S, A) %>%
    mutate(fold = sample(rep(1:K, length.out = n()))) %>%
    ungroup() %>%
    mutate(Q = S * R * (1 - C))

  ## initialize predicted probabilities
  dat$muhat_0 <- NA_real_
  dat$muhat_1 <- NA_real_
  dat$Q_prob  <- NA_real_

  ## loop through K folds
  for (k in 1:K) {

    ## training and testing data
    train_idx <- dat$fold != k
    test_idx  <- dat$fold == k
    train <- dat[train_idx, ]
    test  <- dat[test_idx, ]

    ## train outcome regression function mu
    train_outcome <- train %>%
      filter(S == 1, R == 1, C == 0)

    mu_mod <- nonpar_est(
      x = train_outcome %>% select(all_of(c("A", mu_covariates))),
      y = train_outcome %>% select(Y) %>% unlist() %>% factor(),
      method = method,
      arguments = arguments)

    ## predicted probabilities in test data set with A = 0 and A = 1
    test_mu_idx <- dat$fold == k & (dat$S == 0 | dat$R == 1)
    test_0 <- dat[test_mu_idx,] %>% mutate(A = 0)
    test_1 <- dat[test_mu_idx,] %>% mutate(A = 1)

    dat$muhat_0[test_mu_idx] <- nonpar_pred(
      mod = mu_mod,
      newdata = test_0 %>% select(all_of(c("A", mu_covariates))),
      method = method)

    dat$muhat_1[test_mu_idx] <- nonpar_pred(
      mod = mu_mod,
      newdata = test_1 %>% select(all_of(c("A", mu_covariates))),
      method = method)

    ## train Q predictor
    train_Q0 <- train %>% filter((Q == 1 & A == 0) | S == 0)
    train_Q1 <- train %>% filter((Q == 1 & A == 1) | S == 0)

    Q_reg_0 <- nonpar_est(
      x = train_Q0 %>% select(all_of(pi_covariates)),
      y = train_Q0$Q %>% factor(),
      wts = train_Q0$wt,
      method = method,
      arguments = arguments)

    Q_reg_1 <- nonpar_est(
      x = train_Q1 %>% select(all_of(pi_covariates)),
      y = train_Q1$Q %>% factor(),
      wts = train_Q1$wt,
      method = method,
      arguments = arguments)

    ## predicted Q probabilities among uncensored responders in testing data
    test_Q_0_ind <- dat$fold == k & dat$Q == 1 & dat$A == 0
    test_Q_1_ind <- dat$fold == k & dat$Q == 1 & dat$A == 1

    dat$Q_prob[test_Q_0_ind] <- nonpar_pred(
      mod = Q_reg_0,
      newdata = dat[test_Q_0_ind,] %>% select(all_of(pi_covariates)),
      method = method)

    dat$Q_prob[test_Q_1_ind] <- nonpar_pred(
      mod = Q_reg_1,
      newdata = dat[test_Q_1_ind,] %>% select(all_of(pi_covariates)),
      method = method)

    print(paste0("fold ", k, "/", K, " complete"))

  }

  ## estimated propensity scores in trial data
  Q_ind <- dat$Q == 1
  dat$pihat <- NA_real_
  dat$pihat[Q_ind] <- dat$Q_prob[Q_ind] /
    (1 - dat$Q_prob[Q_ind])

  ## trim probabilities
  dat$pihat <- pmin(1-1E-4, pmax(1E-4, dat$pihat))


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



## proposed AIPW estimator with xgboost and sample splitting
aipw_fit_nss_v2 <- function(dat, mu_covariates, pi_covariates, K,
                            method, arguments = NULL) {

  dat <- dat %>%
    group_by(S, A) %>%
    mutate(fold = sample(rep(1:K, length.out = n()))) %>%
    ungroup() %>%
    mutate(Q = S * R * (1 - C))

  n <- sum(dat$wt)

  ## initialize results
  etahats <- data.frame(
    fold = 1:K,
    etahat_0 = NA, etahat_1 = NA, cov_00 = NA, cov_01 = NA, cov_11 = NA)

  ## loop through K folds
  for (k in 1:K) {

    ## data sets of fold k and complement
    dat_k <- dat %>% filter(fold == k)
    dat_ck <- dat %>% filter(fold != k)

    ## train outcome regression function mu on fold ck
    mu_mod <- nonpar_est(
      x = dat_ck %>% filter(Q == 1) %>% select(all_of(c("A", mu_covariates))),
      y = dat_ck %>% filter(Q == 1) %>% select(Y) %>% unlist() %>% factor(),
      method = method,
      arguments = arguments)

    ## predict probabilities in fold k with with A set to 0 and 1
    mu_k_ind <- dat_k$S == 0 | dat_k$R == 1
    dat_k$muhat_0 <- dat_k$muhat_1 <- NA
    dat_k$muhat_0[mu_k_ind] <- nonpar_pred(
      mod = mu_mod,
      newdata = dat_k[mu_k_ind,] %>%
        mutate(A = 0) %>%
        select(all_of(c("A", mu_covariates))),
      method = method)
    dat_k$muhat_1[mu_k_ind] <- nonpar_pred(
      mod = mu_mod,
      newdata = dat_k[mu_k_ind,] %>%
        mutate(A = 1) %>%
        select(all_of(c("A", mu_covariates))),
      method = method)

    ## train Q predictor on fold ck
    Q0_ck <- dat_ck %>% filter((Q == 1 & A == 0) | S == 0)
    Q1_ck <- dat_ck %>% filter((Q == 1 & A == 1) | S == 0)

    Q_reg_0 <- nonpar_est(
      x = Q0_ck %>% select(all_of(pi_covariates)),
      y = Q0_ck$Q %>% factor(),
      wts = Q0_ck$wt,
      method = method,
      arguments = arguments)

    Q_reg_1 <- nonpar_est(
      x = Q1_ck %>% select(all_of(pi_covariates)),
      y = Q1_ck$Q %>% factor(),
      wts = Q1_ck$wt,
      method = method,
      arguments = arguments)

    ## predicted Q probabilities among uncensored responders in fold k
    Q0_k_ind <- dat_k$Q == 1 & dat_k$A == 0
    Q1_k_ind <- dat_k$Q == 1 & dat_k$A == 1
    dat_k$Q_prob <- NA
    dat_k$Q_prob[Q0_k_ind] <- nonpar_pred(
      mod = Q_reg_0,
      newdata = dat_k[Q0_k_ind,] %>% select(all_of(pi_covariates)),
      method = method)
    dat_k$Q_prob[Q1_k_ind] <- nonpar_pred(
      mod = Q_reg_1,
      newdata = dat_k[Q1_k_ind,] %>% select(all_of(pi_covariates)),
      method = method)

    ## estimated propensity scores in fold k trial data
    Q_k_ind <- dat_k$Q == 1
    dat_k$pihat <- NA_real_
    dat_k$pihat[Q_k_ind] <- dat_k$Q_prob[Q_k_ind] /
      (1 - dat_k$Q_prob[Q_k_ind])

    ## Hajek denominator in fold k
    n_trial_hat_k <- dat_k %>%
      filter(Q == 1) %>%
      mutate(
        term_0 = (1 - A) / pihat,
        term_1 = A / pihat) %>%
      summarise(
        nhat_0 = sum(term_0),
        nhat_1 = sum(term_1)) %>%
      unlist()

    ## auxiliary sample size in fold k
    n_aux_k <- sum(dat_k$wt[dat_k$S == 0])

    ## total sample size in fold k
    n_k <- sum(dat_k$wt)

    ## AIPW estimator in fold k
    etahat_k <- dat_k %>%
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
        if_0 = (S * R * (1 - C) * ipw_0 / n_trial_hat_k[1] +
                  (1 - S) * or_0 / n_aux_k) * n_k,
        if_1 = (S * R * (1 - C) * ipw_1 / n_trial_hat_k[2] +
                  (1 - S) * or_1 / n_aux_k) * n_k) %>%

      summarise(

        ## AIPW estimator
        etahat_0 = mean(if_0),
        etahat_1 = mean(if_1),

        ## covariance estimator
        cov_00 = var(if_0) / n,
        cov_01 = cov(if_0, if_1) / n,
        cov_11 = var(if_1) / n)

    etahats[k, 2:6] <- etahat_k

    print(paste0("fold ", k, "/", K, " complete"))

  }

  ## aggregate across k folds
  etahat <- colMeans(etahats[,2:6])


  # return list of results --------------------------------------------------

  res <- list(

    # causal parameter estimates and covariance
    eta_results =
      data.frame(
        etahat_0 = etahat[1],
        etahat_1 = etahat[2],
        cov_00 = etahat[3],
        cov_01 = etahat[4],
        cov_11 = etahat[5]),

    # data set
    dat = dat)

  return(res)

}
