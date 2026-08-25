#' Crown AIPW NSS1 estimator
#'
#' estimate eta(0) and eta(1) using a crown AIPW estimator, i.e., accounting for
#' nonresponse, with nonparametric sample splitting, option 1
#'
#' @param dat data frame containing the following columns:
#' \itemize{
#'\item `S`: a binary indicator for whether the observation belongs to the trial
#'data (`S`=1) the auxiliary data (`S`=0)
#' \item `R`: a binary indicator for whether the observation is a responder
#' (`R`=1) or not (`R`=0)
#' \item `C`: a binary indicator for whether the outcome is censored (`C`=1) or
#' not (`C`=0)
#' \item `A`: a binary exposure
#' \item `Y`: a binary outcome
#' \item other covariates specified in `mu_fmla`
#' \item `wt`: survey weights
#' }
#'
#' @param mu_covariates a character vector of covariate names in `dat` to use
#' in the outcome regression model
#'
#' @param pi_covariates a character vector of covariate names in `dat` to use
#' in the propensity score regression model
#'
#' @param K a positive integer number of folds for sample splitting
#'
#' @param method a character string, method of nonparametric estimator
#'
#' @param arguments an optional list with arguments to pass to `nonpar_est`
#'
#' @return a list containing the following:
#' \itemize{
#' \item `eta_hat`: a numeric vector estimated mean potential outcomes eta(0)
#' and eta(1)
#' \item `eta_hat_cov`: a numeric matrix, estimated covariance of `eta_hat`
#' \item `dat`: a data frame with estimated propensity score weights
#' }
#'
#' @export
aipw_fit_nss_v1 <- function(dat, mu_covariates, pi_covariates, K,
                            method, arguments) {


  ## check input
  stopifnot(
    "dat must contain columns S, R, C, A, Y, wt" =
      all(c("S", "R", "C", "A", "Y", "wt") %in% names(dat)),
    "S must be binary (0/1)" = all(dat$S %in% c(0, 1)),
    "R must be binary (0/1)" = all(dat$R %in% c(0, 1)),
    "C must be binary (0/1)" = all(dat$C %in% c(0, 1)),
    "A must be binary (0/1)" = all(dat$A %in% c(0, 1)),
    "Y must be binary (0/1)" = all(dat$Y %in% c(0, 1)),
    "mu_covariates must be columns of dat" = all(mu_covariates %in% names(dat)),
    "pi_covariates must be columns of dat" = all(pi_covariates %in% names(dat)),
    "K must be a positive integer" = K > 0 && K == round(K))


  ## perform sample splitting
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

    ## train outcome regression function mu with forced A interactions
    train_outcome <- train %>%
      filter(Q == 1)
    train_outcome_0 <- train_outcome %>% filter(A == 0)
    train_outcome_1 <- train_outcome %>% filter(A == 1)

    mu_mod_0 <- nonpar_est(
      x = train_outcome_0 %>% select(all_of(mu_covariates)),
      y = train_outcome_0 %>% select(Y) %>% unlist(),
      method = method,
      arguments = arguments)

    mu_mod_1 <- nonpar_est(
      x = train_outcome_1 %>% select(all_of(mu_covariates)),
      y = train_outcome_1 %>% select(Y) %>% unlist(),
      method = method,
      arguments = arguments)

    ## predicted probabilities in test data set with A = 0 and A = 1
    test_mu_idx <- dat$fold == k & (dat$S == 0 | dat$R == 1)
    test_0 <- dat[test_mu_idx,] %>% mutate(A = 0)
    test_1 <- dat[test_mu_idx,] %>% mutate(A = 1)

    dat$muhat_0[test_mu_idx] <- nonpar_pred(
      mod = mu_mod_0,
      newdata = test_0 %>% select(all_of(mu_covariates)),
      method = method)

    dat$muhat_1[test_mu_idx] <- nonpar_pred(
      mod = mu_mod_1,
      newdata = test_1 %>% select(all_of(mu_covariates)),
      method = method)

    ## train Q predictor
    train_Q0 <- train %>% filter((Q == 1 & A == 0) | S == 0)
    train_Q1 <- train %>% filter((Q == 1 & A == 1) | S == 0)

    Q_reg_0 <- nonpar_est(
      x = train_Q0 %>% select(all_of(pi_covariates)),
      y = train_Q0$Q,
      wts = train_Q0$wt,
      method = method,
      arguments = arguments)

    Q_reg_1 <- nonpar_est(
      x = train_Q1 %>% select(all_of(pi_covariates)),
      y = train_Q1$Q,
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

  ## Hajek estimator of trial sample size
  n_trial_hat <- dat %>%
    filter(Q == 1) %>%
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
      ipw_0 = ifelse(Q == 1 & A == 0,
                     (Y - muhat_0) / pihat,
                     0),

      ipw_1 = ifelse(Q == 1 & A == 1,
                     (Y - muhat_1) / pihat,
                     0),

      ## outcome regression terms
      or_0 = ifelse(S == 0, muhat_0 * wt, 0),
      or_1 = ifelse(S == 0, muhat_1 * wt, 0),

      ## influence function
      if_0 = (Q * ipw_0 / n_trial_hat[1] +
                (1 - S) * or_0 / n_aux) * nrow(dat),
      if_1 = (Q * ipw_1 / n_trial_hat[2] +
                (1 - S) * or_1 / n_aux) * nrow(dat)) %>%

    summarise(

      ## AIPW estimator
      etahat_0 = mean(if_0),
      etahat_1 = mean(if_1),

      ## covariance estimator
      cov_00 = var(if_0) / nrow(dat),
      cov_01 = cov(if_0, if_1) / nrow(dat),
      cov_11 = var(if_1) / nrow(dat))

  ## estimated variance
  est_var = matrix(
    data = c(etahat[3], etahat[4], etahat[4], etahat[5]),
    nrow = 2, ncol = 2, byrow = 2)


  # return list of results --------------------------------------------------

  res <- list(

    ## causal parameter estimate
    eta_hat = c("etahat_0" = etahat$etahat_0,
                "etahat_1" = etahat$etahat_1),

    ## estimated covariance of eta_hat
    eta_hat_cov = matrix(
      c(etahat$cov_00, etahat$cov_01,
        etahat$cov_01, etahat$cov_11),
      nrow = 2, ncol = 2, byrow = TRUE),

    ## data set
    dat = dat)

  return(res)

}



#' Crown AIPW NSS2 estimator
#'
#' estimate eta(0) and eta(1) using a crown AIPW estimator, i.e., accounting for
#' nonresponse, with nonparametric sample splitting, option 2
#'
#' @param dat data frame containing the following columns:
#' \itemize{
#'\item `S`: a binary indicator for whether the observation belongs to the trial
#'data (`S`=1) the auxiliary data (`S`=0)
#' \item `R`: a binary indicator for whether the observation is a responder
#' (`R`=1) or not (`R`=0)
#' \item `C`: a binary indicator for whether the outcome is censored (`C`=1) or
#' not (`C`=0)
#' \item `A`: a binary exposure
#' \item `Y`: a binary outcome
#' \item other covariates specified in `mu_fmla`
#' \item `wt`: survey weights
#' }
#'
#' @param mu_covariates a formula for the outcome regression model using
#' variables in `dat`
#'
#' @param pi_covariates a formula for the propensity score regression model
#' using variables in `dat`
#'
#' @param K a positive integer number of folds for sample splitting
#'
#' @param method a character string, method of nonparametric estimator
#'
#' @param arguments an optional list with arguments to pass to `nonpar_est`
#'
#' @return a list containing the following:
#' \itemize{
#' \item `eta_hat`: a numeric vector estimated mean potential outcomes eta(0)
#' and eta(1)
#' \item `eta_hat_cov`: a numeric matrix, estimated covariance of `eta_hat`
#' \item `dat`: a data frame with estimated propensity score weights
#' }
#'
#' @export
aipw_fit_nss_v2 <- function(dat, mu_covariates, pi_covariates, K,
                            method, arguments = NULL) {

  ## check input
  stopifnot(
    "dat must contain columns S, R, C, A, Y, wt" =
      all(c("S", "R", "C", "A", "Y", "wt") %in% names(dat)),
    "S must be binary (0/1)" = all(dat$S %in% c(0, 1)),
    "R must be binary (0/1)" = all(dat$R %in% c(0, 1)),
    "C must be binary (0/1)" = all(dat$C %in% c(0, 1)),
    "A must be binary (0/1)" = all(dat$A %in% c(0, 1)),
    "Y must be binary (0/1)" = all(dat$Y %in% c(0, 1)),
    "mu_covariates must be columns of dat" = all(mu_covariates %in% names(dat)),
    "pi_covariates must be columns of dat" = all(pi_covariates %in% names(dat)),
    "K must be a positive integer" = K > 0 && K == round(K))


  ## perform sample splitting
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

    ## train outcome regression function mu on fold ck with forces A interaction
    mu_mod_0 <- nonpar_est(
      x = dat_ck %>% filter(Q == 1, A == 0) %>% select(all_of(mu_covariates)),
      y = dat_ck %>% filter(Q == 1, A == 0) %>% select(Y) %>% unlist(),
      method = method,
      arguments = arguments)

    mu_mod_1 <- nonpar_est(
      x = dat_ck %>% filter(Q == 1, A == 1) %>% select(all_of(mu_covariates)),
      y = dat_ck %>% filter(Q == 1, A == 1) %>% select(Y) %>% unlist(),
      method = method,
      arguments = arguments)

    ## predict probabilities in fold k with with A set to 0 and 1
    mu_k_ind <- dat_k$S == 0 | dat_k$R == 1
    dat_k$muhat_0 <- dat_k$muhat_1 <- NA
    dat_k$muhat_0[mu_k_ind] <- nonpar_pred(
      mod = mu_mod_0,
      newdata = dat_k[mu_k_ind,] %>%
        select(all_of(mu_covariates)),
      method = method)
    dat_k$muhat_1[mu_k_ind] <- nonpar_pred(
      mod = mu_mod_1,
      newdata = dat_k[mu_k_ind,] %>%
        select(all_of(mu_covariates)),
      method = method)

    ## train Q predictor on fold ck
    Q0_ck <- dat_ck %>% filter((Q == 1 & A == 0) | S == 0)
    Q1_ck <- dat_ck %>% filter((Q == 1 & A == 1) | S == 0)

    Q_reg_0 <- nonpar_est(
      x = Q0_ck %>% select(all_of(pi_covariates)),
      y = Q0_ck$Q,
      wts = Q0_ck$wt,
      method = method,
      arguments = arguments)

    Q_reg_1 <- nonpar_est(
      x = Q1_ck %>% select(all_of(pi_covariates)),
      y = Q1_ck$Q,
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

    ## Hajek estimator of sample size in fold k
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
        ipw_0 = ifelse(Q == 1 & A == 0,
                       (Y - muhat_0) / pihat,
                       0),

        ipw_1 = ifelse(Q == 1 & A == 1,
                       (Y - muhat_1) / pihat,
                       0),

        ## outcome regression terms
        or_0 = ifelse(S == 0, muhat_0 * wt, 0),
        or_1 = ifelse(S == 0, muhat_1 * wt, 0),

        ## influence function
        if_0 = (Q * ipw_0 / n_trial_hat_k[1] +
                  (1 - S) * or_0 / n_aux_k) * n_k,
        if_1 = (Q * ipw_1 / n_trial_hat_k[2] +
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

  ## estimated variance
  est_var = matrix(
    c(etahat["cov_00"], etahat["cov_01"],
      etahat["cov_01"], etahat["cov_11"]),
    nrow = 2, ncol = 2, byrow = TRUE)

  # return list of results --------------------------------------------------

  res <- list(

    ## causal parameter estimate
    eta_hat = c("etahat_0" = unname(etahat["etahat_0"]),
                "etahat_1" = unname(etahat["etahat_1"])),

    ## estimated covariance of eta_hat
    eta_hat_cov = est_var,

    ## data set
    dat = dat)

  return(res)

}
