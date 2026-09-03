#' Crown DML estimator
#'
#' estimate eta(0) and eta(1) using a crown DML estimator, i.e., accounting for
#' nonresponse, with nonparametric sample splitting
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
dml_fit <- function(dat, mu_covariates, pi_covariates, K,
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
