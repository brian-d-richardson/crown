
## naive AIPW estimator
aipw_fit_naive_nonpar <- function(dat, mu_covariates, pi_covariates, pA = 0.5,
                                  smoothness_orders = 1,
                                  max_degree = 2,
                                  num_knots = c(25, 10)) {

  # fit outcome regression --------------------------------------------------

  ## random forest among uncensored responders
  outcome_reg <- fit_hal(
    X = filter(dat, R == 1, C == 0) %>%
      select(all_of(c("A", mu_covariates))),
    Y = filter(dat, R == 1, C == 0) %>%
      select(Y) %>% unlist(),
    family = "binomial",
    smoothness_orders = smoothness_orders,
    max_degree = max_degree,
    num_knots = num_knots)

  ## data sets with A set to 0, 1
  dat0 <- dat %>% mutate(A = 0, Y = 0)
  dat1 <- dat %>% mutate(A = 1, Y = 0)

  ## predict outcomes under A set to 0, 1
  dat$muhat_0 <- NA_real_
  dat$muhat_0[dat$S == 0 | dat$R == 1] <- predict(
    outcome_reg,
    new_data = dat0 %>%
      select(all_of(c("A", mu_covariates))))

  dat$muhat_1 <- NA_real_
  dat$muhat_1[dat$S == 0 | dat$R == 1] <- predict(
    outcome_reg,
    new_data = dat1 %>%
      select(all_of(c("A", mu_covariates))))


  # censoring model ---------------------------------------------------------

  ## logistic regression
  censor_reg <- fit_hal(
    X = dat %>%
      select(all_of(pi_covariates)),
    Y = dat %>%
      select(C) %>% unlist(),
    family = "binomial",
    smoothness_orders = smoothness_orders,
    max_degree = max_degree,
    num_knots = num_knots)

  ## predicted values
  dat$piC <- NA_real_
  dat$piC <-
    1 - predict(
      censor_reg,
      new_data = dat %>%
        select(all_of(pi_covariates)))

  ## estimated arm * response probabilities in trial data
  dat$piA <- ifelse(dat$A == 1, pA, 1 - pA)

  ## joint probabilities
  dat$pihat <- dat$piA * dat$piC


  # AIPW estimator ----------------------------------------------------------

  ## Hajek denominator
  n_trial_hat <- dat %>%
    filter(C == 0) %>%
    mutate(
      term_0 = (1 - A) / pihat,
      term_1 = A / pihat) %>%
    summarise(
      nhat_0 = sum(term_0),
      nhat_1 = sum(term_1)) %>%
    unlist()

  ## auxiliary sample size
  n_aux <- sum(dat$wt)

  ## AIPW estimator
  etahat <- dat %>%
    mutate(

      ## IPW terms
      ipw_0 = ifelse(C == 0 & A == 0,
                     (Y - muhat_0) / pihat,
                     0),

      ipw_1 = ifelse(C == 0 & A == 1,
                     (Y - muhat_1) / pihat,
                     0),

      ## outcome regression terms
      or_0 = muhat_0 * wt,
      or_1 = muhat_1 * wt,

      ## influence function
      if_0 = ((1 - C) * ipw_0 / n_trial_hat[1] +
                or_0 / n_aux) * nrow(dat),
      if_1 = ((1 - C) * ipw_1 / n_trial_hat[2] +
                or_1 / n_aux) * nrow(dat)) %>%

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

    # outcome regression model results
    outcome_reg_results = outcome_reg,

    # censoring regression model results
    censor_reg_results = censor_reg,

    # data set with weights
    dat = dat)

  return(res)

}



## proposed AIPW estimator with HAL
aipw_fit_nonpar <- function(dat, mu_covariates, pi_covariates,
                            smoothness_orders = 1,
                            max_degree = 2,
                            num_knots = c(25, 10)) {

  # fit outcome regression --------------------------------------------------

  ## highly adaptive LASSO among uncensored responders
  #tic("fit outcome regression")
  outcome_reg <- fit_hal(
    X = filter(dat, S == 1, R == 1, C == 0) %>%
      select(all_of(c("A", mu_covariates))),
    Y = filter(dat, S == 1, R == 1, C == 0) %>%
      select(Y) %>% unlist(),
    family = "binomial",
    smoothness_orders = smoothness_orders,
    max_degree = max_degree,
    num_knots = num_knots)
  #toc()

  ## data sets with A set to 0, 1
  dat0 <- dat %>% mutate(A = 0, Y = 0)
  dat1 <- dat %>% mutate(A = 1, Y = 0)

  ## predict outcomes under A set to 0, 1
  dat$muhat_0 <- NA_real_
  dat$muhat_0[dat$S == 0 | dat$R == 1] <- predict(
    outcome_reg,
    new_data = dat0 %>% filter(dat$S == 0 | dat$R == 1) %>%
      select(all_of(c("A", mu_covariates))))

  dat$muhat_1 <- NA_real_
  dat$muhat_1[dat$S == 0 | dat$R == 1] <- predict(
    outcome_reg,
    new_data = dat1 %>% filter(dat$S == 0 | dat$R == 1) %>%
      select(all_of(c("A", mu_covariates))))


  # response model ----------------------------------------------------------

  ## add Q labels to data
  dat <- dat %>%
    mutate(Q = S * R * (1 - C))

  ## restrict sample to SR(1-C) = 1 or S = 0, and add Q labels
  rdat0 <- dat %>% filter((Q == 1 & A == 0) | S == 0)
  rdat1 <- dat %>% filter((Q == 1 & A == 1) | S == 0)

  ## random forest models for Q membership in restricted data
  #tic("fit response model 0")
  Q_reg_0 <- fit_hal(
    X = rdat0 %>% select(all_of(pi_covariates)),
    Y = rdat0 %>% select(Q) %>% unlist(),
    weights = rdat0$wt,
    family = "binomial",
    smoothness_orders = smoothness_orders,
    max_degree = max_degree,
    num_knots = num_knots)
  #toc()

  #tic("fit response model 1")
  Q_reg_1 <- fit_hal(
    X = rdat1 %>% select(all_of(pi_covariates)),
    Y = rdat1 %>% select(Q) %>% unlist(),
    weights = rdat1$wt,
    family = "binomial",
    smoothness_orders = smoothness_orders,
    max_degree = max_degree,
    num_knots = num_knots)
  #toc()

  ## estimated Q probabilities among uncensored responders
  Q_ind <- dat$Q == 1
  dat$Q_prob <- NA_real_
  dat$Q_prob[Q_ind & dat$A == 0] <-
    predict(
      Q_reg_0,
      new_data = dat %>% filter(Q_ind & A == 0) %>%
        select(all_of(pi_covariates)))

  dat$Q_prob[Q_ind & dat$A == 1] <-
    predict(
      Q_reg_1,
      new_data = dat %>% filter(Q_ind & A == 1) %>%
        select(all_of(pi_covariates)))

  ## estimated propensity scores in trial data
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

    # outcome regression model results
    outcome_reg_results = outcome_reg,

    # Q regression model results
    Q_reg_0_results = Q_reg_0,
    Q_reg_1_results = Q_reg_1,

    # data set with weights
    dat = dat)

  return(res)

}
