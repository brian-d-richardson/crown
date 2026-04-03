
## naive AIPW estimator
aipw_fit_naive_nonpar <- function(dat, mu_covariates, pi_covariates, pA = 0.5) {

  # fit outcome regression --------------------------------------------------

  ## random forest among uncensored responders
  outcome_reg <- ranger(
    x = filter(dat, R == 1, C == 0) %>%
      select(all_of(c("A", mu_covariates))),
    y = filter(dat, R == 1, C == 0) %>%
      select(Y) %>% unlist() %>% factor(),
    probability = T)
  
  ## data sets with A set to 0, 1
  dat0 <- dat %>% mutate(A = 0, Y = 0)
  dat1 <- dat %>% mutate(A = 1, Y = 0)
  
  ## predict outcomes under A set to 0, 1
  dat$muhat_0 <- NA_real_
  dat$muhat_0[dat$S == 0 | dat$R == 1] <- predict(
    outcome_reg,
    data = dat0 %>% 
      select(all_of(c("A", mu_covariates))))$predictions[,2]
  
  dat$muhat_1 <- NA_real_
  dat$muhat_1[dat$S == 0 | dat$R == 1] <- predict(
    outcome_reg,
    data = dat1 %>% 
      select(all_of(c("A", mu_covariates))))$predictions[,2]


  # censoring model ---------------------------------------------------------

  ## logistic regression
  censor_reg <- ranger(
    x = dat %>%
      select(all_of(pi_covariates)),
    y = dat %>%
      select(C) %>% unlist() %>% factor(),
    probability = T)

  ## predicted values
  dat$piC <- NA_real_
  dat$piC <-
    1 - predict(
      censor_reg,
      data = dat %>% 
        select(all_of(pi_covariates)))$predictions[,2]

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
      or_1 = muhat_1 * wt) %>%

    summarise(

      ## AIPW estimators
      etahat_0 = sum((ipw_0 / n_trial_hat[1]) +
                       (or_0 / n_aux)),
      etahat_1 = sum((ipw_1 / n_trial_hat[2]) +
                       (or_1 / n_aux)))


  # variance estimator ------------------------------------------------------

  ## variance estimator
  est_var <- get.sand.est(
    param = c(etahat$etahat_0, etahat$etahat_1),
    n = nrow(dat),
    get.psi = function(xx) {

      ## extract pieces of combined parameter
      eh0 <- xx[1]
      eh1 <- xx[2]

      ## stacked estimating function
      cbind(

        ## etahat estimating function
        case_when(
          dat$C == 0 & dat$A == 0 ~
            (nrow(dat) / n_trial_hat[1]) * (dat$Y - dat$muhat_0) / dat$pihat,
          .default = 0) +
          (nrow(dat) / n_aux) * dat$wt * dat$muhat_0 -
          eh0,

        case_when(
          dat$C == 0 & dat$A == 1 ~
            (nrow(dat) / n_trial_hat[2]) * (dat$Y - dat$muhat_1) / dat$pihat,
          .default = 0) +
          (nrow(dat) / n_aux) * dat$wt * dat$muhat_1 -
          eh1)

    }
  )


  # return list of results --------------------------------------------------

  res <- list(

    # causal parameter estimates and covariance
    eta_results =
      data.frame(
        etahat_0 = etahat$etahat_0,
        etahat_1 = etahat$etahat_1,
        cov_00 = est_var[1, 1],
        cov_01 = est_var[1, 2],
        cov_11 = est_var[2, 2]),


    # outcome regression model results
    outcome_reg_results = outcome_reg,

    # censoring regression model results
    censor_reg_results = censor_reg,

    # data set with weights
    dat = dat)

  return(res)

}



## proposed AIPW estimator allowing A --> R effect
aipw_fit_nonpar <- function(dat, mu_covariates, pi_covariates) {

  # fit outcome regression --------------------------------------------------

  ## random forest among uncensored responders
  outcome_reg <- ranger(
    x = filter(dat, S == 1, R == 1, C == 0) %>%
      select(all_of(c("A", mu_covariates))),
    y = filter(dat, S == 1, R == 1, C == 0) %>%
      select(Y) %>% unlist() %>% factor(),
    probability = T)

  ## data sets with A set to 0, 1
  dat0 <- dat %>% mutate(A = 0, Y = 0)
  dat1 <- dat %>% mutate(A = 1, Y = 0)

  ## predict outcomes under A set to 0, 1
  dat$muhat_0 <- NA_real_
  dat$muhat_0[dat$S == 0 | dat$R == 1] <- predict(
    outcome_reg,
    data = dat0 %>% filter(dat$S == 0 | dat$R == 1) %>%
      select(all_of(c("A", mu_covariates))))$predictions[,2]

  dat$muhat_1 <- NA_real_
  dat$muhat_1[dat$S == 0 | dat$R == 1] <- predict(
    outcome_reg,
    data = dat1 %>% filter(dat$S == 0 | dat$R == 1) %>%
      select(all_of(c("A", mu_covariates))))$predictions[,2]


  # response model ----------------------------------------------------------

  ## add Q labels to data
  dat <- dat %>%
    mutate(Q = S * R * (1 - C))

  ## restrict sample to SR(1-C) = 1 or S = 0, and add Q labels
  rdat0 <- dat %>% filter((Q == 1 & A == 0) | S == 0)
  rdat1 <- dat %>% filter((Q == 1 & A == 1) | S == 0)

  ## random forest models for Q membership in restricted data
  Q_reg_0 <- ranger(
    x = rdat0 %>% select(all_of(pi_covariates)),
    y = rdat0 %>% select(Q) %>% unlist() %>% factor(),
    case.weights = rdat0$wt,
    probability = T)

  Q_reg_1 <- ranger(
    x = rdat1 %>% select(all_of(pi_covariates)),
    y = rdat1 %>% select(Q) %>% unlist() %>% factor(),
    case.weights = rdat1$wt,
    probability = T)

  ## estimated Q probabilities among uncensored responders
  Q_ind <- dat$Q == 1
  dat$Q_prob <- NA_real_
  dat$Q_prob[Q_ind & dat$A == 0] <-
    predict(
      Q_reg_0,
      data = dat %>% filter(Q_ind & A == 0) %>%
        select(all_of(pi_covariates)))$predictions[,2]

  dat$Q_prob[Q_ind & dat$A == 1] <-
    predict(
      Q_reg_1,
      data = dat %>% filter(Q_ind & A == 1) %>%
        select(all_of(pi_covariates)))$predictions[,2]

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
      or_1 = ifelse(S == 0, muhat_1 * wt, 0)) %>%

    summarise(

      ## AIPW estimators
      etahat_0 = sum((ipw_0 / n_trial_hat[1]) +
                       (or_0 / n_aux)),
      etahat_1 = sum((ipw_1 / n_trial_hat[2]) +
                       (or_1 / n_aux)))


  # variance estimator ------------------------------------------------------

  ## variance estimator
  est_var <- get.sand.est(
    param = c(etahat$etahat_0, etahat$etahat_1),
    n = nrow(dat),
    get.psi = function(xx) {

      ## extract pieces of combined parameter
      eh0 <- xx[1]
      eh1 <- xx[2]

      ## estimating function
      cbind(

        ## etahat estimating function
        case_when(
          dat$Q == 1 & dat$A == 0 ~
            (nrow(dat) / n_trial_hat[1]) * (dat$Y - dat$muhat_0) / dat$pihat,
          dat$S == 0 ~
            (nrow(dat) / n_aux) * dat$wt * dat$muhat_0,
          .default = 0) - eh0,

        case_when(
          dat$Q == 1 & dat$A == 1 ~
            (nrow(dat) / n_trial_hat[2]) * (dat$Y - dat$muhat_1) / dat$pihat,
          dat$S == 0 ~
            (nrow(dat) / n_aux) * dat$wt * dat$muhat_1,
          .default = 0) - eh1)

    }
  )


  # return list of results --------------------------------------------------

  res <- list(

    # causal parameter estimates and covariance
    eta_results =
      data.frame(
        etahat_0 = etahat$etahat_0,
        etahat_1 = etahat$etahat_1,
        cov_00 = est_var[1, 1],
        cov_01 = est_var[1, 2],
        cov_11 = est_var[2, 2]),

    # outcome regression model results
    outcome_reg_results = outcome_reg,

    # Q regression model results
    Q_reg_0_results = Q_reg_0,
    Q_reg_1_results = Q_reg_1,

    # data set with weights
    dat = dat)

  return(res)

}
