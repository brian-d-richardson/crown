#' Naive AIPW estimator
#'
#' estimate eta(0) and eta(1) using a naive AIPW estimator, i.e., ignoring
#' nonresponse
#'
#' @param dat data frame containing the following columns:
#' \itemize{
#' \item `C`: a binary indicator for whether the outcome is censored (`C`=1) or
#' not (`C`=0)
#' \item `A`: a binary exposure
#' \item `Y`: a binary outcome
#' \item other covariates specified in `mu_fmla`
#' \item `wt`: survey weights
#' }
#'
#' @param C_fmla a formula for the censoring mechanism regression model using
#' variables in `dat`
#'
#' @param mu_fmla a formula for the outcome regression model using variables in
#' `dat`
#'
#' @param pA an optional number in (0, 1), the marginal probability of
#' treatment (A = 1), default is 0.5.
#'
#' @return a list containing the following:
#' \itemize{
#' \item `eta_hat`: a numeric vector estimated mean potential outcomes eta(0)
#' and eta(1)
#' \item `eta_hat_cov`: a numeric matrix, estimated covariance of `eta_hat`
#' \item `outcome_reg`: a list, results of outcome regression model
#' \item `censor_reg`: a list, results of censoring regression model
#' \item `dat`: a data frame including a column for propensity score weights
#' }
#'
#' @export
aipw_fit_naive <- function(dat, mu_fmla, C_fmla, pA = 0.5) {


  # check input -------------------------------------------------------------

  ## required columns present
  stopifnot(
    "dat must contain columns C, A, Y, wt" =
      all(c("C", "A", "Y", "wt") %in% names(dat)))

  ## required columns are binary (0/1)
  stopifnot(
    "C must be binary (0/1)" = all(dat$C %in% c(0, 1)),
    "A must be binary (0/1)" = all(dat$A %in% c(0, 1)),
    "Y must be binary (0/1)" = all(dat$Y %in% c(0, 1)))


  # fit outcome regression --------------------------------------------------

  ## logistic regression among uncensored responders
  outcome_reg <- glm(
    formula = mu_fmla,
    family = "binomial",
    data = filter(dat, C == 0),
    weights = wt)

  ## data sets with A set to 0, 1
  dat0 <- dat %>% mutate(A = 0, Y = 0)
  dat1 <- dat %>% mutate(A = 1, Y = 0)

  ## predict outcomes under A set to 0, 1
  dat$muhat_0 <- predict(
    outcome_reg,
    newdata = dat0,
    type = "response")

  dat$muhat_1 <- predict(
    outcome_reg,
    newdata = dat1,
    type = "response")


  # censoring model ---------------------------------------------------------

  ## logistic regression model for censoring mechanism
  censor_reg <- glm(
    formula = C_fmla,
    family = "binomial",
    data = dat,
    weights = wt)

  ## predicted values from censoring model
  dat$piC <-
    1 - predict(
      censor_reg,
      newdata = dat,
      type = "response")

  ## treatment assignment probabilities
  dat$piA <- ifelse(dat$A == 1, pA, 1 - pA)


  # AIPW estimator ----------------------------------------------------------

  ## joint probabilities of treatment assignment and uncensored status
  dat$pihat <- dat$piA * dat$piC

  ## Hajek estimator of trial sample size
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

  ## design matrices for variance estimation
  X0 <- model.matrix(outcome_reg, data = dat0)
  X1 <- model.matrix(outcome_reg, data = dat1)
  X_cens <- model.matrix(censor_reg, data = dat)

  ## variance estimator
  est_var <- get.sand.est(
    param = c(etahat$etahat_0, etahat$etahat_1,
              coef(outcome_reg),
              coef(censor_reg)),
    n = nrow(dat),
    get.psi = function(xx) {

      ## extract pieces of combined parameter
      eh0 <- xx[1]
      eh1 <- xx[2]
      bb_outcome <- xx[2 + 1:length(coef(outcome_reg))]
      bb_cens <- xx[2 + length(coef(outcome_reg)) +
                      1:length(coef(censor_reg))]

      ## re-create predicted outcomes using supplied parameter
      muhat_0_ <- plogis(as.vector(X0 %*% bb_outcome))
      muhat_1_ <- plogis(as.vector(X1 %*% bb_outcome))

      ## re-create weights using supplied parameter
      piC_ <- 1 - plogis(as.vector(X_cens %*% bb_cens))
      pihat_ <- piC_ * dat$piA

      ## recreate Hajek denominator
      n_trial_hat_ <- dat %>%
        mutate(pihat = pihat_) %>%
        filter(C == 0) %>%
        mutate(
          term_0 = (1 - A) / pihat,
          term_1 = A / pihat) %>%
        summarise(
          nhat_0 = sum(term_0),
          nhat_1 = sum(term_1)) %>%
        unlist()

      ## stacked estimating function
      cbind(

        ## outcome regression estimating function
        psi.lr(data = dat,
               beta = bb_outcome,
               formula = mu_fmla) *
          (1 - dat$C),

        ## censoring regression estimating function
        psi.lr(data = dat,
               beta = bb_cens,
               formula = C_fmla),

        ## etahat estimating function
        case_when(
          dat$C == 0 & dat$A == 0 ~
            (nrow(dat) / n_trial_hat_[1]) * (dat$Y - muhat_0_) / pihat_,
          .default = 0) +
          (nrow(dat) / n_aux) * dat$wt * muhat_0_ -
          eh0,

        case_when(
          dat$C == 0 & dat$A == 1 ~
            (nrow(dat) / n_trial_hat_[2]) * (dat$Y - muhat_1_) / pihat_,
          .default = 0) +
          (nrow(dat) / n_aux) * dat$wt * muhat_1_ -
          eh1)

    }
  )


  # return list of results --------------------------------------------------

  res <- list(

    ## causal parameter estimate
    eta_hat = c("etahat_0" = etahat$etahat_0,
                "etahat_1" = etahat$etahat_1),

    ## estimated covariance of eta_hat
    eta_hat_cov = est_var[1:2, 1:2],

    ## outcome regression model results
    outcome_reg = outcome_reg,

    ## censoring regression model results
    censor_reg = censor_reg,

    ## data set including estimated propensity score weights
    dat = dat)

  return(res)

}



#' Crown AIPW estimator
#'
#' estimate eta(0) and eta(1) using a crown AIPW estimator, i.e., accounting for
#' nonresponse
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
#' @param mu_fmla a formula for the outcome regression model using variables in
#' `dat`
#'
#'
#' @param pi_fmla a formula for the propensity score regression model using
#' variables in `dat`
#'
#' @return a list containing the following:
#' \itemize{
#' \item `eta_hat`: a numeric vector estimated mean potential outcomes eta(0)
#' and eta(1)
#' \item `eta_hat_cov`: a numeric matrix, estimated covariance of `eta_hat`
#' \item `outcome_reg`: a list, results of outcome regression model
#' \item `Q_reg_0`: a list, results of Q0 membership regression model
#' \item `Q_reg_1`: a list, results of Q1 membership regression model
#' \item `dat`: a data frame including a column for propensity score weights
#' }
#'
#' @export
aipw_fit <- function(dat, mu_fmla, pi_fmla) {


  # check input -------------------------------------------------------------

  ## required columns present
  stopifnot(
    "dat must contain columns S, R, C, A, Y, wt" =
      all(c("S", "R", "C", "A", "Y", "wt") %in% names(dat)))

  ## required columns are binary (0/1)
  stopifnot(
    "S must be binary (0/1)" = all(dat$S %in% c(0, 1)),
    "R must be binary (0/1)" = all(dat$R %in% c(0, 1)),
    "C must be binary (0/1)" = all(dat$C %in% c(0, 1)),
    "A must be binary (0/1)" = all(dat$A %in% c(0, 1)),
    "Y must be binary (0/1)" = all(dat$Y %in% c(0, 1)))


  # fit outcome regression --------------------------------------------------

  ## logistic regression among uncensored responders
  outcome_reg <- glm(
    formula = mu_fmla,
    family = "binomial",
    data = filter(dat, S == 1, R == 1, C == 0))

  ## data sets with A set to 0, 1
  dat0 <- dat %>% mutate(A = 0, Y = 0)
  dat1 <- dat %>% mutate(A = 1, Y = 0)

  ## predict outcomes under A set to 0, 1
  dat$muhat_0 <- NA_real_
  dat$muhat_0[dat$S == 0 | dat$R == 1] <- predict(
    outcome_reg,
    newdata = dat0 %>% filter(S == 0 | R == 1),
    type = "response")

  dat$muhat_1 <- NA_real_
  dat$muhat_1[dat$S == 0 | dat$R == 1] <- predict(
    outcome_reg,
    newdata = dat1 %>% filter(S == 0 | R == 1),
    type = "response")


  # response model ----------------------------------------------------------

  ## add Q labels to data
  dat <- dat %>%
    mutate(Q = S * R * (1 - C))

  ## restrict sample to Q = 1 or S = 0
  restricted_dat <- dat %>%
    filter(Q == 1 | S == 0)

  ## fit logistic regression models for Q membership in restricted data
  Q_reg_0 <- glm(
    formula = pi_fmla,
    family = "binomial",
    data = filter(restricted_dat, A == 0 | Q == 0),
    weights = wt)

  Q_reg_1 <- glm(
    formula = pi_fmla,
    family = "binomial",
    data = filter(restricted_dat, A == 1 | Q == 0),
    weights = wt)

  ## estimated Q probabilities among uncensored responders
  Q_ind <- dat$Q == 1
  dat$Q_prob <- NA_real_
  dat$Q_prob[Q_ind & dat$A == 0] <-
    predict(
      Q_reg_0,
      newdata = dat[Q_ind & dat$A == 0,],
      type = "response")

  dat$Q_prob[Q_ind & dat$A == 1] <-
    predict(
      Q_reg_1,
      newdata = dat[Q_ind & dat$A == 1,],
      type = "response")

  ## estimated propensity scores in trial data
  dat$pihat <- NA_real_
  dat$pihat[Q_ind] <- dat$Q_prob[Q_ind] /
    (1 - dat$Q_prob[Q_ind])


  # AIPW estimator ----------------------------------------------------------

  ## Hajek estimator of trial sample sizes by A
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
      or_1 = ifelse(S == 0, muhat_1 * wt, 0)) %>%

    summarise(

      ## AIPW estimators
      etahat_0 = sum((ipw_0 / n_trial_hat[1]) +
                       (or_0 / n_aux)),
      etahat_1 = sum((ipw_1 / n_trial_hat[2]) +
                       (or_1 / n_aux)))


  # variance estimator ------------------------------------------------------

  ## design matrices for variance estimation
  X0 <- model.matrix(outcome_reg, data = dat0)
  X1 <- model.matrix(outcome_reg, data = dat1)
  X_Q <- model.matrix(Q_reg_0, data = dat)

  ## variance estimator
  est_var <- get.sand.est(
    param = c(etahat$etahat_0, etahat$etahat_1,
              coef(outcome_reg),
              coef(Q_reg_0), coef(Q_reg_1)),
    n = nrow(dat),
    get.psi = function(xx) {

      ## extract pieces of combined parameter
      eh0 <- xx[1]
      eh1 <- xx[2]
      bb_outcome <- xx[2 + 1:length(coef(outcome_reg))]
      bb_Q_0 <- xx[2 + length(coef(outcome_reg)) +
                      1:length(coef(Q_reg_0))]
      bb_Q_1 <- tail(xx, length(coef(Q_reg_1)))

      ## re-create predicted outcomes using supplied parameter
      muhat_0_ <- plogis(as.vector(X0 %*% bb_outcome))
      muhat_1_ <- plogis(as.vector(X1 %*% bb_outcome))

      ## re-create weights using supplied parameter
      Q_prob_ <- case_when(
        dat$Q == 1 & dat$A == 0 ~
          plogis(as.vector(X_Q %*% bb_Q_0)),
        dat$Q == 1 & dat$A == 1 ~
          plogis(as.vector(X_Q %*% bb_Q_1)),
        .default = 0)
      pihat_ <- Q_prob_ / (1 - Q_prob_)

      ## recreate Hajek denominator
      n_trial_hat_ <- dat %>%
        mutate(pihat = pihat_) %>%
        filter(Q == 1) %>%
        mutate(
          term_0 = (1 - A) / pihat,
          term_1 = A / pihat) %>%
        summarise(
          nhat_0 = sum(term_0),
          nhat_1 = sum(term_1)) %>%
        unlist()

      ## stacked estimating function
      cbind(

        ## outcome regression estimating function
        psi.lr(data = dat,
               beta = bb_outcome,
               formula = mu_fmla) *
          dat$S * dat$R * (1 - dat$C),

        ## propensity score regression estimating functions
        psi.lr(data = dat,
               beta = bb_Q_0,
               formula = pi_fmla) *
          as.numeric((dat$Q == 1 & dat$A == 0) | dat$S == 0),

        psi.lr(data = dat,
               beta = bb_Q_1,
               formula = pi_fmla) *
          as.numeric((dat$Q == 1 & dat$A == 1) | dat$S == 0),

        ## etahat estimating function
        case_when(
          dat$Q == 1 & dat$A == 0 ~
            (nrow(dat) / n_trial_hat_[1]) * (dat$Y - muhat_0_) / pihat_,
          dat$S == 0 ~
            (nrow(dat) / n_aux) * dat$wt * muhat_0_,
          .default = 0) - eh0,

        case_when(
          dat$Q == 1 & dat$A == 1 ~
            (nrow(dat) / n_trial_hat_[2]) * (dat$Y - muhat_1_) / pihat_,
          dat$S == 0 ~
            (nrow(dat) / n_aux) * dat$wt * muhat_1_,
          .default = 0) - eh1)

    }
  )


  # return list of results --------------------------------------------------

  res <- list(

    ## causal parameter estimate
    eta_hat = c("etahat_0" = etahat$etahat_0,
                "etahat_1" = etahat$etahat_1),

    ## estimated covariance of eta_hat
    eta_hat_cov = est_var[1:2, 1:2],

    ## outcome regression model results
    outcome_reg = outcome_reg,

    ## Q membership regression model results
    Q_reg_0 = Q_reg_0,
    Q_reg_1 = Q_reg_1,

    ## data set including estimated propensity score weights
    dat = dat)

  return(res)

}
