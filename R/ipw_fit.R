#' Naive IPW estimator
#'
#' estimate eta(0) and eta(1) using a naive IPW estimator, i.e., ignoring
#' nonresponse
#'
#' @param dat data frame containing the following columns:
#' \itemize{
#' \item `RY`: a binary indicator for whether the outcome is observed (`RY`=1) or
#' not (`RY`=0)
#' \item `A`: a binary exposure
#' \item `Y`: a binary outcome
#' \item other covariates specified in `mu_fmla`
#' \item `wt`: survey weights
#' }
#'
#' @param RY_fmla a formula for the outcome observation mechanism regression
#' model using variables in `dat`
#'
#' @param pA an optional number in (0, 1), the marginal probability of
#' treatment (A = 1), default is 0.5.
#'
#' @return a list containing the following:
#' \itemize{
#' \item `eta_hat`: a numeric vector estimated mean potential outcomes eta(0)
#' and eta(1)
#' \item `eta_hat_cov`: a numeric matrix, estimated covariance of `eta_hat`
#' \item `RY_reg`: a list, results of outcome observation regression model
#' \item `dat`: a data frame including a column for propensity score weights
#' }
#'
#' @export
ipw_fit_naive <- function(dat, RY_fmla, pA = 0.5) {


  # check input -------------------------------------------------------------

  ## required columns present
  stopifnot(
    "dat must contain columns RY, A, Y, wt" =
      all(c("RY", "A", "Y") %in% names(dat)))

  ## required columns are binary (0/1)
  stopifnot(
    "RY must be binary (0/1)" = all(dat$RY %in% c(0, 1)),
    "A must be binary (0/1)" = all(dat$A %in% c(0, 1)),
    "Y must be binary (0/1)" = all(dat$Y %in% c(0, 1)))


  # RY model ---------------------------------------------------------

  ## logistic regression model for outcome observation mechanism
  RY_reg <- glm(
    formula = RY_fmla,
    family = "binomial",
    data = dat,
    weights = wt)

  ## predicted values from RY model
  dat$piRY <- predict(
    RY_reg,
    newdata = dat,
    type = "response")

  ## treatment assignment probabilities
  dat$piA <- ifelse(dat$A == 1, pA, 1 - pA)


  # IPW estimator -----------------------------------------------------------

  ## joint probabilities of having treatment A and observed outcome
  dat$pihat <- dat$piRY * dat$piA

  ## Hajek estimator of trial sample sizes by A
  n_trial_hat <- dat %>%
    filter(RY == 1) %>%
    mutate(
      term_0 = (1 - A) / pihat,
      term_1 = A / pihat) %>%
    summarise(
      nhat_0 = sum(term_0),
      nhat_1 = sum(term_1)) %>%
    unlist()

  ## IPW estimator
  etahat <- dat %>%
    mutate(
      term_0 = ifelse(RY == 1,
                      (1 - A) * Y / pihat,
                      0),
      term_1 = ifelse(RY == 1,
                      A * Y / pihat,
                      0)) %>%
    summarise(
      etahat_0 = sum(term_0) / n_trial_hat[1],
      etahat_1 = sum(term_1) / n_trial_hat[2])


  # variance estimator ------------------------------------------------------

  ## design matrix for variance estimation
  X_RY <- model.matrix(RY_reg, data = dat)

  ## variance estimator
  est_var <- get.sand.est(
    param = c(etahat$etahat_0, etahat$etahat_1,
              n_trial_hat,
              coef(RY_reg)),
    n = nrow(dat),
    get.psi = function(xx) {

      ## extract pieces of combined parameter
      eh0 <- xx[1]
      eh1 <- xx[2]
      nh0 <- xx[3]
      nh1 <- xx[4]
      bb_RY <- xx[4 + 1:length(coef(RY_reg))]

      ## re-create weights using supplied parameter
      piRY_ <- plogis(as.vector(X_RY %*% bb_RY))
      pihat_ <- piRY_ * ifelse(dat$A == 1, pA, 1 - pA)

      ## stacked estimating function
      cbind(

        ## RY regression estimating function
        psi.lr(data = dat,
               beta = bb_RY,
               formula = RY_fmla),

        ## etahat estimating function
        ifelse(
          dat$RY == 1 & dat$A == 0,
          1 / pihat_,
          0) - nh0,

        ifelse(
          dat$RY == 1 & dat$A == 1,
          1 / pihat_,
          0) - nh1,

        ## etahat estimating function
        ifelse(
          dat$RY == 1 & dat$A == 0,
          (nrow(dat) / nh0) * dat$Y / pihat_,
          0) - eh0,

        ifelse(
          dat$RY == 1 & dat$A == 1,
          (nrow(dat) / nh1) * dat$Y / pihat_,
          0) - eh1)
    }
  )


  # return list of results --------------------------------------------------

  res <- list(

    ## causal parameter estimate
    eta_hat = c("etahat_0" = etahat$etahat_0,
                "etahat_1" = etahat$etahat_1),

    ## estimated covariance of eta_hat
    eta_hat_cov = est_var[1:2, 1:2],

    ## RY regression model results
    RY_reg = RY_reg,

    ## data set including estimated propensity score weights
    dat = dat)


  return(res)

}



#' Crown IPW estimator
#'
#' estimate eta(0) and eta(1) using a crown IPW estimator, i.e., accounting for
#' nonresponse
#'
#' @param dat data frame containing the following columns:
#' \itemize{
#'\item `S`: a binary indicator for whether the observation belongs to the trial
#'data (`S`=1) the auxiliary data (`S`=0)
#' \item `RW`: a binary indicator for whether the observation is a responder
#' (`RW`=1) or not (`RW`=0)
#' \item `RY`: a binary indicator for whether the outcome is observed (`RY`=1) or
#' not (`RY`=0)
#' \item `A`: a binary exposure
#' \item `Y`: a binary outcome
#' \item other covariates specified in `mu_fmla`
#' \item `wt`: survey weights
#'
#' }
#' @param pi_fmla a formula for the propensity score regression model using
#' variables in `dat`
#'
#' @return a list containing the following:
#' \itemize{
#' \item `eta_hat`: a numeric vector estimated mean potential outcomes eta(0)
#' and eta(1)
#' \item `eta_hat_cov`: a numeric matrix, estimated covariance of `eta_hat`
#' \item `Q_reg_0`: a list, results of Q0 membership regression model
#' \item `Q_reg_1`: a list, results of Q1 membership regression model
#' \item `dat`: a data frame including a column for propensity score weights
#' }
#'
#' @export
ipw_fit <- function(dat, pi_fmla) {


  # check input -------------------------------------------------------------

  ## required columns present
  stopifnot(
    "dat must contain columns S, RW, RY, A, Y, wt" =
      all(c("S", "RW", "RY", "A", "Y", "wt") %in% names(dat)))

  ## required columns are binary (0/1)
  stopifnot(
    "S must be binary (0/1)" = all(dat$S %in% c(0, 1)),
    "RW must be binary (0/1)" = all(dat$RW %in% c(0, 1)),
    "RY must be binary (0/1)" = all(dat$RY %in% c(0, 1)),
    "A must be binary (0/1)" = all(dat$A %in% c(0, 1)),
    "Y must be binary (0/1)" = all(dat$Y %in% c(0, 1)))


  # response model ----------------------------------------------------------

  ## add Q labels to data
  dat <- dat %>%
    mutate(Q = S * RW * RY)

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

  ## estimated Q probabilities among responders with observed outcome
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


  # IPW estimator -----------------------------------------------------------

  ## Hajek estimator of trial sample size by A
  n_trial_hat <- dat %>%
    filter(Q == 1) %>%
    mutate(
      term_0 = (1 - A) / pihat,
      term_1 = A / pihat) %>%
    summarise(
      nhat_0 = sum(term_0),
      nhat_1 = sum(term_1)) %>%
    unlist()

  ## IPW estimator
  etahat <- dat %>%
    mutate(
      term_0 = ifelse(Q == 1,
                      (1 - A) * Y / pihat,
                      0),
      term_1 = ifelse(Q == 1,
                      A * Y / pihat,
                      0)) %>%
    summarise(
      etahat_0 = sum(term_0) / n_trial_hat[1],
      etahat_1 = sum(term_1) / n_trial_hat[2])


  # variance estimator ------------------------------------------------------

  ## design matrix for variance estimation
  X_Q <- model.matrix(Q_reg_0, data = dat)

  ## variance estimator
  est_var <- get.sand.est(
    param = c(etahat$etahat_0, etahat$etahat_1,
              n_trial_hat,
              coef(Q_reg_0), coef(Q_reg_1)),
    n = nrow(dat),
    get.psi = function(xx) {

      ## extract pieces of combined parameter
      eh0 <- xx[1]
      eh1 <- xx[2]
      nh0 <- xx[3]
      nh1 <- xx[4]
      bb_Q_0 <- xx[4 + 1:length(coef(Q_reg_0))]
      bb_Q_1 <- tail(xx, length(coef(Q_reg_1)))

      ## re-create weights using supplied parameter
      Q_prob_ <- case_when(
        dat$Q == 1 & dat$A == 0 ~
         plogis(as.vector(X_Q %*% bb_Q_0)),
        dat$Q == 1 & dat$A == 1 ~
          plogis(as.vector(X_Q %*% bb_Q_1)),
        .default = 0)
      pihat_ <- Q_prob_ / (1 - Q_prob_)

      ## stacked estimating function
      cbind(

        ## propensity score regression estimating functions
        psi.lr(data = dat,
               beta = bb_Q_0,
               formula = pi_fmla) *
          as.numeric((dat$Q == 1 & dat$A == 0) | dat$S == 0),

        psi.lr(data = dat,
               beta = bb_Q_1,
               formula = pi_fmla) *
          as.numeric((dat$Q == 1 & dat$A == 1) | dat$S == 0),

        ## n hat estimating function
        ifelse(
          dat$Q == 1 & dat$A == 0,
          1 / pihat_,
          0) - nh0,

        ifelse(
          dat$Q == 1 & dat$A == 1,
          1 / pihat_,
          0) - nh1,

        ## eta hat estimating function
        ifelse(
          dat$Q == 1 & dat$A == 0,
          (nrow(dat) / nh0) * dat$Y / pihat_,
          0) - eh0,

        ifelse(
          dat$Q == 1 & dat$A == 1,
          (nrow(dat) / nh1) * dat$Y / pihat_,
          0) - eh1)
    }
  )


  # return list of results --------------------------------------------------

  res <- list(

    ## causal parameter estimate
    eta_hat = c("etahat_0" = etahat$etahat_0,
                "etahat_1" = etahat$etahat_1),

    ## estimated covariance of eta_hat
    eta_hat_cov = est_var[1:2, 1:2],

    ## Q membership regression model results
    Q_reg_0 = Q_reg_0,
    Q_reg_1 = Q_reg_1,

    ## data set including estimated propensity score weights
    dat = dat)

  return(res)

}

