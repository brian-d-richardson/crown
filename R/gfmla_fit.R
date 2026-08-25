#' Naive g-formula estimator
#'
#' estimate eta(0) and eta(1) using a naive g-formula estimator, i.e., ignoring
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
#' @param mu_fmla a formula for the outcome regression model using variables in
#' `dat`
#'
#' @return a list containing the following:
#' \itemize{
#' \item `eta_hat`: a numeric vector estimated mean potential outcomes eta(0)
#' and eta(1)
#' \item `eta_hat_cov`: a numeric matrix, estimated covariance of `eta_hat`
#' \item `outcome_reg`: a list, results of outcome regression model
#' }
#'
#' @export
gfmla_fit_naive <- function(dat, mu_fmla) {


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

  ## logistic regression among uncensored
  outcome_reg <- glm(
    formula = mu_fmla,
    family = "binomial",
    data = filter(dat, C == 0),
    weights = wt)


  # g-formula estimator -----------------------------------------------------

  ## auxiliary data sets with A set to 0, 1
  dat0 <- dat %>% mutate(A = 0, Y = 0)
  dat1 <- dat %>% mutate(A = 1, Y = 0)

  ## predict outcomes under A set to 0, 1
  muhat_0 <- predict(
    outcome_reg,
    newdata = dat0,
    type = "response")

  muhat_1 <- predict(
    outcome_reg,
    newdata = dat1,
    type = "response")

  ## g-formula estimator
  etahat_0 = weighted.mean(x = muhat_0, w = dat$wt)
  etahat_1 = weighted.mean(x = muhat_1, w = dat$wt)


  # variance estimator ------------------------------------------------------

  ## design matrices for variance estimation
  X0 <- model.matrix(outcome_reg, data = dat0)
  X1 <- model.matrix(outcome_reg, data = dat1)

  ## variance estimator
  est_var <- get.sand.est(
    param = c(etahat_0, etahat_1, coef(outcome_reg)),
    n = nrow(dat),
    get.psi = function(xx) {

      ## extract pieces of combined parameter
      eh0 <- xx[1]
      eh1 <- xx[2]
      bb <- tail(xx, -2)

      ## stacked estimating function
      cbind(

        ## outcome regression estimating function
        psi.lr(data = dat,
               beta = bb,
               formula = mu_fmla) *
          (1 - dat$C),

        ## etahat estimating function
        dat$wt * (plogis(as.vector(X0 %*% bb)) - eh0),
        dat$wt * (plogis(as.vector(X1 %*% bb)) - eh1)
      )
    }
  )


  # return list of results --------------------------------------------------

  res <- list(

    ## causal parameter estimate
    eta_hat = c("etahat_0" = etahat_0,
                "etahat_1" = etahat_1),

    ## estimated covariance of eta_hat
    eta_hat_cov = est_var[1:2, 1:2],

    ## outcome regression model results
    outcome_reg = outcome_reg)

  return(res)
}


#' Crown g-formula estimator
#'
#' estimate eta(0) and eta(1) using a crown g-formula estimator, i.e.,
#' accounting for nonresponse
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
#' @return a list containing the following:
#' \itemize{
#' \item `eta_hat`: a numeric vector estimated mean potential outcomes eta(0)
#' and eta(1)
#' \item `eta_hat_cov`: a numeric matrix, estimated covariance of `eta_hat`
#' \item `outcome_reg`: a list, results of outcome regression model
#' }
#'
#' @export
gfmla_fit <- function(dat, mu_fmla) {


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
    data = filter(dat, S == 1, R == 1, C == 0),
    weights = wt)


  # g-formula estimator -----------------------------------------------------

  ## auxiliary data sets with A set to 0, 1
  dat0 <- dat %>% mutate(A = 0, Y = 0)
  dat1 <- dat %>% mutate(A = 1, Y = 0)

  ## predict outcomes under A set to 0, 1 in auxiliary data
  muhat_0 <- predict(
    outcome_reg,
    newdata = dat0 %>% filter(S == 0),
    type = "response")

  muhat_1 <- predict(
    outcome_reg,
    newdata = dat1 %>% filter(S == 0),
    type = "response")

  ## g-formula estimator
  etahat_0 = weighted.mean(x = muhat_0, w = dat$wt[dat$S == 0])
  etahat_1 = weighted.mean(x = muhat_1, w = dat$wt[dat$S == 0])


  # variance estimator ------------------------------------------------------

  ## design matrices for variance estimation
  X0 <- model.matrix(outcome_reg, data = dat0)
  X1 <- model.matrix(outcome_reg, data = dat1)

  ## variance estimator
  est_var <- get.sand.est(
    param = c(etahat_0, etahat_1, coef(outcome_reg)),
    n = nrow(dat),
    get.psi = function(xx) {

      ## extract pieces of combined parameter
      eh0 <- xx[1]
      eh1 <- xx[2]
      bb <- tail(xx, -2)

      ## stacked estimating function
      cbind(

        ## outcome regression estimating function
        psi.lr(data = dat,
               beta = bb,
               formula = mu_fmla) *
          dat$S * dat$R * (1 - dat$C),

        ## etahat estimating function
        (1 - dat$S) * dat$wt * (plogis(as.vector(X0 %*% bb)) - eh0),
        (1 - dat$S) * dat$wt * (plogis(as.vector(X1 %*% bb)) - eh1)
      )
    }
  )


  # return list of results --------------------------------------------------

  res <- list(

    ## causal parameter estimate
    eta_hat = c("etahat_0" = etahat_0,
                "etahat_1" = etahat_1),

    ## estimated covariance of eta_hat
    eta_hat_cov = est_var[1:2, 1:2],

    ## outcome regression model results
    outcome_reg = outcome_reg)

  return(res)
}
