## naive g-formula estimator
gfmla_fit_naive <- function(dat, mu_fmla) {


  # fit outcome regression --------------------------------------------------

  ## logistic regression among uncensored
  outcome_reg <- glm(
    formula = mu_fmla,
    family = "binomial",
    data = filter(dat, C == 0))


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

  # return list of results
  res <- list(

    # causal parameter estimates and covariance
    eta_results =
      data.frame(
        etahat_0 = etahat_0,
        etahat_1 = etahat_1,
        cov_00 = est_var[1, 1],
        cov_01 = est_var[1, 2],
        cov_11 = est_var[2, 2]),

    # outcome regression model results
    outcome_reg_results = outcome_reg)

  return(res)
}



## proposed g-formula estimator
gfmla_fit <- function(dat, mu_fmla) {


  # fit outcome regression --------------------------------------------------

  ## logistic regression among uncensored responders
  outcome_reg <- glm(
    formula = mu_fmla,
    family = "binomial",
    data = filter(dat, S == 1, R == 1, C == 0))


  # g-formula estimator -----------------------------------------------------

  ## auxiliary data sets with A set to 0, 1
  dat0 <- dat %>% mutate(A = 0, Y = 0)
  dat1 <- dat %>% mutate(A = 1, Y = 0)

  ## predict outcomes under A set to 0, 1
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

  # return list of results
  res <- list(

    # causal parameter estimates and covariance
    eta_results =
      data.frame(
        etahat_0 = etahat_0,
        etahat_1 = etahat_1,
        cov_00 = est_var[1, 1],
        cov_01 = est_var[1, 2],
        cov_11 = est_var[2, 2]),

    # outcome regression model results
    outcome_reg_results = outcome_reg)

  return(res)
}
