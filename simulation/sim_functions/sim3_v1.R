###############################################################################
###############################################################################

# Simulation 3 Development: nonparametric sample splitting

# Brian Richardson

# 2026-04-09

###############################################################################
###############################################################################


sim3_fun <- function(m, n_trial, n_aux, p_resp, p_cens,
                     K, seed, run.checks = F) {

  # (for code checking only) ------------------------------------------------

  if (FALSE) {

    ## clear workspace
    rm(list = ls())

    ## load packages
    library(dplyr)
    library(tidyr)
    library(xgboost)
    library(devtools)

    ## set work directory
    setwd("C:/Users/brich/OneDrive - University of North Carolina at Chapel Hill/Desktop/CIRL/PopART/crown")

    ## load crown package
    load_all()

    ## number of clusters
    m <- 20

    ## trial sample size
    n_trial <- 5000

    ## auxiliary sample size
    n_aux <- 5000

    ## response rate
    p_resp <- 0.5

    ## censoring probability among responders
    p_cens <- 0.3

    ## indicator for correctly specified outcome and censoring models
    mu_correct <- T
    pi_correct <- T

    ## seed for random number generation
    seed <- 1

    ## indicator for whether to run checks
    run.checks <- T

    ## probability of treatment assignment
    pA <- 0.5

    ## number of folds
    K <- 5
  }

  # simulate data -----------------------------------------------------------

  ## intercepts for R and C models (justified below)
  R_int <- 0.4449358
  C_int <- -0.8711302
  if (FALSE) {

    ## simulate covariates once
    set.seed(1)
    n_large <- 1e6
    A  <- rbinom(n_large, 1, 0.5)
    W1 <- rbinom(n_large, 1, 0.5)
    W2 <- rnorm(n_large)
    W3 <- rnorm(n_large)

    ## solve for R model intercept: 0.4546268
    R_int <- uniroot(
      f = function(ri) mean(plogis(ri + -2 * A * W1 + W2 * 0.2*W3^2)) - p_resp,
      interval = c(-10, 10))$root

    ## solve for C model intercept: -1.015664
    C_int <- uniroot(
      f = function(ci) mean(plogis(ci - 0.25*A + 0.25*W1 + 0.2*W2^2*W3)) - p_cens,
      interval = c(-10, 10))$root
  }

  ## possibly adjust trial and auxiliary sample sizes to be multiples of m
  n_trial <- m * round(n_trial / m)
  n_aux <- m * round(n_aux / m)

  ## combined sample size
  n_comb <- n_trial + n_aux

  ## set seed
  set.seed(seed)

  ## cluster-level covariate generated once from N(0, 1)
  XX1 <- c(0.92, 0.78, 0.07, -1.99, 0.62, -0.06, -0.16, -1.47, -0.48, 0.42,
           1.36, -0.10, 0.39, -0.05, -1.38, -0.41, -0.39, -0.06, 1.10, 0.76)

  ## cluster-level randomization
  AA <- sample(rep(c(0, 1), times = m / 2), replace = F)

  ## simulate trial data
  trial_dat <- data.frame(
    id = 1:n_trial) %>%
    mutate(

      ## cluster assignment
      clust = rep(1:m, length = n_trial),

      ## one cluster-level covariate
      X1 = XX1[clust],

      ## arm assigned based on cluster randomization
      A = AA[clust],

      ## one binary covariate
      W1 = rbinom(n_trial, 1, 0.5),

      ## two normal covariates
      W2 = rnorm(n_trial, 0, 1),
      W3 = rnorm(n_trial, 0, 1),

      ## response indicator
      pR = plogis(R_int - 2*A*W1 + 0.2*W2*W3^2),
      R = rbinom(n_trial, 1, pR),

      ## censoring indicator
      pC = plogis(C_int - 0.25*A + 0.25*W1 + 0.2*W2^2*W3),
      C = rbinom(n_trial, 1, pC),
      C = ifelse(R == 1, C, 0),

      ## potential outcomes
      mu0 = plogis(-1 + 2*W1),
      mu1 = plogis(0 - 1*W1 + 0.2*sin(W2) - 0.2*W3^2 + 0.25*X1),
      Y0 = rbinom(n_trial, 1, mu0),
      Y1 = rbinom(n_trial, 1, mu1),

      ## outcome
      Y = (1 - A)*Y0 + A*Y1,
      Y = ifelse(R == 1 & C == 0, Y, 0),

      ## survey weight of 1 for trial data
      wt = 1)

  ## simulate auxiliary sample
  aux_dat <- data.frame(
    id = 1:n_aux + n_trial) %>%
    mutate(

      ## cluster assignment
      clust = rep(1:m, length = n_aux),

      ## one cluster-level covariate
      X1 = XX1[clust],

      ## arm assigned based on cluster randomization
      A = AA[clust],

      ## one binary covariate: W1 = 1 over-represented
      W1 = rbinom(n_aux, 1, 0.75),

      ## two normal covariates
      W2 = rnorm(n_aux, 0, 1),
      W3 = rnorm(n_aux, 0, 1),

      ## survey weights to reflect over-representation of W_1
      wt = ifelse(W1 == 0, 0.75, 0.25),

      ## missing response, censoring indicator, and outcome
      R = 0,
      C = 0,
      Y = 0)
  aux_dat$wt <- aux_dat$wt / mean(aux_dat$wt)

  ## combined data
  dat <- bind_rows(
    trial_dat %>% select(!c(pR, pC, mu0, mu1, Y0, Y1)) %>% mutate(S = 1),
    aux_dat %>% mutate(S = 0))

  ## potential outcome means
  eta0 <- mean(trial_dat$mu0)
  eta1 <- mean(trial_dat$mu1)

  # check simulated data ----------------------------------------------------

  if (run.checks) {

    ## check response rate
    print(mean(trial_dat$R)); print(p_resp)
    trial_dat %>% group_by(A, W1) %>% summarise(p_resp = mean(R)) %>% print()
    plot(trial_dat$W2, trial_dat$pR)
    plot(trial_dat$W3, trial_dat$pR)

    ## check censoring rate
    trial_dat %>% filter(R == 1) %>% group_by(A, W1) %>%
      summarise(p_cens = mean(C)) %>% print()
    trial_dat %>% filter(R == 1) %>%
      summarise(p_cens = mean(C)) %>% print()
    print(p_cens)
    plot(trial_dat$W2, trial_dat$pC)
    plot(trial_dat$W3, trial_dat$pC)

    ## check mean potential outcomes in population
    print(eta0); print(eta1); print(eta1 - eta0)
    trial_dat %>% group_by(W1) %>% summarise(
      eta0 = mean(Y0), eta1 = mean(Y1)) %>%
      mutate(rd = eta1 - eta0) %>%
      print()
    plot(trial_dat$W2, trial_dat$mu0)
    plot(trial_dat$W2, trial_dat$mu1)
    plot(trial_dat$W3, trial_dat$mu0)
    plot(trial_dat$W3, trial_dat$mu1)

    ## check mean potential outcomes among responders
    print(c(eta0, eta1))
    trial_dat %>% filter(R == 1) %>% summarise(m0 = mean(Y0), m1 = mean(Y1))
  }

  # analyze data ------------------------------------------------------------

  ## proposed AIPW with (incorrect) parametric models
  mu_fmla <- Y ~ A * (X1 + W1 + W2 + W3)
  pi_fmla <- Q ~ X1 + W1 + W2 + W3
  aipw_par <- aipw_fit(
    dat = dat,
    mu_fmla = mu_fmla,
    pi_fmla = pi_fmla)

  ## proposed AIPW with nonparametric sample splitting
  aipw_nss <- aipw_fit_nss(
    dat = dat,
    mu_covariates = c("X1", "W1", "W2", "W3"),
    pi_covariates = c("X1", "W1", "W2", "W3"),
    K = K,
    num.trees = num.trees)

  # combine results ---------------------------------------------------------

  ## make data frame with results
  res <- bind_rows(

    mutate(aipw_par$eta_results, name = "aipw_par"),
    mutate(aipw_nss$eta_results, name = "aipw_nss")) %>%

    separate(
      name,
      into = c("est", "version"),
      sep = "_") %>%

    mutate(

      rdhat = etahat_1 - etahat_0,
      rrhat = etahat_1 / etahat_0,
      var_rd = cov_11 + cov_00 - 2*cov_01,
      var_rr = (cov_11 / (etahat_0^2)) +
        (cov_00 * (etahat_1^2) / (etahat_0^4)) -
        cov_01 * 2 * etahat_1 / (etahat_0^3),

      Estimator = factor(
        est,
        levels = c("aipw"),
        labels = c("AIPW")),
      Version = factor(
        version,
        levels = c("par", "nss"),
        labels = c("Parametric", "Nonparametric SS"))) %>%

    select(!c(est, version)) %>%

    mutate(

      # true etas
      eta_0 = eta0,
      eta_1 = eta1,

      # true risk differences and ratios
      rd = eta_1 - eta_0,
      rr = eta_1 / eta_0,

      # simulation settings
      seed = seed,
      m = m,
      n_trial = n_trial,
      n_aux = n_aux,
      p_resp = p_resp,
      p_cens = p_cens,
      K = K)

  # Plot Results ------------------------------------------------------------

  if (run.checks) {

    ## plot risk differences
    library(ggplot2)
    res %>%
      mutate(rd_lower = rdhat - qnorm(0.975) * sqrt(var_rd),
             rd_upper = rdhat + qnorm(0.975) * sqrt(var_rd)) %>%
    ggplot(aes(x = Estimator,
               y = rdhat,
               ymin = rd_lower,
               ymax = rd_upper)) +

      geom_point() +
      geom_linerange() +
      geom_hline(aes(yintercept = rd),
                 linetype = "dashed",
                 color = "blue") +

      facet_wrap(~ Version) +
      labs(y = expression(hat(RD)),
           x = NULL) +
      ggtitle("") +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"),
            legend.position = "none")
  }

  # return results ----------------------------------------------------------

  return(res)
}


