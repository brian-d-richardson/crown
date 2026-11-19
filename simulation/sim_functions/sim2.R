###############################################################################
###############################################################################

# Simulation 2 Function: complex DGP

# Brian Richardson

# 2026-11-19

###############################################################################
###############################################################################


sim2_fun <- function(m, n_trial, n_aux, p_resp, p_obsout,
                     K = 5, seed,
                     run.checks = F) {
  
  # (for code checking only) ------------------------------------------------
  
  if (FALSE) {
    
    ## clear workspace
    rm(list = ls())
    
    ## load packages
    library(dplyr)
    library(tidyr)
    library(devtools)
    library(here)
    library(xgboost)
    
    ## set work directory
    #setwd(here::here(()); setwd(dirname(dirname(getwd())))
    setwd("C:/Users/richabr/Desktop/crown")
    
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
    
    ## observed outcome probability among responders
    p_obsout <- 0.7
    
    ## number of folds
    K <- 5
    
    ## seed for random number generation
    seed <- 1
    
    ## probability of treatment assignment
    pA <- 0.5
    
    ## covariates
    mu_covariates <- c("X1", "W1", "W2", "W3")
    pi_covariates <- c("X1", "W1", "W2", "W3")
    
    ## indicator for whether to run checks
    run.checks <- T
  }
  
  # simulate data -----------------------------------------------------------
  
  ## intercepts for R and C models (justified below)
  RW_int <- -1.417151
  RY_int <- 0.8532847
  if (FALSE) {
    
    ## simulate covariates once
    set.seed(1)
    n_large <- 1e6
    A  <- rbinom(n_large, 1, 0.5)
    W1 <- rbinom(n_large, 1, 0.5)
    W2 <- rnorm(n_large)
    W3 <- rnorm(n_large)
    
    ## solve for RW model intercept: -1.417151
    RW_int <- uniroot(
      f = function(ri) mean(plogis(
        ri +  + 2*as.numeric(W2^2 < 1))) - p_resp,
      interval = c(-10, 10))$root
    
    ## solve for C model intercept: 0.8532847
    RY_int <- uniroot(
      f = function(ryi) mean(plogis(
        ryi + 0.25*A - 0.25*W1)) - p_obsout,
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
      pRW = plogis(RW_int + 2*as.numeric(W2^2 < 1)),
      RW = rbinom(n_trial, 1, pRW),
      
      ## observed outcome indicator
      pRY = plogis(RY_int + 0.25*A - 0.25*W1),
      RY = rbinom(n_trial, 1, pRY),
      RY = ifelse(RW == 1, RY, 0),
      
      ## potential outcomes
      mu0 = ifelse(W2^2 < 1, 0.9, plogis(0.25*sin(pi*W3/4))),
      mu1 = ifelse(W2^2 < 1, 0.1, plogis(0.25*sin(pi*W3/4))),
      Y0 = rbinom(n_trial, 1, mu0),
      Y1 = rbinom(n_trial, 1, mu1),
      
      ## outcome
      Y = (1 - A)*Y0 + A*Y1,
      Y = ifelse(RW == 1 & RY == 1, Y, 0),
      
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
      
      ## missing response, observed outcome indicator, and outcome
      RW = 0,
      RY = 0,
      Y = 0)
  
  ## combined data
  dat <- bind_rows(
    trial_dat %>% select(!c(Y0, Y1)) %>% mutate(S = 1),
    aux_dat %>% mutate(S = 0))
  
  ## potential outcome means
  eta0 <- mean(trial_dat$Y0)
  eta1 <- mean(trial_dat$Y1)
  
  # check simulated data ----------------------------------------------------
  
  if (run.checks) {
    
    ## check response rate
    print(mean(trial_dat$RW)); print(p_resp)
    trial_dat %>% group_by(W2^2 < 1) %>% summarise(p_resp = mean(RW)) %>% print()
    plot(trial_dat$W2, trial_dat$pR)
    plot(trial_dat$W3, trial_dat$pR)
    
    ## check observed outcome rate
    trial_dat %>% filter(RW == 1) %>% group_by(A, W1) %>%
      summarise(p_obsout = mean(RY)) %>% print()
    trial_dat %>% filter(RW == 1) %>%
      summarise(p_obsout = mean(RY)) %>% print()
    print(p_obsout)
    plot(trial_dat$W2, trial_dat$pRY)
    plot(trial_dat$W3, trial_dat$pRY)
    
    ## check mean potential outcomes in population
    print(eta0); print(eta1); print(eta1 - eta0)
    trial_dat %>% group_by(W2^2 < 1) %>% summarise(
      eta0 = mean(Y0), eta1 = mean(Y1)) %>%
      mutate(rd = eta1 - eta0) %>%
      print()
    plot(trial_dat$W2, trial_dat$mu0)
    plot(trial_dat$W2, trial_dat$mu1)
    plot(trial_dat$W3, trial_dat$mu0)
    plot(trial_dat$W3, trial_dat$mu1)
    
    ## check mean potential outcomes among responders
    print(c(eta0, eta1))
    trial_dat %>% filter(RW == 1) %>% summarise(m0 = mean(Y0), m1 = mean(Y1))
  }
  
  # define parameters for analysis ------------------------------------------
  
  ## outcome regression and propensity score formulas
  mu_fmla <- Y ~ A * (X1 + W1 + W2 + W3)
  pi_fmla <- Q ~ X1 + W1 + W2 + W3
  
  
  # analyze data ------------------------------------------------------------
  
  ## proposed g-formula
  gfmla_crown <- fit_gformula(
    data = dat,
    outcome_formula = mu_fmla,
    version = "proposed")
  
  ## proposed IPW
  ipw_crown <- fit_ipw(
    data = dat,
    weight_formula = pi_fmla,
    version = "proposed")
  
  ## proposed AIPW
  aipw_crown <- fit_aipw(
    data = dat,
    outcome_formula = mu_fmla,
    weight_formula = pi_fmla,
    version = "proposed")
  
  ## proposed DML
  dml_crown_fit <- fit_dml(
    data = dat,
    mu_covariates = c("X1", "W1", "W2", "W3"),
    pi_covariates = c("X1", "W1", "W2", "W3"))
  dml_crown <- data.frame(
    etahat_0 = dml_crown_fit$eta_hat[1],
    etahat_1 = dml_crown_fit$eta_hat[2],
    cov_00 = dml_crown_fit$eta_hat_cov[1, 1],
    cov_01 = dml_crown_fit$eta_hat_cov[1, 2],
    cov_11 = dml_crown_fit$eta_hat_cov[2, 2])

  # combine results ---------------------------------------------------------
  
  ## helper function to format results
  extract_eta_info <- function(name) {
    
    obj <- get(name)
    
    data.frame(
      etahat_0 = obj$etahat_0,
      etahat_1 = obj$etahat_1,
      cov_00 = obj$cov_00,
      cov_01 = obj$cov_01,
      cov_11 = obj$cov_11,
      name = name,
      stringsAsFactors = FALSE)
  }
  
  ## make data frame with results
  res <- bind_rows(
    
    extract_eta_info("gfmla_crown"),
    extract_eta_info("ipw_crown"),
    extract_eta_info("aipw_crown"),
    extract_eta_info("dml_crown")) %>%
    
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
        levels = c("ipw", "gfmla", "aipw", "dml"),
        labels = c("IPW", "G-Formula", "AIPW", "DML"))) %>%
    
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
      p_obsout = p_obsout)
  
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
                 ymax = rd_upper,
                 fill = Estimator,
                 color = Estimator)) +
      
      geom_point() +
      geom_linerange() +
      geom_hline(aes(yintercept = rd),
                 linetype = "dashed",
                 color = "blue") +
      labs(y = expression(hat(RD)),
           x = NULL) +
      ggtitle("") +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"),
            legend.position = "none")
  }
  
  # Return Results ----------------------------------------------------------
  
  return(res)
}

#sim2_fun(m = 20, n_trial = 4E3, n_aux = 1E3, p_resp = 0.5, p_obsout = 0.7, seed = 1)

