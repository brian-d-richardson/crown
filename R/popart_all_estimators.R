###############################################################################
###############################################################################

# PopART Analysis Function

# Brian Richardson

# 2026-02-10

###############################################################################
###############################################################################

popart_all_estimators <- function(
    dat, outcome_reg_fmla, censor_reg_fmla, response_reg_fmla,
    pA = 0.5, run.checks = F) {
  
  
  # STEP 1: G-formula -------------------------------------------------------
  
  ## naive g-formula
  gfmla_naive <- gfmla_fit_naive(
    dat = filter(dat, S == 1, R ==1),
    outcome_reg_fmla = outcome_reg_fmla)
  
  ## proposed g-formula
  gfmla_prop <- gfmla_fit(
    dat = dat,
    outcome_reg_fmla = outcome_reg_fmla)
  

  # STEP 2: IPW -------------------------------------------------------------

  ## naive IPW
  ipw_naive <- ipw_fit_naive(
    dat = filter(dat, S == 1, R == 1),
    censor_reg_fmla = censor_reg_fmla)
  
  ## IPW (assuming no A --> R)
  ipw_prop1 <- ipw_fit_op1(
    dat = dat,
    censor_reg_fmla = censor_reg_fmla,
    response_reg_fmla = response_reg_fmla)

  ## IPW (assuming A --> R)
  ipw_prop2 <- ipw_fit_op2(
    dat = dat,
    censor_reg_fmla = censor_reg_fmla,
    response_reg_fmla = response_reg_fmla)
  

  # STEP 4: AIPW ------------------------------------------------------------

  ## naive AIPW
  aipw_naive <- aipw_fit_naive(
    dat = filter(dat, S == 1, R == 1),
    outcome_reg_fmla = outcome_reg_fmla,
    censor_reg_fmla = censor_reg_fmla)
  
  ## AIPW (assuming no A --> R)
  aipw_prop1 <- aipw_fit_op1(
    dat = dat,
    outcome_reg_fmla = outcome_reg_fmla,
    censor_reg_fmla = censor_reg_fmla,
    response_reg_fmla = response_reg_fmla)
  
  ## AIPW (assuming A --> R)
  aipw_prop2 <- aipw_fit_op2(
    dat = dat,
    outcome_reg_fmla = outcome_reg_fmla,
    censor_reg_fmla = censor_reg_fmla,
    response_reg_fmla = response_reg_fmla)
  

  # STEP 5: format results --------------------------------------------------
  
  ## make data frame with results
  res <- bind_rows(
    
    mutate(gfmla_naive$eta_results, name = "gfmla_naive"),
    mutate(gfmla_prop$eta_results, name = "gfmla_prop1"),
    mutate(gfmla_prop$eta_results, name = "gfmla_prop2"),
    
    mutate(ipw_naive$eta_results, name = "ipw_naive"),
    mutate(ipw_prop1$eta_results, name = "ipw_prop1"),
    mutate(ipw_prop2$eta_results, name = "ipw_prop2"),
    
    mutate(aipw_naive$eta_results, name = "aipw_naive"),
    mutate(aipw_prop1$eta_results, name = "aipw_prop1"),
    mutate(aipw_prop2$eta_results, name = "aipw_prop2")) %>% 
    
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
        levels = c("ipw", "gfmla", "aipw"),
        labels = c("IPW", "G-Formula", "AIPW")),
      Version = factor(
        version,
        levels = c("naive", "prop1", "prop2"),
        labels = c("Naive", "Proposed (Option 1)", "Proposed (Option 2)"))) %>% 
    
    select(!c(est, version))
  
  ## make data frame with weights
  dat <- dat %>% 
    mutate(
      piC = ipw_prop1$dat$piC,
      piAR_op1 = ipw_prop1$dat$piAR,
      piAR_op2 = ipw_prop2$dat$piAR,
      pihat_naive = ipw_prop2$dat$piC *
        ifelse(ipw_prop2$dat$A == 1, pA, 1 - pA),
      pihat_op1 = ipw_prop1$dat$pihat,
      pihat_op2 = ipw_prop2$dat$pihat)
  
  return(list(
    results = res,
    data = dat))
}
