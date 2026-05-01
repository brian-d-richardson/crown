## nonparametric nuisance estimation wrapper funtion
nonpar_est <- function(x, y, method, arguments = NULL, wts = NULL) {
  
  mod <- NULL
  
  if (method == "SuperLearner") {
    
    mod <- SuperLearner::SuperLearner(
      Y = y,
      X = x,
      family = binomial(),
      obsWeights = wts,
      SL.library = arguments$SL.library,
      cvControl = arguments$cvControl)
    
  } else if (method == "xgboost") {
    
    mod <- xgboost::xgboost(
      x = x,
      y = y,
      objective = "binary:logistic",
      weights = wts)
    
  } else {
    print("unrecognized method")
  }
  return(mod)
}


## nonparametric nuisance prediction wrapper function
nonpar_pred <- function(mod, newdata, method) {
  
  pred <- NULL
  
  if (method == "SuperLearner") {
    
    pred <- predict(
      mod,
      newdata = newdata,
      onlySL = T)$pred
    
  } else if (method == "xgboost") {
    
    pred <- predict(
      mod,
      newdata = newdata)
    
  } else {
    print("unrecognized method")
  }
  
  return(pred)
  
}