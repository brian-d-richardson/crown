## estimating function for logistic regression
psi.lr <- function(data, beta, formula) {
  
  # model frame (ensures same handling of factors, etc.)
  mf <- model.frame(formula, data)
  
  # response
  y <- model.response(mf)
  
  # design matrix
  X <- model.matrix(formula, data)
  
  # linear predictor
  lin.pred <- as.vector(X %*% beta)
  
  # inverse logit of linear predictor
  mu <- plogis(lin.pred)
  
  # score (estimating function)
  psi <- X * as.vector(y - mu) * data$wt
  
  return(psi)
  
}
