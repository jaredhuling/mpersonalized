unique_rule_lasso_method = function(modelYlist, modelXlist, Ybar, Xbar, Xsd, lambda){

  x = do.call(rbind, modelXlist)
  y = unlist(modelYlist)
  total_n = length(y)

  x = sqrt(total_n) * x; y = sqrt(total_n) * y

  lasso_model = glmnet(x = x, y = y, family = "gaussian", standardize = FALSE,
                       intercept  = FALSE, lambda = lambda)

  if (is.null(lambda))
    lambda = lasso_model$lambda

  decreasing_order = order(lambda, decreasing = TRUE)

  nlambda = length(lambda)
  betalist = vector("list", nlambda)
  interceptlist = vector("list", nlambda)

  for (ind in 1:nlambda){

    beta = lasso_model$beta[,ind]
    intercept = Ybar - sum(Xbar * beta)
    beta = beta / Xsd

    betalist[[decreasing_order[ind]]] = beta
    interceptlist[[decreasing_order[ind]]] = intercept
  }

  return(list(interceptlist = interceptlist, betalist = betalist, lambda = lambda))

}
