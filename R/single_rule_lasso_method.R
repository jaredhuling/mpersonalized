single_rule_lasso_method = function(modelYlist, modelXlist, Wlist, Ybar, Xbar, Xsd, lambda, surrogate = c("squared_error", "logistic"))
{

  surrogate <- match.arg(surrogate)
  x = do.call(rbind, modelXlist)
  y = unlist(modelYlist)
  total_n = length(y)

  x = sqrt(total_n) * x; y = sqrt(total_n) * y

  if (surrogate == "squared_error")
  {
    family <- "gaussian"

    lasso_model = glmnet(x = x, y = y, family = family,
                         standardize = FALSE,
                         intercept  = FALSE, lambda = lambda)


    decreasing_order = order(lambda, decreasing = TRUE)

    nlambda = length(lambda)
    betalist = vector("list", nlambda)
    interceptlist = vector("list", nlambda)

    for (ind in 1:nlambda)
    {

      beta = lasso_model$beta[,ind]
      intercept = Ybar - sum(Xbar * beta)
      beta = beta / Xsd

      betalist[[decreasing_order[ind]]] = beta
      interceptlist[[decreasing_order[ind]]] = intercept
    }
  } else
  {
    family <- "binomial"

    lasso_model = glmnet(x = x, y = y, family = family,
                         standardize = FALSE,
                         intercept  = TRUE, lambda = lambda)


    decreasing_order = order(lambda, decreasing = TRUE)

    nlambda = length(lambda)
    betalist = vector("list", nlambda)
    interceptlist = vector("list", nlambda)

    for (ind in 1:nlambda)
    {

      beta = lasso_model$beta[,ind]
      intercept = Ybar - sum(Xbar * beta)
      beta = beta / Xsd

      betalist[[decreasing_order[ind]]] = beta
      interceptlist[[decreasing_order[ind]]] = intercept
    }
  }



  return(list(interceptlist = interceptlist, betalist = betalist, lambda = lambda))

}
